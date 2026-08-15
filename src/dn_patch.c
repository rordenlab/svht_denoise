// patch geometry and the single-patch denoising kernel.
//
// The kernel implements Eq. (4) of Gavish & Donoho 2014 on the Casorati matrix
// of a cubic patch: M patch voxels by N volumes, threshold at
// omega(beta) * median singular value, keep the surviving components.
//
// Two shortcuts make this much cheaper than it looks, and neither changes the
// answer:
//
//  1. Hard thresholding is an orthogonal projection, Xhat = Y V_r V_r', so only
//     the RIGHT singular vectors are needed.  Those are the eigenvectors of the
//     N x N Gram matrix Y'Y, whose eigenvalues are the squared singular values.
//     So a 125 x 102 SVD becomes a 102 x 102 symmetric eigenproblem.
//  2. Only the CENTRE voxel of the patch is written out, so only one row of Xhat
//     is needed.  The projector V_r V_r' is never formed: the centre row is
//     pushed through the retained eigenvectors in two multiplies.
//
// Layout note: the patch is stored TRANSPOSED (volume-major, pt[j*m + i]) rather
// than in the Python's voxel-major order.  That makes each Gram entry a dot
// product of two contiguous length-M rows, and makes the gather read each input
// volume in a compact spatial neighbourhood instead of striding across volumes.

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_coef.h"
#include "dn_patch.h"

// Gram columns computed per pass over a row.  Swept 2, 4 and 8 in the full
// workload on 100x100x58x138 -- 1127, 795 and 821 seconds of CPU respectively,
// against 1221 unblocked.  Do not infer this from clang's vectoriser: at 4 and 8
// the cost model declines to vectorise the inner loop at all and only 2 keeps
// the vectorised form, which would pick the width that recovers 8% over the one
// that recovers 35%.
#define DN_GRAM_BLOCK 4

// Origin of the patch centred on `coord` along an axis of length `dim`, clamped
// so the patch stays inside the image: near an edge the patch is SHIFTED INWARDS
// rather than truncated or padded, so every voxel gets a full-size patch.
static int dn_patch_origin(int coord, int dim, int extent) {
	int o = coord - extent / 2;
	if (o < 0) o = 0;
	if (o > dim - extent) o = dim - extent;
	return o;
}

int dn_auto_extent(int nvol) {
	for (int k = 3; k <= DN_MAX_EXTENT; k += 2) {
		double cube = (double)k * k * k;   // double avoids overflow in the test
		if (cube > (double)nvol) return k;
	}
	return 0;
}

int dn_geom_init(dn_geom *g, int nx, int ny, int nz, int nvol, int extent) {
	if (!g) return 1;
	memset(g, 0, sizeof(*g));

	if (extent < 3 || extent > DN_MAX_EXTENT || (extent % 2) == 0) {
		dn_err("patch extent must be an odd number between 3 and %d (got %d)\n",
		       DN_MAX_EXTENT, extent);
		return 1;
	}
	if (extent > nx || extent > ny || extent > nz) {
		dn_err("patch extent %d exceeds an image dimension (%dx%dx%d)\n",
		       extent, nx, ny, nz);
		return 1;
	}
	const int m = extent * extent * extent;
	if (m <= nvol) {
		dn_err("patch extent %d gives %d voxels for %d volumes; this build requires\n",
		       extent, m, nvol);
		dn_err("  more voxels than volumes (M > N). Use a larger -extent: the default\n");
		dn_err("  rule is the smallest odd k with k^3 > N, which here is %d.\n",
		       dn_auto_extent(nvol));
		return 1;
	}

	g->nx = nx; g->ny = ny; g->nz = nz;
	g->nvol = nvol;
	g->extent = extent;
	g->m = m;
	g->nvox3d = (size_t)nx * (size_t)ny * (size_t)nz;
	g->beta = (double)nvol / (double)m;

	// dn_work's offset table is int32_t, and the largest entry it must hold is
	// (extent-1)*(nx*ny + nx + 1).  Nothing else rejects a slice big enough to
	// overflow that: dn_image_read checks each dimension against INT32_MAX and
	// the total against SIZE_MAX, neither of which bounds nx*ny.  A truncated
	// offset would go NEGATIVE and read out of bounds in the gather, so refuse
	// the geometry instead.  Unreachable with any real image -- it needs a single
	// slice of over a billion voxels -- but the check is two lines and the
	// alternative failure is silent.
	const double max_off = (double)(extent - 1) *
	                       ((double)nx * (double)ny + (double)nx + 1.0);
	if (max_off > 2147483647.0) {
		dn_err("image slice %dx%d is too large for a %d-voxel patch offset table\n",
		       nx, ny, extent);
		return 1;
	}

	if (!dn_beta_valid(g->beta)) {
		dn_err("internal: aspect ratio %g is out of range\n", g->beta);
		return 1;
	}
	const double mu = dn_mp_median(g->beta);
	g->omega = dn_lambda_star(g->beta) / sqrt(mu);
	if (!(mu > 0.0) || !isfinite(mu) || !isfinite(g->omega) || g->omega <= 0.0) {
		dn_err("internal: threshold coefficient failed for beta=%g\n", g->beta);
		return 1;
	}
	// sigma_hat = median_s / sqrt(max(M,N) * mu); M > N is guaranteed above.
	g->sigma_scale = 1.0 / sqrt((double)m * mu);
	return 0;
}

// ---------------------------------------------------------------------------
// worker scratch
// ---------------------------------------------------------------------------

dn_work *dn_work_create(const dn_geom *g) {
	if (!g) return NULL;
	dn_work *w = (dn_work *)dn_calloc(1, sizeof(dn_work));
	if (!w) return NULL;
	const int m = g->m, n = g->nvol;
	w->pt = (float *)dn_malloc((size_t)n * m, sizeof(float));
	w->gram = (double *)dn_malloc((size_t)n * n, sizeof(double));
	w->evecs = (double *)dn_malloc((size_t)n * n, sizeof(double));
	w->evals = (double *)dn_malloc((size_t)n, sizeof(double));
	w->s = (double *)dn_malloc((size_t)n, sizeof(double));
	w->yc = (double *)dn_malloc((size_t)n, sizeof(double));
	w->acc = (double *)dn_malloc((size_t)n, sizeof(double));
	w->offset = (int32_t *)dn_malloc((size_t)m, sizeof(int32_t));
	w->eig = dn_eig_create(n);
	if (!w->pt || !w->gram || !w->evecs || !w->evals || !w->s || !w->yc || !w->acc ||
	    !w->offset || !w->eig) {
		dn_work_free(w);
		return NULL;
	}

	// Flat spatial offsets of the patch voxels from the patch origin.  Constant
	// for the whole run because every patch is full size.
	//
	// NIfTI orders a volume with X FASTEST: index = x + y*nx + z*nx*ny.  That is
	// the opposite of a C-order flattening, and getting it
	// backwards produces a plausible-looking image denoised along the wrong axes.
	// dx is innermost here so the gather walks contiguous memory.
	const int k = g->extent;
	int t = 0;
	for (int dz = 0; dz < k; dz++)
		for (int dy = 0; dy < k; dy++)
			for (int dx = 0; dx < k; dx++)
				w->offset[t++] = (int32_t)((size_t)dz * g->nx * g->ny + (size_t)dy * g->nx + dx);
	return w;
}

void dn_work_free(dn_work *w) {
	if (!w) return;
	free(w->pt);
	free(w->gram);
	free(w->evecs);
	free(w->evals);
	free(w->s);
	free(w->yc);
	free(w->acc);
	free(w->offset);
	dn_eig_free(w->eig);
	free(w);
}

// ---------------------------------------------------------------------------
// thresholding
// ---------------------------------------------------------------------------

// Median of an ASCENDING array, using the mean of the middle pair for even n --
// the same convention as MATLAB's median() and NumPy's np.median().
static double dn_median_ascending(const double *s, int n) {
	if (n < 1) return 0.0;
	if (n % 2) return s[n / 2];
	return 0.5 * (s[n / 2 - 1] + s[n / 2]);
}

// Number of retained components for ASCENDING singular values.
//
// The rule is s >= omega * median(s), with one refinement: a component whose
// singular value is exactly zero is never counted.  When the threshold is
// positive the two are identical, because s >= threshold > 0 already implies
// s > 0.  They differ only when the median is zero, i.e. when at least half the
// spectrum has collapsed -- and there the extra components are null directions
// of Y, which contribute EXACTLY nothing to the projection (y_c . v_k is an
// entry of Y v_k = 0).  So this changes no denoised value anywhere; it only
// stops the rank map reading N over empty background, and it makes the all-zero
// all-zero patch fall out on its own: every s is zero, so the rank is 0,
// the projection is empty, and sigma is 0.
static int dn_retained_count(const double *s, int n, double omega, double *median_out) {
	const double med = dn_median_ascending(s, n);
	if (median_out) *median_out = med;
	const double thresh = omega * med;
	int r = 0;
	for (int i = n - 1; i >= 0; i--) {
		if (s[i] > 0.0 && s[i] >= thresh) r++;
		else break;   // ascending, so everything below also fails
	}
	return r;
}

// ---------------------------------------------------------------------------
// the kernel
// ---------------------------------------------------------------------------

// Gather, and report whether anything in the patch was non-zero.  The flag is
// accumulated here rather than by a second pass over w->pt: the gather already
// touches every element, and at 125x102 a second pass is 102 kB of pointless
// traffic on the majority of patches (these datasets arrive brain-masked, so
// most patches are empty).
static int patch_gather_nonzero(const dn_geom *g, dn_work *w, const float *img,
                                int ix, int iy, int iz, int *centre_out) {
	const int m = g->m, n = g->nvol;
	const size_t nvox3d = g->nvox3d;

	const int ox = dn_patch_origin(ix, g->nx, g->extent);
	const int oy = dn_patch_origin(iy, g->ny, g->extent);
	const int oz = dn_patch_origin(iz, g->nz, g->extent);
	const size_t base = (size_t)ox + (size_t)oy * g->nx + (size_t)oz * g->nx * g->ny;

	// Row of the patch matrix holding the voxel being denoised, in the same
	// (dz, dy, dx) enumeration order used to build w->offset.
	const int cx = ix - ox, cy = iy - oy, cz = iz - oz;
	if (centre_out) *centre_out = (cz * g->extent + cy) * g->extent + cx;

	// Gather volume by volume.  Each volume is read over a compact
	// neighbourhood; the outer stride over volumes is paid m times, not m*n.
	int nonzero = 0;
	for (int j = 0; j < n; j++) {
		const float *vol = img + (size_t)j * nvox3d + base;
		float *dst = w->pt + (size_t)j * m;
		for (int i = 0; i < m; i++) {
			const float v = vol[w->offset[i]];
			dst[i] = v;
			nonzero |= (v != 0.0f);
		}
	}
	return nonzero;
}

int dn_denoise_voxel(const dn_geom *g, dn_work *w, const float *img,
                     int ix, int iy, int iz,
                     float *out, float *sigma, uint16_t *rank) {
	const int m = g->m, n = g->nvol;

	// 1. Gather.  An entirely zero patch has an entirely zero spectrum, so the
	//    median, the threshold, the rank and the projection are all zero --
	//    so the shortcut is behaviour-preserving (verified byte-identical) and
	//    worth doing: these datasets arrive brain-masked, so most of the volume
	//    is empty background, and each such patch was otherwise paying a full
	//    Gram and eigensolve to arrive at zero.
	int centre = 0;
	if (!patch_gather_nonzero(g, w, img, ix, iy, iz, &centre)) {
		for (int j = 0; j < n; j++) out[j] = 0.0f;
		if (sigma) *sigma = 0.0f;
		if (rank) *rank = 0;
		return 0;
	}

	// 2. Gram = Y' Y, upper triangle then mirrored.
	// Four independent accumulators: a single running sum is a serial dependency
	// chain through the FMA latency, which stops the compiler vectorising and
	// leaves most of the unit idle.  The split is a FIXED partition, so the
	// summation order is still identical on every run and at every thread count.
	//
	// Blocked over j on top of that: DN_GRAM_BLOCK columns share one pass over
	// row i, so row i is read once per block instead of once per j.  That is
	// where the time went -- blocking alone is 35% of the whole run at 138x343,
	// while the float rows it feeds on are worth 11% with it and nothing without.
	//
	// Neither is a numerics change.  The rows widen back to double here, and
	// since they were copied from a float32 image the widened value IS the value
	// a double array would have held; blocking only interleaves independent dot
	// products, leaving the four accumulators and the fixed partition inside each
	// one untouched.  So the output is bit-identical -- asserted against a
	// preserved binary on 80 M voxels, not assumed.
	for (int i = 0; i < n; i++) {
		const float *a = w->pt + (size_t)i * m;
		int j = i;
		for (; j + DN_GRAM_BLOCK <= n; j += DN_GRAM_BLOCK) {
			const float *b[DN_GRAM_BLOCK];
			double s[DN_GRAM_BLOCK][4] = {{0.0}};
			for (int q = 0; q < DN_GRAM_BLOCK; q++) b[q] = w->pt + (size_t)(j + q) * m;
			int v = 0;
			for (; v + 3 < m; v += 4) {
				const double x0 = a[v], x1 = a[v + 1], x2 = a[v + 2], x3 = a[v + 3];
				for (int q = 0; q < DN_GRAM_BLOCK; q++) {
					s[q][0] += x0 * (double)b[q][v];
					s[q][1] += x1 * (double)b[q][v + 1];
					s[q][2] += x2 * (double)b[q][v + 2];
					s[q][3] += x3 * (double)b[q][v + 3];
				}
			}
			for (; v < m; v++) {
				const double x = a[v];
				for (int q = 0; q < DN_GRAM_BLOCK; q++) s[q][0] += x * (double)b[q][v];
			}
			for (int q = 0; q < DN_GRAM_BLOCK; q++) {
				const double acc = (s[q][0] + s[q][1]) + (s[q][2] + s[q][3]);
				w->gram[(size_t)i * n + j + q] = acc;
				w->gram[(size_t)(j + q) * n + i] = acc;
			}
		}
		for (; j < n; j++) {
			const float *b = w->pt + (size_t)j * m;
			double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
			int v = 0;
			for (; v + 3 < m; v += 4) {
				a0 += (double)a[v] * (double)b[v];
				a1 += (double)a[v + 1] * (double)b[v + 1];
				a2 += (double)a[v + 2] * (double)b[v + 2];
				a3 += (double)a[v + 3] * (double)b[v + 3];
			}
			for (; v < m; v++) a0 += (double)a[v] * (double)b[v];
			const double acc = (a0 + a1) + (a2 + a3);
			w->gram[(size_t)i * n + j] = acc;
			w->gram[(size_t)j * n + i] = acc;
		}
	}

	// 3. Eigenvalues of the Gram are the squared singular values.
	//
	// Eigenvalues indistinguishable from zero are flushed to zero rather than
	// carried as tiny positives.  Forming Y'Y squares the condition number, so
	// the smallest computed eigenvalues of a rank-deficient patch are pure
	// round-off whose sign and magnitude depend on the solver: keeping them made
	// the RANK MAP differ between two independent eigensolvers at 4% of voxels
	// (all in near-empty background) even though the
	// denoised values agreed to 1e-12, because those directions are null and
	// contribute nothing to the projection.  The cut is the standard numerical
	// rank tolerance, relative to the largest eigenvalue.
	// Only the eigenvalues are needed to reach the threshold decision; the
	// eigenvectors are fetched below, and only the ones that survive it.
	int rc = dn_eig_values(w->eig, w->gram, w->evals);
	if (rc) return rc;
	const double lam_max = w->evals[n - 1];   // ascending
	const double lam_tol = (lam_max > 0.0) ? lam_max * (double)n * DBL_EPSILON : 0.0;
	for (int i = 0; i < n; i++)
		w->s[i] = (w->evals[i] > lam_tol) ? sqrt(w->evals[i]) : 0.0;

	// 4. Threshold.
	double med = 0.0;
	const int r = dn_retained_count(w->s, n, g->omega, &med);

	if (sigma) *sigma = (float)(med * g->sigma_scale);
	if (rank) *rank = (uint16_t)r;

	// 5. Project the centre row onto the retained subspace, without ever
	//    forming the N x N projector.  The accumulator stays in double and is
	//    rounded to float32 exactly once, at the end.
	if (r <= 0) {
		for (int j = 0; j < n; j++) out[j] = 0.0f;
		return 0;
	}

	// Now, and only now, do we know which eigenvectors matter.
	rc = dn_eig_vectors(w->eig, n - r, r, w->evecs);
	if (rc) return rc;

	double *yc = w->yc, *acc = w->acc;
	for (int j = 0; j < n; j++) {
		yc[j] = (double)w->pt[(size_t)j * m + centre];
		acc[j] = 0.0;
	}
	const double *v = w->evecs;
	for (int k = n - r; k < n; k++) {
		double dot = 0.0;
		for (int j = 0; j < n; j++) dot += yc[j] * v[(size_t)j * n + k];
		for (int j = 0; j < n; j++) acc[j] += dot * v[(size_t)j * n + k];
	}
	for (int j = 0; j < n; j++) out[j] = (float)acc[j];
	return 0;
}
