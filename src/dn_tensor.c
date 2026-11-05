// Recursive tensor denoising of one patch against a GIVEN noise variance.  See
// dn_tensor.h for the method; the variance comes from the caller.
//
// LAYOUT.  The patch is a column-major tensor and every step is a reshape of a
// flat buffer, never a permutation: level k views the current tensor as a
// d_k x c matrix (rows fastest), and the reduced tensor it passes on, V*S, is a
// c x r matrix whose flat buffer IS level k+1's d_{k+1} x (c*r/d_{k+1})
// unfolding.  The reverse pass undoes it in the same buffers.
//
// SVD THROUGH THE SMALLER GRAM.  An unfolding A (m x c) is decomposed by
// diagonalising A A^T when m <= c and A^T A otherwise; singular values are the
// square roots of its eigenvalues and the other factor is one product away
// (A^T U = V S, or A V / s = U).  Squaring the condition number costs nothing
// at the noise edge, because the cut is applied to squared values anyway, but
// it does turn roundoff into tiny positive eigenvalues that an exact SVD would
// report as zero.  retained() never keeps a component below a numerical-rank
// floor for that reason; it can only bind at sigma2 == 0.
//
// THE CUT.  Level k keeps the components with s^2 > sigma2*(sqrt(m)+sqrt(c))^2,
// m x c being THAT level's unfolding -- the edge of the noise bulk, which is
// also exactly where the Frobenius shrinker the last level applies reaches
// zero, so every level discards what matrix-mode optshrink would.

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_eig.h"
#include "dn_tensor.h"

// Eigensolver workspaces kept per worker, by size.  A patch decomposes K
// unfoldings, and their sizes depend on the ranks found, so they vary from
// patch to patch; a handful of slots keeps the common sizes resident without
// holding one per size ever seen.
#define DN_TENSOR_EIG_SLOTS 6

// Numerical-rank floor, relative to the largest eigenvalue and per row of the
// Gram: well above roundoff in forming and diagonalising it, and far below any
// component a noise level from real data could keep.
#define DN_TENSOR_RANK_EPS (64.0 * DBL_EPSILON)

typedef struct {
	int n;
	unsigned long used;   // 0: empty
	dn_eig *e;
} eig_slot;

struct dn_tensor_work {
	dn_tensor_opts o;
	size_t n;              // elements in the patch
	int gmax;              // largest Gram any unfolding can need
	double *a, *t;         // n each: a centered copy; a transpose
	double *b[2];          // n each: the forward chain, alternately
	double *gram, *evecs;  // gmax^2
	double *evals, *f;     // gmax
	double *u[DN_TENSOR_MAX_MODES];      // retained left factors, dims[k] x r
	double *mean[DN_TENSOR_MAX_MODES];   // dims[k]
	int mrow[DN_TENSOR_MAX_MODES];       // 1: row means, 0: column means, -1: none
	eig_slot slot[DN_TENSOR_EIG_SLOTS];
	unsigned long tick, fallbacks;       // fallbacks of workspaces already freed
};

// min(d_k, n / d_k): the Gram size of mode k's unfolding of the full patch,
// which bounds every unfolding of that level.
static size_t mode_gram(const dn_tensor_opts *o, size_t n, int k) {
	const size_t d = (size_t)o->dims[k], rest = n / d;
	return d < rest ? d : rest;
}

static int check_opts(const dn_tensor_opts *o, size_t *n_out, int *gmax_out) {
	if (!o || o->nmodes < 3 || o->nmodes > DN_TENSOR_MAX_MODES) {
		dn_err("tensor denoising needs 3 to %d modes\n", DN_TENSOR_MAX_MODES);
		return 1;
	}
	size_t n = 1;
	for (int k = 0; k < o->nmodes; k++) {
		if (o->dims[k] < 2) {
			dn_err("tensor mode %d has size %d; every mode needs at least 2\n", k, o->dims[k]);
			return 1;
		}
		if (dn_mul_size(n, (size_t)o->dims[k], &n)) {
			dn_err("the tensor patch is too large for this build's address space\n");
			return 1;
		}
	}
	size_t gmax = 0;
	for (int k = 0; k < o->nmodes; k++) {
		const size_t g = mode_gram(o, n, k);
		if (g > gmax) gmax = g;
	}
	if (gmax > DN_MAX_VOL) {
		dn_err("a tensor unfolding needs a %zu x %zu eigensolve; the limit is %d\n",
		       gmax, gmax, DN_MAX_VOL);
		return 1;
	}
	*n_out = n;
	*gmax_out = (int)gmax;
	return 0;
}

int dn_tensor_check(const dn_tensor_opts *o) {
	size_t n;
	int gmax;
	return check_opts(o, &n, &gmax);
}

size_t dn_tensor_work_bytes(const dn_tensor_opts *o) {
	size_t n = 1, g2 = 0, uf = 0;
	for (int k = 0; k < o->nmodes; k++) n *= (size_t)o->dims[k];
	for (int k = 0; k < o->nmodes; k++) {
		const size_t g = mode_gram(o, n, k);
		if (g * g > g2) g2 = g * g;
		uf += (size_t)o->dims[k] * g;
	}
	// Four n-blocks, the Gram and its vectors, the left factors, and up to one
	// eigensolver per mode resident at a time (each holds two g^2 blocks).
	return (4 * n + 2 * g2 + uf + 2 * g2 * (size_t)o->nmodes) * sizeof(double);
}

dn_tensor_work *dn_tensor_work_create(const dn_tensor_opts *o) {
	size_t n;
	int gmax;
	if (check_opts(o, &n, &gmax)) return NULL;
	dn_tensor_work *w = (dn_tensor_work *)dn_calloc(1, sizeof(dn_tensor_work));
	if (!w) return NULL;
	w->o = *o;
	w->n = n;
	w->gmax = gmax;
	const size_t g = (size_t)gmax;
	w->a = (double *)dn_malloc(n, sizeof(double));
	w->t = (double *)dn_malloc(n, sizeof(double));
	w->b[0] = (double *)dn_malloc(n, sizeof(double));
	w->b[1] = (double *)dn_malloc(n, sizeof(double));
	w->gram = (double *)dn_malloc(g * g, sizeof(double));
	w->evecs = (double *)dn_malloc(g * g, sizeof(double));
	w->evals = (double *)dn_malloc(g, sizeof(double));
	w->f = (double *)dn_malloc(g, sizeof(double));
	int bad = !w->a || !w->t || !w->b[0] || !w->b[1] || !w->gram || !w->evecs ||
	          !w->evals || !w->f;
	for (int k = 0; !bad && k < o->nmodes; k++) {
		const size_t d = (size_t)o->dims[k];
		w->u[k] = (double *)dn_malloc(d * mode_gram(o, n, k), sizeof(double));
		w->mean[k] = (double *)dn_malloc(d, sizeof(double));
		bad = !w->u[k] || !w->mean[k];
	}
	if (bad) { dn_tensor_work_free(w); return NULL; }
	return w;
}

void dn_tensor_work_free(dn_tensor_work *w) {
	if (!w) return;
	free(w->a); free(w->t); free(w->b[0]); free(w->b[1]);
	free(w->gram); free(w->evecs); free(w->evals); free(w->f);
	for (int k = 0; k < DN_TENSOR_MAX_MODES; k++) { free(w->u[k]); free(w->mean[k]); }
	for (int s = 0; s < DN_TENSOR_EIG_SLOTS; s++) dn_eig_free(w->slot[s].e);
	free(w);
}

unsigned long dn_tensor_fallbacks(const dn_tensor_work *w) {
	if (!w) return 0ul;
	unsigned long fb = w->fallbacks;
	for (int s = 0; s < DN_TENSOR_EIG_SLOTS; s++) fb += dn_eig_fallbacks(w->slot[s].e);
	return fb;
}

// The workspace for a g x g solve: resident, or replacing the least recently
// used.  NULL on allocation failure.
static dn_eig *eig_for(dn_tensor_work *w, int g) {
	w->tick++;
	int lru = 0;
	for (int s = 0; s < DN_TENSOR_EIG_SLOTS; s++) {
		if (w->slot[s].e && w->slot[s].n == g) {
			w->slot[s].used = w->tick;
			return w->slot[s].e;
		}
		if (w->slot[s].used < w->slot[lru].used) lru = s;
	}
	eig_slot *s = &w->slot[lru];
	if (s->e) {
		w->fallbacks += dn_eig_fallbacks(s->e);
		dn_eig_free(s->e);
	}
	s->e = dn_eig_create(g);
	s->n = g;
	s->used = s->e ? w->tick : 0;
	return s->e;
}

// G = B^T B for B (len x g, column-major), into g x g full storage.  Each entry
// is one dot product of two contiguous columns, summed in a fixed order.
static void gram_cols(const double *B, size_t len, int g, double *G) {
	for (int k = 0; k < g; k++) {
		const double *bk = B + (size_t)k * len;
		for (int l = 0; l <= k; l++) {
			const double *bl = B + (size_t)l * len;
			double s = 0.0;
			for (size_t i = 0; i < len; i++) s += bk[i] * bl[i];
			G[(size_t)k * g + l] = G[(size_t)l * g + k] = s;
		}
	}
}

// Decompose the m x c matrix A through its smaller Gram.  On return w->evals
// holds the g = min(m, c) squared singular values, ASCENDING, and the returned
// workspace is staged for dn_eig_vectors.  NULL on failure.
static dn_eig *spectrum(dn_tensor_work *w, const double *A, int m, size_t c, int *g) {
	if ((size_t)m <= c) {
		// A A^T: transpose first so that every Gram entry is a contiguous dot
		// product rather than a stride-m walk.
		for (size_t j = 0; j < c; j++)
			for (int i = 0; i < m; i++) w->t[j + c * (size_t)i] = A[(size_t)i + (size_t)m * j];
		*g = m;
		gram_cols(w->t, c, m, w->gram);
	} else {
		*g = (int)c;
		gram_cols(A, (size_t)m, (int)c, w->gram);
	}
	dn_eig *e = eig_for(w, *g);
	if (!e || dn_eig_values(e, w->gram, w->evals) != 0) return NULL;
	return e;
}

// The top r components of the staged decomposition of A (m x c, Gram size g):
// the left factor U (m x r) into u, and Y = V*S, each column scaled by f[k]
// unless f is NULL, into y (c x r).  Column-major, largest first.  Returns
// non-zero if the eigenvectors could not be found.
//
// SIGNS.  Each component is returned with the sign that makes the largest-
// magnitude entry of its V column positive (the first, on a tie).  Without
// centring the sign is immaterial -- U and V*S flip together and the
// reconstruction is bit-identical either way, since (-a)(-b) == ab exactly --
// but with -tensor_demean it is NOT: the next level takes means across the
// reshaped V*S, where a flipped column flips a block of entries, so the result
// depends on which sign the SVD happened to return.  LAPACK defines no sign at
// all (measured on ME-PCA, which inherits it: one flip moves a 3x7x5 patch by
// 7% relative L2), so a demeaned run is reproducible only under a FIXED
// convention.  It also keeps the Accelerate and portable builds, whose
// eigenvectors may differ in sign, in agreement.  A column SUM would not do:
// centring rows makes every V column sum to zero, leaving the sign to roundoff.
static int factors(dn_tensor_work *w, dn_eig *e, const double *A, int m, size_t c, int g, int r,
                   const double *f, double *u, double *y) {
	if (r == 0) return 0;
	if (dn_eig_vectors(e, g - r, r, w->evecs) != 0) return 1;
	const double *ev = w->evecs;
	if ((size_t)m <= c) {
		// The eigenvectors are U; V*S = A^T U.
		for (int k = 0; k < r; k++) {
			const int col = g - 1 - k;
			double *uk = u + (size_t)m * k;
			for (int i = 0; i < m; i++) uk[i] = ev[(size_t)i * g + col];
			double *yk = y + c * (size_t)k, big = 0.0;   // V*S is V scaled by s > 0
			for (size_t j = 0; j < c; j++) {
				const double *aj = A + (size_t)m * j;
				double s = 0.0;
				for (int i = 0; i < m; i++) s += aj[i] * uk[i];
				yk[j] = s;
				if (fabs(s) > fabs(big)) big = s;
			}
			const double fk = (f ? f[k] : 1.0) * (big < 0.0 ? -1.0 : 1.0);
			if (big < 0.0)
				for (int i = 0; i < m; i++) uk[i] = -uk[i];
			for (size_t j = 0; j < c; j++) yk[j] *= fk;
		}
		return 0;
	}
	// The eigenvectors are V; U = A V / s.  r <= g = c here.
	for (int k = 0; k < r; k++) {
		const int col = g - 1 - k;
		double big = 0.0;
		for (size_t j = 0; j < c; j++)
			if (fabs(ev[j * (size_t)g + col]) > fabs(big)) big = ev[j * (size_t)g + col];
		const double sgn = big < 0.0 ? -1.0 : 1.0;
		const double sk = sqrt(fmax(w->evals[col], 0.0));
		const double fk = f ? f[k] : 1.0;
		double *yk = y + c * (size_t)k;
		for (size_t j = 0; j < c; j++) yk[j] = sgn * ev[j * (size_t)g + col] * sk * fk;
		double *uk = u + (size_t)m * k;
		for (int i = 0; i < m; i++) uk[i] = 0.0;
		for (size_t j = 0; j < c; j++) {
			const double v = sgn * ev[j * (size_t)g + col] / sk;
			const double *aj = A + (size_t)m * j;
			for (int i = 0; i < m; i++) uk[i] += aj[i] * v;
		}
	}
	return 0;
}

// The matrix level k decomposes: X itself, or X centered along its longer side
// into w->a, with the means kept for the reverse pass.
static const double *prepare(dn_tensor_work *w, const double *X, int m, size_t c, int k) {
	if (!w->o.subtract_mean) {
		w->mrow[k] = -1;
		return X;
	}
	double *mu = w->mean[k], *A = w->a;
	if ((size_t)m < c) {
		// More columns than rows: each row's mean over its columns.
		w->mrow[k] = 1;
		for (int i = 0; i < m; i++) mu[i] = 0.0;
		for (size_t j = 0; j < c; j++)
			for (int i = 0; i < m; i++) mu[i] += X[(size_t)i + (size_t)m * j];
		for (int i = 0; i < m; i++) mu[i] /= (double)c;
		for (size_t j = 0; j < c; j++)
			for (int i = 0; i < m; i++) A[(size_t)i + (size_t)m * j] = X[(size_t)i + (size_t)m * j] - mu[i];
	} else {
		w->mrow[k] = 0;
		for (size_t j = 0; j < c; j++) {
			const double *xj = X + (size_t)m * j;
			double s = 0.0;
			for (int i = 0; i < m; i++) s += xj[i];
			mu[j] = s / (double)m;
			for (int i = 0; i < m; i++) A[(size_t)i + (size_t)m * j] = xj[i] - mu[j];
		}
	}
	return A;
}

// Components above the bulk edge `cutoff`, never below the numerical-rank
// floor.  ev is ASCENDING, so the kept ones are a run from the top.
static int retained(const double *ev, int g, double cutoff) {
	const double top = ev[g - 1];
	const double floor = top > 0.0 ? top * g * DN_TENSOR_RANK_EPS : 0.0;
	int r = 0;
	while (r < g && ev[g - 1 - r] > cutoff && ev[g - 1 - r] > floor) r++;
	return r;
}

int dn_tensor_denoise(const dn_tensor_opts *o, dn_tensor_work *w, double *x,
                      double sigma2, int *ranks) {
	if (!(sigma2 >= 0.0) || !isfinite(sigma2)) return 1;
	const int K = o->nmodes;

	// Forward: cut each level at its own bulk edge, pass V*S on.  A level that
	// keeps nothing ends the recursion.
	int lm[DN_TENSOR_MAX_MODES], lr[DN_TENSOR_MAX_MODES];
	size_t lc[DN_TENSOR_MAX_MODES];
	const double *X = x;
	double *Y = NULL;
	int m = o->dims[0], levels = 0;
	size_t c = w->n / (size_t)m;
	for (int lev = 0; lev < K; lev++) {
		const double *A = prepare(w, X, m, c, lev);
		int g;
		dn_eig *e = spectrum(w, A, m, c, &g);
		if (!e) return 1;
		const double root = sqrt((double)m) + sqrt((double)c);
		const int r = retained(w->evals, g, sigma2 * (root * root));
		const double *f = NULL;
		if (lev == K - 1 && o->shrink_last && r > 0) {
			// Gavish & Donoho 2017's Frobenius-optimal shrinkage, in unnormalized
			// units, against the residual sizes left once r are kept.
			const double rows = fmax((double)(m - r), 1.0), cols = fmax((double)c - r, 1.0);
			for (int k = 0; k < r; k++) {
				const double s2 = w->evals[g - 1 - k];
				const double v = s2 - 2.0 * (cols + rows) * sigma2 +
				                 (cols - rows) * (cols - rows) * (sigma2 * sigma2) / s2;
				w->f[k] = sqrt(fmax(v, 0.0)) / sqrt(s2);
			}
			f = w->f;
		}
		Y = w->b[lev & 1];
		if (factors(w, e, A, m, c, g, r, f, w->u[lev], Y)) return 1;
		lm[lev] = m;
		lc[lev] = c;
		lr[lev] = r;
		levels = lev + 1;
		if (r == 0 || lev == K - 1) break;
		X = Y;
		c = c * (size_t)r / (size_t)o->dims[lev + 1];
		m = o->dims[lev + 1];
	}

	// Reverse: X_k = U_k Y^T + mean, where Y is level k+1's reconstruction
	// viewed as c_k x r_k.  A level that kept nothing is its mean alone.
	const double *Yc = Y;
	double *spare[2] = {w->a, w->t};
	for (int lev = levels - 1; lev >= 0; lev--) {
		double *out = lev == 0 ? x : spare[lev & 1];
		const int mm = lm[lev], rr = lr[lev];
		const size_t cc = lc[lev];
		const double *U = w->u[lev], *mu = w->mean[lev];
		const int mrow = w->mrow[lev];
		for (size_t j = 0; j < cc; j++) {
			double *oj = out + (size_t)mm * j;
			for (int i = 0; i < mm; i++) oj[i] = 0.0;
			for (int k = 0; k < rr; k++) {
				const double yk = Yc[j + cc * (size_t)k];
				const double *uk = U + (size_t)mm * k;
				for (int i = 0; i < mm; i++) oj[i] += uk[i] * yk;
			}
			if (mrow == 1)
				for (int i = 0; i < mm; i++) oj[i] += mu[i];
			else if (mrow == 0)
				for (int i = 0; i < mm; i++) oj[i] += mu[j];
		}
		Yc = out;
	}

	for (int k = 0; k < K; k++) ranks[k] = k < levels ? lr[k] : 0;
	return 0;
}
