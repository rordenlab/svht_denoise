// Sliding-window tensor denoising over a multi-echo series.
//
// The window rules are ME-PCA's (its denoise_recursive_tensor), so that, given
// the same noise variance, the result can be checked against that package with
// tools/tensor_parity.py.  They differ from the matrix denoiser's on purpose:
//   - a window sits at every corner where it fits INSIDE the image, on a stride
//     grid from the origin.  Nothing is shifted inwards at an edge, so a stride
//     that does not land on the far edge leaves voxels no window covers;
//   - a window is denoised when it touches the mask, and then writes every voxel
//     it covers, in the mask or not;
//   - a voxel is the unweighted mean of every estimate written to it.  One that
//     no window covered keeps its input value and reads NaN in every map.
//
// Determinism, as dn_agg.c: window corners are split into tiles one window wide
// (rounded up to whole strides) and coloured 2x2x2.  A window reaches w-1
// voxels past its corner and same-colour tiles start two tiles apart, so they
// never write the same voxel; colours run in a fixed order and each tile's
// windows in raster order, so every voxel receives its contributions in the
// same order at any thread count.

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_run.h"
#include "dn_tensor.h"

typedef struct {
	const dn_trun *r;
	pthread_mutex_t lock;
	const size_t *tiles;   // this colour's tiles
	size_t ntiles, next;
	int T[3], nt[3], npos[3];
	size_t st[3];          // patch-tensor stride of DN_TMODE_SPACE, _TIME, _ECHO
	int m;                 // voxels per window
	size_t *off;           // m voxel offsets from the window corner, x fastest
	int *rows;             // window rows a window WRITES: all m, or its centre
	int nrows;
	int *crow;             // window rows of the centre voxel, or column
	int ncrow;
	uint32_t *count;       // nvox3d: estimates received so far
	double *sacc, *racc;   // running means of sigma2, and of each mode's rank
	int failed;
	unsigned long eig_fail, fallbacks;
} tpool;

static int tile_colour(size_t t, const int nt[3]) {
	const size_t tx = t % (size_t)nt[0], ty = (t / (size_t)nt[0]) % (size_t)nt[1];
	const size_t tz = t / ((size_t)nt[0] * nt[1]);
	return (int)((tx & 1) | (ty & 1) << 1 | (tz & 1) << 2);
}

// Corner positions, tile width (whole strides, at least one window) and tile
// count per axis.
static void tile_grid(const dn_trun *r, int T[3], int nt[3], int npos[3]) {
	const int dims[3] = {r->nx, r->ny, r->nz};
	for (int a = 0; a < 3; a++) {
		npos[a] = dims[a] - r->win[a] + 1;
		T[a] = (r->win[a] + r->stride[a] - 1) / r->stride[a] * r->stride[a];
		nt[a] = (npos[a] + T[a] - 1) / T[a];
	}
}

int dn_tensor_setup(dn_trun *r) {
	const int dims[3] = {r->nx, r->ny, r->nz};
	for (int a = 0; a < 3; a++) {
		if (r->win[a] < 1 || r->win[a] > dims[a]) {
			dn_err("the tensor window (%d x %d x %d) does not fit in the image (%d x %d x %d)\n",
			       r->win[0], r->win[1], r->win[2], r->nx, r->ny, r->nz);
			return 1;
		}
		if (r->stride[a] < 1) {
			dn_err("the tensor stride must be positive\n");
			return 1;
		}
		// A stride past the last corner leaves corner 0 alone, as the stride to
		// that corner does.  Clamping keeps the tile and corner arithmetic in
		// range for any stride the parser accepts (up to INT32_MAX).
		const int npos = dims[a] - r->win[a] + 1;
		if (r->stride[a] > npos) r->stride[a] = npos;
	}
	if (r->noise_var == NULL && !(r->sigma2 >= 0.0 && isfinite(r->sigma2))) {
		dn_err("the noise variance must be finite and non-negative (got %g)\n", r->sigma2);
		return 1;
	}
	const int size[3] = {r->win[0] * r->win[1] * r->win[2], r->nt, r->ne};
	r->k.nmodes = 3;
	for (int k = 0; k < 3; k++) r->k.dims[k] = size[r->order[k]];
	return dn_tensor_check(&r->k);
}

static size_t patch_size(const dn_trun *r) {
	size_t n = 1;
	for (int k = 0; k < r->k.nmodes; k++) n *= (size_t)r->k.dims[k];
	return n;
}

int dn_tensor_threads(const dn_trun *r, int requested) {
	int T[3], nt[3], npos[3], count[8] = {0}, maxper = 0;
	tile_grid(r, T, nt, npos);
	const size_t ntile = (size_t)nt[0] * nt[1] * nt[2];
	for (size_t t = 0; t < ntile; t++) count[tile_colour(t, nt)]++;
	for (int c = 0; c < 8; c++) if (count[c] > maxper) maxper = count[c];
	const size_t per = dn_tensor_work_bytes(&r->k) + patch_size(r) * sizeof(double);
	return dn_cap_threads(requested, (size_t)maxper, per);
}

// Does any of the window rows `idx` (all m when NULL) lie in the mask?
static int touches(const tpool *p, size_t base, const int *idx, int n) {
	for (int i = 0; i < n; i++)
		if (p->r->mask[base + (size_t)p->off[idx ? idx[i] : i]]) return 1;
	return 0;
}

static void *tworker(void *arg) {
	tpool *p = (tpool *)arg;
	const dn_trun *r = p->r;
	const int nx = r->nx, ny = r->ny, nt = r->nt, ne = r->ne, K = r->k.nmodes;
	const size_t nv = (size_t)nx * ny * r->nz;
	const size_t sS = p->st[DN_TMODE_SPACE], sT = p->st[DN_TMODE_TIME], sE = p->st[DN_TMODE_ECHO];

	dn_tensor_work *w = dn_tensor_work_create(&r->k);
	double *x = (double *)dn_malloc(patch_size(r), sizeof(double));
	double *q = (double *)dn_malloc((size_t)p->nrows, sizeof(double));
	if (!w || !x || !q) {
		pthread_mutex_lock(&p->lock);
		p->failed = 1;
		pthread_mutex_unlock(&p->lock);
		dn_tensor_work_free(w); free(x); free(q);
		return NULL;
	}
	unsigned long fail = 0;

	for (;;) {
		pthread_mutex_lock(&p->lock);
		if (p->failed || p->next >= p->ntiles) {
			pthread_mutex_unlock(&p->lock);
			break;
		}
		const size_t t = p->tiles[p->next++];
		pthread_mutex_unlock(&p->lock);

		const int tc[3] = {(int)(t % (size_t)p->nt[0]), (int)((t / (size_t)p->nt[0]) % (size_t)p->nt[1]),
		                   (int)(t / ((size_t)p->nt[0] * p->nt[1]))};
		int lo[3], hi[3];
		for (int a = 0; a < 3; a++) {
			lo[a] = tc[a] * p->T[a];
			hi[a] = lo[a] + p->T[a] < p->npos[a] ? lo[a] + p->T[a] : p->npos[a];
		}
		for (int cz = lo[2]; cz < hi[2]; cz += r->stride[2])
			for (int cy = lo[1]; cy < hi[1]; cy += r->stride[1])
				for (int cx = lo[0]; cx < hi[0]; cx += r->stride[0]) {
					const size_t base = (size_t)cx + (size_t)nx * ((size_t)cy + (size_t)ny * cz);
					if (r->mask) {
						if (!touches(p, base, NULL, p->m)) continue;
						if (r->centre_assign && !touches(p, base, p->crow, p->ncrow)) continue;
					}
					for (int e = 0; e < ne; e++)
						for (int j = 0; j < nt; j++) {
							const float *src = r->echo[e] + nv * (size_t)j + base;
							double *dst = x + sT * (size_t)j + sE * (size_t)e;
							for (int i = 0; i < p->m; i++) dst[sS * (size_t)i] = (double)src[p->off[i]];
						}
					// The window's noise variance: the map's mean over the voxels it
					// covers, in a fixed order, or the one known value.
					double sigma2 = r->sigma2;
					if (r->noise_var) {
						double s = 0.0;
						for (int i = 0; i < p->m; i++) s += (double)r->noise_var[base + p->off[i]];
						sigma2 = s / (double)p->m;
					}
					int ranks[DN_TENSOR_MAX_MODES];
					if (dn_tensor_denoise(&r->k, w, x, sigma2, ranks)) {
						fail++;
						continue;
					}
					// Running means: an estimate enters with coefficient 1/count, the
					// first with exactly 1 (see AGENTS.md on averaging).
					for (int i = 0; i < p->nrows; i++) {
						const size_t v = base + (size_t)p->off[p->rows[i]];
						q[i] = 1.0 / (double)++p->count[v];
						if (p->sacc) p->sacc[v] += q[i] * (sigma2 - p->sacc[v]);
						if (p->racc)
							for (int k = 0; k < K; k++) {
								double *rk = p->racc + (size_t)k * nv + v;
								*rk += q[i] * ((double)ranks[k] - *rk);
							}
					}
					for (int e = 0; e < ne; e++)
						for (int j = 0; j < nt; j++) {
							float *o = r->out[e] + nv * (size_t)j + base;
							const double *xr = x + sT * (size_t)j + sE * (size_t)e;
							for (int i = 0; i < p->nrows; i++) {
								const int row = p->rows[i];
								float *oi = o + p->off[row];
								*oi = (float)((double)*oi + q[i] * (xr[sS * (size_t)row] - (double)*oi));
							}
						}
				}
	}
	const unsigned long fb = dn_tensor_fallbacks(w);
	dn_tensor_work_free(w);
	free(x);
	free(q);
	if (fail || fb) {
		pthread_mutex_lock(&p->lock);
		p->eig_fail += fail;
		p->fallbacks += fb;
		pthread_mutex_unlock(&p->lock);
	}
	return NULL;
}

int dn_tensor_execute(const dn_trun *r) {
	const int nx = r->nx, ny = r->ny, K = r->k.nmodes;
	const size_t nv = (size_t)nx * ny * r->nz;
	int status = 1;
	size_t *tiles = NULL;
	tpool p;
	memset(&p, 0, sizeof(p));
	p.r = r;
	tile_grid(r, p.T, p.nt, p.npos);

	// Patch-tensor strides: mode k of the tensor is DN_TMODE_* order[k], and the
	// tensor is column-major in that order.
	size_t s = 1;
	for (int k = 0; k < K; k++) {
		p.st[r->order[k]] = s;
		s *= (size_t)r->k.dims[k];
	}

	// Window rows, x fastest then y then z -- the reference's spatial index.  The
	// centre is ceil(w/2) - 1 per axis; a two-element window spans every slice,
	// so its centre is the whole column.
	p.m = r->win[0] * r->win[1] * r->win[2];
	p.off = (size_t *)dn_malloc((size_t)p.m, sizeof(size_t));
	p.rows = (int *)dn_malloc((size_t)p.m, sizeof(int));
	p.crow = (int *)dn_malloc((size_t)r->win[2], sizeof(int));
	p.count = (uint32_t *)dn_calloc(nv, sizeof(uint32_t));
	if (r->sigma2_out || r->gain) p.sacc = (double *)dn_calloc(nv, sizeof(double));
	if (r->rank || r->gain) p.racc = (double *)dn_calloc(nv * (size_t)K, sizeof(double));
	if (!p.off || !p.rows || !p.crow || !p.count || ((r->sigma2_out || r->gain) && !p.sacc) ||
	    ((r->rank || r->gain) && !p.racc)) goto done;
	const int ccx = (r->win[0] + 1) / 2 - 1, ccy = (r->win[1] + 1) / 2 - 1, ccz = (r->win[2] + 1) / 2 - 1;
	for (int c = 0, i = 0; c < r->win[2]; c++)
		for (int b = 0; b < r->win[1]; b++)
			for (int a = 0; a < r->win[0]; a++, i++) {
				p.off[i] = (size_t)a + (size_t)nx * ((size_t)b + (size_t)ny * (size_t)c);
				if (a == ccx && b == ccy && (r->full_z || c == ccz)) p.crow[p.ncrow++] = i;
			}
	if (r->centre_assign) {
		for (int i = 0; i < p.ncrow; i++) p.rows[i] = p.crow[i];
		p.nrows = p.ncrow;
	} else {
		for (int i = 0; i < p.m; i++) p.rows[i] = i;
		p.nrows = p.m;
	}

	const size_t ntile = (size_t)p.nt[0] * p.nt[1] * p.nt[2];
	tiles = (size_t *)dn_malloc(ntile, sizeof(size_t));
	if (!tiles) goto done;
	if (pthread_mutex_init(&p.lock, NULL) != 0) {
		dn_err("cannot create the worker mutex\n");
		goto done;
	}
	const int nthreads = dn_tensor_threads(r, r->nthreads);
	for (int colour = 0; colour < 8 && !p.failed; colour++) {
		size_t k = 0;
		for (size_t t = 0; t < ntile; t++)
			if (tile_colour(t, p.nt) == colour) tiles[k++] = t;
		p.tiles = tiles;
		p.ntiles = k;
		p.next = 0;
		if (k) dn_thread_run(tworker, &p, (size_t)nthreads < k ? nthreads : (int)k);
	}
	pthread_mutex_destroy(&p.lock);
	if (dn_pool_report(p.failed, p.eig_fail, p.fallbacks)) goto done;

	// Uncovered voxels keep their input and read NaN in every map.  The SNR gain
	// is the reference's heuristic, computed from the AVERAGED ranks, which is
	// not the mean of each window's own gain; with nothing kept anywhere it has
	// no meaning and reads NaN (the reference divides by zero).
	double nel = 1.0;
	for (int k = 0; k < K; k++) nel *= r->k.dims[k];
	for (size_t v = 0; v < nv; v++) {
		if (!p.count[v]) {
			for (int e = 0; e < r->ne; e++)
				for (int j = 0; j < r->nt; j++) r->out[e][v + nv * (size_t)j] = r->echo[e][v + nv * (size_t)j];
			if (r->sigma2_out) r->sigma2_out[v] = NAN;
			if (r->rank)
				for (int k = 0; k < K; k++) r->rank[(size_t)k * nv + v] = NAN;
			if (r->gain) r->gain[v] = NAN;
			continue;
		}
		if (r->sigma2_out) r->sigma2_out[v] = (float)p.sacc[v];
		if (r->rank)
			for (int k = 0; k < K; k++) r->rank[(size_t)k * nv + v] = (float)p.racc[(size_t)k * nv + v];
		if (r->gain) {
			double prod = 1.0, sum = 0.0;
			for (int k = 0; k < K; k++) {
				const double pk = p.racc[(size_t)k * nv + v];
				prod *= pk;
				sum += (r->k.dims[k] - pk) * pk;
			}
			r->gain[v] = prod + sum > 0.0 ? (float)sqrt(nel / (prod + sum)) : NAN;
		}
	}
	status = 0;

done:
	free(tiles);
	free(p.off);
	free(p.rows);
	free(p.crow);
	free(p.count);
	free(p.sacc);
	free(p.racc);
	return status;
}
