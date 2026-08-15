/* Gibbs ringing removal by the local subvoxel-shift method of
 *   Kellner E, Dhital B, Kiselev VG, Reisert M.  Gibbs-ringing artifact removal
 *   based on local subvoxel-shifts.  Magn Reson Med 2016; 76:1574-1581.
 *
 * Adapted from MRtrix3's cmd/mrdegibbs.cpp, by Ben Jeurissen and J-Donald
 * Tournier.  The arithmetic is a deliberate transcription -- same accumulation
 * order in the TV window, the same iteratively-built ramp -- so that a
 * divergence from the reference means something.  Two groupings do differ:
 * the spectral split multiplies before dividing where upstream scales in place,
 * and the inverse transform's 1/n rides on the ramp rather than being applied
 * after.  Both were rebuilt with upstream's exact grouping and compared: float32
 * output is IDENTICAL at 9x9, 13x11, 104x104, 128x128 and 140x92.
 *
 * What is NOT carried over is Eigen itself, or the -axes option; see the notes
 * below.  The transform IS carried over, from the kissfft backend Eigen::FFT
 * uses -- adapted into dg_fft.c beside this file, which see.
 *
 * Copyright (c) 2008-2025 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

/* WITHIN-SLICE AXES ARE FIXED AT x,y, and for NIfTI input that is not a
 * restriction.  mrdegibbs picks them from the "SliceEncodingDirection" keyval,
 * which MRtrix populates only from DICOM or an imported JSON -- its NIfTI reader
 * never sets it -- so mrdegibbs on a .nii/.nii.gz always uses {0,1} too.  A BIDS
 * sidecar sitting next to the file changes nothing: nothing reads it.
 *
 * THIS FILE DOES NOT USE ACCELERATE, and that is a measured decision rather than
 * a portability one -- the rest of the tool does use it, for the eigensolver.
 *
 * vDSP is not an option either: its DFT accepts only lengths f*2^n with
 * f in {1,3,5,15}, and 104, 140, 92, 100, 58, 138 and 297 all return NULL from
 * vDSP_DFT_zop_CreateSetupD.  Hence dg_fft.c, which carries its own transform.
 * A BLAS-based alternative was built and measured and lost; AGENTS.md records
 * why, so that it does not get rebuilt.
 *
 * The 41 subvoxel shifts are what the reference makes them: a phase ramp applied
 * to the line's spectrum before each inverse transform.
 */

#include <complex.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_run.h"       // dn_default_threads, for the thread cap
#include "dg.h"
#include "mrdegibbs/dg_fft.h"

// mrdegibbs's defaults, and the only values anyone uses.  They are constants
// rather than options because every one of them would have to be threaded
// through the axis setup, where the ramps are built.
#define DG_NSH   20
#define DG_MINW  1
#define DG_MAXW  3
#define DG_NSHIFT (2*DG_NSH + 1)

// Smallest transform length the TV window arithmetic is defined for: the
// reference indexes (l - maxW - 1 + n) and (n - maxW - 1), which need n >= 5 to
// stay non-negative before the modulo.  mrdegibbs itself guards nothing, so this
// is the arithmetic bound and not a judgement about what is clinically useful.
#define DG_MIN_DIM (DG_MAXW + 2)

// Written out rather than taken from M_PI, which is not in C99's <math.h>.
// Same value as MRtrix's Math::pi.
#define DG_PI 3.14159265358979323846

typedef double _Complex dcx;

// shifts[] from the reference: 0, then +1..+nsh, then -1..-nsh.
static int dg_shift(int j) { return (j <= DG_NSH) ? j : -(j - DG_NSH); }

// Everything that depends only on a transform LENGTH, built once and shared by
// every worker.  All of it is small: the ramps are 41*n complex doubles, 92 KiB
// at n = 140.
typedef struct {
	int n;
	double *cw;    // n,      (1 + cos(2 pi L/n)) / 2
	dg_fft fft, ifft;
	dcx *ramp;     // 41 * n, with the inverse transform's 1/n folded in
} dg_axis;

static void axis_free(dg_axis *a) {
	free(a->cw);
	free(a->ramp);
	dg_fft_free(&a->fft);
	dg_fft_free(&a->ifft);
	memset(a, 0, sizeof(*a));
}

static int axis_init(dg_axis *a, int n) {
	memset(a, 0, sizeof(*a));
	a->n  = n;
	a->cw = (double *)dn_calloc((size_t)n, sizeof(double));
	dcx *ramp = (dcx *)dn_calloc((size_t)DG_NSHIFT*n, sizeof(dcx));
	if (!a->cw || !ramp) { free(ramp); axis_free(a); return 1; }

	for (int L = 0; L < n; L++)
		a->cw[L] = (1.0 + cos(2.0*DG_PI*(double)L/(double)n))*0.5;

	if (dg_fft_init(&a->fft, n, 0) || dg_fft_init(&a->ifft, n, 1)) {
		free(ramp); axis_free(a); return 1;
	}
	const int maxn = (n & 1) ? (n-1)/2 : n/2 - 1;
	for (int j = 0; j < DG_NSHIFT; j++) {
		dcx *r = ramp + (size_t)j*n;
		if (j == 0) {
			// Column 0 of `shifted` is the line's spectrum, untouched: the
			// reference's ramp loop starts at j = 1.  So no ramp here -- and, in
			// particular, the Nyquist bin is NOT zeroed, though it is for every
			// other j.
			for (int L = 0; L < n; L++) r[L] = 1.0;
		} else {
			// Built by repeated multiplication, exactly as the reference does,
			// rather than from cos/sin of the accumulated angle.  The two drift
			// apart by ~1e-14 over n/2 steps; matching the drift removes a source
			// of divergence from mrdegibbs for nothing.
			const double phi = DG_PI*(double)dg_shift(j)/((double)n*(double)DG_NSH);
			const dcx u = cos(phi) + sin(phi)*I;
			dcx e = 1.0;
			r[0] = 1.0;
			if (!(n & 1)) r[n/2] = 0.0;
			for (int l = 0; l < maxn; l++) {
				e = u*e;
				r[l+1]   = e;
				r[n-1-l] = conj(e);
			}
		}
		// The FFT is unscaled, so the inverse's 1/n rides on the ramp instead.
		for (int L = 0; L < n; L++) r[L] /= (double)n;
	}
	a->ramp = ramp;
	return 0;
}

// One line's worth of the reference's unring_1d inner body: pick the subvoxel
// shift that minimises total variation either side of each sample, and
// interpolate to it.  `sh` is the 41 shifted copies, sh[j*n + l]; the result
// goes back into `eig` at stride `es`.
//
// The interleaving of the argmin scan and the sliding-window update inside the
// j loop is the reference's, and it is load-bearing: at position l, candidate j
// is compared using the TV as it stood BEFORE j's own update at that l.  Hoist
// the scan out of the loop and the answer changes.
static void tv_select(const dcx *sh, int n, dcx *eig, size_t es) {
	double tv1[DG_NSHIFT], tv2[DG_NSHIFT];

	for (int j = 0; j < DG_NSHIFT; j++) {
		const dcx *cj = sh + (size_t)j*n;
		double t1 = 0.0, t2 = 0.0;
		for (int t = DG_MINW; t <= DG_MAXW; t++) {
			t1 += fabs(creal(cj[(n-t)%n]) - creal(cj[(n-t-1)%n]));
			t1 += fabs(cimag(cj[(n-t)%n]) - cimag(cj[(n-t-1)%n]));
			t2 += fabs(creal(cj[(n+t)%n]) - creal(cj[(n+t+1)%n]));
			t2 += fabs(cimag(cj[(n+t)%n]) - cimag(cj[(n+t+1)%n]));
		}
		tv1[j] = t1;
		tv2[j] = t2;
	}

	for (int l = 0; l < n; l++) {
		double mintv = DBL_MAX;
		int minidx = 0;
		for (int j = 0; j < DG_NSHIFT; j++) {
			if (tv1[j] < mintv) { mintv = tv1[j]; minidx = j; }
			if (tv2[j] < mintv) { mintv = tv2[j]; minidx = j; }

			const dcx *cj = sh + (size_t)j*n;
			tv1[j] += fabs(creal(cj[(l-DG_MINW+1+n)%n]) - creal(cj[(l-DG_MINW+n)%n]));
			tv1[j] -= fabs(creal(cj[(l-DG_MAXW+n)%n])   - creal(cj[(l-DG_MAXW-1+n)%n]));
			tv2[j] += fabs(creal(cj[(l+DG_MAXW+1+n)%n]) - creal(cj[(l+DG_MAXW+2+n)%n]));
			tv2[j] -= fabs(creal(cj[(l+DG_MINW+n)%n])   - creal(cj[(l+DG_MINW+1+n)%n]));

			tv1[j] += fabs(cimag(cj[(l-DG_MINW+1+n)%n]) - cimag(cj[(l-DG_MINW+n)%n]));
			tv1[j] -= fabs(cimag(cj[(l-DG_MAXW+n)%n])   - cimag(cj[(l-DG_MAXW-1+n)%n]));
			tv2[j] += fabs(cimag(cj[(l+DG_MAXW+1+n)%n]) - cimag(cj[(l+DG_MAXW+2+n)%n]));
			tv2[j] -= fabs(cimag(cj[(l+DG_MINW+n)%n])   - cimag(cj[(l+DG_MINW+1+n)%n]));
		}

		const dcx *cm = sh + (size_t)minidx*n;
		const double a0r = creal(cm[(l-1+n)%n]), a1r = creal(cm[l]), a2r = creal(cm[(l+1+n)%n]);
		const double a0i = cimag(cm[(l-1+n)%n]), a1i = cimag(cm[l]), a2i = cimag(cm[(l+1+n)%n]);
		const double s = (double)dg_shift(minidx)/(2.0*(double)DG_NSH);

		if (s > 0.0) eig[(size_t)l*es] = (a1r*(1.0-s) + a0r*s) + (a1i*(1.0-s) + a0i*s)*I;
		else         eig[(size_t)l*es] = (a1r*(1.0+s) - a2r*s) + (a1i*(1.0+s) - a2i*s)*I;
	}
}

typedef struct {
	dcx *a, *b;                        // nx*ny each
	dcx *s;                            // the 41 shifted copies: 41 * max(nx,ny)
	dcx *line, *spec, *tmp, *fftscr;   // max(nx,ny) each
} dg_scratch;

// Unring `im` along one axis.  `im` is a column-major nx-by-ny plane with
// leading dimension `ld`; `transposed` says the lines run across it (the y
// pass), which is the reference's unring_1d(im2.transpose()) writing back
// through the view.
static void unring(const dg_axis *ax, const dg_scratch *w, dcx *im, int ld,
                   int numlines, int transposed) {
	const int n = ax->n;
	// One line at a time: apply each ramp to the spectrum and inverse transform
	// it, which is what the reference does.
	for (int k = 0; k < numlines; k++) {
		dcx *eig = transposed ? im + (size_t)k : im + (size_t)k*ld;
		const size_t es = transposed ? (size_t)ld : (size_t)1;
		for (int L = 0; L < n; L++) w->spec[L] = eig[(size_t)L*es];
		for (int j = 0; j < DG_NSHIFT; j++) {
			const dcx *r = ax->ramp + (size_t)j*n;
			for (int L = 0; L < n; L++) w->tmp[L] = w->spec[L]*r[L];
			dg_fft_exec(&ax->ifft, w->s + (size_t)j*n, w->tmp, 1, w->fftscr);
		}
		tv_select(w->s, n, eig, es);
	}
}

// Transform every line of `im` along one axis, in place.  `transposed` says the
// lines run across the plane rather than down it.
static void plane_lines(const dg_axis *ax, const dg_scratch *w, dcx *im, int ld,
                        int nlines, int transposed, int inverse) {
	const int n = ax->n;
	const dg_fft *p = inverse ? &ax->ifft : &ax->fft;
	const double sc = inverse ? 1.0/(double)n : 1.0;
	for (int k = 0; k < nlines; k++) {
		dcx *base = transposed ? im + (size_t)k : im + (size_t)k*ld;
		const size_t st = transposed ? (size_t)ld : (size_t)1;
		dg_fft_exec(p, w->line, base, (int)st, w->fftscr);
		for (int i = 0; i < n; i++) base[(size_t)i*st] = w->line[i]*sc;
	}
}

static void scratch_free(dg_scratch *w) {
	free(w->a); free(w->b); free(w->s);
	free(w->line); free(w->spec); free(w->tmp); free(w->fftscr);
	memset(w, 0, sizeof(*w));
}

static int scratch_init(dg_scratch *w, int nx, int ny) {
	memset(w, 0, sizeof(*w));
	const size_t np = (size_t)nx*(size_t)ny;
	// The x pass runs ny lines of length nx and the y pass nx of length ny, so
	// one buffer sized on the larger of each serves both.
	const int nmax = (nx > ny) ? nx : ny;
	w->a = (dcx *)dn_malloc(np, sizeof(dcx));
	w->b = (dcx *)dn_malloc(np, sizeof(dcx));
	// One line's 41 shifted copies: 92 KiB at n = 140.
	w->s      = (dcx *)dn_malloc((size_t)DG_NSHIFT*(size_t)nmax, sizeof(dcx));
	w->line   = (dcx *)dn_malloc((size_t)nmax, sizeof(dcx));
	w->spec   = (dcx *)dn_malloc((size_t)nmax, sizeof(dcx));
	w->tmp    = (dcx *)dn_malloc((size_t)nmax, sizeof(dcx));
	// maxradix <= n always, so sizing on nmax is a safe bound for either axis.
	w->fftscr = (dcx *)dn_malloc((size_t)nmax, sizeof(dcx));
	if (!w->a || !w->b || !w->s || !w->line || !w->spec || !w->tmp || !w->fftscr) {
		scratch_free(w); return 1;
	}
	return 0;
}

// One (z, volume) plane, in place: `p` is both source and destination, which is
// safe because the plane is copied into scratch before anything is written back.
static void plane_run(const dg_axis *ax, const dg_axis *ay, dg_scratch *w,
                      float *p, int nx, int ny) {
	const size_t np = (size_t)nx*(size_t)ny;

	for (size_t i = 0; i < np; i++) w->a[i] = (double)p[i];

	plane_lines(ay, w, w->a, nx, nx, 1, 0);   // row_FFT: nx lines of length ny
	plane_lines(ax, w, w->a, nx, ny, 0, 0);   // col_FFT: ny lines of length nx

	// Raised-cosine partition of unity, splitting the spectrum into the part
	// unrung along x and the part unrung along y.  im2 must be taken from the
	// ORIGINAL im1 before im1 is scaled.
	for (int k = 0; k < ny; k++) {
		const double ck = ay->cw[k];
		for (int j = 0; j < nx; j++) {
			const double cj = ax->cw[j];
			const size_t q = (size_t)j + (size_t)k*nx;
			if (ck + cj != 0.0) {
				w->b[q] = w->a[q]*cj/(ck+cj);
				w->a[q] = w->a[q]*ck/(ck+cj);
			} else {
				w->a[q] = 0.0;
				w->b[q] = 0.0;
			}
		}
	}

	// Each half is left in the frequency domain along the axis it is about to be
	// unrung, and spatial along the other.
	plane_lines(ay, w, w->a, nx, nx, 1, 1);   // row_iFFT(im1), in place
	plane_lines(ax, w, w->b, nx, ny, 0, 1);   // col_iFFT(im2), in place
	dcx *const im1 = w->a, *const im2 = w->b;

	unring(ax, w, im1, nx, ny, 0);
	unring(ay, w, im2, nx, nx, 1);

	for (size_t i = 0; i < np; i++) p[i] = (float)(creal(im1[i]) + creal(im2[i]));
}

typedef struct {
	const dg_axis *ax, *ay;
	float *data;
	int nx, ny;
	size_t nplane, npix;
	pthread_mutex_t lock;
	size_t next;
	int failed;
} dg_pool;

static void *dg_worker(void *arg) {
	dg_pool *pool = (dg_pool *)arg;
	dg_scratch w;

	if (scratch_init(&w, pool->nx, pool->ny)) {
		pthread_mutex_lock(&pool->lock);
		pool->failed = 1;
		pthread_mutex_unlock(&pool->lock);
		return NULL;
	}

	// One plane per hand-out: a plane is milliseconds of work, so the mutex is
	// noise, and no two planes share an output element.
	for (;;) {
		pthread_mutex_lock(&pool->lock);
		const size_t p = (pool->failed || pool->next >= pool->nplane) ? pool->nplane : pool->next++;
		pthread_mutex_unlock(&pool->lock);
		if (p >= pool->nplane) break;
		plane_run(pool->ax, pool->ay, &w, pool->data + p*pool->npix, pool->nx, pool->ny);
	}

	scratch_free(&w);
	return NULL;
}

int dn_degibbs_check(int nx, int ny) {
	if (nx < DG_MIN_DIM || ny < DG_MIN_DIM) {
		dn_err("-degibbs needs in-plane dimensions of at least %d (got %d x %d)\n",
		       DG_MIN_DIM, nx, ny);
		return 1;
	}
	return 0;
}

// Scratch bytes per worker, to the dominant term.  Two planes, the 41 shifted
// copies of one line, and four line-length temporaries.
static size_t dg_worker_bytes(int nx, int ny) {
	const size_t nmax = (size_t)((nx > ny) ? nx : ny);
	return (2*(size_t)nx*(size_t)ny + (DG_NSHIFT + 4)*nmax) * sizeof(dcx);
}

int dn_degibbs_threads(int nx, int ny, int nz, int nvol, int requested) {
	if (requested < 1) requested = 1;
	const int hw = dn_default_threads();
	if (requested > hw) requested = hw;
	// More workers than planes is pure overhead: the extras allocate a full
	// scratch set, find the queue empty and exit.
	const size_t planes = (size_t)nz*(size_t)nvol;
	if (planes > 0 && (size_t)requested > planes) requested = (int)planes;
	// Same 1 GiB total budget the denoiser uses.  Unbounded, this is 34 MB per
	// worker at 1024 in plane and 134 MB at 2048 -- an audit found it capped by
	// nothing but the core count.
	const size_t per = dg_worker_bytes(nx, ny);
	const size_t budget = (size_t)1 << 30;
	if (per > 0 && (size_t)requested * per > budget) {
		const int cap = (int)(budget / per);
		requested = (cap < 1) ? 1 : cap;
	}
	return requested;
}

int dn_degibbs(float *data, int nx, int ny, int nz, int nvol, int nthreads) {
	if (!data || nz < 1 || nvol < 1) return 1;
	if (dn_degibbs_check(nx, ny)) return 1;
	// Capped here as well as at the CLI, so the budget holds however it is called.
	nthreads = dn_degibbs_threads(nx, ny, nz, nvol, nthreads);

	dg_axis ax, ay;
	if (axis_init(&ax, nx)) { dn_err("cannot allocate the %d-point transform\n", nx); return 1; }
	const int share = (ny == nx);
	if (!share && axis_init(&ay, ny)) {
		dn_err("cannot allocate the %d-point transform\n", ny);
		axis_free(&ax);
		return 1;
	}

	dg_pool pool;
	memset(&pool, 0, sizeof(pool));
	pool.ax     = &ax;
	pool.ay     = share ? &ax : &ay;
	pool.data   = data;
	pool.nx     = nx;
	pool.ny     = ny;
	pool.npix   = (size_t)nx*(size_t)ny;
	pool.nplane = (size_t)nz*(size_t)nvol;

	int rc = 0;
	if (pthread_mutex_init(&pool.lock, NULL) != 0) {
		dn_err("cannot create the degibbs worker mutex\n");
		rc = 1;
	} else {
		dn_thread_run(dg_worker, &pool, nthreads);
		pthread_mutex_destroy(&pool.lock);
	}

	if (pool.failed) {
		dn_err("a degibbs worker could not allocate its scratch space\n");
		rc = 1;
	}
	axis_free(&ax);
	if (!share) axis_free(&ay);
	return rc;
}
