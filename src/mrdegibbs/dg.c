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
 * THE PARTIAL-FOURIER (-pF) CODE BELOW IS NOT MRtrix3's, and the copyright line
 * that follows does not reach it: mrdegibbs has no such option.  It is an
 * independent implementation of the RPG method of
 *   Lee H-H, Novikov DS, Fieremans E.  Removal of partial Fourier-induced Gibbs
 *   (RPG) ringing artifacts in MRI.  Magn Reson Med 2021; 86:2733-2750,
 * written from the paper and a procedural description; no code from the authors'
 * toolbox or from TORTOISE was read.  Its reasoning is inline at the two plane
 * functions rather than here, because it is per-factor: 6/8 reduces to the
 * ordinary method on odd and even columns, 7/8 is a resampling construction.
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
#include <limits.h>
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

// Distinct transform lengths dn_degibbs can ask for: nx, ny, and one for 6/8 or
// two for 7/8.  Sized by the number of axis_get CALLS rather than by how many
// can fire, so the bound is something a reader can count.
#define DG_MAX_AXES 5

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

// The axis for length `n`, built once and then reused, since an axis depends on
// nothing else.  NULL, having reported, if it cannot be allocated -- or if `cap`
// is exhausted, which is a bug rather than a runtime condition but is checked
// here so the array bound does not rest on counting call sites correctly.
static const dg_axis *axis_get(dg_axis *v, int *nv, int cap, int n) {
	for (int i = 0; i < *nv; i++) if (v[i].n == n) return &v[i];
	if (*nv >= cap) {
		dn_err("degibbs asked for more than %d transform lengths\n", cap);
		return NULL;
	}
	if (axis_init(&v[*nv], n)) {
		dn_err("cannot allocate the %d-point transform\n", n);
		return NULL;
	}
	return &v[(*nv)++];
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
	float *up;                         // 7/8: the 3x upsampled plane
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
	free(w->up);
	memset(w, 0, sizeof(*w));
}

static int scratch_init(dg_scratch *w, int nx, int ny, double pf) {
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
	// 7/8 RPG expands y threefold; the interleave it also needs lives in the
	// caller's plane, which is dead for exactly that span.  Tested separately
	// because a NULL here is legitimate off the 7/8 path.
	if (pf == DG_PF_7_8 && !(w->up = (float *)dn_malloc(3*np, sizeof(float)))) {
		scratch_free(w); return 1;
	}
	return 0;
}

// One (z, volume) plane, in place: `p` is both source and destination, which is
// safe because the plane is copied into scratch before anything is written back.
static void plane_run(const dg_axis *ax, const dg_axis *ay, const dg_axis *ay2,
                      dg_scratch *w, float *p, int nx, int ny) {
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

	// PARTIAL FOURIER.  Asymmetric truncation at -v*Kmax and +Kmax, v = 2*pf-1,
	// makes the point spread function a superposition of TWO sincs (Lee et al.
	// 2021, Eq 5), of intervals dy and dy/v.  They superpose, so they come out
	// independently -- the conventional method applied twice.  At pf = 6/8,
	// v = 1/2 and the second interval is 2*dy, which on the odd and even
	// y-columns taken separately IS dy: the ordinary method applies unchanged.
	//
	// This needs no new kernel.  `unring` and `plane_lines` already take the
	// element stride as `ld`, so a sub-image is `im2 + par*nx` with `ld = 2*nx`
	// and a half-length axis.  The forward transform is needed because unring
	// consumes a spectrum and `im2` is spatial by this point.
	//
	// Only the y half is treated.  The x half keeps the wide ringing, as the
	// paper notes: no filter can suppress both y_r and y_v at once.
	if (ay2) {
		for (int par = 0; par < 2; par++) {
			dcx *sub = im2 + (size_t)par*nx;
			plane_lines(ay2, w, sub, 2*nx, nx, 1, 0);
			unring(ay2, w, sub, 2*nx, nx, 1);
		}
	}

	for (size_t i = 0; i < np; i++) p[i] = (float)(creal(im1[i]) + creal(im2[i]));
}

// RPG at 7/8 follows Lee et al.'s resampling construction: expand the
// phase-encode grid by three, SuShi each of four interleaves, contract to the
// original grid, then remove the remaining ordinary y ringing with one 1-D
// pass.  The repeated values are nearest-neighbour samples; positions 3k+1 are
// exactly the original sample centres, so that is the matching downsample.
static void plane_run_7_8(const dg_axis *ax, const dg_axis *ay,
                           const dg_axis *ay78, const dg_axis *ay78short,
                           dg_scratch *w,
                           float *p, int nx, int ny) {
	const size_t np = (size_t)nx*(size_t)ny;
	const int ny3 = 3*ny;
	float *const up = w->up;
	// `p` itself is the interleave buffer: it is dead from the expansion below
	// until the downsample at the end rewrites every row of it, and an interleave
	// is ceil(3ny/4) <= ny rows, so it fits.  A fourth plane of scratch bought
	// nothing but a smaller worker cap.
	float *const sub = p;

	for (int k = 0; k < ny; k++) {
		for (int r = 0; r < 3; r++)
			memcpy(up + (size_t)(3*k + r)*nx, p + (size_t)k*nx,
			       (size_t)nx*sizeof(*p));
	}
	for (int par = 0; par < 4; par++) {
		const int nsub = (ny3 - par + 3)/4;
		const dg_axis *const asub = nsub == ay78->n ? ay78 : ay78short;
		for (int k = 0; k < nsub; k++)
			memcpy(sub + (size_t)k*nx, up + (size_t)(par + 4*k)*nx,
			       (size_t)nx*sizeof(*p));
		plane_run(ax, asub, NULL, w, sub, nx, nsub);
		for (int k = 0; k < nsub; k++)
			memcpy(up + (size_t)(par + 4*k)*nx, sub + (size_t)k*nx,
			       (size_t)nx*sizeof(*p));
	}
	for (int k = 0; k < ny; k++)
		memcpy(p + (size_t)k*nx, up + (size_t)(3*k + 1)*nx,
		       (size_t)nx*sizeof(*p));

	// `unring` takes a y-spectrum and leaves a spatial line.  This is the final
	// ordinary-ringing pass after the four interleaves have removed the wider PF
	// component.
	for (size_t i = 0; i < np; i++) w->a[i] = (double)p[i];
	plane_lines(ay, w, w->a, nx, nx, 1, 0);
	unring(ay, w, w->a, nx, nx, 1);

	// Clamped because THE REFERENCE CLAMPS: dropping this moves the phantom output
	// by relative L2 1.3e-4, where the clamped output agrees with the reference to
	// 1.6e-6.  Note these negatives are ringing UNDERSHOOT created by the
	// transform, not noise about zero -- the input was already truncated in
	// dn_degibbs, so the argument made there does not apply here.  Neither the
	// conventional path nor 6/8 has an equivalent step, and giving them one would
	// cost byte-identity with mrdegibbs, which does not clamp.
	for (size_t i = 0; i < np; i++) {
		const double v = creal(w->a[i]);
		p[i] = (float)(v > 0.0 ? v : 0.0);
	}
}

typedef struct {
	const dg_axis *ax, *ay;
	const dg_axis *ay2;      // half-length y axis for partial Fourier, or NULL
	const dg_axis *ay78;     // 3/4-length y axis for 7/8 partial Fourier, or NULL
	const dg_axis *ay78short;// shorter 7/8 interleave; == ay78 when the two coincide
	float *data;
	int nx, ny;
	double pf;
	size_t nplane, npix;
	pthread_mutex_t lock;
	size_t next;
	int failed;
} dg_pool;

static void *dg_worker(void *arg) {
	dg_pool *pool = (dg_pool *)arg;
	dg_scratch w;

	if (scratch_init(&w, pool->nx, pool->ny, pool->pf)) {
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
		if (pool->pf == DG_PF_7_8)
			plane_run_7_8(pool->ax, pool->ay, pool->ay78, pool->ay78short, &w,
			              pool->data + p*pool->npix, pool->nx, pool->ny);
		else
			plane_run(pool->ax, pool->ay, pool->ay2, &w, pool->data + p*pool->npix,
			          pool->nx, pool->ny);
	}

	scratch_free(&w);
	return NULL;
}

int dn_degibbs_check(int nx, int ny, double pf) {
	if (nx < DG_MIN_DIM || ny < DG_MIN_DIM) {
		dn_err("-degibbs needs in-plane dimensions of at least %d (got %d x %d)\n",
		       DG_MIN_DIM, nx, ny);
		return 1;
	}
	if (pf == DG_PF_FULL) return 0;
	if (pf != DG_PF_6_8 && pf != DG_PF_7_8) {
		// The pipeline differs per factor rather than varying continuously, so
		// the nearest implemented one is NOT a safe substitute.
		// %.9g on the value: %g rounds to 6 significant digits, so -pF 0.750000001
		// was refused with a message naming 0.75 as both the input and a supported
		// factor.  The two constants are exact, so %g is right for them.
		dn_err("-pF %.9g is not implemented; this build supports %g (7/8) and %g (6/8)\n",
		       pf, DG_PF_7_8, DG_PF_6_8);
		return 1;
	}
	if (pf == DG_PF_6_8 && ny % 2) {
		// The 6/8 pipeline splits y into odd and even columns; an odd count
		// makes the two halves different lengths, which nothing here handles.
		dn_err("-pF %g needs an even y dimension (got %d)\n", pf, ny);
		return 1;
	}
	// 7/8 forms 3*ny + 3 as int, and dn_image_read admits a dimension all the way
	// to INT32_MAX.  Reaching this needs a ~16 GB file, so nothing real comes
	// close -- but the overflow would be undefined behaviour rather than a clean
	// allocation failure, so it is refused here beside the other arithmetic bound.
	if (pf == DG_PF_7_8 && ny > (INT_MAX - 3)/3) {
		dn_err("-pF %g cannot handle a y dimension of %d\n", pf, ny);
		return 1;
	}
	// The 6/8 sibling: plane_run strides its sub-image by 2*nx, also in int.
	// Needs a 4 GB plane to reach, but leaving one of the pair guarded and the
	// other not is worse than guarding neither.
	if (pf == DG_PF_6_8 && nx > INT_MAX/2) {
		dn_err("-pF %g cannot handle an x dimension of %d\n", pf, nx);
		return 1;
	}
	const int ysub = pf == DG_PF_6_8 ? ny/2 : 3*ny/4;
	if (ysub < DG_MIN_DIM) {
		dn_err("-pF %g makes its interleave y dimension fall below %d (got %d)\n",
		       pf, DG_MIN_DIM, ysub);
		return 1;
	}
	return 0;
}

// Scratch bytes per worker, to the dominant term.  Two planes, the 41 shifted
// copies of one line, four line-length temporaries, and at 7/8 the threefold
// y expansion.  Must track scratch_init exactly: pricing a worker low lets the
// 1 GiB cap admit more of them than the budget allows.
static size_t dg_worker_bytes(int nx, int ny, double pf) {
	const size_t nmax = (size_t)((nx > ny) ? nx : ny);
	size_t bytes = (2*(size_t)nx*(size_t)ny + (DG_NSHIFT + 4)*nmax) * sizeof(dcx);
	if (pf == DG_PF_7_8) bytes += 3*(size_t)nx*(size_t)ny*sizeof(float);
	return bytes;
}

int dn_degibbs_threads(int nx, int ny, int nz, int nvol, int requested, double pf) {
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
	const size_t per = dg_worker_bytes(nx, ny, pf);
	const size_t budget = (size_t)1 << 30;
	if (per > 0 && (size_t)requested * per > budget) {
		const int cap = (int)(budget / per);
		requested = (cap < 1) ? 1 : cap;
	}
	return requested;
}

int dn_degibbs(float *data, int nx, int ny, int nz, int nvol, int nthreads, double pf) {
	if (!data || nz < 1 || nvol < 1) return 1;
	if (dn_degibbs_check(nx, ny, pf)) return 1;

	// THE METHOD ASSUMES MAGNITUDE DATA, so the input is made to satisfy that
	// rather than the assumption being left implicit.  A negative here is a noise
	// fluctuation about a true zero; it is TRUNCATED, not folded to +|v|, which
	// would turn zero-mean noise into a Rician-like positive floor and hand the
	// TV window spurious steps to fit.  A flat zero baseline is what the subvoxel
	// shift search wants in background regions.
	//
	// No-op on magnitude input, which is why byte-identity with mrdegibbs is
	// unaffected: only a caller that fed signed data -- a phase-rotated series --
	// sees any change, and that caller was violating the assumption already.
	//
	// fmaxf, NOT `if (data[i] < 0.0f) data[i] = 0.0f`, and the difference is
	// worth stating because the second form looks equivalent and is not:
	//   - It DOES NOT VECTORISE: turning a conditional store into an unconditional
	//     one would invent a write C11 forbids, so clang emits a scalar
	//     compare-and-branch per element.  2.1x, and AGENTS.md carries the figures.
	//   - fmaxf maps NaN and -Inf to 0 where the comparison preserves them, which
	//     is wanted: a NaN here would poison every downstream analysis.  +Inf is
	//     NOT touched.  The CLI refuses non-finite input at read time, so this
	//     binds a direct caller -- and `-degibbs y` hands over the denoiser's
	//     output buffer, not the one that check ran on.
	{
		const size_t n = (size_t)nx*(size_t)ny*(size_t)nz*(size_t)nvol;
		for (size_t i = 0; i < n; i++) data[i] = fmaxf(data[i], 0.0f);
	}
	// Capped here as well as at the CLI, so the budget holds however it is called.
	nthreads = dn_degibbs_threads(nx, ny, nz, nvol, nthreads, pf);

	// Every transform length needed, deduplicated by length.  A square plane
	// wants one axis rather than two, and so do the 7/8 interleaves when their
	// two lengths coincide -- one rule covers both, and also the case ny/2 == nx,
	// which used to build a second identical axis.  At most four distinct lengths
	// are ever asked for: nx, ny, plus one for 6/8 or two for 7/8.
	// At most four of these five fire -- 6/8 and 7/8 are mutually exclusive -- and
	// axis_get refuses past DG_MAX_AXES either way.  nav counts successes, so it
	// bounds both the NULL check and the teardown.
	dg_axis av[DG_MAX_AXES];
	int nav = 0;
	const dg_axis *const ax  = axis_get(av, &nav, DG_MAX_AXES, nx);
	const dg_axis *const ay  = axis_get(av, &nav, DG_MAX_AXES, ny);
	const dg_axis *const ay2 = pf == DG_PF_6_8 ? axis_get(av, &nav, DG_MAX_AXES, ny/2) : NULL;
	const dg_axis *const ay78  = pf == DG_PF_7_8 ? axis_get(av, &nav, DG_MAX_AXES, (3*ny + 3)/4) : NULL;
	const dg_axis *const ay78s = pf == DG_PF_7_8 ? axis_get(av, &nav, DG_MAX_AXES, 3*ny/4) : NULL;
	if (!ax || !ay || (pf == DG_PF_6_8 && !ay2) ||
	    (pf == DG_PF_7_8 && (!ay78 || !ay78s))) {
		for (int i = 0; i < nav; i++) axis_free(&av[i]);
		return 1;
	}

	dg_pool pool;
	memset(&pool, 0, sizeof(pool));
	pool.ax     = ax;
	pool.ay     = ay;
	pool.ay2    = ay2;
	pool.ay78   = ay78;
	pool.ay78short = ay78s;
	pool.data   = data;
	pool.nx     = nx;
	pool.ny     = ny;
	pool.pf     = pf;
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
	for (int i = 0; i < nav; i++) axis_free(&av[i]);
	return rc;
}
