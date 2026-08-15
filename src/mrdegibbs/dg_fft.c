/* Mixed-radix complex FFT: radix 2, 3, 4 and 5 with a generic butterfly for any
 * other prime factor, so it transforms ANY length.  That is what the degibbs
 * kernel needs -- 104 = 2^3 * 13 and 140 = 2^2 * 5 * 7 both have a prime factor
 * no power-of-two FFT will take.
 *
 * Adapted to C from Eigen's unsupported/Eigen/src/FFT/ei_kissfft_impl.h,
 * Copyright (c) 2009 Mark Borgerding, itself derived from kissfft
 * (http://sourceforge.net/projects/kissfft), Copyright 2003-2009 Mark
 * Borgerding.  That file is MPL-2.0, as is this project.
 *
 * This is the same implementation MRtrix3's mrdegibbs reaches through
 * Eigen::FFT -- it links no FFTW and no BLAS -- which is deliberate: matching
 * the reference's algorithm keeps the portable path's numbers close to it as
 * well as its speed.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "mrdegibbs/dg_fft.h"

#define DG_FFT_PI 3.14159265358979323846
#define DG_FFT_MAXSTAGE 32   // 2^32 exceeds any image dimension by a wide margin

// Four quadrants built from small angles rather than one running phase, which
// is what keeps the twiddles accurate at large n.  Transcribed from the
// reference; the symmetry tw[n-i] = conj(tw[i]) is exact by construction.
static void make_twiddles(dg_cx *tw, int n, int inverse) {
	const double phinc = 0.25*DG_FFT_PI/(double)n;
	const double flip = inverse ? 1.0 : -1.0;
	tw[0] = 1.0;
	if ((n & 1) == 0) tw[n/2] = -1.0;
	int i = 1;
	for (; i*8 < n; ++i) {
		const double c = cos(i*8*phinc), s = sin(i*8*phinc);
		tw[i]   = c + (s*flip)*I;
		tw[n-i] = c - (s*flip)*I;
	}
	for (; i*4 < n; ++i) {
		const double c = cos((2*n - 8*i)*phinc), s = sin((2*n - 8*i)*phinc);
		tw[i]   = s + (c*flip)*I;
		tw[n-i] = s - (c*flip)*I;
	}
	for (; i*8 < 3*n; ++i) {
		const double c = cos((8*i - 2*n)*phinc), s = sin((8*i - 2*n)*phinc);
		tw[i]   = -s + (c*flip)*I;
		tw[n-i] = -s - (c*flip)*I;
	}
	for (; i*2 < n; ++i) {
		const double c = cos((4*n - 8*i)*phinc), s = sin((4*n - 8*i)*phinc);
		tw[i]   = -c + (s*flip)*I;
		tw[n-i] = -c - (s*flip)*I;
	}
}

// Factor out 4s first, then 2s, then odd numbers.  Radix 4 ahead of 2 is worth
// it: it halves the stage count and does the same work with fewer twiddles.
static int factorize(dg_fft *p, int n) {
	int rem = n, radix = 4, ns = 0;
	do {
		while (rem % radix) {
			switch (radix) {
				case 4:  radix = 2; break;
				case 2:  radix = 3; break;
				default: radix += 2; break;
			}
			if (radix*radix > rem) radix = rem;   // no factor can exceed sqrt(rem)
		}
		if (ns >= DG_FFT_MAXSTAGE) return 1;
		rem /= radix;
		p->radix[ns] = radix;
		p->remainder[ns] = rem;
		if (radix > p->maxradix) p->maxradix = radix;
		ns++;
	} while (rem > 1);
	p->nstage = ns;
	return 0;
}

static void bfly2(dg_cx *F, size_t fstride, const dg_cx *tw, int m) {
	for (int k = 0; k < m; ++k) {
		const dg_cx t = F[m+k]*tw[(size_t)k*fstride];
		F[m+k] = F[k] - t;
		F[k]  += t;
	}
}

static void bfly4(dg_cx *F, size_t fstride, const dg_cx *tw, int m, int inverse) {
	const double neg = inverse ? -1.0 : 1.0;
	for (int k = 0; k < m; ++k) {
		dg_cx s0 = F[k+m]*tw[(size_t)k*fstride];
		dg_cx s1 = F[k+2*m]*tw[(size_t)k*fstride*2];
		dg_cx s2 = F[k+3*m]*tw[(size_t)k*fstride*3];
		const dg_cx s5 = F[k] - s1;
		F[k] += s1;
		const dg_cx s3 = s0 + s2;
		dg_cx s4 = s0 - s2;
		s4 = (cimag(s4)*neg) + (-creal(s4)*neg)*I;
		F[k+2*m] = F[k] - s3;
		F[k]    += s3;
		F[k+m]   = s5 + s4;
		F[k+3*m] = s5 - s4;
	}
}

static void bfly3(dg_cx *F, size_t fstride, const dg_cx *tw, int m) {
	const int m2 = 2*m;
	const dg_cx epi3 = tw[fstride*(size_t)m];
	const dg_cx *tw1 = tw, *tw2 = tw;
	int k = m;
	do {
		const dg_cx s1 = F[m]*(*tw1);
		const dg_cx s2 = F[m2]*(*tw2);
		const dg_cx s3 = s1 + s2;
		dg_cx s0 = s1 - s2;
		tw1 += fstride;
		tw2 += fstride*2;
		F[m] = (creal(F[0]) - 0.5*creal(s3)) + (cimag(F[0]) - 0.5*cimag(s3))*I;
		s0 *= cimag(epi3);
		F[0] += s3;
		F[m2] = (creal(F[m]) + cimag(s0)) + (cimag(F[m]) - creal(s0))*I;
		F[m] += (-cimag(s0)) + creal(s0)*I;
		++F;
	} while (--k);
}

static void bfly5(dg_cx *F, size_t fstride, const dg_cx *tw, int m) {
	const dg_cx ya = tw[fstride*(size_t)m];
	const dg_cx yb = tw[fstride*2*(size_t)m];
	dg_cx *F0 = F, *F1 = F+m, *F2 = F+2*m, *F3 = F+3*m, *F4 = F+4*m;

	for (int u = 0; u < m; ++u) {
		const dg_cx s0 = *F0;
		const dg_cx s1 = *F1*tw[(size_t)u*fstride];
		const dg_cx s2 = *F2*tw[2*(size_t)u*fstride];
		const dg_cx s3 = *F3*tw[3*(size_t)u*fstride];
		const dg_cx s4 = *F4*tw[4*(size_t)u*fstride];

		const dg_cx s7  = s1 + s4, s10 = s1 - s4;
		const dg_cx s8  = s2 + s3, s9  = s2 - s3;

		*F0 += s7;
		*F0 += s8;

		const dg_cx s5 = s0 + ((creal(s7)*creal(ya) + creal(s8)*creal(yb))
		                    + (cimag(s7)*creal(ya) + cimag(s8)*creal(yb))*I);
		const dg_cx s6 = (cimag(s10)*cimag(ya) + cimag(s9)*cimag(yb))
		               + (-creal(s10)*cimag(ya) - creal(s9)*cimag(yb))*I;
		*F1 = s5 - s6;
		*F4 = s5 + s6;

		const dg_cx s11 = s0 + ((creal(s7)*creal(yb) + creal(s8)*creal(ya))
		                     + (cimag(s7)*creal(yb) + cimag(s8)*creal(ya))*I);
		const dg_cx s12 = (-cimag(s10)*cimag(yb) + cimag(s9)*cimag(ya))
		                + (creal(s10)*cimag(yb) - creal(s9)*cimag(ya))*I;
		*F2 = s11 + s12;
		*F3 = s11 - s12;

		++F0; ++F1; ++F2; ++F3; ++F4;
	}
}

// Any prime radix above 5.  O(p^2) per group, which is why the small radices
// are written out: 13 here costs 169 multiplies where a radix-13 kernel would
// cost far less, but a factor of 13 is rare and a kernel for every prime is not
// a trade worth making.
static void bfly_generic(dg_cx *F, size_t fstride, const dg_cx *tw, int norig,
                         int m, int p, dg_cx *scr) {
	for (int u = 0; u < m; ++u) {
		int k = u;
		for (int q1 = 0; q1 < p; ++q1) { scr[q1] = F[k]; k += m; }
		k = u;
		for (int q1 = 0; q1 < p; ++q1) {
			int twidx = 0;
			F[k] = scr[0];
			for (int q = 1; q < p; ++q) {
				twidx += (int)fstride*k;
				if (twidx >= norig) twidx -= norig;
				F[k] += scr[q]*tw[twidx];
			}
			k += m;
		}
	}
}

static void kf_work(const dg_fft *p, int stage, dg_cx *xout, const dg_cx *xin,
                    size_t fstride, int in_stride, dg_cx *scr) {
	const int radix = p->radix[stage];
	const int m = p->remainder[stage];
	dg_cx *const beg = xout;
	dg_cx *const end = xout + (size_t)radix*m;

	if (m > 1) {
		do {
			kf_work(p, stage+1, xout, xin, fstride*(size_t)radix, in_stride, scr);
			xin += fstride*(size_t)in_stride;
		} while ((xout += m) != end);
	} else {
		do {
			*xout = *xin;
			xin += fstride*(size_t)in_stride;
		} while (++xout != end);
	}
	xout = beg;

	switch (radix) {
		case 2:  bfly2(xout, fstride, p->tw, m); break;
		case 3:  bfly3(xout, fstride, p->tw, m); break;
		// The twiddle sign already encodes the direction everywhere else; radix 4
		// is the one butterfly that also needs it explicitly, for its i multiply.
		case 4:  bfly4(xout, fstride, p->tw, m, p->inverse); break;
		case 5:  bfly5(xout, fstride, p->tw, m); break;
		default: bfly_generic(xout, fstride, p->tw, p->n, m, radix, scr); break;
	}
}

void dg_fft_exec(const dg_fft *p, dg_cx *out, const dg_cx *in, int in_stride,
                 dg_cx *scratch) {
	kf_work(p, 0, out, in, 1, in_stride, scratch);
}

int dg_fft_init(dg_fft *p, int n, int inverse) {
	memset(p, 0, sizeof(*p));
	// dn_calloc rounds a zero-element request up to one byte, so without this
	// make_twiddles would write tw[0] past the end.  Unreachable through
	// axis_init, which is guarded by DG_MIN_DIM; here so it stays unreachable.
	if (n < 1) return 1;
	p->n = n;
	p->inverse = inverse ? 1 : 0;
	p->tw        = (dg_cx *)dn_calloc((size_t)n, sizeof(dg_cx));
	p->radix     = (int *)dn_calloc(DG_FFT_MAXSTAGE, sizeof(int));
	p->remainder = (int *)dn_calloc(DG_FFT_MAXSTAGE, sizeof(int));
	if (!p->tw || !p->radix || !p->remainder) { dg_fft_free(p); return 1; }
	if (factorize(p, n)) { dg_fft_free(p); return 1; }
	make_twiddles(p->tw, n, inverse);
	return 0;
}

void dg_fft_free(dg_fft *p) {
	if (!p) return;
	free(p->tw); free(p->radix); free(p->remainder);
	memset(p, 0, sizeof(*p));
}
