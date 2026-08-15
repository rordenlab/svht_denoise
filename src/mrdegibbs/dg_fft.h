// Mixed-radix complex FFT: the degibbs transform, on every platform.
// See dg_fft.c for provenance.

#ifndef DG_FFT_H
#define DG_FFT_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double _Complex dg_cx;

// A plan for one length and one direction.  Read-only once built, so one plan
// is shared by every worker; the per-call scratch is passed in instead.
typedef struct {
	int n;
	int inverse;
	int nstage;
	int maxradix;      // largest stage radix, i.e. the scratch dg_fft_exec needs
	dg_cx *tw;         // n twiddles
	int *radix;        // per stage
	int *remainder;    // per stage
} dg_fft;

// `inverse` selects the twiddle sign only: the result is UNSCALED, so an
// inverse transform still has to be divided by n by the caller.
int  dg_fft_init(dg_fft *p, int n, int inverse);
void dg_fft_free(dg_fft *p);

// out[0..n-1] = transform of in[0], in[in_stride], ...  `out` must not alias
// `in`.  `scratch` must hold at least p->maxradix elements.
void dg_fft_exec(const dg_fft *p, dg_cx *out, const dg_cx *in, int in_stride,
                 dg_cx *scratch);

#ifdef __cplusplus
}
#endif

#endif // DG_FFT_H
