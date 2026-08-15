// Householder reduction to tridiagonal form (tred2) followed by
// the implicit-shift symmetric QL algorithm with eigenvector accumulation
// (tql2).  Both are derived from the Algol procedures of Bowdler, Martin,
// Reinsch and Wilkinson (Handbook for Automatic Computation, Vol. II - Linear
// Algebra) and the corresponding EISPACK Fortran, via the public-domain JAMA
// translation.  The immediate source is niimath's src/tensor.c (BSD-2-Clause,
// rordenlab/niimath), where the same routines are specialised to n == 3; here
// they are generalised to a runtime n on flat row-major storage, and an
// iteration cap plus a status return are added -- the original loops until
// convergence with no bound, which is not acceptable in a kernel that runs
// millions of times unattended.
//
// WHY NOT CYCLIC JACOBI: niimath's md_jacobi_eigh (medic.c) is the obvious
// candidate and is a good deal simpler, but its cost is ~6*n^3 per sweep and it
// needs several sweeps.  At n = 102 that is roughly an order of magnitude more
// arithmetic than tred2/tql2 for the same answer, in a loop that runs millions
// of times per image.

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_eig.h"

#ifdef DN_USE_ACCELERATE
#include <Accelerate/Accelerate.h>
// Evaluation build (`make ACCELERATE=1`): LAPACK's blocked dsytrd replaces the
// reduction, and dormtr the back-transform, on the SPLIT path only.  The full
// path (n < DN_SPLIT_MIN) is untouched, because the reduction is not where its
// time goes.  `hh` carries LAPACK's tau instead of our h values -- same length,
// same lifetime, so it is reused rather than duplicated.
#endif

struct dn_eig {
	int n;
	double *d;   // n : diagonal / eigenvalues
	double *e;   // n : off-diagonal workspace
	double *w;   // n*n : eigenvector scratch for the split path and its fallback

	// State carried from dn_eig_values() to dn_eig_vectors().
	double *refl;    // n*n : Householder reflectors from the reduction
	double *hh;      // n   : their h values
	double *td, *te; // n   : the tridiagonal, preserved (tql2 destroys its copy)
	double *ev;      // n   : the eigenvalues, ascending
	double *z;       // n   : one eigenvector under construction
	double *gt;      // 6*n : tridiagonal LU workspace (4n doubles + an int block)
	double tnorm;    // ||T||, for the inverse-iteration tolerances
#ifdef DN_USE_ACCELERATE
	double *lwk;     // LAPACK workspace, sized once at create
	int lwork;
#endif
	unsigned long fallbacks;  // times inverse iteration lost and the full solve ran
	int staged;      // 1 if the split path has valid state
	int full_valid;  // 1 if `refl` already holds a complete eigenvector set
};

// ---------------------------------------------------------------------------
// tred2 + tql2
// ---------------------------------------------------------------------------

// The Householder reduction, shared by both callers below.  They differ only in
// the TAIL, and only because d is left holding the per-step h values here:
// tred_reduce puts the diagonal back and keeps the reflectors, while tred2
// consumes those same h values to expand the transform.
static void tred_reduce_raw(double *V, double *d, double *e, double *hh, int n) {
	for (int j = 0; j < n; j++) d[j] = V[(size_t)(n - 1) * n + j];
	for (int i = 0; i < n; i++) hh[i] = 0.0;

	for (int i = n - 1; i > 0; i--) {
		double scale = 0.0, h = 0.0;
		for (int k = 0; k < i; k++) scale += fabs(d[k]);
		if (scale == 0.0) {
			e[i] = d[i - 1];
			for (int j = 0; j < i; j++) {
				d[j] = V[(size_t)(i - 1) * n + j];
				V[(size_t)i * n + j] = 0.0;
				V[(size_t)j * n + i] = 0.0;
			}
		} else {
			for (int k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
			double f = d[i - 1];
			double g = sqrt(h);
			if (f > 0) g = -g;
			e[i] = scale * g;
			h -= f * g;
			d[i - 1] = f - g;
			for (int j = 0; j < i; j++) e[j] = 0.0;

			for (int j = 0; j < i; j++) {
				f = d[j];
				V[(size_t)j * n + i] = f;          // the reflector, kept
				g = e[j] + V[(size_t)j * n + j] * f;
				for (int k = j + 1; k <= i - 1; k++) {
					g += V[(size_t)k * n + j] * d[k];
					e[k] += V[(size_t)k * n + j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
			double hhh = f / (h + h);
			for (int j = 0; j < i; j++) e[j] -= hhh * d[j];
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i - 1; k++)
					V[(size_t)k * n + j] -= (f * e[k] + g * d[k]);
				d[j] = V[(size_t)(i - 1) * n + j];
				V[(size_t)i * n + j] = 0.0;
			}
		}
		hh[i] = h;
		d[i] = h;
	}
}

// Householder reduction of the symmetric matrix V (n*n row-major, destroyed) to
// tridiagonal form.  On exit V holds the accumulated orthogonal transform, d the
// diagonal and e the sub-diagonal.  hh is scratch: this caller wants the
// transform expanded rather than the reflectors kept, so it never reads it back.
static void tred2(double *V, double *d, double *e, double *hh, int n) {
	tred_reduce_raw(V, d, e, hh, n);   // leaves d holding the per-step h values

	for (int i = 0; i < n - 1; i++) {
		V[(size_t)(n - 1) * n + i] = V[(size_t)i * n + i];
		V[(size_t)i * n + i] = 1.0;
		double h = d[i + 1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) d[k] = V[(size_t)k * n + i + 1] / h;
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) g += V[(size_t)k * n + i + 1] * V[(size_t)k * n + j];
				for (int k = 0; k <= i; k++) V[(size_t)k * n + j] -= g * d[k];
			}
		}
		for (int k = 0; k <= i; k++) V[(size_t)k * n + i + 1] = 0.0;
	}
	for (int j = 0; j < n; j++) {
		d[j] = V[(size_t)(n - 1) * n + j];
		V[(size_t)(n - 1) * n + j] = 0.0;
	}
	V[(size_t)(n - 1) * n + n - 1] = 1.0;
	e[0] = 0.0;
}


// sqrt(a*a + b*b), replacing libm's hypot in tql2's rotation loop -- the two
// call sites there are 7% of the whole run, and appeared in the profile as
// DYLD-STUB$$hypot, so libm's was not even being inlined.  Worth 8.7% of CPU
// time on 100x100x58x138.
//
// THE FAST PATH MUST NOT DIVIDE, and that is the whole point rather than a
// micro-optimisation.  The obvious scaling form -- hi * sqrt(1 + (lo/hi)^2) for
// every input -- was written first and measured: it is worth 1.3%, because the
// division costs about what the call it replaced did.  Only the version below,
// which reaches sqrt() with no division at all in the ordinary case, converts
// the 7% the profile promised.
//
// The guarded range is nowhere near the exponent limits (a double reaches
// ~1.8e308), so squaring inside it cannot overflow or flush to zero, and the
// scaled form still catches anything outside.
//
// THIS IS NOT A DROP-IN FOR libm's hypot.  `a*a + b*b` rounds before the sqrt
// does, where hypot is correctly rounded once; over 4 million pairs spanning
// 2^[-300, 300], 0.56% disagree, always by exactly 1 ulp.  Swapping back to
// hypot therefore does NOT reproduce a previous output bit for bit -- end to end
// it moves the denoised series by ~2e-10 relative L2, well inside float32, with
// the rank map unchanged.  What the tool promises is determinism, which holds
// because dn_hypot is used consistently, not because it matches libm.
//
// Measure that with FULL-ENTROPY inputs: 24-bit mantissas make a*a and b*b
// exactly representable, and duly report agreement on every one of 4 M pairs.
//
// Non-finite behaviour is NOT simply "NaN in, NaN out".  (+-Inf, finite) gives
// Inf, as hypot does.  Two cases deviate, both because a comparison against a
// NaN is false and so `hi` takes the other operand:
//
//   (Inf, NaN) and (NaN, Inf) give NaN where hypot must give Inf.
//   (NaN, +-0.0) gives 0.0 -- it SWALLOWS the NaN.  (0.0, NaN) gives NaN, so
//   the asymmetry is real rather than a typo.
//
// Only the second could hide anything, and only from `dn_hypot(p, e[i])` with a
// NaN p and an exactly-zero e[i]; the other call site passes a literal 1.0.
// Reaching it needs a non-finite Gram, and every other rotation in the sweep
// still carries the NaN into the spectrum that dn_eig_values checks.  Left
// alone deliberately: guarding it means a NaN test on the fast path, which is
// the one thing this function exists to keep clear.
static inline double dn_hypot(double a, double b) {
	const double s = a * a + b * b;
	if (s > 1e-250 && s < 1e250) return sqrt(s);   // no division in the common case
	a = fabs(a);
	b = fabs(b);
	const double hi = (a > b) ? a : b;
	const double lo = (a > b) ? b : a;
	if (hi == 0.0) return 0.0;   // the only case where lo/hi would be 0/0
	const double r = lo / hi;
	return hi * sqrt(1.0 + r * r);
}

#define DN_TQL2_MAXITER 50

// Symmetric tridiagonal QL with implicit shifts.  Returns 0 on success, or the
// 1-based index of the eigenvalue that failed to converge.
//
// V == NULL requests EIGENVALUES ONLY, which is what the split path needs: the
// median wants the whole spectrum but the projection wants only a few vectors.
// The guard is hoisted out of the length-n inner loop, so the values-only path
// costs one predictable branch per rotation rather than a duplicated algorithm.
static int tql2(double *V, double *d, double *e, int n) {
	for (int i = 1; i < n; i++) e[i - 1] = e[i];
	e[n - 1] = 0.0;

	double f = 0.0, tst1 = 0.0;
	const double eps = DBL_EPSILON;   // 2^-52

	for (int l = 0; l < n; l++) {
		double t = fabs(d[l]) + fabs(e[l]);
		if (t > tst1) tst1 = t;

		// Bound the deflation search at n-1, as EISPACK's own tql2 does.  The
		// JAMA translation this came from searches to n and relies on e[n-1]
		// being exactly zero to terminate.  That invariant holds for finite
		// input but not for a NaN: ANY comparison against a NaN is false, so
		// `fabs(e[m]) <= eps*tst1` never breaks and the loop leaves m == n.  The
		// QL step that follows then indexes d[l+1] and e[i+1] at l == n-1,
		// writing one double PAST THE END of both arrays -- a heap overflow into
		// whatever dn_malloc handed out next.  Bounding the search removes the
		// dependence on that invariant entirely.
		int m = l;
		while (m < n - 1) {
			if (fabs(e[m]) <= eps * tst1) break;
			m++;
		}
		if (m > l) {
			int iter = 0;
			do {
				if (++iter > DN_TQL2_MAXITER) return l + 1;

				double g = d[l];
				double p = (d[l + 1] - g) / (2.0 * e[l]);
				double r = dn_hypot(p, 1.0);
				if (p < 0) r = -r;
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				double dl1 = d[l + 1];
				double h = g - d[l];
				for (int i = l + 2; i < n; i++) d[i] -= h;
				f += h;

				p = d[m];
				double c = 1.0, c2 = c, c3 = c;
				double el1 = e[l + 1];
				double s = 0.0, s2 = 0.0;
				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = dn_hypot(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);
					if (V) {
						for (int k = 0; k < n; k++) {
							h = V[(size_t)k * n + i + 1];
							V[(size_t)k * n + i + 1] = s * V[(size_t)k * n + i] + c * h;
							V[(size_t)k * n + i] = c * V[(size_t)k * n + i] - s * h;
						}
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;
			} while (fabs(e[l]) > eps * tst1);
		}
		d[l] += f;
		e[l] = 0.0;
	}

	// Sort ascending, permuting eigenvector columns with the eigenvalues.
	for (int i = 0; i < n - 1; i++) {
		int k = i;
		double p = d[i];
		for (int j = i + 1; j < n; j++) {
			if (d[j] < p) { k = j; p = d[j]; }
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			if (V) {
				for (int j = 0; j < n; j++) {
					p = V[(size_t)j * n + i];
					V[(size_t)j * n + i] = V[(size_t)j * n + k];
					V[(size_t)j * n + k] = p;
				}
			}
		}
	}
	return 0;
}

// ---------------------------------------------------------------------------
// split path: eigenvalues first, then only the eigenvectors that are wanted
// ---------------------------------------------------------------------------
//
// tred2 above both reduces to tridiagonal form AND accumulates the orthogonal
// transform, and tql2 then accumulates every eigenvector through its rotations.
// Together those two accumulations are roughly 7/8 of the arithmetic, and this
// tool throws almost all of it away: the median needs the eigenvalues, the
// projection needs about 8 eigenvectors out of 102.
//
// So: reduce while KEEPING the Householder reflectors instead of expanding them
// (tred_reduce), take the eigenvalues off the tridiagonal with no vectors at all
// (tql2 with V == NULL), then build only the wanted eigenvectors by inverse
// iteration on
// the tridiagonal -- O(n) per solve -- back-transformed through the reflectors.
// This is the shape of LAPACK's dsytrd / dsterf / dstein / dormtr, written small.

#ifndef DN_USE_ACCELERATE
// Reduction half of tred2.  On exit V holds the Householder vector for step p in
// V[0..p-1][p], hh[p] its h, d the diagonal and e the sub-diagonal.
static void tred_reduce(double *V, double *d, double *e, double *hh, int n) {
	tred_reduce_raw(V, d, e, hh, n);
	// Recover the diagonal, which the shared loop overwrote with the h values.
	for (int j = 0; j < n; j++) d[j] = V[(size_t)j * n + j];
	e[0] = 0.0;
}
#endif

#ifndef DN_USE_ACCELERATE
// z <- Q z, where Q is the product of the stored reflectors.  The accumulation
// loop in tred2 applies reflector p to column j for every p > j in increasing p,
// so Q = H_{n-1} ... H_1 and the same order applies here.
static void tred_apply_q(const double *V, const double *hh, double *z, int n) {
	for (int p = 1; p <= n - 1; p++) {
		const double h = hh[p];
		if (h == 0.0) continue;
		double g = 0.0;
		for (int k = 0; k < p; k++) g += V[(size_t)k * n + p] * z[k];
		g /= h;
		for (int k = 0; k < p; k++) z[k] -= g * V[(size_t)k * n + p];
	}
}
#else
// dsytrd wants a symmetric matrix in column-major order; ours is row-major with
// BOTH triangles filled, and a full symmetric matrix is its own transpose, so it
// is passed straight through with no copy and no repacking.
//
// The sub-diagonal goes into e + 1, not e.  LAPACK indexes it e[i] = T(i, i+1),
// while everything here indexes it e[i] = T(i-1, i) with e[0] unused; offsetting
// the pointer by one converts between them for free.
static void tred_reduce_accel(dn_eig *eg, double *V, double *d, double *e,
                              double *tau, int n) {
	__LAPACK_int N = n, LDA = n, LW = eg->lwork, INFO = 0;
	dsytrd_("U", &N, V, &LDA, d, e + 1, tau, eg->lwk, &LW, &INFO);
	e[0] = 0.0;
}

// Q applied to `nc` vectors at once.  They are stored one per contiguous run of
// n doubles, which is exactly a column-major n x nc block with ldc = n, so the
// per-vector loop this replaces collapses into one call.
static void tred_apply_q_accel(dn_eig *eg, const double *V, const double *tau,
                               double *C, int nc, int n) {
	__LAPACK_int N = n, NC = nc, LDA = n, LDC = n, LW = eg->lwork, INFO = 0;
	dormtr_("L", "U", "N", &N, &NC, (double *)V, &LDA, (double *)tau,
	        C, &LDC, eg->lwk, &LW, &INFO);
}
#endif

// Solve (T - lam*I) x = b for a symmetric tridiagonal T, in place on x.
// Gaussian elimination with partial pivoting, LAPACK-style: a pivot that
// underflows is replaced by eps*||T|| so an exact eigenvalue shift, which makes
// the matrix singular by construction, still yields a usable direction.
static void tridiag_solve(const double *d, const double *e, int n, double lam,
                          double *x, double *work, double tnorm) {
	double *dl = work;             // n : sub-diagonal
	double *dd = work + n;         // n : diagonal
	double *du = work + 2 * n;     // n : super-diagonal
	double *du2 = work + 3 * n;    // n : second super-diagonal (fill-in)
	int *piv = (int *)(work + 4 * n);

	const double tiny = DBL_EPSILON * (tnorm > 0.0 ? tnorm : 1.0);

	for (int i = 0; i < n; i++) {
		dd[i] = d[i] - lam;
		dl[i] = (i + 1 < n) ? e[i + 1] : 0.0;   // e[i+1] couples i and i+1
		du[i] = (i + 1 < n) ? e[i + 1] : 0.0;
		du2[i] = 0.0;
		piv[i] = 0;
	}

	for (int i = 0; i < n - 1; i++) {
		if (fabs(dd[i]) >= fabs(dl[i])) {
			if (fabs(dd[i]) < tiny) dd[i] = (dd[i] < 0.0) ? -tiny : tiny;
			const double m = dl[i] / dd[i];
			dl[i] = m;
			dd[i + 1] -= m * du[i];
			if (i + 2 < n) du2[i] = 0.0;
		} else {
			// Pivot: swap rows i and i+1.
			//
			// Deliberately NO tiny-pivot clamp here, unlike the two branches that
			// bound dd[].  A review suggested adding one for symmetry; measured,
			// it changed 34% of voxels on a test volume and moved 85794
			// ranks.  The reason is that this branch is only taken when
			// |dd[i]| < |dl[i]|, so a "tiny" dl[i] means the whole row is
			// negligible -- and in inverse iteration (T - lam*I) is deliberately
			// near-singular, so those rows are exactly where the eigenvector's
			// information lives.  Clamping them destroys it.  A genuinely
			// unusable factorisation still cannot escape: tridiag_eigvec checks
			// isfinite() on the result and falls back to the full solve.
			const double m = dd[i] / dl[i];
			dd[i] = dl[i];
			const double t = du[i];
			du[i] = dd[i + 1];
			dd[i + 1] = t - m * dd[i + 1];
			if (i + 2 < n) {
				du2[i] = du[i + 1];
				du[i + 1] = -m * du[i + 1];
			}
			dl[i] = m;
			piv[i] = 1;
		}
	}
	if (fabs(dd[n - 1]) < tiny) dd[n - 1] = (dd[n - 1] < 0.0) ? -tiny : tiny;

	// Forward substitution with the recorded row swaps.
	for (int i = 0; i < n - 1; i++) {
		if (piv[i]) {
			const double t = x[i];
			x[i] = x[i + 1];
			x[i + 1] = t - dl[i] * x[i];
		} else {
			x[i + 1] -= dl[i] * x[i];
		}
	}
	// Back substitution.
	x[n - 1] /= dd[n - 1];
	if (n > 1) x[n - 2] = (x[n - 2] - du[n - 2] * x[n - 1]) / dd[n - 2];
	for (int i = n - 3; i >= 0; i--)
		x[i] = (x[i] - du[i] * x[i + 1] - du2[i] * x[i + 2]) / dd[i];
}

static double vec_norm(const double *v, int n) {
	double s = 0.0;
	for (int i = 0; i < n; i++) s += v[i] * v[i];
	return sqrt(s);
}

// One eigenvector of the tridiagonal for eigenvalue `lam`, by inverse iteration,
// fully re-orthogonalised against the vectors already produced in this batch
// (`done` columns of `out`, stride n).  Deterministic: the start vector comes
// from a fixed sequence seeded by the eigenvalue's index, never from the data or
// from rand(), so the result cannot vary between runs or thread counts.
static int tridiag_eigvec(const double *d, const double *e, int n, double lam,
                          double tnorm, int idx, double *z, double *work,
                          double *out, int done, int outstride) {
	unsigned seed = 2463534242u + (unsigned)idx * 2654435761u;
	for (int i = 0; i < n; i++) {
		seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
		z[i] = ((double)(seed >> 8) / 8388608.0) - 1.0;   // in [-1, 1)
	}
	double nz = vec_norm(z, n);
	if (!(nz > 0.0)) { for (int i = 0; i < n; i++) z[i] = 1.0; nz = sqrt((double)n); }
	for (int i = 0; i < n; i++) z[i] /= nz;

	for (int it = 0; it < 5; it++) {
		tridiag_solve(d, e, n, lam, z, work, tnorm);
		nz = vec_norm(z, n);
		if (!(nz > 0.0) || !isfinite(nz)) return 1;
		for (int i = 0; i < n; i++) z[i] /= nz;

		for (int p = 0; p < done; p++) {
			const double *v = out + (size_t)p * outstride;
			double g = 0.0;
			for (int i = 0; i < n; i++) g += v[i] * z[i];
			for (int i = 0; i < n; i++) z[i] -= g * v[i];
		}
		nz = vec_norm(z, n);
		if (!(nz > 0.0) || !isfinite(nz)) return 1;
		for (int i = 0; i < n; i++) z[i] /= nz;

		// Residual ||T z - lam z||, which for a converged vector is O(eps*||T||).
		double res = 0.0;
		for (int i = 0; i < n; i++) {
			double t = (d[i] - lam) * z[i];
			if (i > 0) t += e[i] * z[i - 1];
			if (i + 1 < n) t += e[i + 1] * z[i + 1];
			res += t * t;
		}
		if (sqrt(res) <= 1e-8 * (tnorm > 0.0 ? tnorm : 1.0)) return 0;
	}
	// Not converged.  Report it: returning success here would hand back a
	// plausible-looking but wrong direction, and a wrong projection is invisible
	// downstream -- the output is still finite, still smooth, still the right
	// shape.  The caller falls back to the full solve rather than using this.
	return 1;
}

// ---------------------------------------------------------------------------
// wrapper
// ---------------------------------------------------------------------------

// The header promises that a non-finite result is reported rather than returned.
// One O(n) sweep at the wrapper covers that, including the n == 1 fast path;
// putting it inside tql2 would miss the latter.
static int spectrum_finite(const double *ev, int n) {
	for (int i = 0; i < n; i++)
		if (!isfinite(ev[i])) return 0;
	return 1;
}

dn_eig *dn_eig_create(int n) {
	if (n < 1 || n > DN_MAX_VOL) return NULL;
	dn_eig *e = (dn_eig *)dn_calloc(1, sizeof(dn_eig));
	if (!e) return NULL;
	e->n = n;
	e->d = (double *)dn_malloc((size_t)n, sizeof(double));
	e->e = (double *)dn_malloc((size_t)n, sizeof(double));
	e->w = (double *)dn_malloc((size_t)n * n, sizeof(double));
	e->refl = (double *)dn_malloc((size_t)n * n, sizeof(double));
	e->hh = (double *)dn_malloc((size_t)n, sizeof(double));
	e->td = (double *)dn_malloc((size_t)n, sizeof(double));
	e->te = (double *)dn_malloc((size_t)n, sizeof(double));
	e->ev = (double *)dn_malloc((size_t)n, sizeof(double));
	e->z = (double *)dn_malloc((size_t)n, sizeof(double));
	// 4 doubles plus an int array per column, kept in one allocation; the int
	// block is sized generously so alignment is never an issue.
	e->gt = (double *)dn_malloc((size_t)n * 6, sizeof(double));
	if (!e->d || !e->e || !e->w || !e->refl || !e->hh || !e->td || !e->te ||
	    !e->ev || !e->z || !e->gt) { dn_eig_free(e); return NULL; }
#ifdef DN_USE_ACCELERATE
	// Size the LAPACK workspace ONCE per worker, not per voxel: both routines
	// report their optimum when called with lwork = -1, and this arena is reused
	// for millions of solves.  Take the larger of the two, and never less than
	// the documented minimum each guarantees to work with.
	{
		__LAPACK_int N = n, LDA = n, NC = n, LDC = n, QUERY = -1, INFO = 0;
		double wq1 = 0.0, wq2 = 0.0;
		dsytrd_("U", &N, e->refl, &LDA, e->td, e->te, e->hh, &wq1, &QUERY, &INFO);
		dormtr_("L", "U", "N", &N, &NC, e->refl, &LDA, e->hh, e->w, &LDC,
		        &wq2, &QUERY, &INFO);
		int want = (int)wq1 > (int)wq2 ? (int)wq1 : (int)wq2;
		if (want < 2 * n) want = 2 * n;
		e->lwork = want;
		e->lwk = (double *)dn_malloc((size_t)want, sizeof(double));
		if (!e->lwk) { dn_eig_free(e); return NULL; }
	}
#endif
	return e;
}

void dn_eig_free(dn_eig *e) {
	if (!e) return;
	free(e->d);
	free(e->e);
	free(e->w);
	free(e->refl);
	free(e->hh);
	free(e->td);
	free(e->te);
	free(e->ev);
	free(e->z);
	free(e->gt);
#ifdef DN_USE_ACCELERATE
	free(e->lwk);
#endif
	free(e);
}

// Full solve: tred2 to tridiagonal form, then tql2 with eigenvector accumulation.
static int dn_eig_solve(dn_eig *e, const double *a, double *evals, double *evecs) {
	if (!e || !a || !evals || !evecs) return -1;
	const int n = e->n;

	if (n == 1) {
		evals[0] = a[0];
		evecs[0] = 1.0;
		return spectrum_finite(evals, 1) ? 0 : -2;
	}

	memcpy(evecs, a, (size_t)n * n * sizeof(double));
	tred2(evecs, e->d, e->e, e->hh, n);
	const int rc = tql2(evecs, e->d, e->e, n);
	if (rc) return rc;
	memcpy(evals, e->d, (size_t)n * sizeof(double));
	return spectrum_finite(evals, n) ? 0 : -2;
}

int dn_eig_values(dn_eig *e, const double *a, double *evals) {
	if (!e || !a || !evals) return -1;
	const int n = e->n;
	e->staged = 0;
	e->full_valid = 0;

	if (n == 1) {
		if (!spectrum_finite(a, 1)) return -2;
		evals[0] = a[0];
		e->staged = 1;
		return 0;
	}

	// Below this size the split path loses: setting up the inverse iteration and
	// back-transforming costs more than simply accumulating every eigenvector.
	// Measured on this machine, 27x21 patches take 15.7 us split against 14.9 us
	// for a full solve, while 125x102 patches take 357 us against 872 us.  Two
	// data points, so the crossover is bracketed rather than pinpointed; 48 sits
	// inside the bracket and the small-n case is not where any time is spent.
	// dn_testgen's `big48` fixture is sized to exactly this number so that
	// `make test` reaches the split path at all -- it is the ONLY fixture that
	// does.  Raise this and resize the fixture in the same change: the two paths
	// agree to ~1e-12 while the fixture's pinned RMS tolerance is 1e-5, so a
	// fixture that quietly fell back to the full solver would keep passing while
	// covering nothing.
	#define DN_SPLIT_MIN 48

	if (n < DN_SPLIT_MIN) {
		// Below the crossover; do the whole thing and cache the vectors in `refl`.
		const int rc = dn_eig_solve(e, a, evals, e->refl);
		if (rc) return rc;   // dn_eig_solve already checked finiteness
		e->full_valid = 1;
		e->staged = 1;
		return 0;
	}

	memcpy(e->refl, a, (size_t)n * n * sizeof(double));
#ifdef DN_USE_ACCELERATE
	tred_reduce_accel(e, e->refl, e->td, e->te, e->hh, n);
#else
	tred_reduce(e->refl, e->td, e->te, e->hh, n);
#endif

	// tql2 destroys its inputs, so it works on the copies.
	memcpy(e->d, e->td, (size_t)n * sizeof(double));
	memcpy(e->e, e->te, (size_t)n * sizeof(double));
	const int rc = tql2(NULL, e->d, e->e, n);
	if (rc) return rc;

	if (!spectrum_finite(e->d, n)) return -2;
	memcpy(evals, e->d, (size_t)n * sizeof(double));
	memcpy(e->ev, e->d, (size_t)n * sizeof(double));

	e->tnorm = 0.0;
	for (int i = 0; i < n; i++) {
		double t = fabs(e->td[i]);
		if (i > 0) t += fabs(e->te[i]);
		if (i + 1 < n) t += fabs(e->te[i + 1]);
		if (t > e->tnorm) e->tnorm = t;
	}
	e->staged = 1;
	return 0;
}

// When inverse iteration will not converge, fall back to the path we already
// trust: run the full tridiagonal QL from an identity basis to get every
// eigenvector of T, then back-transform only the wanted columns.  This needs no
// state beyond what dn_eig_values() already kept -- in particular it does NOT
// need the Gram matrix back.
static int eig_vectors_full_fallback(dn_eig *e, int first, int count, double *evecs) {
	const int n = e->n;
	double *Z = e->w;

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Z[(size_t)i * n + j] = (i == j) ? 1.0 : 0.0;
	memcpy(e->d, e->td, (size_t)n * sizeof(double));
	memcpy(e->e, e->te, (size_t)n * sizeof(double));
	if (tql2(Z, e->d, e->e, n) != 0) return 1;

	for (int c = 0; c < count; c++) {
		const int k = first + c;
		for (int i = 0; i < n; i++) e->z[i] = Z[(size_t)i * n + k];
#ifdef DN_USE_ACCELERATE
		tred_apply_q_accel(e, e->refl, e->hh, e->z, 1, n);
#else
		tred_apply_q(e->refl, e->hh, e->z, n);
#endif
		for (int i = 0; i < n; i++) evecs[(size_t)i * n + k] = e->z[i];
	}
	e->fallbacks++;
	return 0;
}

int dn_eig_vectors(dn_eig *e, int first, int count, double *evecs) {
	if (!e || !evecs || !e->staged) return -1;
	const int n = e->n;
	if (first < 0 || count < 0 || first > n || count > n - first) return -1;
	if (count == 0) return 0;

	if (e->full_valid) {
		for (int k = first; k < first + count; k++)
			for (int i = 0; i < n; i++)
				evecs[(size_t)i * n + k] = e->refl[(size_t)i * n + k];
		return 0;
	}
	if (n == 1) { evecs[0] = 1.0; return 0; }

	// Build into e->w as contiguous rows so re-orthogonalisation is cheap, then
	// transpose into the caller's column layout.
	for (int c = 0; c < count; c++) {
		const int k = first + c;
		if (tridiag_eigvec(e->td, e->te, n, e->ev[k], e->tnorm, k,
		                   e->z, e->gt, e->w, c, n) != 0)
			return eig_vectors_full_fallback(e, first, count, evecs);
		memcpy(e->w + (size_t)c * n, e->z, (size_t)n * sizeof(double));
	}
#ifdef DN_USE_ACCELERATE
	tred_apply_q_accel(e, e->refl, e->hh, e->w, count, n);
#else
	for (int c = 0; c < count; c++)
		tred_apply_q(e->refl, e->hh, e->w + (size_t)c * n, n);
#endif
	for (int c = 0; c < count; c++)
		for (int i = 0; i < n; i++)
			evecs[(size_t)i * n + first + c] = e->w[(size_t)c * n + i];
	return 0;
}

unsigned long dn_eig_fallbacks(const dn_eig *e) {
	return e ? e->fallbacks : 0ul;
}
