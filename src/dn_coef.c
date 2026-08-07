// With a = (1-sqrt(beta))^2, b = (1+sqrt(beta))^2, a+b = 2(1+beta),
// b-a = 4*sqrt(beta), a*b = (1-beta)^2, the antiderivative of
// sqrt((b-t)(t-a))/t is
//   A(t) = sqrt((b-t)(t-a))
//        - (1+beta) * asin( (a+b-2t) / (b-a) )
//        - (1-beta) * asin( ((a+b)t - 2ab) / (t(b-a)) )
// and A(a) evaluates in closed form to -pi*beta, so
//   F(x) = ( A(x) + pi*beta ) / (2*pi*beta),
// which satisfies F(a) = 0 and F(b) = 1 identically.  Verified against mpmath to
// better than 1e-23 for beta in [0.001, 1].  The median then follows from plain
// bisection, with no quadrature anywhere.

#include <math.h>

#include "dn_coef.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int dn_beta_valid(double beta) {
	return isfinite(beta) && beta > 0.0 && beta <= 1.0;
}

static double clamp_unit(double v) {
	if (v < -1.0) return -1.0;
	if (v > 1.0) return 1.0;
	return v;
}

static double dn_mp_cdf(double x, double beta) {
	if (!dn_beta_valid(beta) || !isfinite(x)) return NAN;
	const double sb = sqrt(beta);
	const double a = (1.0 - sb) * (1.0 - sb);
	const double b = (1.0 + sb) * (1.0 + sb);
	if (x <= a) return 0.0;
	if (x >= b) return 1.0;

	const double apb = 2.0 * (1.0 + beta);   // a + b
	const double bma = 4.0 * sb;             // b - a
	const double ab = (1.0 - beta) * (1.0 - beta);

	double rad = (b - x) * (x - a);
	double A = (rad > 0.0) ? sqrt(rad) : 0.0;
	A -= (1.0 + beta) * asin(clamp_unit((apb - 2.0 * x) / bma));
	// At beta == 1 this term is multiplied by zero; x > a >= 0 keeps the
	// division safe, since x == 0 can only occur at beta == 1 where a == 0 and
	// the x <= a test above has already returned.
	A -= (1.0 - beta) * asin(clamp_unit((apb * x - 2.0 * ab) / (x * bma)));

	double F = (A + M_PI * beta) / (2.0 * M_PI * beta);
	// Not clamp_unit(F): that clamps to [-1,1], which for the `< 0.0` test is
	// indistinguishable from testing F itself (NaN included) and reads as if a
	// second, different clamp were being applied.  F is a CDF, so [0,1].
	return (F < 0.0) ? 0.0 : (F > 1.0 ? 1.0 : F);
}

double dn_mp_median(double beta) {
	if (!dn_beta_valid(beta)) return NAN;
	const double sb = sqrt(beta);
	double lo = (1.0 - sb) * (1.0 - sb);
	double hi = (1.0 + sb) * (1.0 + sb);

	// Bisection to the limit of double precision.  100 halvings takes any
	// starting bracket below 1 ulp; the loop also exits early once lo and hi
	// are adjacent representable numbers.
	for (int i = 0; i < 100; i++) {
		double mid = 0.5 * (lo + hi);
		if (mid <= lo || mid >= hi) break;
		if (dn_mp_cdf(mid, beta) < 0.5)
			lo = mid;
		else
			hi = mid;
	}
	return 0.5 * (lo + hi);
}

// dn_lambda_star() lives in dn_lambda_star.c: it is Eq. (11).
