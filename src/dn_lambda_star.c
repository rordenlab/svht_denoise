// lambda_star(beta) is the optimal hard threshold coefficient for a KNOWN noise
// level, expressed as a multiple of sigma*sqrt(n), for an m-by-n matrix with
// beta = m/n.  It is Equation (11) of the paper included in this package:
//
//   M. Gavish and D. L. Donoho, "The Optimal Hard Threshold for Singular Values
//   is 4/sqrt(3)", IEEE Trans. Inf. Theory 60(8):5040-5053, 2014.
//
// Contract:
//   - domain is 0 < beta <= 1;
//   - return NAN outside that domain, INCLUDING for a NaN input.  Note that
//     `beta > 0.0 && beta <= 1.0` written naively lets a NaN through, because
//     every comparison against a NaN is false.  dn_beta_valid() in dn_coef.h
//     already handles this;
//   - accurate to ~1e-12 absolute across the domain.
//
// It is called once per run, not per voxel.  Do not optimise it.

#include <math.h>

#include "dn_coef.h"

double dn_lambda_star(double beta) {
	if (!dn_beta_valid(beta)) return NAN;

	/* Equation (11), Gavish & Donoho (2014):
	 *
	 *   lambda_*(beta)^2 = 2(beta+1) + 8*beta / (beta + 1 + sqrt(beta^2 + 14*beta + 1))
	 *
	 * At beta = 1 this reduces to sqrt(16/3) = 4/sqrt(3), the square-matrix
	 * case the paper is titled for.  Every term is a sum of positive
	 * quantities, so there is no cancellation anywhere on 0 < beta <= 1 and a
	 * plain double evaluation lands at machine epsilon, far inside the
	 * required 1e-12.
	 */
	const double bp1 = beta + 1.0;
	const double disc = sqrt(beta * beta + 14.0 * beta + 1.0);
	return sqrt(2.0 * bp1 + 8.0 * beta / (bp1 + disc));
}
