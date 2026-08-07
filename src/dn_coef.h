#ifndef DN_COEF_H
#define DN_COEF_H

#ifdef __cplusplus
extern "C" {
#endif

// Median of the Marcenko-Pastur distribution with ratio beta.
double dn_mp_median(double beta);

// lambda_*(beta), Eq. (11): the optimal threshold for a KNOWN noise level, as a
// multiple of sigma*sqrt(n).  omega(beta) = lambda_*(beta) / sqrt(mu_beta) is the
// Eq. (4) coefficient this tool actually thresholds on; dn_patch.c computes it
// once per run, reusing mu_beta for the sigma estimate.
double dn_lambda_star(double beta);

// 1 if beta is a usable aspect ratio (finite, 0 < beta <= 1).
int dn_beta_valid(double beta);

#ifdef __cplusplus
}
#endif

#endif // DN_COEF_H
