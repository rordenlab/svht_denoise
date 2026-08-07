// symmetric eigensolver, behind a backend-neutral wrapper.

#ifndef DN_EIG_H
#define DN_EIG_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct dn_eig dn_eig;

// Workspace for repeated n x n solves.  One per worker thread.
dn_eig *dn_eig_create(int n);
void dn_eig_free(dn_eig *e);

// Split form, for callers that need the whole spectrum but only a few
// eigenvectors -- which is this tool's entire workload: the median needs every
// eigenvalue, the projection needs only the handful above the threshold.
//
// Computing all N eigenvectors when 8 of 102 are wanted is most of the cost of a
// full solve, so dn_eig_values() stops after the eigenvalues and dn_eig_vectors()
// then produces just the requested contiguous range by inverse iteration on the
// tridiagonal form, back-transformed through the stored Householder reflectors.
//
// dn_eig_vectors(e, first, count, evecs) fills evecs columns [first, first+count)
// and leaves the others untouched; it must be called after dn_eig_values() with
// the same workspace and no intervening solve.
//
// `a` is n*n, row-major, full symmetric storage, and ITS CONTENTS ARE DESTROYED.
// evals comes back ASCENDING; evecs[i*n+j] is element i of the eigenvector for
// evals[j], so eigenvectors are the COLUMNS.  Both return 0 on success and
// non-zero on failure to converge -- callers must check, because a silently
// unconverged basis would corrupt the projection with no other symptom.
int dn_eig_values(dn_eig *e, double *a, double *evals);
int dn_eig_vectors(dn_eig *e, int first, int count, double *evecs);

// How many times dn_eig_vectors() had to abandon inverse iteration and fall back
// to the full solve.  Expected to be zero; a non-zero count is not a wrong
// answer -- the fallback is the trusted path -- but it is worth knowing about,
// so the driver reports it.
unsigned long dn_eig_fallbacks(const dn_eig *e);

#ifdef __cplusplus
}
#endif

#endif // DN_EIG_H
