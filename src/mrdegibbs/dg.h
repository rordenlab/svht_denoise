// Gibbs ringing removal by local subvoxel shifts.

#ifndef DG_H
#define DG_H

#ifdef __cplusplus
extern "C" {
#endif

// Partial-Fourier factors this build implements.  DG_PF_FULL means symmetric
// k-space, i.e. the conventional method and the default.
#define DG_PF_FULL 1.0
#define DG_PF_6_8  0.75
#define DG_PF_7_8  0.875

// Can these in-plane dimensions be unrung at all, at this partial-Fourier
// factor?  Returns 0 if so, non-zero having reported why.  Exposed so the CLI
// can refuse BEFORE denoising: the check needs only the header, and discovering
// it afterwards throws away the whole run.
int dn_degibbs_check(int nx, int ny, double pf);

// The team size that will ACTUALLY run: `requested`, capped by the core count,
// by the number of planes, and by a scratch-memory budget.  Exposed so the CLI
// reports the effective count rather than the requested one -- the same contract
// dn_effective_threads has for the denoiser, and for the same reason.
int dn_degibbs_threads(int nx, int ny, int nz, int nvol, int requested, double pf);

// Unring every (z, volume) plane of `data` in place.  `data` is the volume-major
// float32 buffer everything else here uses: data[voxel + volume*nx*ny*nz].
// Within-slice axes are x and y, fixed; see the note in dg.c.
//
// `pf` is the partial-Fourier factor: DG_PF_FULL for symmetric k-space,
// DG_PF_6_8 or DG_PF_7_8, each with its corresponding RPG correction for the
// wider ringing asymmetric acquisition adds along y.
//
// The method assumes MAGNITUDE data, and negative input is truncated to zero to
// make that hold -- a no-op on a magnitude image, and lossy on a phase-rotated
// one, which is signed.  A NaN or -Inf input becomes zero too; +Inf is NOT
// touched, because the truncation is fmaxf and that is maxNum.  Destroying a NaN
// beats propagating it, and the CLI refuses non-finite input at read time
// anyway -- though note `-degibbs y` is handed the DENOISER's output buffer, not
// the one that check ran on.  Ringing correction still undershoots, so the
// OUTPUT can be negative on every path except 7/8, which clamps as its
// reference does.
// Returns 0 on success, non-zero on failure (already reported).
int dn_degibbs(float *data, int nx, int ny, int nz, int nvol, int nthreads, double pf);

#ifdef __cplusplus
}
#endif

#endif // DG_H
