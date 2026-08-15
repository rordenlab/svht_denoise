// Gibbs ringing removal by local subvoxel shifts.

#ifndef DG_H
#define DG_H

#ifdef __cplusplus
extern "C" {
#endif

// Can these in-plane dimensions be unrung at all?  Returns 0 if so, non-zero
// having reported why.  Exposed so the CLI can refuse BEFORE denoising: the
// check needs only the header, and discovering it afterwards throws away the
// whole run.
int dn_degibbs_check(int nx, int ny);

// The team size that will ACTUALLY run: `requested`, capped by the core count,
// by the number of planes, and by a scratch-memory budget.  Exposed so the CLI
// reports the effective count rather than the requested one -- the same contract
// dn_effective_threads has for the denoiser, and for the same reason.
int dn_degibbs_threads(int nx, int ny, int nz, int nvol, int requested);

// Unring every (z, volume) plane of `data` in place.  `data` is the volume-major
// float32 buffer everything else here uses: data[voxel + volume*nx*ny*nz].
// Within-slice axes are x and y, fixed; see the note in dg.c.
// Returns 0 on success, non-zero on failure (already reported).
int dn_degibbs(float *data, int nx, int ny, int nz, int nvol, int nthreads);

#ifdef __cplusplus
}
#endif

#endif // DG_H
