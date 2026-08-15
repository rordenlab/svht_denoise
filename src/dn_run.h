// whole-volume driver.

#ifndef DN_RUN_H
#define DN_RUN_H

#include <stdint.h>

#include "dn_patch.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	const dn_geom *g;
	const float *img;       // nvox3d * nvol, volume-major
	const uint8_t *mask;    // nvox3d, or NULL for "everywhere"
	float *out;             // nvox3d * nvol, zero-initialised by the caller
	float *noise;           // nvox3d, or NULL
	uint16_t *rank;         // nvox3d, or NULL
	size_t n_work;      // voxels the mask actually selects (nvox3d when unmasked)
	int nthreads;
} dn_run;

// Denoise every requested voxel.  Returns 0 on success.
//
// Output is byte-identical for any nthreads: each voxel's arithmetic depends
// only on its own patch, in a fixed order, and no two workers ever touch the
// same output element.  Scheduling therefore cannot influence the result.
int dn_run_execute(const dn_run *r);

// Number of hardware threads to use by default.
int dn_default_threads(void);

// The team size that will ACTUALLY be used for this geometry: `requested`,
// capped by the core count, by the number of work chunks, and by a scratch-memory
// budget.  Exposed
// so the CLI reports the effective count rather than the requested one.
int dn_effective_threads(const dn_geom *g, size_t n_work, int requested);

#ifdef __cplusplus
}
#endif

#endif // DN_RUN_H
