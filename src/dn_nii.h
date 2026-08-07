// NIfTI input/output.  This is the only file that knows about
// nifti_io.h; everything else works on plain float32 buffers.

#ifndef DN_NII_H
#define DN_NII_H

#include <stdint.h>
#include <stddef.h>

#include "nifti_io.h"

#ifdef __cplusplus
extern "C" {
#endif

// A loaded image, converted to float32 with scl_slope/scl_inter already applied.
// Voxel order is NIfTI's own: VOLUME-MAJOR, data[voxel + volume*nvox3d].
typedef struct {
	nifti_image *nim;   // header, retained for geometry and for writing outputs
	float *data;        // nvox3d * nvol floats
	int nx, ny, nz;
	int nvol;
	size_t nvox3d;
} dn_image;

// Read `fname` and convert to float32.  `role` appears in error messages ("input",
// "mask").  If require4d, the image must have nvol >= 2.  Rejects complex, RGB and
// unsupported datatypes, and any non-finite sample.  Returns 0 on success.
int dn_image_read(const char *fname, const char *role, int require4d, dn_image *img);

void dn_image_free(dn_image *img);

// Build a byte mask from a mask image, validated against `ref`.  Mask truth is
// "finite and > 0".  The mask must be 3D, or 4D with a single volume, and must
// match all three spatial dimensions of `ref`.  Returns a freshly allocated
// nvox3d-length buffer, or NULL on failure (already reported).
uint8_t *dn_mask_build(const dn_image *ref, const char *fname, size_t *n_in_mask);

// Resolve a user-supplied output prefix to the header and image files nifti will
// really create ("out" -> "out.nii"; "x.hdr" -> "x.hdr" + "x.img").  Needed
// because a collision check on the typed strings misses names that normalise
// onto the same file.  Caller frees both.  Returns 0 on success.
int dn_resolve_names(const char *prefix, int nifti_type, char **hdr, char **img);

// Write outputs derived from `ref`'s geometry.  `nvol` volumes of float32, or a
// single volume of uint16.  Returns 0 on success.
// NOTE: these MUTATE *ref (datatype, dims, scaling, filenames, data pointer)
// and do not restore it. See the comment on write_image in dn_nii.c.
int dn_write_f32(dn_image *ref, const char *fname, const float *data, int nvol);
int dn_write_u16(dn_image *ref, const char *fname, const uint16_t *data);

#ifdef __cplusplus
}
#endif

#endif // DN_NII_H
