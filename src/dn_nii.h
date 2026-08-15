// NIfTI input/output.  Almost the only file that knows about nifti_io.h;
// everything else works on plain float32 buffers, with one exception --
// dn_phase.c reads nim->sform_code / qform_code directly to decide whether a
// world transform is present at all.

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
	// The nifti_type the INPUT was read with, kept separately because nim's own
	// copy does not survive the first write: nifti_set_filenames() rewrites it to
	// match whatever names it just built.  Since all outputs share one nim, a
	// single .hdr/.img output would otherwise silently switch every LATER
	// extensionless output from .nii to .hdr/.img -- names the collision check
	// never modelled, which is how an output came to land on an input at exit 0.
	int nifti_type;
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

// Resolve an INPUT path to the header and image files the reader will really
// OPEN.  This is not the same question dn_resolve_names answers: that one
// derives where an output would be WRITTEN, and the two diverge exactly where it
// matters.  Given "x" with an existing x.hdr/x.img pair on disk, the reader
// opens the pair while the output derivation says "x.nii", so a collision check
// built on the latter declares "x" and "x.hdr" to be different files and lets an
// output land on an input.  Delegated to the reader itself, header only, so the
// answer cannot drift from what the real read will do.
// Caller frees both.  Returns 0 on success, non-zero if the image cannot be
// opened at all -- which is worth failing on early, before anything is written.
int dn_resolve_input_names(const char *path, char **hdr, char **img);

// How far apart, in millimetres, do two headers place the same voxel?
//
// Matching dimensions do not establish that two images sample the same physical
// grid: they can be the same size and still differ in voxel size, orientation or
// origin, and combining them then mixes unrelated tissue into a plausible-looking
// result.  Measured over the eight volume corners, because that is the only
// comparison that means anything -- element-wise equality of the affines both
// over-rejects (two encodings of one transform) and under-rejects (a small
// matrix change is a large displacement far from the origin).  Both headers are
// normalised to millimetres via xyz_units first.
//
// Follows niimath's md_same_grid (BSD-2-Clause, rordenlab/niimath), including
// its gate.  Corners are taken from `a`; compare dimensions first.
//
// Compare against DN_GRID_TOL_MM, and write the test so that a non-finite offset
// fails closed.  Real magnitude/phase pairs out of one reconstruction agree
// EXACTLY (measured at 0.0 mm on the EDDEN dataset), so the tolerance is there
// for headers that have been round-tripped, not for genuine disagreement.
#define DN_GRID_TOL_MM 0.001f
float dn_grid_offset_mm(const dn_image *a, const dn_image *b);

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
