// NIfTI input/output.
//
// The float32 conversion follows the pattern of niimath's dtifit.c / core.c
// (BSD-2-Clause, rordenlab/niimath): read in the file's own datatype, then apply
// scl_slope/scl_inter once and clear them, so downstream code never has to think
// about scaling again.

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_nii.h"

// Names only the datatypes we REJECT: this is called from the default arm of the
// conversion switch, so every type that switch handles is unreachable here.
static const char *dt_name(int dt) {
	switch (dt) {
	case DT_COMPLEX64:  return "complex64";
	case DT_COMPLEX128: return "complex128";
	case DT_COMPLEX256: return "complex256";
	case DT_RGB24:      return "rgb24";
	case DT_RGBA32:     return "rgba32";
	case DT_FLOAT128:   return "float128";
	default:            return "unknown";
	}
}

// Convert nim->data into a fresh float32 buffer, applying scl_slope/scl_inter.
static float *to_f32(const nifti_image *nim, const char *role, size_t nvox) {
	float *out = (float *)dn_malloc(nvox, sizeof(float));
	if (!out) return NULL;

	double slope = nim->scl_slope;
	double inter = nim->scl_inter;
	// A zero (or non-finite) slope means "no scaling" per the NIfTI standard.
	if (!isfinite(slope) || slope == 0.0) { slope = 1.0; inter = 0.0; }
	if (!isfinite(inter)) inter = 0.0;

	const void *src = nim->data;
	#define DN_CONVERT(TYPE) do { \
		const TYPE *p = (const TYPE *)src; \
		for (size_t i = 0; i < nvox; i++) out[i] = (float)(slope * (double)p[i] + inter); \
	} while (0)

	switch (nim->datatype) {
	case DT_UINT8:   DN_CONVERT(uint8_t);  break;
	case DT_INT8:    DN_CONVERT(int8_t);   break;
	case DT_INT16:   DN_CONVERT(int16_t);  break;
	case DT_UINT16:  DN_CONVERT(uint16_t); break;
	case DT_INT32:   DN_CONVERT(int32_t);  break;
	case DT_UINT32:  DN_CONVERT(uint32_t); break;
	case DT_INT64:   DN_CONVERT(int64_t);  break;
	case DT_UINT64:  DN_CONVERT(uint64_t); break;
	case DT_FLOAT32: DN_CONVERT(float);    break;
	case DT_FLOAT64: DN_CONVERT(double);   break;
	default:
		dn_err("%s image has datatype %s, which is not supported.\n", role, dt_name(nim->datatype));
		dn_err("  Supported: uint8/int8, int16/uint16, int32/uint32, int64/uint64, float32, float64.\n");
		dn_err("  Complex and RGB data are out of scope for this tool.\n");
		free(out);
		return NULL;
	}
	#undef DN_CONVERT
	return out;
}

int dn_image_read(const char *fname, const char *role, int require4d, dn_image *img) {
	if (!img) return 1;
	memset(img, 0, sizeof(*img));

	nifti_image *nim = nifti_image_read(fname, 1);
	if (!nim) {
		dn_err("unable to read %s image '%s'\n", role, fname);
		return 1;
	}
	if (!nim->data) {
		dn_err("%s image '%s' has no voxel data\n", role, fname);
		nifti_image_free(nim);
		return 1;
	}

	// Guard the geometry before trusting any product of the dimensions.
	if (nim->nx < 1 || nim->ny < 1 || nim->nz < 1 || nim->nt < 1) {
		dn_err("%s image '%s' has a degenerate shape (%lldx%lldx%lldx%lld)\n", role, fname,
		       (long long)nim->nx, (long long)nim->ny, (long long)nim->nz, (long long)nim->nt);
		nifti_image_free(nim);
		return 1;
	}
	if (nim->ndim > 4 && (nim->nu > 1 || nim->nv > 1 || nim->nw > 1)) {
		dn_err("%s image '%s' has more than four dimensions; this tool handles 3D and 4D only\n",
		       role, fname);
		nifti_image_free(nim);
		return 1;
	}
	if (nim->nx > INT32_MAX || nim->ny > INT32_MAX || nim->nz > INT32_MAX || nim->nt > INT32_MAX) {
		dn_err("%s image '%s' has a dimension beyond this build's limits\n", role, fname);
		nifti_image_free(nim);
		return 1;
	}

	size_t nvox3d, nvox;
	if (dn_mul_size((size_t)nim->nx, (size_t)nim->ny, &nvox3d) ||
	    dn_mul_size(nvox3d, (size_t)nim->nz, &nvox3d) ||
	    dn_mul_size(nvox3d, (size_t)nim->nt, &nvox)) {
		dn_err("%s image '%s' is too large for this build's address space\n", role, fname);
		nifti_image_free(nim);
		return 1;
	}
	if ((int64_t)nvox != nim->nvox) {
		dn_err("%s image '%s' has an inconsistent header (nvox %lld, dims give %zu)\n",
		       role, fname, (long long)nim->nvox, nvox);
		nifti_image_free(nim);
		return 1;
	}

	const int nvol = (int)nim->nt;
	if (require4d && nvol < 2) {
		dn_err("%s image '%s' has %d volume%s; denoising needs a 4D series with at least 2\n",
		       role, fname, nvol, nvol == 1 ? "" : "s");
		nifti_image_free(nim);
		return 1;
	}
	if (require4d && nvol > DN_MAX_VOL) {
		dn_err("%s image '%s' has %d volumes; this build's limit is %d\n",
		       role, fname, nvol, DN_MAX_VOL);
		nifti_image_free(nim);
		return 1;
	}

	float *data = to_f32(nim, role, nvox);
	if (!data) { nifti_image_free(nim); return 1; }

	// Every input voxel can enter a patch, so a single non-finite sample would
	// silently poison a whole neighbourhood of outputs.  Reject globally.
	for (size_t i = 0; i < nvox; i++) {
		if (!isfinite(data[i])) {
			dn_err("%s image '%s' contains a non-finite value at index %zu.\n", role, fname, i);
			dn_err("  Every input voxel can enter a patch, so this is rejected rather than\n");
			dn_err("  propagated into its neighbourhood.\n");
			free(data);
			nifti_image_free(nim);
			return 1;
		}
	}

	// The scaling has been folded in; make sure nothing downstream re-applies it.
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;

	// The image is now held as float32 in img->data, so the raw buffer that
	// nifti_image_read allocated is dead weight.  Release it here rather than
	// leaving it for nifti_image_free: the write path swaps an output buffer
	// into nim->data and clears the pointer afterwards, so whatever was there
	// would otherwise be orphaned (measured: 8.1 MB leaked on a test
	// fixture).  nifti_image_free tolerates a NULL data pointer.
	free(nim->data);
	nim->data = NULL;

	img->nim = nim;
	img->data = data;
	img->nx = (int)nim->nx;
	img->ny = (int)nim->ny;
	img->nz = (int)nim->nz;
	img->nvol = nvol;
	img->nvox3d = nvox3d;
	return 0;
}

void dn_image_free(dn_image *img) {
	if (!img) return;
	free(img->data);
	if (img->nim) nifti_image_free(img->nim);
	memset(img, 0, sizeof(*img));
}

uint8_t *dn_mask_build(const dn_image *ref, const char *fname, size_t *n_in_mask) {
	dn_image m;
	if (dn_image_read(fname, "mask", 0, &m)) return NULL;

	if (m.nvol != 1) {
		dn_err("mask '%s' has %d volumes; a mask must be 3D (or 4D with one volume)\n",
		       fname, m.nvol);
		dn_image_free(&m);
		return NULL;
	}
	if (m.nx != ref->nx || m.ny != ref->ny || m.nz != ref->nz) {
		dn_err("mask '%s' is %dx%dx%d but the input is %dx%dx%d\n",
		       fname, m.nx, m.ny, m.nz, ref->nx, ref->ny, ref->nz);
		dn_image_free(&m);
		return NULL;
	}

	uint8_t *out = (uint8_t *)dn_malloc(ref->nvox3d, sizeof(uint8_t));
	if (!out) { dn_image_free(&m); return NULL; }

	// dn_image_read has already rejected every non-finite sample, so "> 0" is the
	// whole of the truth test.  (NaN > 0 is false in any case.)
	size_t count = 0;
	for (size_t i = 0; i < ref->nvox3d; i++) {
		const int on = m.data[i] > 0.0f;
		out[i] = (uint8_t)on;
		count += (size_t)on;
	}
	dn_image_free(&m);

	if (count == 0) {
		dn_err("mask '%s' selects no voxels\n", fname);
		free(out);
		return NULL;
	}
	if (n_in_mask) *n_in_mask = count;
	return out;
}

// Resolve a user-supplied output prefix to the two files nifti will actually
// create.  This exists because the names the user types are NOT the names that
// get written: nifti_set_filenames() appends an extension when it does not
// recognise one ("out" -> "out.nii"), and rewrites .img <-> .hdr so a two-file
// pair always resolves to both members.  A collision check that compares the
// typed strings therefore misses "in.nii" vs "in", which land on the same file.
//
// The derivation is delegated to nifti_set_filenames() itself, on a throwaway
// header, so this answer cannot drift from what write_image() does later.  It
// depends on nifti_type because that is what decides .nii versus .hdr/.img for
// a prefix with no extension, so pass the type of the image being written.
//
// Both names are heap-allocated and belong to the caller.  Returns 0 on success.
int dn_resolve_names(const char *prefix, int nifti_type, char **hdr, char **img) {
	if (hdr) *hdr = NULL;
	if (img) *img = NULL;
	if (!prefix || !hdr || !img) return 1;

	// Only nifti_type is read; fname/iname must start NULL because
	// nifti_set_filenames() frees whatever it finds there.
	nifti_image probe;
	memset(&probe, 0, sizeof probe);
	probe.nifti_type = nifti_type;

	if (nifti_set_filenames(&probe, prefix, 0, 0) != 0) {
		free(probe.fname);
		free(probe.iname);
		return 1;
	}
	*hdr = probe.fname;
	*img = probe.iname;
	return 0;
}

// Write `data` using ref's geometry.
//
// This MUTATES ref->nim -- datatype, dimensions, scaling, filenames and the data
// pointer -- and does not restore it.  The parameter is deliberately non-const
// to say so: an earlier version took `const dn_image *` and mutated it anyway,
// which made the declaration a lie.
//
// Call ORDER does not matter: every field this reads is either untouched (the
// spatial dims, qform/sform) or re-established on entry, and nvox is recomputed
// from ref->nvox3d rather than the mutated nim->nvox, so each call is
// self-contained.  What is not safe is reading geometry back out of `ref`
// afterwards and expecting the input's values.
static int write_image(dn_image *ref, const char *fname,
                       void *data, int datatype, int nbyper, int nvol) {
	nifti_image *nim = ref->nim;

	nim->datatype = datatype;
	nim->nbyper = nbyper;
	nim->data = data;

	nim->nt = nvol;
	nim->dim[4] = nvol;
	nim->ndim = (nvol > 1) ? 4 : 3;
	nim->dim[0] = nim->ndim;
	nim->nu = nim->nv = nim->nw = 1;
	nim->dim[5] = nim->dim[6] = nim->dim[7] = 1;
	nim->nvox = (int64_t)ref->nvox3d * nvol;

	// Scaling was folded into the values; leaving a stale slope would rescale
	// them again on the next read.
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;
	// cal_min/cal_max describe the INPUT's display range and are meaningless for
	// a noise map or a rank map, and wrong for a rescaled denoised image.
	nim->cal_min = 0.0;
	nim->cal_max = 0.0;

	if (nifti_set_filenames(nim, fname, 0, 1) != 0) {
		dn_err("cannot use '%s' as an output name\n", fname);
		nim->data = NULL;
		return 1;
	}
	const int rc = nifti_image_write_status(nim);
	nim->data = NULL;   // the buffer belongs to the caller, not to nifti_image_free
	if (rc != 0) {
		dn_err("failed to write '%s'\n", fname);
		return 1;
	}
	return 0;
}

int dn_write_f32(dn_image *ref, const char *fname, const float *data, int nvol) {
	return write_image(ref, fname, (void *)(uintptr_t)data, DT_FLOAT32, 4, nvol);
}

int dn_write_u16(dn_image *ref, const char *fname, const uint16_t *data) {
	return write_image(ref, fname, (void *)(uintptr_t)data, DT_UINT16, 2, 1);
}
