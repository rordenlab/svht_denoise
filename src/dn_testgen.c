// Fixture generator and comparator for the test suite.  Not part of svht_denoise.
//
// This writes and reads NIfTI-1 by hand rather than through nifti_io, on purpose:
// a test whose inputs are produced by the library under test cannot detect a
// fault in that library.  Uncompressed .nii only, which is why it needs no zlib
// and why test.sh names every file .nii.
//
//   dn_testgen mk  <kind> <out.nii>            write a fixture ("kind@" -> .hdr/.img pair)
//   dn_testgen cmp <a.nii> <b.nii> <tol>       exit 0 if max|a-b| <= tol
//   dn_testgen rms <file.nii> <want> <reltol>  exit 0 if the RMS still matches
//
// Fixtures are 9x9x3 x 8 volumes (the "mask" kind is a single volume),
// deterministic, and every phase kind bar phase-siemens encodes the SAME
// underlying phase, so any two of them must denoise to the same answer.

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NX 9
#define NY 9
#define NZ 3
#define NT 8
#define N3 (NX * NY * NZ)
#define NVOX (N3 * NT)

typedef struct {
	int32_t sizeof_hdr; char data_type[10]; char db_name[18];
	int32_t extents; int16_t session_error; char regular; char dim_info;
	int16_t dim[8]; float intent_p1, intent_p2, intent_p3; int16_t intent_code;
	int16_t datatype, bitpix, slice_start; float pixdim[8];
	float vox_offset, scl_slope, scl_inter; int16_t slice_end;
	char slice_code, xyzt_units; float cal_max, cal_min, slice_duration, toffset;
	int32_t glmax, glmin; char descrip[80]; char aux_file[24];
	int16_t qform_code, sform_code;
	float quatern_b, quatern_c, quatern_d, qoffset_x, qoffset_y, qoffset_z;
	float srow_x[4], srow_y[4], srow_z[4];
	char intent_name[16]; char magic[4];
} nhdr1;

// A cheap deterministic generator: no rand(), so fixtures are identical
// everywhere and a failure always reproduces.
static double rng(uint32_t i) {
	uint32_t x = i * 2654435761u + 1013904223u;
	x ^= x >> 16; x *= 2246822519u; x ^= x >> 13; x *= 3266489917u; x ^= x >> 16;
	return (double)x / 4294967296.0;   // [0, 1)
}

// The shared underlying phase, in radians.  Two voxels are pinned to -pi and
// +pi so the field spans EXACTLY one turn.
//
// That is not cosmetic.  Real acquired phase covers the full circle, because its
// background is noise with uniform phase, and `auto` relies on it.  A fixture
// spanning 0.1% less would take the "already radians" branch in one encoding and
// the "normalise the span to one turn" branch in the others, leaving the two
// answers ~1% apart -- a property of the fixture, not of the code, but one that
// would look exactly like a bug in the comparison below.
static double phase_rad(int v) {
	if (v == 0) return -M_PI;
	if (v == 1) return M_PI;
	return (rng((uint32_t)v + 7777u) - 0.5) * 2.0 * M_PI;
}

// Two-file NIfTI-1: <prefix>.hdr plus <prefix>.img, magic "ni1", data in the
// .img.  Needed because the input/output name-resolution bug only exists for a
// prefix that the READER expands to an existing pair while the OUTPUT rule would
// expand it to a .nii; a fixture set of single .nii files cannot reach it.
static int write_pair(const char *prefix, const float *data, int nt, double vox);

static int write_nii_dim(const char *path, const float *data, int nx, int ny, int nz, int nt, double vox) {
	nhdr1 h;
	if (sizeof h != 348) { fprintf(stderr, "header is %zu bytes, not 348\n", sizeof h); return 1; }
	memset(&h, 0, sizeof h);
	h.sizeof_hdr = 348;
	h.regular = 'r';
	h.dim[0] = (int16_t)(nt > 1 ? 4 : 3);
	h.dim[1] = (int16_t)nx; h.dim[2] = (int16_t)ny; h.dim[3] = (int16_t)nz; h.dim[4] = (int16_t)nt;
	h.dim[5] = h.dim[6] = h.dim[7] = 1;
	h.datatype = 16;   // DT_FLOAT32
	h.bitpix = 32;
	h.pixdim[0] = 1.0f;
	h.pixdim[1] = h.pixdim[2] = h.pixdim[3] = (float)vox;
	h.pixdim[4] = 1.0f;
	h.vox_offset = 352.0f;
	h.scl_slope = 0.0f;         // "no scaling" per the standard
	h.xyzt_units = 2;           // NIFTI_UNITS_MM
	h.sform_code = 1;
	h.srow_x[0] = (float)vox; h.srow_y[1] = (float)vox; h.srow_z[2] = (float)vox;
	memcpy(h.magic, "n+1", 4);

	FILE *f = fopen(path, "wb");
	if (!f) { perror(path); return 1; }
	const char pad[4] = {0};
	const size_t n = (size_t)nx * (size_t)ny * (size_t)nz * (size_t)nt;
	int ok = (fwrite(&h, sizeof h, 1, f) == 1) && (fwrite(pad, 4, 1, f) == 1) &&
	         (fwrite(data, sizeof(float), n, f) == n);
	if (fclose(f) != 0) ok = 0;
	if (!ok) { fprintf(stderr, "failed writing %s\n", path); return 1; }
	return 0;
}

static int write_nii(const char *path, const float *data, int nt, double vox) {
	return write_nii_dim(path, data, NX, NY, NZ, nt, vox);
}

static int write_pair(const char *prefix, const float *data, int nt, double vox) {
	char hdr[1024], img[1024];
	if ((size_t)snprintf(hdr, sizeof hdr, "%s.hdr", prefix) >= sizeof hdr) return 1;
	if ((size_t)snprintf(img, sizeof img, "%s.img", prefix) >= sizeof img) return 1;

	// Reuse write_nii's header construction by writing a .nii, then splitting it.
	if (write_nii(hdr, data, nt, vox)) return 1;
	FILE *f = fopen(hdr, "r+b");
	if (!f) { perror(hdr); return 1; }
	nhdr1 h;
	int ok = (fread(&h, sizeof h, 1, f) == 1);
	if (ok) {
		h.vox_offset = 0.0f;            // data lives in the .img, from byte 0
		memcpy(h.magic, "ni1", 4);
		ok = (fseek(f, 0, SEEK_SET) == 0) && (fwrite(&h, sizeof h, 1, f) == 1);
	}
	if (fclose(f) != 0) ok = 0;
	// Truncate the .hdr to exactly the header, and put the voxels in the .img.
	if (ok) { FILE *t = fopen(hdr, "rb"); char buf[352]; size_t n = 0;
		if (t) { n = fread(buf, 1, sizeof buf, t); fclose(t); }
		FILE *w = fopen(hdr, "wb");
		ok = w && n >= sizeof h && fwrite(buf, sizeof h, 1, w) == 1;
		if (w && fclose(w) != 0) ok = 0; }
	if (ok) { FILE *w = fopen(img, "wb");
		ok = w && fwrite(data, sizeof(float), (size_t)N3 * nt, w) == (size_t)N3 * nt;
		if (w && fclose(w) != 0) ok = 0; }
	if (!ok) fprintf(stderr, "failed writing the %s pair\n", prefix);
	return ok ? 0 : 1;
}

static float *read_nii(const char *path, int *nvox) {
	FILE *f = fopen(path, "rb");
	if (!f) { perror(path); return NULL; }
	nhdr1 h;
	if (fread(&h, sizeof h, 1, f) != 1 || h.sizeof_hdr != 348 || h.datatype != 16) {
		fprintf(stderr, "%s: not a float32 NIfTI-1 single file\n", path);
		fclose(f); return NULL;
	}
	long n = 1;
	for (int i = 1; i <= (h.dim[0] < 4 ? h.dim[0] : 4); i++) n *= (h.dim[i] > 0 ? h.dim[i] : 1);
	if (fseek(f, (long)h.vox_offset, SEEK_SET) != 0) { fclose(f); return NULL; }
	float *d = (float *)malloc((size_t)n * sizeof(float));
	if (!d) { fclose(f); return NULL; }
	if (fread(d, sizeof(float), (size_t)n, f) != (size_t)n) {
		fprintf(stderr, "%s: short read\n", path); free(d); fclose(f); return NULL;
	}
	fclose(f);
	*nvox = (int)n;
	return d;
}

int main(int argc, char **argv) {
	if (argc == 5 && !strcmp(argv[1], "cmp")) {
		int na = 0, nb = 0;
		float *a = read_nii(argv[2], &na), *b = read_nii(argv[3], &nb);
		if (!a || !b) { free(a); free(b); return 2; }
		if (na != nb) { fprintf(stderr, "different sizes: %d vs %d\n", na, nb); free(a); free(b); return 2; }
		double worst = 0.0;
		for (int i = 0; i < na; i++) {
			const double d = fabs((double)a[i] - (double)b[i]);
			if (d > worst) worst = d;
		}
		free(a); free(b);
		const double tol = atof(argv[4]);
		printf("%.6g\n", worst);
		return (worst <= tol) ? 0 : 1;
	}

	// A numerical anchor.  Every other comparison here pits the tool against
	// itself under a different phase encoding, so a change to the ARITHMETIC --
	// the binomial kernel, say -- shows up on both sides and cancels.  Pinning
	// one reduced statistic is what catches that.  RMS, not a byte comparison,
	// so the check survives last-bit differences between platforms.
	if (argc == 5 && !strcmp(argv[1], "rms")) {
		int n = 0;
		float *a = read_nii(argv[2], &n);
		if (!a) return 2;
		double s = 0.0;
		for (int i = 0; i < n; i++) s += (double)a[i] * (double)a[i];
		free(a);
		const double rms = sqrt(s / (n ? n : 1));
		const double want = atof(argv[3]), reltol = atof(argv[4]);
		printf("%.9g\n", rms);
		return (want != 0.0 && fabs(rms - want) / fabs(want) <= reltol) ? 0 : 1;
	}

	if (argc != 4 || strcmp(argv[1], "mk")) {
		fprintf(stderr, "usage: dn_testgen mk   <kind> <out.nii>\n"
		                "       dn_testgen cmp  <a.nii> <b.nii> <tol>\n"
		                "       dn_testgen rms  <file.nii> <expected> <reltol>\n");
		return 2;
	}

	// A trailing "@" on the kind means "write a .hdr/.img pair at this prefix"
	// instead of a single .nii.  Stripped here so the kind names below stay plain.
	char kbuf[64];
	if ((size_t)snprintf(kbuf, sizeof kbuf, "%s", argv[2]) >= sizeof kbuf) return 2;
	size_t klen = strlen(kbuf);
	const int as_pair = (klen && kbuf[klen - 1] == '@');
	if (as_pair) kbuf[--klen] = '\0';
	const char *kind = kbuf;

	float *d = (float *)malloc(NVOX * sizeof(float));
	if (!d) return 2;
	double vox = 1.0;
	int nt = NT, rc = 0;

	if (!strcmp(kind, "mag")) {
		// Volume 0 is brightest, so it is the static-phase reference.
		for (int t = 0; t < NT; t++)
			for (int v = 0; v < N3; v++)
				d[t * N3 + v] = (float)((t == 0 ? 900.0 : 500.0) + 100.0 * rng((uint32_t)(t * N3 + v)));
	} else if (!strcmp(kind, "mask")) {
		nt = 1;
		for (int v = 0; v < N3; v++) d[v] = (v % 3) ? 1.0f : 0.0f;
	} else if (!strcmp(kind, "phase-rad") || !strcmp(kind, "phase-int") ||
	           !strcmp(kind, "phase-turns") || !strcmp(kind, "phase-deg") ||
	           !strcmp(kind, "phase-affine")) {
		// One phase field, four encodings of it.  scale converts radians to the
		// stored unit; svht_denoise must undo whichever it is given.
		double scale = 1.0;
		if (!strcmp(kind, "phase-int"))   scale = 4096.0 / M_PI;
		if (!strcmp(kind, "phase-turns")) scale = 1.0 / (2.0 * M_PI);
		if (!strcmp(kind, "phase-deg"))   scale = 180.0 / M_PI;
		if (!strcmp(kind, "phase-affine")) vox = 2.0;   // same values, wrong grid
		for (int t = 0; t < NT; t++)
			for (int v = 0; v < N3; v++)
				d[t * N3 + v] = (float)(phase_rad(t * N3 + v) * scale);
	} else if (!strcmp(kind, "phase-siemens")) {
		// Integer-valued, as a scanner really writes it.  Kept separate from
		// phase-int (which is an exact rescaling of phase-rad, so the encodings can
		// be compared against each other) because rounding breaks that equality.
		// Its purpose is the degenerate case: with -phaseunits turns every value is
		// a whole number of turns, exp(i*2*pi*n) == 1, and the rotation vanishes.
		for (int t2 = 0; t2 < NT; t2++)
			for (int v = 0; v < N3; v++)
				d[t2 * N3 + v] = (float)((int)(phase_rad(t2 * N3 + v) * (4096.0 / M_PI)));
	} else if (!strcmp(kind, "mask-wrongdim")) {
		// One volume, so it passes the volume-count guard, but a different grid --
		// which is what reaches dn_mask_build's DIMENSION check.  The obvious
		// candidate (phase-affine) is the same shape and trips the volume guard
		// first, so that branch had no coverage at all.
		for (int v = 0; v < N3; v++) d[v] = 1.0f;
		rc = write_nii_dim(argv[3], d, NX - 2, NY, NZ, 1, vox);
		free(d);
		return rc;
	} else if (!strcmp(kind, "phase-halfturns")) {
		// Half-integer turns.  Under -phaseunits turns this gives exp(i*k*pi),
		// i.e. a per-voxel factor of +1 or -1 -- a REAL sign modulation, not the
		// identity.  It exists to pin the no-op warning's condition: a sin-only
		// test passes here and would then claim the output is the magnitude
		// unchanged, which is false by ~100%.
		for (int t2 = 0; t2 < NT; t2++)
			for (int v = 0; v < N3; v++)
				d[t2 * N3 + v] = (float)(0.5 * (double)((int)(rng((uint32_t)(t2 * N3 + v)) * 17.0) - 8));
	} else if (!strcmp(kind, "phase-nan")) {
		// Valid data, but one NaN in the sform: dn_grid_offset_mm must fail CLOSED.
		// nifti_dmat44_spatial_ok only inspects the 3x3 block, so a NaN origin
		// survives the reader and reaches the grid comparison.
		for (int t2 = 0; t2 < NT; t2++)
			for (int v = 0; v < N3; v++)
				d[t2 * N3 + v] = (float)phase_rad(t2 * N3 + v);
		if (write_nii(argv[3], d, nt, vox)) { free(d); return 1; }
		FILE *f = fopen(argv[3], "r+b");
		if (!f) { free(d); return 1; }
		const float nan_v = (float)NAN;
		const long sv_off = (long)(offsetof(nhdr1, srow_x) + 3 * sizeof(float));
		int ok2 = (fseek(f, sv_off, SEEK_SET) == 0) && (fwrite(&nan_v, 4, 1, f) == 1);
		if (fclose(f) != 0) ok2 = 0;
		free(d);
		// This branch writes its own file, so the "@" pair form cannot apply.
		// Say so rather than exiting 0 having written the wrong shape.
		if (as_pair) { fprintf(stderr, "phase-nan does not support the '@' pair form\n"); return 2; }
		return ok2 ? 0 : 1;
	} else if (!strcmp(kind, "phase-const")) {
		for (int i = 0; i < NVOX; i++) d[i] = 1234.0f;
	} else {
		fprintf(stderr, "unknown fixture kind '%s'\n", kind);
		free(d);
		return 2;
	}

	rc = as_pair ? write_pair(argv[3], d, nt, vox) : write_nii(argv[3], d, nt, vox);
	free(d);
	return rc;
}
