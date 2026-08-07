// DWI denoising by optimal singular value hard thresholding.
//
// Implements the threshold of
//   M. Gavish and D. L. Donoho, "The Optimal Hard Threshold for Singular Values
//   is 4/sqrt(3)", IEEE Trans. Inf. Theory 60(8):5040-5053, 2014,
// applied patch-wise to a 4D diffusion series.  This is NOT Marchenko-Pastur PCA
// and is not an emulation of MRtrix3 dwidenoise; it is a different estimator
// that happens to share the patch geometry.

#include <limits.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dn.h"
#include "dn_nii.h"
#include "dn_patch.h"
#include "dn_run.h"

#define DN_VERSION "1.0.0"

static void usage(void) {
	printf("svht_denoise %s -- DWI denoising by optimal singular value hard thresholding\n", DN_VERSION);
	printf("\n");
	printf("USAGE\n");
	printf("  svht_denoise <input> <output> [options]\n");
	printf("\n");
	printf("  <input>   4D NIfTI diffusion series, at least 2 volumes.\n");
	printf("            .nii and .nii.gz are supported; the reader also accepts\n");
	printf("            whatever else this NIfTI build supports.\n");
	printf("  <output>  denoised series, always written as float32\n");
	printf("\n");
	printf("OPTIONS\n");
	printf("  -mask <image>     only denoise voxels where the mask is > 0.\n");
	printf("                    Patches still read the whole image, so a masked run is\n");
	printf("                    identical to an unmasked one inside the mask; voxels\n");
	printf("                    outside the mask are written as zero.\n");
	printf("  -noise <image>    write the estimated noise level (float32)\n");
	printf("  -rank <image>     write the number of retained components (uint16)\n");
	printf("  -extent <k>       patch side length, odd, k^3 > number of volumes.\n");
	printf("                    Default: the smallest odd k that satisfies this.\n");
	printf("  -nthreads <n>     worker threads (default: all cores). -p is an alias.\n");
	printf("  -quiet            suppress the run summary on stderr\n");
	printf("  -version, --version   print the version and exit\n");
	printf("  -help, --help, -h     print this message and exit\n");
	printf("\n");
	printf("NOTES\n");
	printf("  The threshold is the unknown-noise-level form, tau = omega(beta) * median\n");
	printf("  singular value of the patch, so no noise level need be supplied.\n");
	printf("  .bval/.bvec are never read: this algorithm does not use gradient directions.\n");
	printf("  Non-finite input is rejected rather than propagated, because any input\n");
	printf("  voxel can enter a patch.\n");
}

// Every input and output must be a different file.  Each output is written
// separately from its own buffer, so two sharing a name leave only the last; an
// output naming an INPUT destroys it.  Both at exit 0 without this.
//
// Existing files are compared by device+inode, which catches hard links and
// symlinks that no amount of string normalisation would.  Files that do not
// exist yet -- the usual case for outputs -- fall back to comparing the names,
// after realpath() on the directory so "out.nii" and "./out.nii" match.
static int same_file(const char *a, const char *b) {
	struct stat sa, sb;
	const int ea = (stat(a, &sa) == 0), eb = (stat(b, &sb) == 0);
	if (ea && eb) return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
	if (ea != eb) return 0;   // one exists, one does not: different files

	// Neither exists: compare directory-resolved names.
	char ra[PATH_MAX], rb[PATH_MAX];
	const char *p[2] = {a, b};
	char *out[2] = {ra, rb};
	for (int i = 0; i < 2; i++) {
		char dir[PATH_MAX], rdir[PATH_MAX];
		const char *slash = strrchr(p[i], '/');
		size_t dlen = slash ? (size_t)(slash - p[i]) : 0;
		if (slash && dlen == 0) dlen = 1;          // "/name" -> "/"
		if (dlen >= sizeof(dir)) dlen = sizeof(dir) - 1;
		if (slash) { memcpy(dir, p[i], dlen); dir[dlen] = '\0'; }
		else snprintf(dir, sizeof(dir), ".");
		if (!realpath(dir, rdir)) snprintf(rdir, sizeof(rdir), "%s", dir);
		// Truncation would make two distinct long paths compare equal and cause a
		// bogus refusal, so a name that does not fit is reported as "different"
		// rather than silently shortened.
		if (snprintf(out[i], PATH_MAX, "%s/%s", rdir, slash ? slash + 1 : p[i]) >= PATH_MAX)
			return 0;
	}
	return strcmp(ra, rb) == 0;
}

// Every input and output must be a different file.
//
// This compares the names nifti will REALLY write, not the ones the user typed.
// The two differ whenever an extension is missing or is half of a .hdr/.img
// pair, and comparing the typed strings let `svht_denoise in.nii in` normalise
// the output onto the input and destroy it, at exit 0.  So each path is first
// expanded through dn_resolve_names() into its header and image files, and every
// resulting file is compared against every other.
//
// `in_hdr`/`in_img` are the files the input was actually READ from, which nifti
// already resolved by searching the disk, so they are authoritative and are used
// as-is.  `nifti_type` decides .nii versus .hdr/.img for an extensionless prefix
// and is taken from the loaded input, because that is what the outputs inherit.
static int check_path_collisions(const char *in_hdr, const char *in_img,
                                 int nifti_type, const char *fmask,
                                 const char *fout, const char *fnoise,
                                 const char *frank) {
	const char *label[5] = {"input", "-mask", "output", "-noise", "-rank"};
	const char *prefix[5] = {NULL, fmask, fout, fnoise, frank};
	// Slot 0 borrows the reader's strings; slots 1-4 own what dn_resolve_names
	// allocated.  Keeping the two apart means only the latter get freed.
	char *owned_hdr[5] = {NULL}, *owned_img[5] = {NULL};
	const char *hdr[5] = {in_hdr}, *img[5] = {in_img};
	int rc = 0;

	for (int i = 1; i < 5; i++) {
		if (!prefix[i]) continue;
		if (dn_resolve_names(prefix[i], nifti_type, &owned_hdr[i], &owned_img[i]) != 0) {
			dn_err("cannot use '%s' as a %s name\n", prefix[i], label[i]);
			rc = 1;
			goto done;
		}
		hdr[i] = owned_hdr[i];
		img[i] = owned_img[i];
	}

	for (int i = 0; i < 5; i++) {
		if (!hdr[i]) continue;
		for (int j = i + 1; j < 5; j++) {
			if (!hdr[j]) continue;
			// Both members of each pair, since .hdr and .img can collide crosswise.
			const char *a[2] = {hdr[i], img[i]};
			const char *b[2] = {hdr[j], img[j]};
			for (int u = 0; u < 2; u++) {
				for (int w = 0; w < 2; w++) {
					if (!same_file(a[u], b[w])) continue;
					dn_err("the %s and %s paths both resolve to '%s'.\n",
					       label[i], label[j], a[u]);
					dn_err("  Refusing to run: outputs are written separately, so this would\n");
					dn_err("  discard one of them%s.\n", (i < 2) ? " and destroy an input" : "");
					rc = 1;
					goto done;
				}
			}
		}
	}

done:
	for (int i = 1; i < 5; i++) { free(owned_hdr[i]); free(owned_img[i]); }
	return rc;
}

static int parse_int(const char *s, int *out) {
	char *end = NULL;
	long v = strtol(s, &end, 10);
	if (!end || *end != '\0' || end == s) return 1;
	if (v < INT32_MIN || v > INT32_MAX) return 1;
	*out = (int)v;
	return 0;
}

int main(int argc, char *argv[]) {
	const char *fin = NULL, *fout = NULL;
	const char *fmask = NULL, *fnoise = NULL, *frank = NULL;
	int extent = 0, nthreads = 0, quiet = 0;

	for (int i = 1; i < argc; i++) {
		const char *a = argv[i];
		if (!strcmp(a, "-help") || !strcmp(a, "--help") || !strcmp(a, "-h")) {
			usage();
			return EXIT_SUCCESS;
		}
		if (!strcmp(a, "-version") || !strcmp(a, "--version")) {
			printf("svht_denoise %s\n", DN_VERSION);
			return EXIT_SUCCESS;
		}
		if (!strcmp(a, "-quiet")) { quiet = 1; continue; }

		if (a[0] == '-' && a[1] != '\0') {
			// Every remaining flag takes exactly one argument.
			if (i + 1 >= argc) {
				dn_err("option '%s' needs a value\n", a);
				return EXIT_FAILURE;
			}
			const char *v = argv[++i];
			if (!strcmp(a, "-mask")) fmask = v;
			else if (!strcmp(a, "-noise")) fnoise = v;
			else if (!strcmp(a, "-rank")) frank = v;
			else if (!strcmp(a, "-extent")) {
				if (parse_int(v, &extent)) { dn_err("-extent needs an integer (got '%s')\n", v); return EXIT_FAILURE; }
			} else if (!strcmp(a, "-nthreads") || !strcmp(a, "-p")) {
				if (parse_int(v, &nthreads) || nthreads < 1) {
					dn_err("-nthreads needs a positive integer (got '%s')\n", v);
					return EXIT_FAILURE;
				}
			} else {
				dn_err("unknown option '%s'\n", a);
				dn_err("  run '%s -help' for the supported options\n", DN_NAME);
				return EXIT_FAILURE;
			}
			continue;
		}

		if (!fin) fin = a;
		else if (!fout) fout = a;
		else {
			dn_err("unexpected extra argument '%s'\n", a);
			return EXIT_FAILURE;
		}
	}

	if (!fin || !fout) {
		dn_err("missing required arguments.\n");
		dn_err("  usage: svht_denoise <input> <output> [options];  -help for details\n");
		return EXIT_FAILURE;
	}
	int status = EXIT_FAILURE;
	dn_image img;
	uint8_t *mask = NULL;
	float *out = NULL, *noise = NULL;
	uint16_t *rank = NULL;

	if (dn_image_read(fin, "input", 1, &img)) return EXIT_FAILURE;

	// After the read, not before: resolving an output prefix needs the input's
	// nifti_type, and img.nim->fname/iname are the files the reader actually
	// opened.  Still ahead of every write, and ahead of the denoising itself.
	if (check_path_collisions(img.nim->fname, img.nim->iname, img.nim->nifti_type,
	                          fmask, fout, fnoise, frank)) goto done;

	const int auto_extent = dn_auto_extent(img.nvol);
	const int extent_was_auto = (extent == 0);
	if (extent == 0) {
		extent = auto_extent;
		if (extent == 0) {
			dn_err("no patch size up to %d fits %d volumes (needs k^3 > N)\n",
			       DN_MAX_EXTENT, img.nvol);
			goto done;
		}
	}

	dn_geom g;
	if (dn_geom_init(&g, img.nx, img.ny, img.nz, img.nvol, extent)) goto done;

	size_t n_in_mask = img.nvox3d;
	if (fmask) {
		mask = dn_mask_build(&img, fmask, &n_in_mask);
		if (!mask) goto done;
	}

	if (nthreads < 1) nthreads = dn_default_threads();
	// Report what will actually run, not what was asked for.
	nthreads = dn_effective_threads(&g, nthreads);

	if (!quiet) {
		dn_err("input        : %s (%dx%dx%d, %d volumes)\n",
		        fin, img.nx, img.ny, img.nz, img.nvol);
		dn_err("patch        : %dx%dx%d = %d voxels x %d volumes%s\n",
		        g.extent, g.extent, g.extent, g.m, g.nvol,
		        extent_was_auto ? " (auto)" : "");
		dn_err("beta         : %.9f\n", g.beta);
		dn_err("omega(beta)  : %.9f\n", g.omega);
		dn_err("threads      : %d\n", nthreads);
		if (mask)
			dn_err("mask         : %s (%zu of %zu voxels)\n", fmask, n_in_mask, img.nvox3d);
	}

	// Zero-initialised: voxels outside the mask are never visited and must read
	// as zero in every output.
	out = (float *)dn_calloc(img.nvox3d * (size_t)img.nvol, sizeof(float));
	if (!out) goto done;
	if (fnoise) {
		noise = (float *)dn_calloc(img.nvox3d, sizeof(float));
		if (!noise) goto done;
	}
	if (frank) {
		rank = (uint16_t *)dn_calloc(img.nvox3d, sizeof(uint16_t));
		if (!rank) goto done;
	}

	dn_run r;
	memset(&r, 0, sizeof(r));
	r.g = &g;
	r.img = img.data;
	r.mask = mask;
	r.out = out;
	r.noise = noise;
	r.rank = rank;
	r.nthreads = nthreads;

	if (dn_run_execute(&r)) goto done;

	// Write the auxiliary maps first: if the disk fills, failing before the main
	// output is written is less confusing than leaving a complete-looking result
	// beside a missing noise map.
	if (noise && dn_write_f32(&img, fnoise, noise, 1)) goto done;
	if (rank && dn_write_u16(&img, frank, rank)) goto done;
	if (dn_write_f32(&img, fout, out, img.nvol)) goto done;

	status = EXIT_SUCCESS;

done:
	free(out);
	free(noise);
	free(rank);
	free(mask);
	dn_image_free(&img);
	return status;
}
