// DWI denoising by optimal singular value hard thresholding.
//
// Implements the threshold of
//   M. Gavish and D. L. Donoho, "The Optimal Hard Threshold for Singular Values
//   is 4/sqrt(3)", IEEE Trans. Inf. Theory 60(8):5040-5053, 2014,
// applied patch-wise to a 4D diffusion series.  This is NOT Marchenko-Pastur PCA
// and is not an emulation of MRtrix3 dwidenoise; it is a different estimator
// that happens to share the patch geometry.

#include <limits.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>    // readlink, unlink

#include "dn.h"
#include "dn_nii.h"
#include "dn_patch.h"
#include "dn_phase.h"
#include "dn_run.h"
#ifdef DN_DEGIBBS
#include "mrdegibbs/dg.h"
#endif

#define DN_VERSION "0.1.20260813"

// -degibbs: no, yes (after denoising), or only (instead of denoising).
enum { DN_DG_NO = 0, DN_DG_YES, DN_DG_ONLY };

static void usage(void) {
	printf("svht_denoise %s -- DWI denoising by optimal singular value hard thresholding\n", DN_VERSION);
	printf("\n");
	printf("USAGE\n");
	printf("  svht_denoise <input> <output> [options]\n");
	printf("\n");
	printf("  <input>   4D NIfTI diffusion series, at least 2 volumes.  Real datatypes\n");
	printf("            only: a complex-valued NIfTI is rejected.  For complex data,\n");
	printf("            pass the magnitude here and the phase via -phase.\n");
	printf("            .nii and .nii.gz are supported; the reader also accepts\n");
	printf("            whatever else this NIfTI build supports.\n");
	printf("  <output>  denoised series, always written as float32\n");
	printf("\n");
	printf("OPTIONS\n");
	printf("  -phase <image>    treat <input> as the MAGNITUDE of a complex series and\n");
	printf("                    <image> as its phase, and rotate the two onto the real\n");
	printf("                    axis before denoising.  Noise in the rotated data is\n");
	printf("                    the original zero-mean Gaussian rather than the Rician\n");
	printf("                    of a magnitude image, so the low-SNR noise floor is\n");
	printf("                    avoided instead of denoised.  The phase must match the\n");
	printf("                    input in dimensions and world transform.\n");
	printf("                    The OUTPUT IS THEN SIGNED.  If the phase turns out to\n");
	printf("                    encode no rotation at all, that is reported rather than\n");
	printf("                    silently returning the magnitude.\n");
	printf("  -phaseunits <u>   how to read the phase values: radians, degrees, turns\n");
	printf("                    (or its synonym cycles), or auto.  Default auto, which\n");
	printf("                    maps the observed range onto one full turn -- except\n");
	printf("                    that a range already within 0.1 of 2*pi is kept as\n");
	printf("                    radians unchanged.  That is correct for any encoding\n");
	printf("                    that COVERS a turn, which every acquired\n");
	printf("                    phase image does -- its background is noise with\n");
	printf("                    uniform phase.  Say the unit explicitly for data that\n");
	printf("                    does not cover one: a range of -0.5 to 0.5 is half a\n");
	printf("                    radian or one whole cycle, and the file cannot say.\n");
	printf("                    The convention used is reported unless -quiet, so\n");
	printf("                    check it.\n");
	printf("  -real <image>     write the real-axis-rotated input, before denoising\n");
	printf("                    (float32).  Needs -phase.\n");
	printf("  -mask <image>     only denoise voxels where the mask is > 0.\n");
	printf("                    Patches still read the whole image, so a masked run is\n");
	printf("                    identical to an unmasked one inside the mask; voxels\n");
	printf("                    outside the mask are written as zero.\n");
	printf("  -noise <image>    write the estimated noise level (float32)\n");
	printf("  -rank <image>     write the number of retained components (uint16)\n");
	printf("  -extent <k>       patch side length, odd, k^3 > number of volumes.\n");
	printf("                    Default: the smallest odd k that satisfies this.\n");
#ifdef DN_DEGIBBS
	printf("  -degibbs <y/n/o>  remove Gibbs ringing by local subvoxel shifts: yes, no\n");
	printf("                    (default) or only.  \"yes\" denoises first and degibbses\n");
	printf("                    the result, which is the order these belong in; \"only\"\n");
	printf("                    skips the denoising entirely and accepts a 3D image.\n");
	printf("                    Do NOT use on partial Fourier acquisitions.\n");
	printf("                    Within-slice axes are x and y.  Adapted from MRtrix3\n");
	printf("                    mrdegibbs (MPL-2.0); method of Kellner et al. 2016.\n");
#endif
	printf("  -nthreads <n>     worker threads (default: all cores). -p is an alias.\n");
	printf("  -quiet            suppress the run summary on stderr\n");
	printf("  -version, --version   print the version and exit\n");
	printf("  -help, --help, -h     print this message and exit\n");
	printf("\n");
	printf("NOTES\n");
	printf("  An input of \"-\" is read from standard input.  It has no name on disk,\n");
	printf("  so it is exempt from the collision rule below.  \"-\" is NOT accepted as\n");
	printf("  an output: <output>, -noise, -rank and -real each need a filename.\n");
	printf("  Every input and output must be a distinct file, compared as the files\n");
	printf("  actually opened and written rather than as typed.  A collision is\n");
	printf("  refused rather than silently overwriting one of them.\n");
	printf("  -nthreads is a request: it is capped by the core count, by the amount of\n");
	printf("  work available, and by a scratch-memory budget.  The summary reports the\n");
	printf("  number that will really run.\n");
	printf("  The threshold is the unknown-noise-level form, tau = omega(beta) * median\n");
	printf("  singular value of the patch, so no noise level need be supplied.\n");
	printf("  .bval/.bvec are never read: this algorithm does not use gradient directions.\n");
	printf("  Non-finite input is rejected rather than propagated, because any input\n");
	printf("  voxel can enter a patch.\n");
}

// Follow a chain of symlinks by hand, into `out`.
//
// realpath() cannot do this job: it fails outright when the final component does
// not exist, which is exactly the case its caller is handling.  A DANGLING
// symlink -- `ln -s noise.nii out.nii` before noise.nii has been written --
// makes stat() fail on BOTH sides, so without this the link's own name is
// compared against its target's, the two look distinct, and one output silently
// overwrote the other at exit 0.
//
// Returns 0 when the chain was followed to a non-link, non-zero when it could
// not be: too many hops, an unreadable link, or a name that will not fit.  That
// distinction matters -- leaving `out` at a partly-resolved name and calling it
// resolved is NOT conservative, because the caller would then compare two names
// that the kernel will still resolve further at write time.
//
// The hop limit is above both platforms' own: macOS SYMLOOP_MAX is 32 and Linux
// follows 40, so a 33-40 link alias would have gone unresolved here and resolved
// during the write -- exactly the window this function exists to close.
#define DN_MAX_SYMLINK_HOPS 64
static int resolve_links(const char *in, char *out, size_t outsz) {
	if (snprintf(out, outsz, "%s", in) >= (int)outsz) return 1;
	for (int hops = 0; hops < DN_MAX_SYMLINK_HOPS; hops++) {
		struct stat st;
		if (lstat(out, &st) != 0 || !S_ISLNK(st.st_mode)) return 0;   // not a link: done
		char target[PATH_MAX];
		const ssize_t n = readlink(out, target, sizeof target - 1);
		if (n <= 0) return 1;
		target[n] = '\0';

		char next[PATH_MAX];
		if (target[0] == '/') {
			if (snprintf(next, sizeof next, "%s", target) >= (int)sizeof next) return 1;
		} else {
			// Relative targets are relative to the LINK's directory, not the cwd.
			char dir[PATH_MAX];
			if (snprintf(dir, sizeof dir, "%s", out) >= (int)sizeof dir) return 1;
			char *slash = strrchr(dir, '/');
			if (slash) {
				*slash = '\0';
				if (snprintf(next, sizeof next, "%s/%s", dir[0] ? dir : "/", target) >= (int)sizeof next) return 1;
			} else if (snprintf(next, sizeof next, "%s", target) >= (int)sizeof next) return 1;
		}
		if (snprintf(out, outsz, "%s", next) >= (int)outsz) return 1;
	}
	return 1;   // hop limit reached: a cycle, or a chain longer than any kernel follows
}

// Are these two paths the same file?  (The POLICY that they must not be lives
// on check_path_collisions below; this answers only the mechanical question.)
//
// By inode, always.  When both already exist that is a plain stat; when one or
// both do not -- the usual case for outputs -- they are CREATED empty, compared,
// and removed again.
//
// Creating them is the point.  Every attempt to answer this from the NAMES was
// wrong in a way that lost a file, because the rule belongs to the filesystem:
//   - case-insensitive mounts exist on both platforms, so strcmp under-matches
//     and strcasecmp over-matches, and which applies is per-mount;
//   - Apple filesystems are normalization-insensitive -- NFC and NFD spellings
//     of one name are one file -- INCLUDING when formatted case-sensitive;
//   - the case fold maps non-ASCII onto ASCII, so U+212A KELVIN SIGN and "k"
//     are one file, which defeats any "pure ASCII must be distinct" shortcut.
// Doing that properly by hand needs CoreFoundation, which this tool does not
// link and package_macos.sh actively verifies it does not.  The filesystem
// answers all of it, exactly as it will at write time, for about ten lines.
static int same_file(const char *a, const char *b) {
	// Resolve first, so that creating through a dangling symlink -- and removing
	// what we created -- both act on the target rather than on the link.
	char la[PATH_MAX], lb[PATH_MAX];
	if (resolve_links(a, la, sizeof la) || resolve_links(b, lb, sizeof lb)) return 1;

	struct stat sa, sb;
	const int had_a = (stat(la, &sa) == 0);
	const int had_b = (stat(lb, &sb) == 0);
	if (had_a && had_b) return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;

	const int fa = open(la, O_CREAT | O_RDWR, 0600);
	const int fb = (fa >= 0) ? open(lb, O_CREAT | O_RDWR, 0600) : -1;
	if (fa >= 0) close(fa);
	if (fb >= 0) close(fb);

	int same;
	if (fa < 0 || fb < 0) {
		// Could not create: a missing or unwritable directory, most likely.
		// Nothing can be written there either, so nothing can be lost.  Fall back
		// to an exact match and let the write report the real problem, which is a
		// far better message than a collision would be.
		same = (strcmp(la, lb) == 0);
	} else {
		struct stat xa, xb;
		same = (stat(la, &xa) == 0 && stat(lb, &xb) == 0 &&
		        xa.st_dev == xb.st_dev && xa.st_ino == xb.st_ino);
	}
	// ONLY what did not exist a moment ago.  Keying the cleanup off "the open
	// succeeded" instead deletes an existing INPUT, which is the opposite of this
	// function's purpose -- the suite caught it on the first run.
	if (!had_b && fb >= 0) unlink(lb);
	if (!had_a && fa >= 0) unlink(la);
	return same;
}

// Every input and output must be a different file.  Each output is written
// separately from its own buffer, so two sharing a name leave only the last; an
// output naming an INPUT destroys it.  Both at exit 0 without this.
//
// The comparison is over the names nifti will REALLY use, not the ones the user
// typed: the two differ whenever an extension is missing or is half of a
// .hdr/.img pair, and comparing the typed strings let `svht_denoise in.nii in`
// normalise the output onto the input and destroy it.  `nifti_type` decides .nii
// versus .hdr/.img for an extensionless OUTPUT prefix, and comes from the loaded
// input because that is what the outputs inherit.
//
// Each path carries its own label and its own input/output classification, in
// one struct.  An earlier version kept those in parallel arrays indexed by slot,
// where the index silently decided three separate things -- which resolver ran,
// which read/write error was printed, and which refusal text was used -- under
// an ordering invariant that lived only in a comment.  Nothing here may depend
// on the order of the table.
typedef struct {
	const char *label;    // how this path is named in messages
	const char *prefix;   // what the user typed, or NULL if unused
	int is_input;         // read from, therefore must survive the run
	char *hdr, *img;      // resolved; owned by this struct
} dn_path;

static int check_path_collisions(dn_path *p, int n, int nifti_type) {
	int rc = 0;

	// Inputs and outputs are resolved by DIFFERENT rules, and using the output
	// rule on an input is not a cosmetic slip: it is what let `-phase ph -noise
	// ph.hdr` overwrite an existing ph.hdr/ph.img pair, because "ph" was modelled
	// as the "ph.nii" an output would create and so never matched.  See
	// dn_resolve_input_names.
	//
	// The main input is resolved here like any other, rather than borrowing the
	// strings the reader already produced.  It is the same code path and gives
	// the same answer, for one extra header open measured at ~0.1 ms; an
	// ownership exception for a single row cost more than that in reader
	// attention.
	for (int i = 0; i < n; i++) {
		if (!p[i].prefix) continue;
		// "-" is stdin, which nifti accepts as an input.  There is no file on disk
		// for an output to collide with, and a header-only reopen of a stream the
		// reader has already drained cannot succeed -- so skip the row rather than
		// fail the run.  (Resolving it unconditionally is what briefly broke
		// `svht_denoise - out.nii < in.nii`.)
		if (p[i].is_input && p[i].prefix[0] == '-' && p[i].prefix[1] == '\0') continue;
		const int bad = p[i].is_input
		              ? dn_resolve_input_names(p[i].prefix, &p[i].hdr, &p[i].img)
		              : dn_resolve_names(p[i].prefix, nifti_type, &p[i].hdr, &p[i].img);
		if (bad) {
			if (p[i].is_input) dn_err("unable to read %s image '%s'\n", p[i].label, p[i].prefix);
			else dn_err("cannot use '%s' as the %s name\n", p[i].prefix, p[i].label);
			rc = 1;
			goto done;
		}
	}

	for (int i = 0; i < n; i++) {
		if (!p[i].hdr) continue;
		for (int j = i + 1; j < n; j++) {
			if (!p[j].hdr) continue;
			// Both members of each pair, since .hdr and .img can collide crosswise.
			const char *a[2] = {p[i].hdr, p[i].img};
			const char *b[2] = {p[j].hdr, p[j].img};
			for (int u = 0; u < 2; u++) {
				for (int w = 0; w < 2; w++) {
					if (!same_file(a[u], b[w])) continue;
					dn_err("the %s and %s paths both resolve to '%s'.\n",
					       p[i].label, p[j].label, a[u]);
					// Two INPUTS colliding costs nothing on disk, so the reason to
					// refuse is different: they are distinct images by definition
					// -- a magnitude is not its own phase -- and one file named for
					// both is a mistake that would otherwise produce quiet nonsense.
					if (p[i].is_input && p[j].is_input) {
						dn_err("  Refusing to run: these are two different images, so one\n");
						dn_err("  file cannot be both.\n");
					} else {
						dn_err("  Refusing to run: outputs are written separately, so this would\n");
						// Either side may be the input: nothing here may assume the
						// inputs come first in the table.
						dn_err("  discard one of them%s.\n",
						       (p[i].is_input || p[j].is_input) ? " and destroy an input" : "");
					}
					rc = 1;
					goto done;
				}
			}
		}
	}

done:
	for (int i = 0; i < n; i++) { free(p[i].hdr); free(p[i].img); }
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
	const char *fphase = NULL, *freal = NULL;
	dn_phase_units punits = DN_PHASE_AUTO;
	// Presence, not value: "-phaseunits auto" leaves punits at its default, so
	// testing the enum would let that one spelling slip past the check below
	// while every other spelling was rejected.
	int punits_given = 0;
	// extent_given, rather than `extent != 0`, is what separates "not supplied"
	// from "supplied as 0".  Sharing 0 for both meant an explicit `-extent 0`
	// silently selected the automatic size and reported it as "(auto)", while
	// every other value the help text calls invalid -- 2, -1 -- was rejected.
	int extent = 0, extent_given = 0, nthreads = 0, quiet = 0;
	int degibbs = DN_DG_NO;

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
			else if (!strcmp(a, "-phase")) fphase = v;
			else if (!strcmp(a, "-real")) freal = v;
			else if (!strcmp(a, "-phaseunits")) {
				punits_given = 1;
				if (!strcmp(v, "auto")) punits = DN_PHASE_AUTO;
				else if (!strcmp(v, "radians")) punits = DN_PHASE_RADIANS;
				else if (!strcmp(v, "degrees")) punits = DN_PHASE_DEGREES;
				else if (!strcmp(v, "turns") || !strcmp(v, "cycles")) punits = DN_PHASE_TURNS;
				else {
					dn_err("-phaseunits must be radians, degrees, turns/cycles or auto (got '%s')\n", v);
					return EXIT_FAILURE;
				}
			}
			else if (!strcmp(a, "-degibbs")) {
#ifdef DN_DEGIBBS
				if (!strcmp(v, "y")) degibbs = DN_DG_YES;
				else if (!strcmp(v, "n")) degibbs = DN_DG_NO;
				else if (!strcmp(v, "o")) degibbs = DN_DG_ONLY;
				else {
					dn_err("-degibbs must be y, n or o (got '%s')\n", v);
					return EXIT_FAILURE;
				}
#else
				(void)v;
				dn_err("-degibbs is not compiled into this build.\n");
				dn_err("  Rebuild with 'make DEGIBBS=1'.\n");
				return EXIT_FAILURE;
#endif
			}
			else if (!strcmp(a, "-extent")) {
				if (parse_int(v, &extent)) { dn_err("-extent needs an integer (got '%s')\n", v); return EXIT_FAILURE; }
				extent_given = 1;
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
	// -real names the rotated input, which only exists when there is a rotation.
	// Silently writing nothing, or writing an unrotated copy, would both be worse
	// than saying so.
	if (freal && !fphase) {
		dn_err("-real writes the real-axis-rotated input, which needs -phase.\n");
		return EXIT_FAILURE;
	}
	// "-" means stdin for an INPUT.  As an output it made nifti write the image
	// to stdout, which contradicts dn.h's stated invariant that nothing here
	// writes an image to stdout, and the collision check models it as the file
	// "-.nii" -- so the checker and the writer disagreed about where it goes.
	// Undocumented either way; rejecting keeps the invariant true.
	const char *sv_out[4] = {fout, fnoise, frank, freal};
	const char *sv_lbl[4] = {"output", "-noise", "-rank", "-real"};
	for (int i = 0; i < 4; i++) {
		if (sv_out[i] && sv_out[i][0] == '-' && sv_out[i][1] == '\0') {
			dn_err("\"-\" is an input (stdin), not an output; %s needs a filename.\n", sv_lbl[i]);
			return EXIT_FAILURE;
		}
	}
	// -degibbs o does no denoising, so everything describing the denoiser -- its
	// mask, its by-products, its patch size, and the phase rotation that exists
	// only to change the noise the denoiser sees -- has nothing to act on.
	// Ignoring them silently would give a run that looked fine and wrote fewer
	// files than were asked for.
	if (degibbs == DN_DG_ONLY && (fmask || fnoise || frank || fphase || freal || extent_given)) {
		dn_err("-degibbs o does no denoising, so -mask, -noise, -rank, -phase, -real\n");
		dn_err("  and -extent have nothing to act on.\n");
		return EXIT_FAILURE;
	}
	// A masked run leaves hard zeros outside the mask, and a hard zero edge is
	// precisely the discontinuity the Kellner method rings on.  Degibbsing that
	// both breaks the "outside the mask reads as zero" promise -- measured, 647
	// of 648 outside voxels came back non-zero, peaking at a quarter of the data
	// range -- and rings the mask boundary INWARDS, corrupting voxels that are
	// inside it.  Neither is recoverable by re-zeroing afterwards, so refuse the
	// combination rather than quietly returning something worse than either
	// stage alone.  `-degibbs o` already refuses -mask, for its own reason.
	if (degibbs == DN_DG_YES && fmask) {
		dn_err("-mask cannot be combined with -degibbs y.\n");
		dn_err("  Masking leaves hard zeros outside the mask, and ringing removal would\n");
		dn_err("  treat that edge as signal -- corrupting voxels inside the mask too.\n");
		dn_err("  Denoise with -mask, then run -degibbs o on the result if you need both.\n");
		return EXIT_FAILURE;
	}
	if (punits_given && !fphase) {
		dn_err("-phaseunits describes the -phase image, which was not given.\n");
		return EXIT_FAILURE;
	}
	int status = EXIT_FAILURE;
	dn_image img;
	uint8_t *mask = NULL;
	float *out = NULL, *noise = NULL;
	uint16_t *rank = NULL;

	// Not require4d in "only" mode: mrdegibbs accepts a 3D image, and the
	// >= 2 volume rule belongs to the denoiser, which is not running.
	if (dn_image_read(fin, "input", degibbs != DN_DG_ONLY, &img)) return EXIT_FAILURE;

	// After the read, not before: resolving an output prefix needs the input's
	// nifti_type.  Still ahead of every write, and ahead of the denoising itself.
	dn_path paths[] = {
		{"input",  fin,    1, NULL, NULL},
		{"-mask",  fmask,  1, NULL, NULL},
		{"-phase", fphase, 1, NULL, NULL},
		{"output", fout,   0, NULL, NULL},
		{"-noise", fnoise, 0, NULL, NULL},
		{"-rank",  frank,  0, NULL, NULL},
		{"-real",  freal,  0, NULL, NULL},
	};
	if (check_path_collisions(paths, (int)(sizeof paths / sizeof paths[0]),
	                          img.nifti_type)) goto done;

#ifdef DN_DEGIBBS
	// Before any work, not after: this needs only the header, and reaching it
	// from the -degibbs y call site meant a full denoise ran and was then thrown
	// away when the geometry turned out to be unusable.
	if (degibbs != DN_DG_NO && dn_degibbs_check(img.nx, img.ny)) goto done;
#endif

#ifdef DN_DEGIBBS
	// "Only" short-circuits everything the denoiser needs -- geometry, mask,
	// scratch budget -- and rewrites the input buffer in place, so there is no
	// second copy of the series.
	if (degibbs == DN_DG_ONLY) {
		if (nthreads < 1) nthreads = dn_default_threads();
		// Report what will really run, not what was asked for -- same contract the
		// denoising path has, and it was missing here.
		nthreads = dn_degibbs_threads(img.nx, img.ny, img.nz, img.nvol, nthreads);
		if (!quiet) {
			dn_err("input        : %s (%dx%dx%d, %d volume%s)\n", fin,
			        img.nx, img.ny, img.nz, img.nvol, img.nvol == 1 ? "" : "s");
			dn_err("degibbs      : only, x-y planes (no denoising)\n");
			dn_err("threads      : %d\n", nthreads);
		}
		if (dn_degibbs(img.data, img.nx, img.ny, img.nz, img.nvol, nthreads)) goto done;
		if (dn_write_f32(&img, fout, img.data, img.nvol)) goto done;
		status = EXIT_SUCCESS;
		goto done;
	}
#endif

	const int auto_extent = dn_auto_extent(img.nvol);
	const int extent_was_auto = !extent_given;
	if (!extent_given) {
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
	nthreads = dn_effective_threads(&g, n_in_mask, nthreads);

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
#ifdef DN_DEGIBBS
		// Its OWN effective count, not the denoiser's: degibbs caps by planes and
		// by its own scratch budget, so the two stages can legitimately run
		// different team sizes and one number would be wrong for one of them.
		if (degibbs == DN_DG_YES)
			dn_err("degibbs      : yes, x-y planes (%d threads)\n",
			        dn_degibbs_threads(img.nx, img.ny, img.nz, img.nvol, nthreads));
#endif
	}

	// Before anything is allocated for the denoising, and before the first write:
	// on failure this leaves img exactly as it was read, and the run stops with
	// nothing on disk.  Afterwards img.data is the real-valued rotated series and
	// every stage below is unchanged -- the denoiser neither knows nor cares that
	// its input is now signed.
	if (fphase && dn_rotate_to_real(&img, fphase, punits, quiet)) goto done;

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
	r.n_work = n_in_mask;
	r.nthreads = nthreads;

	if (dn_run_execute(&r)) goto done;

#ifdef DN_DEGIBBS
	// After the denoiser, never before: degibbsing first would alter the noise
	// structure the denoiser models.  -noise and -rank describe the denoising and
	// are untouched, as is -real, which is the input before either stage.
	if (degibbs == DN_DG_YES &&
	    dn_degibbs(out, img.nx, img.ny, img.nz, img.nvol, nthreads)) goto done;
#endif

	// Write the auxiliary maps first: if the disk fills, failing before the main
	// output is written is less confusing than leaving a complete-looking result
	// beside a missing noise map.
	if (noise && dn_write_f32(&img, fnoise, noise, 1)) goto done;
	if (rank && dn_write_u16(&img, frank, rank)) goto done;
	// img.data is the ROTATED input: dn_run_execute only reads it.
	if (freal && dn_write_f32(&img, freal, img.data, img.nvol)) goto done;
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
