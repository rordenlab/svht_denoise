// DWI denoising by local PCA with the singular value rules of Gavish & Donoho:
// the optimal hard threshold (IEEE Trans. Inf. Theory 60(8):5040-5053, 2014) and
// the Frobenius-optimal shrinkage (63(4):2137-2152, 2017, the default), applied
// patch-wise to a 4D diffusion series.  This is NOT Marchenko-Pastur PCA and is
// not an emulation of MRtrix3 dwidenoise; the noise level comes from their median
// estimator, not a Marchenko-Pastur fit.

#include <limits.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>    // readlink, unlink

#include "dn.h"
#include "dn_bval.h"
#include "dn_nii.h"
#include "dn_patch.h"
#include "dn_phase.h"
#include "dn_run.h"
#include "dn_sigma.h"
#include "dn_tensor.h"
#include "dn_vst.h"
#ifdef DN_DEGIBBS
#include "mrdegibbs/dg.h"

// The factor as the fraction everyone writes, for the run summary: -help and the
// README both say 7/8 and 6/8, not 0.875 and 0.75.  Empty for full k-space.
static const char *pf_note(double pf) {
	return pf == DG_PF_6_8 ? ", partial Fourier 6/8"
	     : pf == DG_PF_7_8 ? ", partial Fourier 7/8" : "";
}
#endif

// Noise-map estimates under -vst: the first from raw magnitudes, each later one
// corrected in the stabilised domain.  Two, measured on Rician synthetic truth
// (benchmark/): RMSE at sigma 30 / 60 was 13.77 / 30.24 with one, 13.49 / 25.20
// with two, 13.75 / 25.30 with three, and worse beyond -- the correction does not
// converge on the true sigma but overshoots it (64 against 60 after six).
#define DN_VST_ROUNDS 2

#define DN_VERSION "0.1.20261001"

// -degibbs: no, yes (after denoising), or only (instead of denoising).
enum { DN_DG_NO = 0, DN_DG_YES, DN_DG_ONLY };

static void usage(void) {
	printf("svht_denoise %s -- DWI denoising by local PCA with optimal singular value shrinkage\n", DN_VERSION);
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
#ifdef DN_DEGIBBS
	printf("                    The OUTPUT IS THEN SIGNED, unless -degibbs is also\n");
	printf("                    given, which gives that up.  If the phase turns out to\n");
#else
	printf("                    The OUTPUT IS THEN SIGNED.  If the phase turns out to\n");
#endif
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
	printf("                    check it.  Needs -phase.\n");
	printf("  -real <image>     write the real-axis-rotated input, before denoising\n");
	printf("                    (float32).  Needs -phase.\n");
	printf("  -mask <image>     only denoise voxels where the mask is > 0.\n");
	printf("                    Patches still read the whole image, so a masked run is\n");
	printf("                    identical to an unmasked one inside the mask; voxels\n");
	printf("                    outside the mask are written as zero.\n");
	printf("  -noise <image>    write the estimated noise level (float32)\n");
	printf("  -rank <image>     write the number of retained components (uint16)\n");
	printf("  -bval <file>      FSL b-values, one per volume.  Default: the input's name\n");
	printf("                    with .nii/.nii.gz replaced by .bval.  Needs -demean y.\n");
	printf("  -demean <y/n>     remove each voxel's mean over each b-value shell before\n");
	printf("                    PCA and restore it after, which lowers the signal rank.\n");
	printf("                    Default n: measured, it does not help.  Shells are\n");
	printf("                    clustered as MRtrix3 does (b <= 10 is b=0; neighbours\n");
	printf("                    less than 80 s/mm^2 apart share a shell); a shell of one\n");
	printf("                    volume is not demeaned.\n");
	printf("  -filter <f>       how components are kept: optshrink (default), the\n");
	printf("                    Frobenius-optimal shrinkage of Gavish & Donoho 2017;\n");
	printf("                    optthresh, their 2014 unknown-noise hard threshold; or\n");
	printf("                    truncate, a hard cut at the noise bulk edge.  optshrink\n");
	printf("                    and truncate read sigma from a noise map estimated\n");
	printf("                    first, from larger patches on every %dth voxel.\n", DN_SIGMA_STRIDE);
	printf("  -shape <s>        patch shape: sphere (default) or cube.  The sphere is\n");
	printf("                    measured in mm, so it stays round on anisotropic voxels,\n");
	printf("                    and holds the fewest whole distance shells with at\n");
	printf("                    least k^3 voxels (k from -extent).  At the image edge it\n");
	printf("                    is the same number of voxels nearest the one denoised.\n");
	printf("  -aggregator <a>   how patches make the output: gaussian (default), every\n");
	printf("                    patch adds its estimate to every voxel it covers,\n");
	printf("                    weighted by distance from its centre (Manjon et al.\n");
	printf("                    2013); uniform, the same unweighted; or exclusive, each\n");
	printf("                    voxel from the one patch centred on it.\n");
	printf("  -aggregator_fwhm <f>  Gaussian width, in units of the patch-centre spacing\n");
	printf("                    (default 2).  Needs -aggregator gaussian.\n");
	printf("  -stride <s>       centre a patch on every s-th voxel per axis, 1 or 2,\n");
	printf("                    which runs ~1/8 as many.  Default 2, or 1 under\n");
	printf("                    -aggregator exclusive, which needs a patch per voxel.\n");
	printf("  -vst <y/n>        variance-stabilise MAGNITUDE data before denoising: map\n");
	printf("                    its Rician / non-central chi noise to unit-variance\n");
	printf("                    Gaussian (Foi 2011), and map back after.  Default y,\n");
	printf("                    or n with -phase, whose rotated data are already\n");
	printf("                    Gaussian.  By default the way back is the exact-unbiased\n");
	printf("                    inverse, which also removes the noise-floor bias.  Also\n");
	printf("                    n by default for input with any negative value, which\n");
	printf("                    cannot be a magnitude image.\n");
	printf("  -noise_dof <L>    receive channels behind each magnitude (sum of squares),\n");
	printf("                    so the noise is non-central chi with 2L degrees of\n");
	printf("                    freedom.  Default 1, i.e. Rician.  Needs -vst y.\n");
	printf("  -preserve_noise_bias  map back algebraically instead, keeping the noise-floor\n");
	printf("                    bias of an ordinary magnitude image.  Needs -vst y.\n");
	printf("  -extent <k>       patch size: cube side length, odd, k^3 > number of\n");
	printf("                    volumes; a sphere holds at least k^3 voxels.\n");
	printf("                    Default: the smallest odd k that satisfies this.\n");
#ifdef DN_DEGIBBS
	printf("  -degibbs <y/n/o>  remove Gibbs ringing by local subvoxel shifts: yes, no\n");
	printf("                    (default) or only.  \"yes\" denoises first and degibbses\n");
	printf("                    the result, which is the order these belong in; \"only\"\n");
	printf("                    skips the denoising entirely and accepts a 3D image.\n");
	printf("                    Use -pF 0.875 or 0.75 for supported partial-Fourier\n");
	printf("                    acquisitions; without it this assumes full k-space.\n");
	printf("                    Assumes MAGNITUDE data: negative input is truncated to\n");
	printf("                    zero, which is a no-op on a magnitude image and lossy on\n");
	printf("                    a -phase rotated one.  The output can still go negative,\n");
	printf("                    because ringing correction undershoots -- except under\n");
	printf("                    -pF 0.875, which clamps its output as its reference does.\n");
	printf("                    In-plane dimensions must be at least 5, and at least 5\n");
	printf("                    again after a -pF factor splits them.\n");
	printf("                    \"only\" also refuses -noise, -rank, -phase, -real and\n");
	printf("                    -extent, which all describe a denoiser it is not running.\n");
	printf("                    Cannot be combined with -mask: a mask leaves a hard\n");
	printf("                    zero edge, which is exactly what this method rings on.\n");
	printf("                    Within-slice axes are x and y.  Adapted from MRtrix3\n");
	printf("                    mrdegibbs (MPL-2.0); method of Kellner et al. 2016.\n");
	printf("  -pF <factor>      partial-Fourier factor: 0.875 (7/8) or 0.75 (6/8).\n");
	printf("                    Asymmetric k-space truncation adds a SECOND ringing along\n");
	printf("                    y that -degibbs alone leaves in place; this removes it too\n");
	printf("                    (Lee et al. 2021).  Needs -degibbs, assumes y is the\n");
	printf("                    phase-encode direction, and 6/8 needs an even y dimension.\n");
	printf("                    1.0 is accepted and means full k-space, the default.\n");
#endif
	printf("\n");
	printf("TENSOR MODE (multi-echo)\n");
	printf("  svht_denoise <echo1> <output> -tensor gd -echo <echo2> [-echo ...] [options]\n");
	printf("\n");
	printf("  Each patch is a space x time x echo TENSOR, denoised mode by mode in the\n");
	printf("  recursive scheme of Olesen et al. (Magn Reson Med 2022): each unfolding is\n");
	printf("  cut at the edge of its noise bulk and the last mode's kept values are\n");
	printf("  shrunk, unless -tensor_filter says otherwise.  The noise level\n");
	printf("  of each echo comes from the same Gavish & Donoho median-estimator map the\n");
	printf("  matrix denoiser reads, and each echo is put in units of its OWN noise\n");
	printf("  before the tensor is built (echoes rarely share one), then mapped back; no\n");
	printf("  Marchenko-Pastur fit is used.  Windows sit wherever they fit inside the\n");
	printf("  image, overlapping estimates are averaged uniformly, and a voxel no window\n");
	printf("  covers keeps its input value (to float rounding) and reads NaN in every map.\n");
	printf("  Echoes are 4D MAGNITUDE series on one grid with one TR, at least 3 volumes;\n");
	printf("  negative values are refused.  Echo k is written to <output> with\n");
	printf("  \"_echo-<k>\" before its extension: out.nii.gz -> out_echo-1.nii.gz, ...\n");
	printf("  -tensor gd            select tensor mode; the only tensor method\n");
	printf("  -echo <image>         the next echo, in echo order; repeat for each\n");
	printf("  -tensor_phase <image> one echo's phase, in echo order; repeat for each.\n");
	printf("                        Each echo is rotated onto the real axis with its\n");
	printf("                        own phase, as -phase does in matrix mode, so the\n");
	printf("                        OUTPUT IS SIGNED.  -phaseunits applies to all of them.\n");
	printf("  -tensor_window <w>    x,y,z for a local window (default 5,5,5), or x,y for\n");
	printf("                        an in-plane window spanning every slice\n");
	printf("  -tensor_stride <s>    window step: one value, or one per window size\n");
	printf("                        (default 1).  A step that misses the far edge leaves\n");
	printf("                        voxels no window covers.\n");
	printf("  -tensor_order <o>     order the modes are denoised in (default\n");
	printf("                        space,time,echo); shrinkage applies to the last\n");
	printf("                        mode\n");
	printf("  -tensor_sigma2 <v>    known noise VARIANCE of every echo, instead of the maps\n");
	printf("  -tensor_noise <n>     echo (default): each echo in units of its OWN noise\n");
	printf("                        map; or pooled: one map, the mean variance over the\n");
	printf("                        echoes, for all of them.  Pooled is the physical\n");
	printf("                        expectation when every echo is reconstructed alike;\n");
	printf("                        echo-specific parallel-imaging reconstructions can\n");
	printf("                        give each echo its own noise level, which the summary's\n");
	printf("                        per-echo medians show.\n");
	printf("  -tensor_vst <y/n>     variance-stabilize each MAGNITUDE echo against its\n");
	printf("                        map (Rician; exact-unbiased inverse), as -vst does in\n");
	printf("                        matrix mode, rather than only scaling it by the map.\n");
	printf("                        Default y for magnitude; refused with -tensor_phase\n");
	printf("                        (already Gaussian) and -tensor_sigma2 (no map).\n");
	printf("  -tensor_filter <f>    how each mode keeps components: optshrink (default),\n");
	printf("                        cut at the bulk edge with the last mode's kept values\n");
	printf("                        shrunk (Gavish & Donoho 2017); truncate, the same cut\n");
	printf("                        with nothing shrunk; or optthresh, Gavish & Donoho\n");
	printf("                        2014's KNOWN-noise optimal hard threshold at every\n");
	printf("                        mode, lambda*(beta) sigma sqrt(n), nothing shrunk.\n");
	printf("                        All three read the same per-echo noise; matrix\n");
	printf("                        optthresh's unknown-noise rule needs a median over\n");
	printf("                        many noise values, which the echo mode lacks.\n");
	printf("  -tensor_demean <y/n>  center each unfolding along its longer side (default n)\n");
	printf("  -tensor_center <y/n>  write only each window's center voxel (its center\n");
	printf("                        column for an x,y window) instead of averaging\n");
	printf("                        (default n).  With -mask, a window whose center is\n");
	printf("                        outside the mask is skipped, as in ME-PCA.\n");
	printf("  -mask <image>         denoise only windows that touch the mask.  Unlike\n");
	printf("                        matrix mode, such a window writes every voxel it\n");
	printf("                        covers, and voxels outside the mask are not zeroed.\n");
	printf("  -noise <image>        noise standard deviation, one volume per echo: that\n");
	printf("                        echo's map where a window reached (float32, NaN where\n");
	printf("                        uncovered)\n");
	printf("  -rank <image>         mean kept rank per mode, one volume per mode in\n");
	printf("                        -tensor_order (float32; the header names them)\n");
	printf("  -snr_gain <image>     Olesen et al.'s SNR-gain heuristic, from the mean ranks;\n");
	printf("                        NaN where nothing was kept (float32)\n");
	printf("  None of the matrix denoiser's options apply, and each is refused.\n");
	printf("\n");
	printf("COMMON OPTIONS\n");
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
	printf("  No noise level need be supplied: it is estimated from the data.\n");
	printf("  '-filter optthresh -aggregator exclusive -shape cube -vst n' reproduces\n");
	printf("  the original pipeline exactly.\n");
	printf("  .bval is read only to group volumes into shells for -demean; gradient\n");
	printf("  directions (.bvec) are never used.\n");
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
	int is_input;         // read from, so must survive the run; 2: a raw file, not NIfTI
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
		// A raw file (the .bval) is opened exactly as named: both "halves" are it.
		const int bad = p[i].is_input == 2
		              ? !(p[i].hdr = strdup(p[i].prefix)) || !(p[i].img = strdup(p[i].prefix))
		              : p[i].is_input
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

// One stabilising pass over the live voxels: magnitude / sigma -> f.  `prev` is
// the sigma each voxel is currently in units of, ignored when `first`; later
// passes recover the raw magnitude through f's exact algebraic inverse instead of
// keeping a second copy of the series.  `prev` is updated to `sigma`.
static void vst_stabilise(const dn_vst *vt, float *data, size_t nvox3d, int nvol, const uint8_t *live,
                          const float *sigma, float *prev, int first) {
	for (size_t v = 0; v < nvox3d; v++) {
		if (!live[v]) continue;
		const double s = sigma[v], p = first ? 0.0 : prev[v];
		for (int j = 0; j < nvol; j++) {
			float *x = data + v + (size_t)j * nvox3d;
			const double raw = first ? *x : p * dn_vst_inverse(vt, *x, 0);
			*x = (float)dn_vst_forward(vt, raw / s);
		}
		prev[v] = (float)s;
	}
}

// Map a denoised, stabilised series back to magnitude: every sample of a live,
// in-mask voxel -- a denoised 0 is a value too, and the algebraic inverse maps it
// onto the noise floor.  The denoiser saw sigma = 1, so the noise map it wrote is
// replaced by the data's own.
static void vst_restore(const dn_vst *vt, float *out, float *noise, size_t nvox3d, int nvol,
                        const uint8_t *live, const uint8_t *mask, const float *sigma, int unbiased) {
	for (size_t v = 0; v < nvox3d; v++) {
		const double s = sigma[v];
		if (live[v] && (!mask || mask[v]))
			for (int j = 0; j < nvol; j++) {
				float *o = out + v + (size_t)j * nvox3d;
				*o = (float)(s * dn_vst_inverse(vt, *o, unbiased));
			}
		if (noise && noise[v] != 0.0f) noise[v] = (float)s;
	}
}

// Index of `v` in the NULL-terminated `names`, into *out; the index is the option's
// enum value.  Returns non-zero, reported, for anything else.
static int parse_choice(const char *flag, const char *v, const char *const *names, int *out) {
	for (int i = 0; names[i]; i++)
		if (!strcmp(v, names[i])) { *out = i; return 0; }
	dn_err("%s must be one of:", flag);
	for (int i = 0; names[i]; i++) fprintf(stderr, " %s", names[i]);
	fprintf(stderr, " (got '%s')\n", v);
	return 1;
}

// The header's voxel spacing, or isotropic where it cannot vouch for one:
// non-positive, non-finite, or axes more than 1000x apart, which no acquisition
// has and which would size the sphere's neighbour box from a nonsense ratio.
static void header_spacing(const dn_image *img, double sp[3]) {
	sp[0] = img->nim->dx;
	sp[1] = img->nim->dy;
	sp[2] = img->nim->dz;
	int sp_ok = 1;
	for (int a = 0; a < 3; a++)
		if (!(sp[a] > 0.0) || !isfinite(sp[a])) sp_ok = 0;
	if (sp_ok && fmax(fmax(sp[0], sp[1]), sp[2]) > 1000.0 * fmin(fmin(sp[0], sp[1]), sp[2])) {
		dn_err("note: voxel spacing %g x %g x %g is implausible; treating it as isotropic\n",
		       sp[0], sp[1], sp[2]);
		sp_ok = 0;
	}
	if (!sp_ok) sp[0] = sp[1] = sp[2] = 1.0;
}

// Largest noise-map patch size k for this image: a cube must fit the image; a
// sphere need not, but its voxels must, and half the image leaves room for the
// whole-shell overshoot.  dn_sigma_extent never goes below the main patch.
static int sigma_cap(const dn_image *img, int shape) {
	int cap = img->nx < img->ny ? img->nx : img->ny;
	if (img->nz < cap) cap = img->nz;
	if (shape == DN_SHAPE_SPHERE)
		for (cap = 3; (double)(cap + 2) * (cap + 2) * (cap + 2) <= 0.5 * (double)img->nvox3d; cap += 2) ;
	return cap;
}

// ---------------------------------------------------------------------------
// Tensor mode (-tensor gd)
// ---------------------------------------------------------------------------

#define DN_TENSOR_MAX_ECHO 32

typedef struct {
	const char *method;                          // NULL: matrix mode
	const char *fin, *fout, *fmask, *fnoise, *frank, *fgain;
	const char *echoes[DN_TENSOR_MAX_ECHO - 1];  // echoes 2.., in order
	int nechos;
	const char *phases[DN_TENSOR_MAX_ECHO];      // -tensor_phase, one per echo, in order
	int nphases;
	dn_phase_units punits;
	const char *window, *stride, *order;         // as typed, or NULL for the default
	int filter;  // -tensor_filter: DN_TFILTER_*
	int demean, center, known;
	int vst;     // -tensor_vst: 0, 1, or -1 for the default
	int pooled;  // -tensor_noise pooled: one noise map for every echo
	double sigma2;
	int nthreads, quiet;
	int given;   // tensor-only options seen, so matrix mode can refuse them
} tensor_cli;

// "a,b" or "a,b,c": positive integers, into v.  Returns the count, or 0 if malformed.
static int parse_list(const char *s, int *v, int maxn) {
	int n = 0;
	const char *p = s;
	for (;;) {
		char *end = NULL;
		const long x = strtol(p, &end, 10);
		if (end == p || x < 1 || x > INT32_MAX || n == maxn) return 0;
		v[n++] = (int)x;
		if (*end == '\0') return n;
		if (*end != ',') return 0;
		p = end + 1;
	}
}

// Echo k's output, "<stem>_echo-<k><ext>", where ext is whichever of .nii.gz,
// .nii, .hdr or .img the output name ends in (any case), or nothing -- nifti
// then appends what the input's type implies, as for any extensionless output.
static char *echo_name(const char *out, int k) {
	static const char *const exts[] = {".nii.gz", ".nii", ".hdr", ".img", NULL};
	const size_t len = strlen(out);
	size_t stem = len;
	for (int i = 0; exts[i]; i++) {
		const size_t el = strlen(exts[i]);
		if (len <= el) continue;
		int match = 1;
		for (size_t j = 0; j < el && match; j++) {
			char c = out[len - el + j];
			if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
			match = c == exts[i][j];
		}
		if (match) { stem = len - el; break; }
	}
	const size_t sz = len + 32;
	char *name = (char *)dn_malloc(sz, 1);
	if (name) snprintf(name, sz, "%.*s_echo-%d%s", (int)stem, out, k, out + stem);
	return name;
}

static int cmp_float_asc(const void *pa, const void *pb) {
	const float a = *(const float *)pa, b = *(const float *)pb;
	return (a > b) - (a < b);
}

// Put every echo in units of its noise, so that the one noise variance the
// tensor kernel assumes (it then reads sigma = 1) is true of each echo.
//
// Each echo's OWN noise by default, because the reconstructed images need not
// share one even though the thermal noise in k-space does.  On ds006131
// (Siemens, MB6 with in-plane GRAPPA, identical bandwidth on every echo) the
// noise-only scans put echo 5's noise at 0.51 of echo 1's, varying from 0.11 to
// 0.77 across the brain -- echo-specific reconstruction, not physics -- and a
// level pooled over the echoes read 29% low for echo 1 and 39% high for echo 5,
// under-denoising the first and over-cutting the last.  Each echo's own map
// was within 1-6% of those scans on every echo.  `pooled` keeps one level for
// all echoes, for reconstructions where they do share one.
//
// Each echo's map is matrix mode's Gavish & Donoho median map of that echo's
// own series (dn_sigma_map), unmasked as in matrix mode.  Then, per echo:
//   vt != NULL  magnitude: variance-stabilized with matrix mode's transform
//               and its DN_VST_ROUNDS re-estimation, so Rician noise becomes
//               unit Gaussian too and the noise-floor bias can be undone;
//   vt == NULL  signed data (rotated to real), or -tensor_vst n: divided by
//               the map, which is all Gaussian noise needs.
// Every map is estimated before any echo is transformed, so noiseless data
// (*none: every echo) are left untouched, and a noiseless echo among noisy
// ones is refused before anything changes.  sigma[e] and live[e] (vt only) are kept
// for tensor_restore; med[e] is echo e's median sigma for the summary.
// Returns 0, or non-zero on failure (reported).
static int tensor_normalize(dn_image *img, int ne, const dn_vst *vt, int pooled, int nthreads,
                            float **sigma, uint8_t **live, double *med, int *none) {
	const int nt = img[0].nvol;
	const size_t nv = img[0].nvox3d;
	double sp[3];
	header_spacing(&img[0], sp);
	dn_geom g1;
	const int k1 = dn_sigma_extent(nt, dn_auto_extent(nt), sigma_cap(&img[0], DN_SHAPE_SPHERE));
	if (dn_geom_init(&g1, DN_SHAPE_SPHERE, img[0].nx, img[0].ny, img[0].nz, nt, sp[0], sp[1], sp[2], k1, 0))
		return 1;
	*none = 0;
	int empty = 0, first_empty = 0;
	for (int e = 0; e < ne; e++) {
		if (!(sigma[e] = dn_sigma_map(&g1, img[e].data, nthreads))) return 1;
		const int filled = dn_sigma_fill(sigma[e], nv);
		if (filled < 0) return 1;
		if (filled && !empty++) first_empty = e;
	}
	// Every echo noiseless: nothing to normalize against, and nothing to remove.
	// SOME of them: an echo holding no noise at all (constant, or blank) beside
	// noisy ones is broken input rather than data, and its noise cannot be
	// borrowed from the others -- treating the whole run as noiseless instead
	// returned every echo unchanged at exit 0.
	if (empty == ne) {
		*none = 1;
		return 0;
	}
	if (empty) {
		dn_err("echo %d%s no noise level at all (it is constant or blank), while %d other\n",
		       first_empty + 1, empty > 1 ? " and others give" : " gives", ne - empty);
		dn_err("  echo%s do.  Check %s, or give one variance for every echo with\n",
		       ne - empty == 1 ? "" : "es", empty > 1 ? "those echoes" : "that echo");
		dn_err("  -tensor_sigma2.\n");
		return 1;
	}

	float *prev = (float *)dn_malloc(nv, sizeof(float));
	float *sorted = (float *)dn_malloc(nv, sizeof(float));
	if (!prev || !sorted) { free(prev); free(sorted); return 1; }
	int rc = 0;
	for (int e = 0; e < ne && !rc; e++) {
		float *x = img[e].data, *s = sigma[e];
		if (vt) {
			// As matrix mode: a voxel zero in every volume is background, and stays
			// zero so an all-zero patch still short-circuits.
			if (!(live[e] = (uint8_t *)dn_calloc(nv, 1))) { rc = 1; break; }
			for (int j = 0; j < nt; j++)
				for (size_t v = 0; v < nv; v++)
					if (x[v + (size_t)j * nv] != 0.0f) live[e][v] = 1;
			for (int it = 0; it < DN_VST_ROUNDS; it++) {
				vst_stabilise(vt, x, nv, nt, live[e], s, prev, it == 0);
				if (it == DN_VST_ROUNDS - 1) break;
				float *ratio = dn_sigma_map(&g1, x, nthreads);
				if (!ratio) { rc = 1; break; }
				for (size_t v = 0; v < nv; v++)
					if (ratio[v] > 0.0f) s[v] *= ratio[v];
				free(ratio);
			}
		} else {
			for (int j = 0; j < nt; j++)
				for (size_t v = 0; v < nv; v++) x[v + (size_t)j * nv] /= s[v];
		}
		for (size_t v = 0; v < nv; v++) sorted[v] = s[v];
		qsort(sorted, nv, sizeof(float), cmp_float_asc);
		med[e] = nv % 2 ? sorted[nv / 2] : 0.5 * ((double)sorted[nv / 2 - 1] + sorted[nv / 2]);
	}
	// -tensor_noise pooled: the echoes share ONE noise level, the mean of their
	// (corrected) variances, and each is re-expressed against it: back to raw
	// through f's exact algebraic inverse and forward again (vst_stabilise with
	// prev = the echo's own map), or rescaled.  For data whose echoes really do
	// share a noise level, which is the physical expectation when every echo is
	// reconstructed the same way, this trades the per-echo maps' estimation noise
	// for one map estimated from all of them.
	if (!rc && pooled) {
		float *pool = (float *)dn_calloc(nv, sizeof(float));
		if (!pool) { free(prev); free(sorted); return 1; }
		for (int e = 0; e < ne; e++)
			for (size_t v = 0; v < nv; v++) pool[v] += sigma[e][v] * sigma[e][v] / (float)ne;
		for (size_t v = 0; v < nv; v++) pool[v] = sqrtf(pool[v]);
		for (int e = 0; e < ne; e++) {
			float *x = img[e].data, *s = sigma[e];
			if (vt) {
				for (size_t v = 0; v < nv; v++) prev[v] = s[v];
				vst_stabilise(vt, x, nv, nt, live[e], pool, prev, 0);
			} else {
				for (int j = 0; j < nt; j++)
					for (size_t v = 0; v < nv; v++) x[v + (size_t)j * nv] *= s[v] / pool[v];
			}
			for (size_t v = 0; v < nv; v++) s[v] = pool[v];
		}
		free(pool);
	}
	free(prev);
	free(sorted);
	return rc;
}

// Map one denoised echo back to its own units.  A covered voxel is an estimate:
// sigma times the exact-unbiased inverse (magnitude, removing the noise floor)
// or sigma times the value.  A voxel no window covered (NaN in `cover`) still
// holds its normalized INPUT, so it goes back through the algebraic inverse,
// which is f's exact inverse, and keeps its input value to float rounding.
static void tensor_restore(float *out, size_t nv, int nt, const float *sigma, const uint8_t *live,
                           const dn_vst *vt, const float *cover) {
	for (size_t v = 0; v < nv; v++) {
		const int covered = !isnan(cover[v]);
		const double s = sigma[v];
		for (int j = 0; j < nt; j++) {
			float *o = out + v + (size_t)j * nv;
			if (!vt) *o = (float)(s * *o);
			else if (covered) *o = (float)(s * dn_vst_inverse(vt, *o, 1));
			else if (live[v]) *o = (float)(s * dn_vst_inverse(vt, *o, 0));
		}
	}
}

static int run_tensor(const tensor_cli *t) {
	static const char *const mode_names[] = {"space", "time", "echo"};
	const int ne = t->nechos + 1;
	int status = EXIT_FAILURE;
	dn_image img[DN_TENSOR_MAX_ECHO];
	float *out[DN_TENSOR_MAX_ECHO];
	char *oname[DN_TENSOR_MAX_ECHO];
	dn_path *paths = NULL;
	uint8_t *mask = NULL;
	float *sig = NULL, *rank = NULL, *gain = NULL;
	float *sigmas[DN_TENSOR_MAX_ECHO] = {NULL};
	uint8_t *lives[DN_TENSOR_MAX_ECHO] = {NULL};
	dn_vst *vt = NULL;
	memset(img, 0, sizeof img);
	memset(out, 0, sizeof out);
	memset(oname, 0, sizeof oname);

	if (ne < 2) {
		dn_err("-tensor gd needs at least two echoes: <input> and one or more -echo.\n");
		return EXIT_FAILURE;
	}
	const char *fin[DN_TENSOR_MAX_ECHO];
	fin[0] = t->fin;
	for (int e = 1; e < ne; e++) fin[e] = t->echoes[e - 1];
	for (int e = 0; e < ne; e++)
		if (!strcmp(fin[e], "-")) {
			dn_err("tensor mode reads every echo from a named file, not stdin.\n");
			return EXIT_FAILURE;
		}

	dn_trun r;
	memset(&r, 0, sizeof r);
	r.k.filter = t->filter;
	r.k.subtract_mean = t->demean;
	r.sigma2 = t->sigma2;
	r.center_assign = t->center;
	r.ne = ne;

	// A two-element window slides over x and y and spans every slice, as
	// ME-PCA's run_TensorMPPCA does; three elements make it a local 3D window.
	// The default is local.
	int win[3] = {5, 5, 5}, nwin = 3;
	if (t->window && !(nwin = parse_list(t->window, win, 3))) {
		dn_err("-tensor_window needs 2 or 3 positive integers, comma separated (got '%s')\n", t->window);
		return EXIT_FAILURE;
	}
	if (nwin < 2) {
		dn_err("-tensor_window needs 2 or 3 sizes: x,y (spanning every slice) or x,y,z.\n");
		return EXIT_FAILURE;
	}
	int str[3] = {1, 1, 1}, nstr = 1;
	if (t->stride && !(nstr = parse_list(t->stride, str, 3))) {
		dn_err("-tensor_stride needs positive integers, comma separated (got '%s')\n", t->stride);
		return EXIT_FAILURE;
	}
	if (nstr != 1 && nstr != nwin) {
		dn_err("-tensor_stride needs one value, or one per -tensor_window size.\n");
		return EXIT_FAILURE;
	}
	if (nstr == 1) str[1] = str[2] = str[0];
	r.order[0] = DN_TMODE_SPACE;
	r.order[1] = DN_TMODE_TIME;
	r.order[2] = DN_TMODE_ECHO;
	if (t->order) {
		char buf[64];
		int seen = 0, n = 0;
		if (snprintf(buf, sizeof buf, "%s", t->order) >= (int)sizeof buf) n = -1;
		for (char *tok = n < 0 ? NULL : strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
			int mode = -1;
			for (int k = 0; k < 3; k++) if (!strcmp(tok, mode_names[k])) mode = k;
			if (mode < 0 || (seen & (1 << mode)) || n == 3) { n = -1; break; }
			seen |= 1 << mode;
			r.order[n++] = mode;
		}
		if (n != 3) {
			dn_err("-tensor_order needs space, time and echo, each once, comma separated (got '%s')\n",
			       t->order);
			return EXIT_FAILURE;
		}
	}

	for (int e = 0; e < ne; e++) {
		if (dn_image_read(fin[e], "echo", 1, &img[e])) goto done;
		if (e == 0) continue;
		// Every echo must be the same acquisition: dimensions, grid, and time axis.
		// A same-sized echo on a shifted grid, or with another TR, would otherwise
		// be stacked into one tensor and denoised as though it were.
		if (img[e].nx != img[0].nx || img[e].ny != img[0].ny || img[e].nz != img[0].nz ||
		    img[e].nvol != img[0].nvol) {
			dn_err("echo '%s' is %dx%dx%d x %d but '%s' is %dx%dx%d x %d\n", fin[e],
			       img[e].nx, img[e].ny, img[e].nz, img[e].nvol, fin[0],
			       img[0].nx, img[0].ny, img[0].nz, img[0].nvol);
			goto done;
		}
		const float off = dn_grid_offset_mm(&img[0], &img[e]);
		if (!(off <= DN_GRID_TOL_MM)) {
			dn_err("echo '%s' is not on the grid of '%s' (%g mm apart)\n", fin[e], fin[0], off);
			goto done;
		}
		const double dt0 = img[0].nim->dt, dte = img[e].nim->dt;
		if (img[e].nim->time_units != img[0].nim->time_units ||
		    !(fabs(dte - dt0) <= 1e-4 * fmax(fabs(dt0), fabs(dte)))) {
			dn_err("echo '%s' has a different repetition time from '%s' (%g against %g, in\n",
			       fin[e], fin[0], dte, dt0);
			dn_err("  its own time units)\n");
			goto done;
		}
	}
	// Magnitude in: ME-PCA, whose window rules these are, takes abs() of its
	// input, which silently destroys a signed series, and -tensor_phase rotates
	// a magnitude.  A series that is already signed is refused either way.
	for (int e = 0; e < ne; e++)
		for (size_t i = 0; i < img[e].nvox3d * (size_t)img[e].nvol; i++)
			if (img[e].data[i] < 0.0f) {
				dn_err("echo '%s' has negative values; tensor mode expects MAGNITUDE images\n", fin[e]);
				dn_err("  (with -tensor_phase for complex data).\n");
				goto done;
			}
	if (t->nphases && t->nphases != ne) {
		dn_err("%d -tensor_phase images for %d echoes; give one per echo, in echo order.\n",
		       t->nphases, ne);
		goto done;
	}
	// -tensor_vst: on for magnitude by default, as -vst is in matrix mode.  It
	// stabilizes against each echo's estimated map, so a known variance leaves
	// it nothing to stabilize against, and rotated data are Gaussian already.
	const int vst = t->vst < 0 ? (!t->nphases && !t->known) : t->vst;
	if (t->pooled && t->known) {
		dn_err("-tensor_noise chooses how the noise maps are used; -tensor_sigma2 replaces them.\n");
		goto done;
	}
	if (vst && (t->nphases || t->known)) {
		dn_err("-tensor_vst stabilizes magnitude against each echo's estimated noise map;\n");
		dn_err("  %s.\n", t->nphases ? "-tensor_phase data are rotated to real and already Gaussian"
		                             : "-tensor_sigma2 gives no map to stabilize against");
		goto done;
	}

	r.nx = img[0].nx;
	r.ny = img[0].ny;
	r.nz = img[0].nz;
	r.nt = img[0].nvol;
	r.win[0] = win[0];
	r.win[1] = win[1];
	r.stride[0] = str[0];
	r.stride[1] = str[1];
	if (nwin == 2) {
		r.full_z = 1;
		r.win[2] = r.nz;
		r.stride[2] = 1;
	} else {
		r.win[2] = win[2];
		r.stride[2] = str[2];
	}
	if (dn_tensor_setup(&r)) goto done;

	for (int e = 0; e < ne; e++)
		if (!(oname[e] = echo_name(t->fout, e + 1))) goto done;
	const int npath = 3 * ne + 4;
	paths = (dn_path *)dn_calloc((size_t)npath, sizeof(dn_path));
	if (!paths) goto done;
	for (int e = 0; e < ne; e++) {
		paths[e] = (dn_path){"echo", fin[e], 1, NULL, NULL};
		paths[ne + e] = (dn_path){"echo output", oname[e], 0, NULL, NULL};
		paths[2 * ne + e] = (dn_path){"-tensor_phase", t->nphases ? t->phases[e] : NULL, 1, NULL, NULL};
	}
	paths[3 * ne] = (dn_path){"-mask", t->fmask, 1, NULL, NULL};
	paths[3 * ne + 1] = (dn_path){"-noise", t->fnoise, 0, NULL, NULL};
	paths[3 * ne + 2] = (dn_path){"-rank", t->frank, 0, NULL, NULL};
	paths[3 * ne + 3] = (dn_path){"-snr_gain", t->fgain, 0, NULL, NULL};
	// Every output is written in echo 1's file type, because that is the one
	// type the check below resolves extensionless names with.  Each echo is
	// written from its OWN header, and left to its own type an Analyze-pair echo
	// turned "out_echo-2" into out_echo-2.hdr/.img while the check modeled
	// out_echo-2.nii -- so an output could land on that echo's input at exit 0.
	for (int e = 1; e < ne; e++) img[e].nifti_type = img[0].nifti_type;
	if (check_path_collisions(paths, npath, img[0].nifti_type)) goto done;

	size_t n_in_mask = img[0].nvox3d;
	if (t->fmask && !(mask = dn_mask_build(&img[0], t->fmask, &n_in_mask))) goto done;

	// Complex data: each echo onto the real axis with its OWN phase, which
	// evolves with TE, exactly as matrix mode's -phase does for one series.  The
	// noise is then zero-mean Gaussian rather than Rician, and the denoised
	// echoes are signed.  After every refusal, so nothing is computed for a run
	// that would not happen.
	for (int e = 0; e < t->nphases; e++)
		if (dn_rotate_to_real(&img[e], t->phases[e], t->punits, t->quiet)) goto done;

	const int nthreads_req = t->nthreads < 1 ? dn_default_threads() : t->nthreads;
	double med[DN_TENSOR_MAX_ECHO];
	int noiseless = 0, normalized = 0;
	if (!t->known) {
		// The median estimator needs a noise bulk: with fewer than 3 columns its
		// "median" averages signal and noise (see the matrix mode's fallback).
		if (r.nt < 3) {
			dn_err("estimating the noise level needs at least 3 volumes per echo (got %d);\n", r.nt);
			dn_err("  give it with -tensor_sigma2 instead.\n");
			goto done;
		}
		if (vst && !(vt = dn_vst_create(1))) goto done;
		if (tensor_normalize(img, ne, vt, t->pooled, nthreads_req, sigmas, lives, med, &noiseless)) goto done;
		if (noiseless) {
			dn_err("note: no noise level could be estimated, so every component is kept\n");
			r.sigma2 = 0.0;
		} else {
			normalized = 1;
			r.sigma2 = 1.0;   // every echo is now in units of its own noise
		}
	}

	const int nthreads = dn_tensor_threads(&r, nthreads_req);
	r.nthreads = nthreads;
	if (!t->quiet) {
		dn_err("input        : %d echoes, %dx%dx%d, %d volumes\n", ne, r.nx, r.ny, r.nz, r.nt);
		for (int e = 0; e < ne; e++) dn_err("  echo %-2d     : %s\n", e + 1, fin[e]);
		if (t->nphases) dn_err("complex      : each echo rotated onto the real axis; output is signed\n");
		dn_err("method       : tensor, recursive mode-wise cut (Olesen et al. 2022 scheme)\n");
		if (r.full_z)
			dn_err("window       : %dx%d in-plane, spanning all %d slices = %d voxels\n",
			       r.win[0], r.win[1], r.nz, r.win[0] * r.win[1] * r.win[2]);
		else
			dn_err("window       : %dx%dx%d = %d voxels\n",
			       r.win[0], r.win[1], r.win[2], r.win[0] * r.win[1] * r.win[2]);
		dn_err("tensor       : %s %d x %s %d x %s %d\n", mode_names[r.order[0]], r.k.dims[0],
		       mode_names[r.order[1]], r.k.dims[1], mode_names[r.order[2]], r.k.dims[2]);
		if (t->known) dn_err("noise        : known, variance %g\n", t->sigma2);
		else if (!noiseless) {
			dn_err("noise        : Gavish-Donoho median map %s; each echo %s\n",
			       t->pooled ? "pooled over the echoes" : "of each echo",
			       vst ? "variance-stabilized against it (Rician, exact-unbiased inverse)" : "scaled by it");
			dn_err("  median sigma by echo%s:", t->pooled ? " before pooling" : "");
			for (int e = 0; e < ne; e++) fprintf(stderr, " %.4g", med[e]);
			fprintf(stderr, "\n");
		}
		dn_err("filter       : %s\n",
		       t->filter == DN_TFILTER_OPTTHRESH ? "optthresh: known-noise optimal hard threshold, every mode"
		       : t->filter == DN_TFILTER_TRUNCATE ? "truncate: bulk edge, nothing shrunk"
		       : "optshrink: bulk edge, Frobenius-optimal shrinkage on the last mode");
		if (t->demean) dn_err("demean       : every unfolding, along its longer side\n");
		dn_err("aggregation  : %s, stride %d,%d,%d\n",
		       t->center ? "window center only" : "uniform mean of every window", r.stride[0], r.stride[1],
		       r.stride[2]);
		dn_err("threads      : %d\n", nthreads);
		if (mask) dn_err("mask         : %s (%zu of %zu voxels)\n", t->fmask, n_in_mask, img[0].nvox3d);
	}

	const size_t nv = img[0].nvox3d, nel = nv * (size_t)r.nt;
	const float *echo_in[DN_TENSOR_MAX_ECHO];
	for (int e = 0; e < ne; e++) {
		if (!(out[e] = (float *)dn_calloc(nel, sizeof(float)))) goto done;
		echo_in[e] = img[e].data;
	}
	// Always kept: it is NaN exactly where no window reached, which the
	// restore needs to know, as well as being the -noise map.
	if (!(sig = (float *)dn_malloc(nv, sizeof(float)))) goto done;
	if (t->frank && !(rank = (float *)dn_malloc(nv * 3, sizeof(float)))) goto done;
	if (t->fgain && !(gain = (float *)dn_malloc(nv, sizeof(float)))) goto done;
	r.echo = echo_in;
	r.out = out;
	r.mask = mask;
	r.sigma2_out = sig;
	r.rank = rank;
	r.gain = gain;
	if (dn_tensor_execute(&r)) goto done;
	if (normalized)
		for (int e = 0; e < ne; e++) tensor_restore(out[e], nv, r.nt, sigmas[e], lives[e], vt, sig);

	// Maps first, as in matrix mode.  -noise is the standard deviation, the unit
	// matrix mode writes too, one volume per echo: the square root of the mean
	// window variance, in that echo's own units.
	if (t->fnoise) {
		float *noise = (float *)dn_malloc(nv * (size_t)ne, sizeof(float));
		if (!noise) goto done;
		for (int e = 0; e < ne; e++)
			for (size_t v = 0; v < nv; v++)
				noise[v + (size_t)e * nv] = sqrtf(sig[v]) * (normalized ? sigmas[e][v] : 1.0f);
		const int rc = dn_write_f32(&img[0], t->fnoise, noise, ne);
		free(noise);
		if (rc) goto done;
	}
	if (rank) {
		// One volume per tensor mode, in mode order, and the header says which.
		char saved[80];
		memcpy(saved, img[0].nim->descrip, sizeof saved);
		snprintf(img[0].nim->descrip, sizeof saved, "tensor ranks: %s,%s,%s",
		         mode_names[r.order[0]], mode_names[r.order[1]], mode_names[r.order[2]]);
		const int rc = dn_write_f32(&img[0], t->frank, rank, 3);
		memcpy(img[0].nim->descrip, saved, sizeof saved);
		if (rc) goto done;
	}
	if (gain && dn_write_f32(&img[0], t->fgain, gain, 1)) goto done;
	for (int e = 0; e < ne; e++)
		if (dn_write_f32(&img[e], oname[e], out[e], r.nt)) goto done;
	status = EXIT_SUCCESS;

done:
	for (int e = 0; e < ne; e++) {
		dn_image_free(&img[e]);
		free(out[e]);
		free(oname[e]);
	}
	for (int e = 0; e < DN_TENSOR_MAX_ECHO; e++) {
		free(sigmas[e]);
		free(lives[e]);
	}
	dn_vst_free(vt);
	free(paths);
	free(mask);
	free(sig);
	free(rank);
	free(gain);
	return status;
}

int main(int argc, char *argv[]) {
	const char *fin = NULL, *fout = NULL;
	const char *fmask = NULL, *fnoise = NULL, *frank = NULL;
	const char *fphase = NULL, *freal = NULL, *fbval = NULL;
	// Off by default: measured neutral under the hard threshold (RMSE -0.4% at
	// sigma 10, +1.3% at sigma 30, synthetic truth on benchmark/).
	int demean = 0;
	// Defaults chosen by measurement against a synthetic truth (AGENTS.md dead ends):
	// each beat the original exclusive / optthresh / cube pipeline, which
	// `-filter optthresh -aggregator exclusive -shape cube -vst n` still reproduces
	// byte for byte.
	int filter = DN_FILTER_OPTSHRINK;
	int shape = DN_SHAPE_SPHERE;
	int aggregator = DN_AGG_GAUSSIAN, stride = 0;   // 0: auto
	double fwhm = 2.0;
	int fwhm_given = 0;
	int vst = -1, ncoil = 1, ncoil_given = 0, keep_bias = 0;   // -1: auto
	// Options that only the denoiser reads, counted so -degibbs o can refuse them.
	int denoise_opts = 0;
	static const char *const yes_no[] = {"n", "y", NULL};
	static const char *const filters[] = {"optthresh", "optshrink", "truncate", NULL};   // DN_FILTER_*
	static const char *const shapes[] = {"cube", "sphere", NULL};                        // DN_SHAPE_*
	static const char *const aggregators[] = {"exclusive", "gaussian", "uniform", NULL}; // DN_AGG_*
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
#ifdef DN_DEGIBBS
	// pf_given, not `pf != DG_PF_FULL`, is what separates "not supplied" from
	// "supplied as 1.0" -- the same distinction -extent and -phaseunits draw, and
	// without it an explicit `-pF 1.0` was the one factor accepted with no
	// -degibbs to act on it.
	double pf = DG_PF_FULL;   // symmetric k-space
	int pf_given = 0;
#endif
	// Tensor mode's own options, kept apart from the matrix denoiser's so that
	// neither mode can silently read the other's.
	tensor_cli tc;
	memset(&tc, 0, sizeof tc);
	tc.vst = -1;

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
		if (!strcmp(a, "-preserve_noise_bias")) { keep_bias = 1; denoise_opts++; continue; }

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
			else if (!strcmp(a, "-bval")) { fbval = v; denoise_opts++; }
			else if (!strcmp(a, "-vst") || !strcmp(a, "-demean") || !strcmp(a, "-filter") ||
			         !strcmp(a, "-shape") || !strcmp(a, "-aggregator")) {
				denoise_opts++;
				int bad;
				if (!strcmp(a, "-vst")) bad = parse_choice(a, v, yes_no, &vst);
				else if (!strcmp(a, "-demean")) bad = parse_choice(a, v, yes_no, &demean);
				else if (!strcmp(a, "-filter")) bad = parse_choice(a, v, filters, &filter);
				else if (!strcmp(a, "-shape")) bad = parse_choice(a, v, shapes, &shape);
				else bad = parse_choice(a, v, aggregators, &aggregator);
				if (bad) return EXIT_FAILURE;
			}
			else if (!strcmp(a, "-noise_dof")) {
				if (parse_int(v, &ncoil) || ncoil < 1 || ncoil > 64) {
					dn_err("-noise_dof needs a channel count from 1 to 64 (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				ncoil_given = 1;
				denoise_opts++;
			}
			else if (!strcmp(a, "-aggregator_fwhm")) {
				char *pend = NULL;
				fwhm = strtod(v, &pend);
				if (*pend != '\0' || pend == v || !(fwhm > 0.0) || !isfinite(fwhm)) {
					dn_err("-aggregator_fwhm needs a positive number (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				fwhm_given = 1;
				denoise_opts++;
			}
			else if (!strcmp(a, "-stride")) {
				if (parse_int(v, &stride) || stride < 1 || stride > 2) {
					dn_err("-stride must be 1 or 2 (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				denoise_opts++;
			}
			else if (!strcmp(a, "-tensor")) {
				if (strcmp(v, "gd")) {
					dn_err("-tensor must be gd, the only tensor method (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				tc.method = v;
			}
			else if (!strcmp(a, "-echo")) {
				if (tc.nechos >= DN_TENSOR_MAX_ECHO - 1) {
					dn_err("at most %d echoes\n", DN_TENSOR_MAX_ECHO);
					return EXIT_FAILURE;
				}
				tc.echoes[tc.nechos++] = v;
				tc.given++;
			}
			else if (!strcmp(a, "-snr_gain")) { tc.fgain = v; tc.given++; }
			else if (!strcmp(a, "-tensor_phase")) {
				if (tc.nphases >= DN_TENSOR_MAX_ECHO) {
					dn_err("at most %d -tensor_phase images\n", DN_TENSOR_MAX_ECHO);
					return EXIT_FAILURE;
				}
				tc.phases[tc.nphases++] = v;
				tc.given++;
			}
			else if (!strcmp(a, "-tensor_window")) { tc.window = v; tc.given++; }
			else if (!strcmp(a, "-tensor_stride")) { tc.stride = v; tc.given++; }
			else if (!strcmp(a, "-tensor_order")) { tc.order = v; tc.given++; }
			else if (!strcmp(a, "-tensor_noise")) {
				static const char *const noise_modes[] = {"echo", "pooled", NULL};
				if (parse_choice(a, v, noise_modes, &tc.pooled)) return EXIT_FAILURE;
				tc.given++;
			}
			else if (!strcmp(a, "-tensor_filter")) {
				static const char *const tfilters[] = {"optshrink", "truncate", "optthresh", NULL};  // DN_TFILTER_*
				if (parse_choice(a, v, tfilters, &tc.filter)) return EXIT_FAILURE;
				tc.given++;
			}
			else if (!strcmp(a, "-tensor_demean") || !strcmp(a, "-tensor_center") ||
			         !strcmp(a, "-tensor_vst")) {
				int *dst = !strcmp(a, "-tensor_demean") ? &tc.demean
				         : !strcmp(a, "-tensor_vst") ? &tc.vst : &tc.center;
				if (parse_choice(a, v, yes_no, dst)) return EXIT_FAILURE;
				tc.given++;
			}
			else if (!strcmp(a, "-tensor_sigma2")) {
				char *pend = NULL;
				tc.sigma2 = strtod(v, &pend);
				if (*pend != '\0' || pend == v || !(tc.sigma2 >= 0.0) || !isfinite(tc.sigma2)) {
					dn_err("-tensor_sigma2 needs a finite, non-negative noise VARIANCE (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				tc.known = 1;
				tc.given++;
			}
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
#ifdef DN_DEGIBBS
			else if (!strcmp(a, "-degibbs")) {
				if (!strcmp(v, "y")) degibbs = DN_DG_YES;
				else if (!strcmp(v, "n")) degibbs = DN_DG_NO;
				else if (!strcmp(v, "o")) degibbs = DN_DG_ONLY;
				else {
					dn_err("-degibbs must be y, n or o (got '%s')\n", v);
					return EXIT_FAILURE;
				}
			}
			else if (!strcmp(a, "-pF")) {
				char *pend = NULL;
				pf = strtod(v, &pend);
				// Only the PARSE is checked here.  Which factors exist is
				// dn_degibbs_check's to say, and it says it precisely; a range
				// test here told a user that 0.9 was acceptable, then a second
				// message refused it.
				if (*pend != '\0' || pend == v) {
					dn_err("-pF needs a partial-Fourier factor (got '%s')\n", v);
					return EXIT_FAILURE;
				}
				pf_given = 1;
			}
#else
			// One refusal for both, since the message is the same and `a` names
			// whichever was given.
			else if (!strcmp(a, "-degibbs") || !strcmp(a, "-pF")) {
				dn_err("%s is not compiled into this build.\n", a);
				dn_err("  Rebuild with 'make DEGIBBS=1'.\n");
				return EXIT_FAILURE;
			}
#endif
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
	const char *sv_out[5] = {fout, fnoise, frank, freal, tc.fgain};
	const char *sv_lbl[5] = {"output", "-noise", "-rank", "-real", "-snr_gain"};
	for (int i = 0; i < 5; i++) {
		if (sv_out[i] && sv_out[i][0] == '-' && sv_out[i][1] == '\0') {
			dn_err("\"-\" is an input (stdin), not an output; %s needs a filename.\n", sv_lbl[i]);
			return EXIT_FAILURE;
		}
	}
	// Tensor mode is a different estimator with its own options, so each mode
	// refuses the other's rather than quietly ignoring them.
	if (!tc.method && tc.given) {
		dn_err("-echo, -snr_gain and the -tensor_* options describe tensor denoising,\n");
		dn_err("  which needs -tensor gd.\n");
		return EXIT_FAILURE;
	}
	if (tc.method) {
		int matrix_opts = denoise_opts + (fphase != NULL) + (freal != NULL) + extent_given +
		                  (degibbs != DN_DG_NO);
#ifdef DN_DEGIBBS
		matrix_opts += pf_given;
#endif
		if (matrix_opts) {
			dn_err("-tensor gd takes none of the matrix denoiser's options (-phase, -real,\n");
			dn_err("  -bval, -demean, -filter, -shape, -aggregator(_fwhm), -stride, -vst,\n");
			dn_err("  -noise_dof, -preserve_noise_bias, -extent, -degibbs, -pF); its own are\n");
			dn_err("  the -tensor_* options, and -tensor_phase takes each echo's phase.\n");
			return EXIT_FAILURE;
		}
		// -phaseunits is shared: it describes whichever phase images were given.
		if (punits_given && !tc.nphases) {
			dn_err("-phaseunits describes the -tensor_phase images, which were not given.\n");
			return EXIT_FAILURE;
		}
		tc.punits = punits;
		tc.fin = fin;
		tc.fout = fout;
		tc.fmask = fmask;
		tc.fnoise = fnoise;
		tc.frank = frank;
		tc.nthreads = nthreads;
		tc.quiet = quiet;
		return run_tensor(&tc);
	}
	// -degibbs o does no denoising, so everything describing the denoiser -- its
	// mask, its by-products, its patch size, and the phase rotation that exists
	// only to change the noise the denoiser sees -- has nothing to act on.
	// Ignoring them silently would give a run that looked fine and wrote fewer
	// files than were asked for.
	if (degibbs == DN_DG_ONLY && (fmask || fnoise || frank || fphase || freal || extent_given ||
	                              denoise_opts)) {
		dn_err("-degibbs o does no denoising, so none of the denoiser's options have\n");
		dn_err("  anything to act on.\n");
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
#ifdef DN_DEGIBBS
	// -pF changes how ringing is removed, so it means nothing without a stage
	// that removes it.  Silently ignoring it would leave PF ringing in place
	// while the command line says otherwise.
	if (pf_given && degibbs == DN_DG_NO) {
		dn_err("-pF describes how -degibbs should work, but -degibbs is not enabled.\n");
		return EXIT_FAILURE;
	}
#endif
	// Auto: stride 2 under averaging, 1 when only centres are kept.  An automatic
	// VST is decided once the input is read (see below).
	if (!stride) stride = aggregator == DN_AGG_EXCLUSIVE ? 1 : 2;
	if (vst == 1 && fphase) {
		dn_err("-vst stabilises magnitude noise; -phase data are rotated to real and are\n");
		dn_err("  already Gaussian.\n");
		return EXIT_FAILURE;
	}
	if (fwhm_given && aggregator != DN_AGG_GAUSSIAN) {
		dn_err("-aggregator_fwhm sets the Gaussian width, which needs -aggregator gaussian.\n");
		return EXIT_FAILURE;
	}
	if (stride > 1 && aggregator == DN_AGG_EXCLUSIVE) {
		dn_err("-stride %d leaves voxels with no patch centred on them; it needs\n", stride);
		dn_err("  -aggregator gaussian or uniform.\n");
		return EXIT_FAILURE;
	}
	if (fbval && !demean) {
		dn_err("-bval is only used by -demean y.\n");
		return EXIT_FAILURE;
	}
	if (punits_given && !fphase) {
		dn_err("-phaseunits describes the -phase image, which was not given.\n");
		return EXIT_FAILURE;
	}
	int status = EXIT_FAILURE;
	dn_image img;
	uint8_t *mask = NULL, *live = NULL;
	float *out = NULL, *noise = NULL, *means = NULL, *sigma_map = NULL, *unit_map = NULL;
	dn_vst *vt = NULL;
	int need_map = 0;
	int *group = NULL;
	double *shell_b = NULL;
	int *shell_n = NULL;
	char *bval_auto = NULL;
	uint16_t *rank = NULL;

	// Not require4d in "only" mode: mrdegibbs accepts a 3D image, and the
	// >= 2 volume rule belongs to the denoiser, which is not running.
	if (dn_image_read(fin, "input", degibbs != DN_DG_ONLY, &img)) return EXIT_FAILURE;

	// After the read, not before: resolving an output prefix needs the input's
	// nifti_type.  Still ahead of every write, and ahead of the denoising itself.
	// The b-values -demean y will read, explicit or the input's sidecar, so the
	// collision check below covers them: an output named over them would replace
	// the text it read at exit 0.
	if (demean && !fbval) {
		if (!strcmp(fin, "-")) {
			dn_err("-demean y on stdin needs -bval: there is no input name to find it by.\n");
			goto done;
		}
		bval_auto = dn_bval_sidecar(fin);
		if (!bval_auto) goto done;
		fbval = bval_auto;
	}
	dn_path paths[] = {
		{"input",  fin,    1, NULL, NULL},
		{"-bval",  fbval,  2, NULL, NULL},
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
	if (degibbs != DN_DG_NO && dn_degibbs_check(img.nx, img.ny, pf)) goto done;
#endif

#ifdef DN_DEGIBBS
	// "Only" short-circuits everything the denoiser needs -- geometry, mask,
	// scratch budget -- and rewrites the input buffer in place, so there is no
	// second copy of the series.
	if (degibbs == DN_DG_ONLY) {
		if (nthreads < 1) nthreads = dn_default_threads();
		// Report what will really run, not what was asked for -- same contract the
		// denoising path has, and it was missing here.
		nthreads = dn_degibbs_threads(img.nx, img.ny, img.nz, img.nvol, nthreads, pf);
		if (!quiet) {
			dn_err("input        : %s (%dx%dx%d, %d volume%s)\n", fin,
			        img.nx, img.ny, img.nz, img.nvol, img.nvol == 1 ? "" : "s");
			dn_err("degibbs      : only, x-y planes (no denoising)%s\n", pf_note(pf));
			dn_err("threads      : %d\n", nthreads);
		}
		if (dn_degibbs(img.data, img.nx, img.ny, img.nz, img.nvol, nthreads, pf)) goto done;
		if (dn_write_f32(&img, fout, img.data, img.nvol)) goto done;
		status = EXIT_SUCCESS;
		goto done;
	}
#endif

	// A magnitude image is never negative, so a negative voxel means signed data
	// -- typically a series already rotated to real elsewhere -- whose noise is
	// Gaussian, and stabilising it as Rician would bias it.  Only the automatic
	// choice yields to this; an explicit -vst y is obeyed.
	// -phase also says the input is not magnitude.
	int vst_off_signed = 0;
	if (vst < 0) {
		vst = !fphase;
		for (size_t i = 0; vst && i < img.nvox3d * (size_t)img.nvol; i++)
			if (img.data[i] < 0.0f) { vst = 0; vst_off_signed = 1; }
	}
	if ((ncoil_given || keep_bias) && !vst) {
		dn_err("-noise_dof and -preserve_noise_bias describe -vst y, which is %s.\n",
		       vst_off_signed ? "off because the input has negative values" : "not enabled");
		goto done;
	}

	const int extent_was_auto = !extent_given;
	if (!extent_given) extent = dn_auto_extent(img.nvol);

	int ngroups = 0;
	if (demean) {
		group = (int *)dn_malloc((size_t)img.nvol, sizeof(int));
		shell_b = (double *)dn_malloc((size_t)img.nvol, sizeof(double));
		shell_n = (int *)dn_malloc((size_t)img.nvol, sizeof(int));
		if (!group || !shell_b || !shell_n) goto done;
		ngroups = dn_bval_groups(fbval, img.nvol, group, shell_b, shell_n);
		if (ngroups < 0) goto done;
	}
	// The noise level is a MEDIAN singular value, which needs a noise bulk: with
	// two columns left after demeaning, the "median" averages signal and noise and
	// read 5.8 against a true 1, and the VST then zeroed 80% of the output.  A
	// noise map from such data is meaningless, so -noise is refused; otherwise the
	// hard threshold, which uses no sigma, takes over rather than refuse a run the
	// minimum input (2 volumes) is allowed.
	need_map = filter != DN_FILTER_OPTTHRESH || vst;
	if (img.nvol - ngroups < 3) {
		if (fnoise) {
			dn_err("-noise needs at least 3 noise columns (volumes minus demeaned shells);\n");
			dn_err("  this run has %d, from which no noise level can be estimated.\n", img.nvol - ngroups);
			goto done;
		}
		if (need_map) {
			dn_err("note: %d noise column%s cannot give a noise level, so this run uses\n",
			       img.nvol - ngroups, img.nvol - ngroups == 1 ? "" : "s");
			dn_err("  -filter optthresh without -vst.\n");
			filter = DN_FILTER_OPTTHRESH;
			vst = 0;
			need_map = 0;
		}
	}

	dn_geom g, g1;
	double sp[3];
	header_spacing(&img, sp);
	if (dn_geom_init(&g, shape, img.nx, img.ny, img.nz, img.nvol, sp[0], sp[1], sp[2], extent, ngroups))
		goto done;
	if (need_map) {
		const int k1 = dn_sigma_extent(img.nvol, extent, sigma_cap(&img, shape));
		if (dn_geom_init(&g1, shape, img.nx, img.ny, img.nz, img.nvol, sp[0], sp[1], sp[2], k1, ngroups))
			goto done;
	}
	g.filter = filter;

	size_t n_in_mask = img.nvox3d;
	if (fmask) {
		mask = dn_mask_build(&img, fmask, &n_in_mask);
		if (!mask) goto done;
	}

	if (nthreads < 1) nthreads = dn_default_threads();
	// The REQUEST, kept separate from the denoiser's effective count below.
	// -degibbs has its own cap and none of the denoiser's constraints -- no mask,
	// no patch arena -- so handing it the throttled number would let a small mask
	// or the eigensolver's scratch budget starve a stage limited by neither.  The
	// sigma pass has its own geometry and visits no mask, so the same applies.
	const int nthreads_req = nthreads;
	// Report what will actually run, not what was asked for.  Averaging visits
	// patch centres in tiles, not masked voxels, so it sizes its own team.
	nthreads = aggregator == DN_AGG_EXCLUSIVE ? dn_effective_threads(&g, n_in_mask, nthreads)
	                                          : dn_agg_threads(&g, stride, nthreads);

	if (!quiet) {
		dn_err("input        : %s (%dx%dx%d, %d volumes)\n",
		        fin, img.nx, img.ny, img.nz, img.nvol);
		if (shape == DN_SHAPE_SPHERE)
			dn_err("patch        : sphere, radius %.3g mm = %d voxels x %d volumes\n",
			        sqrt(g.r2), g.m, g.nvol);
		else
			dn_err("patch        : %dx%dx%d = %d voxels x %d volumes%s\n",
			        g.extent, g.extent, g.extent, g.m, g.nvol,
			        extent_was_auto ? " (auto)" : "");
		if (ngroups) {
			dn_err("demean       : %d shell%s from %s:", ngroups, ngroups == 1 ? "" : "s", fbval);
			for (int s = 0; s < ngroups; s++)
				fprintf(stderr, " b=%.0f x%d", shell_b[s], shell_n[s]);
			fprintf(stderr, "\n");
		} else if (demean)
			dn_err("demean       : none (no shell of 2+ volumes in %s)\n", fbval);
		dn_err("filter       : %s\n", filters[filter]);
		if (need_map)
			dn_err("sigma map    : from %d-voxel patches on every %dth voxel\n", g1.m, DN_SIGMA_STRIDE);
		if (vst_off_signed)
			dn_err("vst          : off, the input has negative values so is not magnitude\n");
		if (vst)
			dn_err("vst          : non-central chi, %d channel%s%s, %s inverse\n", ncoil,
			        ncoil == 1 ? "" : "s", ncoil == 1 ? " (Rician)" : "",
			        keep_bias ? "algebraic" : "exact-unbiased");
		if (aggregator != DN_AGG_EXCLUSIVE)
			dn_err("aggregator   : %s, patch centred on every %svoxel\n",
			        aggregators[aggregator],
			        stride == 1 ? "" : "2nd ");
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
			dn_err("degibbs      : yes, x-y planes%s (%d threads)\n", pf_note(pf),
			        dn_degibbs_threads(img.nx, img.ny, img.nz, img.nvol, nthreads_req, pf));
#endif
	}


	// Before anything is allocated for the denoising, and before the first write:
	// on failure this leaves img exactly as it was read, and the run stops with
	// nothing on disk.  Afterwards img.data is the real-valued rotated series and
	// every stage below is unchanged -- the denoiser neither knows nor cares that
	// its input is now signed.
	if (fphase && dn_rotate_to_real(&img, fphase, punits, quiet)) goto done;

#ifdef DN_DEGIBBS
	// AFTER the rotation, so it describes a run that got this far, and OUTSIDE
	// any -quiet test: -quiet suppresses the run SUMMARY, and dn_phase.c's
	// whole-turn warning already establishes that a warning about the data is not
	// a summary line.  The two options pull opposite ways -- -phase keeps
	// background noise zero-mean and degibbs truncates it -- and a scripted
	// -quiet run is the case that most needs telling.  One dn_err per line: it
	// prefixes per CALL.
	// Any signed input, not only -phase: a series rotated to real elsewhere is
	// truncated just the same, and once did so silently.
	int signed_input = fphase || vst_off_signed;
	for (size_t i = 0; degibbs == DN_DG_YES && !signed_input && i < img.nvox3d * (size_t)img.nvol; i++)
		signed_input = img.data[i] < 0.0f;
	if (degibbs == DN_DG_YES && signed_input) {
		dn_err("warning: %s, and -degibbs truncates it to zero,\n",
		       fphase ? "-phase output is signed" : "the input has negative values");
		dn_err("  reinstating a positive bias in background voxels.  Use -degibbs n if\n");
		dn_err("  the denoised series needs to stay signed.\n");
	}
#endif

	// After the rotation: the means are of the data the denoiser actually sees.
	if (ngroups) {
		means = dn_bval_means(img.data, img.nvox3d, img.nvol, group, ngroups);
		if (!means) goto done;
		g.group = g1.group = group;
		g.mean = g1.mean = means;
	}
	if (need_map) {
		sigma_map = dn_sigma_map(&g1, img.data, nthreads_req);
		if (!sigma_map) goto done;
		g.sigma_map = sigma_map;
	}
	// Stabilise in place, in units of the local sigma, so the denoiser then sees
	// unit-variance Gaussian noise everywhere and reads sigma = 1.
	//
	// The map above was estimated from RAW magnitudes, whose spread shrinks near
	// the noise floor, so it reads low there: 47 against a true 60 on benchmark/
	// at that noise level.  So sigma is re-estimated in the stabilised domain,
	// where it should come out as 1, and corrected by that factor (dwidenoise2
	// iterates for the same reason).  Each round recovers the raw value through
	// f's exact algebraic inverse rather than keeping a second copy of the series.
	if (vst) {
		const int filled = dn_sigma_fill(sigma_map, img.nvox3d);
		if (filled < 0) goto done;
		if (filled) {
			// Noiseless input: nothing to stabilise against, and nothing to denoise.
			dn_err("note: no noise level could be estimated, so the data are not stabilised\n");
			vst = 0;
		}
	}
	if (vst) {
		vt = dn_vst_create(ncoil);
		unit_map = (float *)dn_malloc(img.nvox3d, sizeof(float));
		live = (uint8_t *)dn_calloc(img.nvox3d, 1);
		if (!vt || !unit_map || !live) goto done;
		// A voxel zero in EVERY volume is masked-out background, not magnitude
		// samples, and stays zero so all-zero patches still short-circuit.  A lone
		// zero sample in a live voxel is real data (quantised magnitude can hit 0)
		// and is transformed with the rest of its vector.
		for (int j = 0; j < img.nvol; j++)
			for (size_t v = 0; v < img.nvox3d; v++)
				if (img.data[v + (size_t)j * img.nvox3d] != 0.0f) live[v] = 1;
		for (int it = 0; it < DN_VST_ROUNDS; it++) {
			// unit_map holds the sigma each voxel is in units of, until the end.
			vst_stabilise(vt, img.data, img.nvox3d, img.nvol, live, sigma_map, unit_map, it == 0);
			if (ngroups) {
				free(means);
				means = dn_bval_means(img.data, img.nvox3d, img.nvol, group, ngroups);
				if (!means) goto done;
				g.mean = g1.mean = means;
			}
			if (it == DN_VST_ROUNDS - 1) break;
			float *ratio = dn_sigma_map(&g1, img.data, nthreads_req);
			if (!ratio) goto done;
			for (size_t v = 0; v < img.nvox3d; v++)
				if (ratio[v] > 0.0f) sigma_map[v] *= ratio[v];
			free(ratio);
		}
		for (size_t v = 0; v < img.nvox3d; v++) unit_map[v] = 1.0f;
		g.sigma_map = unit_map;
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
	r.n_work = n_in_mask;
	r.nthreads = nthreads;
	r.aggregator = aggregator;
	r.fwhm = fwhm;
	r.stride = stride;

	if (aggregator == DN_AGG_EXCLUSIVE ? dn_run_execute(&r) : dn_agg_execute(&r)) goto done;

	if (vst) {
		vst_restore(vt, out, noise, img.nvox3d, img.nvol, live, mask, sigma_map, !keep_bias);
		// The unbiased inverse returns 0 for anything at or below the noise floor,
		// which is right for the model given -- and a mis-stated -noise_dof puts all
		// of the data there (measured: -noise_dof 4 on Rician data zeroed every
		// sample, at exit 0).  Say so rather than hand back an empty image.
		size_t nlive = 0, nzero = 0;
		for (size_t v = 0; v < img.nvox3d; v++) {
			if (!live[v] || (mask && !mask[v])) continue;
			for (int j = 0; j < img.nvol; j++) {
				nlive++;
				nzero += out[v + (size_t)j * img.nvox3d] == 0.0f;
			}
		}
		if (nzero * 2 > nlive) {
			dn_err("warning: %.0f%% of the output is zero, i.e. at or below the noise floor.\n",
			       100.0 * nzero / nlive);
			dn_err("  Check -noise_dof (%d), or use -preserve_noise_bias or -vst n.\n", ncoil);
		}
	}

#ifdef DN_DEGIBBS
	// After the denoiser, never before: degibbsing first would alter the noise
	// structure the denoiser models.  -noise and -rank describe the denoising and
	// are untouched, as is -real, which is the input before either stage.
	if (degibbs == DN_DG_YES &&
	    dn_degibbs(out, img.nx, img.ny, img.nz, img.nvol, nthreads_req, pf)) goto done;
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
	free(means);
	free(sigma_map);
	free(unit_map);
	dn_vst_free(vt);
	free(group);
	free(shell_b);
	free(shell_n);
	free(bval_auto);
	free(rank);
	free(mask);
	free(live);
	dn_image_free(&img);
	return status;
}
