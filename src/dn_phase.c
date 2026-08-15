// Rotation of complex (magnitude + phase) data onto the real axis.
// See dn_phase.h for why this is worth doing at all.
//
// THE METHOD, in the order the code performs it:
//
//   1. Rescale the phase image onto one full turn and form z = |mag| * exp(i*phase).
//   2. Divide out the STATIC phase -- the per-voxel background phase common to the
//      whole series, read off a single high-signal volume.  This removes coil and
//      B0 phase, which is by far the largest term and is not diffusion related.
//   3. Divide out the residual DYNAMIC phase -- the per-volume, per-slice phase
//      left over after step 2, which for diffusion data is dominated by bulk
//      motion during the diffusion gradients.  It is estimated by low-pass
//      filtering the complex image in plane, on the argument that true phase
//      structure is smooth while noise is not.
//   4. Keep the real part.
//
// Steps 2 and 3 both divide by a unit complex number, computed as conj(z)/|z|
// rather than exp(-i*atan2(...)): same value, no transcendental, and no wrapped
// angle ever exists.  Where |z| is exactly zero the rotation is taken to be 1,
// matching MATLAB's angle(0) == 0, which is the convention the reference
// implementation inherits.
//
// RELATION TO THE REFERENCE.  Written from NIFTI_COMP_to_REAL.m as a
// specification; see dn_phase.h for the provenance.  Two deliberate divergences:
//
//   a) The reference normalises the whole series by the smallest non-zero
//      magnitude of the first volume, then multiplies it back at the end.  Every
//      operation in between is scale invariant -- they are all divisions by
//      conj(z)/|z| -- so this is the identity up to float rounding, and is
//      omitted.
//
//   b) The reference estimates the dynamic phase with a per-slice, per-volume
//      2D FFT pair, multiplying k-space by a separable tukeywin(N,1)^3 window.
//      That is done here as a 7-tap convolution instead, with NO FFT anywhere in
//      this tool.  The two are the same operation:
//
//        tukeywin(N,1) is a Hann window.  Written against the SIGNED frequency
//        index k, with DC at k = 0 where the shifted array puts it, the periodic
//        Hann window is
//            w[k] = 0.5*(1 - cos(2*pi*(k + N/2)/N))
//                 = 0.5*(1 + cos(2*pi*k/N))
//                 = 0.5 + 0.25*exp(i*2*pi*k/N) + 0.25*exp(-i*2*pi*k/N).
//        Multiplying k-space by exp(-+ i*2*pi*k/N) is a circular shift of one
//        voxel in image space, so multiplying by w is exactly a circular
//        convolution with [1 2 1]/4.  Cubing the window cubes the kernel:
//            [1 2 1]^3 / 64 = [1 6 15 20 15 6 1]/64,
//        the 7-tap binomial below, applied along x and then y.  Verified
//        numerically against the FFT form at 5e-16 relative error.
//
//      The reference uses MATLAB's SYMMETRIC Hann, whose denominator is N-1
//      rather than N, which puts the window's peak half a sample away from DC
//      and makes it not quite a shift-invariant kernel.  The periodic window
//      used here is the one that is actually centred on DC.  On the EDDEN test
//      data the two agree over the brain to a correlation of 0.9999985, with an
//      RMS difference of 0.105% of the reference RMS.
//
// COST.  The reference is O(N log N) per slice with a large constant and needs
// an FFT library; this is 14 multiply-adds per voxel with none.  The convolution
// is the cheap part of the rotation anyway -- on the EDDEN 104x104x72x21 series
// it is ~11% of it, against ~36% for the one cos/sin pair per sample and ~40%
// for reading the phase image.

#include <math.h>
#include <stdlib.h>

#include "dn.h"
#include "dn_nii.h"
#include "dn_phase.h"

// M_PI is XSI, not C99, and this builds under a bare _POSIX_C_SOURCE.
#define DN_TWO_PI 6.283185307179586476925286766559

// [1 2 1]^3 / 64 -- the image-space equivalent of the k-space Hann^3 low-pass.
#define DN_TAPS 7
#define DN_HALF 3
static const float DN_B7[DN_TAPS] = {
	1.0f / 64.0f, 6.0f / 64.0f, 15.0f / 64.0f, 20.0f / 64.0f,
	15.0f / 64.0f, 6.0f / 64.0f, 1.0f / 64.0f
};

// Circular convolution along x, over nrow contiguous rows of length nx.
//
// The boundary is circular because that is what the FFT form it replaces does;
// it is also the only choice that keeps the operation shift invariant.  It is
// harmless here in any case: this filter only ever supplies a phase estimate,
// and the wrap-around contamination is confined to three voxels of background
// air at each edge.
static void conv_x(const float *src, float *dst, int nx, size_t nrow) {
	for (size_t r = 0; r < nrow; r++) {
		const float *s = src + r * (size_t)nx;
		float *d = dst + r * (size_t)nx;
		// Interior first, where no index can leave the row.
		const int lo = (nx > DN_TAPS) ? DN_HALF : nx;
		const int hi = (nx > DN_TAPS) ? nx - DN_HALF : nx;
		for (int i = lo; i < hi; i++) {
			float acc = 0.0f;
			for (int t = 0; t < DN_TAPS; t++) acc += DN_B7[t] * s[i + t - DN_HALF];
			d[i] = acc;
		}
		// Edges, and the whole row when it is too short to have an interior.
		for (int i = 0; i < nx; i++) {
			if (i >= lo && i < hi) continue;
			float acc = 0.0f;
			for (int t = 0; t < DN_TAPS; t++) {
				int j = i + t - DN_HALF;
				// nx can be as small as 1, so wrap by looping rather than by a
				// single modulo of a possibly still-negative value.
				while (j < 0) j += nx;
				while (j >= nx) j -= nx;
				acc += DN_B7[t] * s[j];
			}
			d[i] = acc;
		}
	}
}

// Circular convolution along y, slice by slice.  Accumulating a whole output row
// at a time keeps both operands contiguous in x and reduces the wrap arithmetic
// to once per tap per row.
static void conv_y(const float *src, float *dst, int nx, int ny, int nz) {
	const size_t slice = (size_t)nx * (size_t)ny;
	for (int z = 0; z < nz; z++) {
		const float *sp = src + (size_t)z * slice;
		float *dp = dst + (size_t)z * slice;
		for (int y = 0; y < ny; y++) {
			float *d = dp + (size_t)y * (size_t)nx;
			for (int i = 0; i < nx; i++) d[i] = 0.0f;
			for (int t = 0; t < DN_TAPS; t++) {
				int j = y + t - DN_HALF;
				while (j < 0) j += ny;
				while (j >= ny) j -= ny;
				const float *s = sp + (size_t)j * (size_t)nx;
				const float b = DN_B7[t];
				for (int i = 0; i < nx; i++) d[i] += b * s[i];
			}
		}
	}
}

static double volume_mean_abs(const float *mag, size_t n3, int t) {
	double sum = 0.0;
	const float *m = mag + (size_t)t * n3;
	for (size_t v = 0; v < n3; v++) sum += fabs((double)m[v]);
	return sum / (double)n3;
}

// The volume whose static phase is used as the reference: the FIRST volume whose
// mean magnitude reaches 95% of the largest.  For a diffusion series that is a
// b = 0 image, which is what you want -- it has the best SNR, so its phase
// estimate is the least noisy, and it is the least affected by the motion-induced
// phase that step 3 then mops up.
//
// If the magnitude is zero everywhere the test can never fire and volume 0 is
// used; the output is all zeros either way.
static int reference_volume(const float *mag, size_t n3, int nvol) {
	double best = 0.0;
	for (int t = 0; t < nvol; t++) {
		const double mean = volume_mean_abs(mag, n3, t);
		if (mean > best) best = mean;
	}
	// Rescanning rather than storing the means keeps this allocation-free, and so
	// failure-free -- there is no half-answer to report.  The second pass stops at
	// the first hit, which for a diffusion series is normally volume 0.
	for (int t = 0; t < nvol; t++)
		if (volume_mean_abs(mag, n3, t) > 0.95 * best) return t;
	return 0;
}

int dn_rotate_to_real(dn_image *img, const char *phase_fname,
                      dn_phase_units units, int quiet) {
	if (!img || !img->data || !phase_fname) return 1;

	dn_image ph;
	if (dn_image_read(phase_fname, "phase", 0, &ph)) return 1;

	int rc = 1;
	float *re = NULL, *im = NULL, *sre = NULL, *sim = NULL, *tmp = NULL;
	float *refc = NULL, *refs = NULL;

	if (ph.nx != img->nx || ph.ny != img->ny || ph.nz != img->nz || ph.nvol != img->nvol) {
		dn_err("phase image '%s' is %dx%dx%dx%d but the input is %dx%dx%dx%d.\n",
		       phase_fname, ph.nx, ph.ny, ph.nz, ph.nvol,
		       img->nx, img->ny, img->nz, img->nvol);
		dn_err("  The two must be the magnitude and phase of the SAME series, one\n");
		dn_err("  phase volume per magnitude volume.\n");
		goto done;
	}
	// Equal dimensions do not make two images the same picture.  Written as
	// !(<=) so a non-finite offset fails closed.
	const float off = dn_grid_offset_mm(img, &ph);
	if (!(off <= DN_GRID_TOL_MM)) {
		dn_err("phase image '%s' does not share the input's grid.\n", phase_fname);
		// With no qform and no sform, nifti synthesises a transform from pixdim
		// with a ZERO origin.  Comparing that against a partner that does carry an
		// origin yields a large displacement that says nothing about the data, so
		// quoting it would send the user looking for the wrong problem.
		const int a_pos = (img->nim->sform_code > 0 || img->nim->qform_code > 0);
		const int b_pos = (ph.nim->sform_code > 0 || ph.nim->qform_code > 0);
		if (a_pos != b_pos) {
			dn_err("  One of the two carries no qform and no sform, so its position is\n");
			dn_err("  unknown and was taken as the origin.  Give both images a spatial\n");
			dn_err("  transform, or copy one header's geometry onto the other.\n");
		} else {
			dn_err("  The two headers place the same voxel up to %g mm apart; voxel size,\n", (double)off);
			dn_err("  orientation and origin must match, or unrelated voxels get combined.\n");
		}
		goto done;
	}

	const size_t n3 = img->nvox3d;
	const int nvol = img->nvol;
	const size_t nvox = n3 * (size_t)nvol;   // dn_image_read already proved this fits
	// Rows for the x pass: every product in this file is widened before it is
	// formed, and ny*nz in plain int would be the one exception.
	const size_t nrow = (size_t)img->ny * (size_t)img->nz;

	// Work out what the phase is measured in.  Only the SCALE has to be right:
	// any constant offset -- including the half turn between a [0, 2pi) and a
	// [-pi, pi) convention -- is a global unit factor that step 2 divides
	// straight back out.  That is why the explicit cases below set no intercept.
	//
	// AUTO maps the OBSERVED RANGE onto one full turn, as the reference
	// implementation and niimath's md_rescale_phase both do.  This is right for
	// EVERY encoding that covers a whole turn, whatever its unit -- scanner
	// integers, degrees, cycles, [-1, 1] -- and every real acquired phase image
	// covers one, because its background is noise whose phase is uniform on
	// [-pi, pi].
	//
	// It is wrong for data that does not span a turn, and no heuristic can fix
	// that: a range of [-0.5, 0.5] is half a radian or one full cycle, and the
	// file does not say which.  An earlier attempt to resolve it by the value
	// MAGNITUDES made that guess silently, and got cycles wrong -- the tool then
	// no-opped, returning the magnitude and reporting success.  Guessing was
	// replaced by asking: `units`.
	double lo = (double)ph.data[0], hi = lo;
	for (size_t i = 1; i < nvox; i++) {
		const double v = (double)ph.data[i];
		if (v < lo) lo = v;
		if (v > hi) hi = v;
	}
	const double span = hi - lo;

	double slope = 1.0, inter = 0.0;
	const char *unit;
	switch (units) {
	case DN_PHASE_RADIANS:
		unit = "radians (given)";
		break;
	case DN_PHASE_DEGREES:
		unit = "degrees (given)";
		slope = DN_TWO_PI / 360.0;
		break;
	case DN_PHASE_TURNS:
		unit = "turns (given)";
		slope = DN_TWO_PI;
		break;
	case DN_PHASE_AUTO:
	default:
		if (fabs(span - DN_TWO_PI) <= 0.1) {
			// Already a full turn in radians; rescaling would only add rounding.
			unit = "radians (a full turn observed)";
		} else if (span > 0.0) {
			unit = "observed range mapped to one turn";
			slope = DN_TWO_PI / span;
			inter = -0.5 * (lo + hi) * slope;
		} else {
			// Constant.  Not an error: a constant phase is a constant rotation and
			// step 2 removes it whatever it is, so the output is the magnitude.
			unit = "constant";
		}
		break;
	}

	// Allocated up front, so that a failure here still leaves img untouched.
	re   = (float *)dn_malloc(n3, sizeof(float));
	im   = (float *)dn_malloc(n3, sizeof(float));
	sre  = (float *)dn_malloc(n3, sizeof(float));
	sim  = (float *)dn_malloc(n3, sizeof(float));
	tmp  = (float *)dn_malloc(n3, sizeof(float));
	refc = (float *)dn_malloc(n3, sizeof(float));
	refs = (float *)dn_malloc(n3, sizeof(float));
	if (!re || !im || !sre || !sim || !tmp || !refc || !refs) goto done;

	const int ref = reference_volume(img->data, n3, nvol);
	// exp(i*rad) == 1 at every voxel iff sin(rad) is everywhere 0 AND cos(rad) is
	// everywhere +1.  Both are needed: sin alone admits rad = k*pi, where the
	// factor is -1 at odd k, which is a real per-voxel sign modulation and NOT the
	// identity -- a half-turn fixture changes the output by ~100% while passing a
	// sin-only test, so the warning below would have stated something false.
	double max_sin = 0.0, min_cos = 1.0;

	// Step 2 setup: the unit conjugate of the reference volume.  Because
	// z = a*exp(i*p) with a >= 0, |z| == a and conj(z)/|z| is just
	// (cos p, -sin p) -- the magnitude cancels and never has to be divided by.
	{
		const float *m = img->data + (size_t)ref * n3;
		const float *p = ph.data + (size_t)ref * n3;
		for (size_t v = 0; v < n3; v++) {
			if (!(fabs((double)m[v]) > 0.0)) { refc[v] = 1.0f; refs[v] = 0.0f; continue; }
			const double rad = (double)p[v] * slope + inter;
			refc[v] = (float)cos(rad);
			refs[v] = (float)-sin(rad);
		}
	}

	if (!quiet) {
		dn_err("phase        : %s (range %g to %g, read as %s)\n",
		       phase_fname, lo, hi, unit);
		dn_err("static phase : volume %d of %d\n", ref, nvol);
	}

	// One volume at a time, so the scratch is a fixed number of 3D buffers rather
	// than two full complex copies of the series.
	for (int t = 0; t < nvol; t++) {
		float *dst = img->data + (size_t)t * n3;      // the magnitude, overwritten below
		const float *p = ph.data + (size_t)t * n3;

		// Form the complex sample and divide out the static phase.  The magnitude
		// is taken as |mag| because a magnitude image has no meaningful sign and
		// the reference does the same.
		for (size_t v = 0; v < n3; v++) {
			const double a = fabs((double)dst[v]);
			const double rad = (double)p[v] * slope + inter;
			const double sr = sin(rad), cr = cos(rad);
			const double asr = fabs(sr);
			if (asr > max_sin) max_sin = asr;
			if (cr < min_cos) min_cos = cr;
			const double c = a * cr, d = a * sr;
			re[v] = (float)(c * (double)refc[v] - d * (double)refs[v]);
			im[v] = (float)(c * (double)refs[v] + d * (double)refc[v]);
		}

		// Step 3: the low-pass estimate of what phase is left.
		conv_x(re, tmp, img->nx, nrow);
		conv_y(tmp, sre, img->nx, img->ny, img->nz);
		conv_x(im, tmp, img->nx, nrow);
		conv_y(tmp, sim, img->nx, img->ny, img->nz);

		// Divide it out and keep the real part.  Only the real part is wanted, so
		// the imaginary half of the product is never formed:
		//   Re[ w * conj(s)/|s| ] = (Re w * Re s + Im w * Im s) / |s|.
		// Reading re[v]/im[v] and writing dst[v] in the same pass is safe: this
		// volume's magnitude was consumed by the loop above and is not read again.
		for (size_t v = 0; v < n3; v++) {
			const double a = (double)sre[v], b = (double)sim[v];
			const double q = a * a + b * b;
			dst[v] = (q > 0.0)
			       ? (float)(((double)re[v] * a + (double)im[v] * b) / sqrt(q))
			       : re[v];
		}
	}

	// An earlier version warned whenever the output contained no negative value.
	// That is not the same question: a high-SNR or noiseless series can be
	// rotated correctly and still come back entirely positive, so the test
	// produced false alarms on exactly the data it was least useful for.
	//
	// This is the real condition, and it is the one that was actually observed:
	// a unit that makes every phase a whole number of turns -- `-phaseunits
	// turns` on integer-valued phase -- gives exp(i*2*pi*n) == 1 at every voxel,
	// so the rotation is the identity and the run reports success having handed
	// back the magnitude unchanged.  Well clear of any rounding: the residual on
	// a genuine whole-turn phase is ~1e-11, and real data reaches ~1.
	// Only a false NEGATIVE remains possible: a phase a hair off a whole turn
	// escapes the threshold while still rounding to no change.  A missed warning
	// is acceptable; a warning that misdescribes the data is not, which is why
	// this tests the phase exactly rather than guessing from the output.
	if (max_sin < 1e-6 && min_cos > 0.0) {
		dn_err("warning: every phase value is a whole number of turns, so the rotation\n");
		dn_err("  is the identity and the output is the magnitude unchanged.  The phase\n");
		dn_err("  was read as %s; check that this is the right unit.\n", unit);
	}

	rc = 0;

done:
	free(re); free(im); free(sre); free(sim); free(tmp);
	free(refc); free(refs);
	dn_image_free(&ph);
	return rc;
}
