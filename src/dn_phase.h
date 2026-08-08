// Rotation of complex (magnitude + phase) data onto the real axis.
//
// WHY.  Everything downstream of this file assumes additive noise that is
// Gaussian and zero-mean, because that is the model both the Marchenko-Pastur
// law and the Gavish-Donoho threshold are derived under.  Magnitude
// reconstruction breaks that assumption: taking |z| of a complex sample folds
// the noise distribution onto the positive half-line, so it becomes Rician (or
// non-central chi for a multi-coil reconstruction), no longer zero-mean, and no
// longer additive.  The bias this leaves -- the "noise floor" -- is negligible
// at high SNR and dominant at low SNR, which is precisely the regime where
// denoising is asked to do the most work.  A patch that straddles that floor has
// its small singular values inflated by a rectified DC term that no threshold
// can distinguish from signal.
//
// The fix is not to estimate the floor but to avoid creating it.  Given the
// phase, the complex samples can be rotated so that the object's signal lies
// along the real axis; the real part is then kept and the imaginary part
// discarded.  Noise in that real channel is still the original zero-mean
// Gaussian -- rotation is unitary, so it cannot change a circularly symmetric
// distribution -- and the values are signed, so nothing is rectified.
//
// This is the "MPPCA*" preprocessing of
//   J.P. Manzano Patron, S. Moeller, J.L.R. Andersson et al., "Denoising
//   diffusion MRI: considerations and implications for analysis", Imaging
//   Neuroscience 2, 2024,
// who report it as both more stable and better performing than denoising the
// complex data directly.  The rotation itself is the phase-removal front end of
// NORDIC (Moeller et al., NeuroImage 226:117539, 2021), distributed as
// NIFTI_COMP_to_REAL.m -- which is NIFTI_NORDIC.m with the denoising removed.
// That code was read as a specification, not translated: see dn_phase.c for the
// two places this implementation deliberately differs, and for the derivation
// that replaces its per-slice FFTs with a 7-tap convolution.
//
// No phase UNWRAPPING is involved, here or anywhere in this tool.  Every step
// works on the complex samples themselves, and a rotation by exp(-i*theta) is
// computed as conj(z)/|z| directly from the complex value, so wraps never arise
// and there is nothing for them to corrupt.

#ifndef DN_PHASE_H
#define DN_PHASE_H

#include "dn_nii.h"

#ifdef __cplusplus
extern "C" {
#endif

// How to read the phase values.  AUTO maps the observed range onto one turn,
// which is right for any encoding that covers a turn; data that does not cover
// one has to name its unit, because no heuristic can recover it.  The reasoning,
// and the failed attempt at guessing, are in dn_phase.c.
typedef enum {
	DN_PHASE_AUTO = 0,
	DN_PHASE_RADIANS,
	DN_PHASE_DEGREES,
	DN_PHASE_TURNS,     // cycles: one unit is a full turn
} dn_phase_units;

// Rotate `img` (a magnitude series) onto the real axis using the phase series in
// `phase_fname`, replacing img->data with the real part.  The result is SIGNED:
// background voxels scatter about zero instead of piling up on the noise floor.
//
// The phase image must match the magnitude in all four dimensions AND share its
// world transform; equal dimensions alone do not make two images the same
// picture.  Whichever unit convention is used is reported unless quiet.
//
// Returns 0 on success; on failure the reason has been reported and img is
// unmodified.
int dn_rotate_to_real(dn_image *img, const char *phase_fname,
                      dn_phase_units units, int quiet);

#ifdef __cplusplus
}
#endif

#endif // DN_PHASE_H
