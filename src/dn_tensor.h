// Tensor denoising: the single-patch recursive kernel and the sliding-window driver.
//
// Each patch is a space x time x echo TENSOR rather than a matrix, decomposed
// one mode at a time in the recursive scheme of Olesen et al. (Tensor denoising
// of multidimensional MRI data, Magn Reson Med 2022, doi:10.1002/mrm.29478):
// unfold, keep the components above the noise bulk edge, pass the reduced
// V*S on to the next mode, and rebuild on the way back.  The noise level is
// NOT estimated here.  The kernel takes it as given, and the CLI supplies the
// same Gavish & Donoho median-estimator map the matrix denoiser reads, so no
// Marchenko-Pastur fit is used anywhere -- the one place this departs from
// Olesen et al., whose MP-PCA pooled fit is deliberately not implemented.  By
// default the last mode's kept values get Gavish & Donoho's 2017 Frobenius
// shrinkage; DN_TFILTER_OPTTHRESH replaces the edge cut at every level with
// their 2014 known-noise hard threshold instead.

#ifndef DN_TENSOR_H
#define DN_TENSOR_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DN_TENSOR_MAX_MODES 8

// How each level chooses, and the last level weights, its components.  The
// order matches -tensor_filter's choices.
enum {
	DN_TFILTER_OPTSHRINK = 0, // cut at the bulk edge; GD 2017 shrinkage on the last mode
	DN_TFILTER_TRUNCATE,      // cut at the bulk edge; nothing shrunk
	DN_TFILTER_OPTTHRESH      // GD 2014 known-sigma hard threshold at every level
};

// Kernel settings.  Everything here is constant for a whole run.
typedef struct {
	int nmodes;                       // K >= 3; a matrix is not denoised here
	int dims[DN_TENSOR_MAX_MODES];    // every one >= 2
	int filter;         // DN_TFILTER_*
	int shrink_observed; // optshrink: GD's observed m x c, not Olesen's (m-r) x (c-r)
	int subtract_mean;  // center each unfolding along its longer side first
} dn_tensor_opts;

typedef struct dn_tensor_work dn_tensor_work;

// Can the kernel take this geometry?  0, or non-zero with the reason reported.
int dn_tensor_check(const dn_tensor_opts *o);

// Per-worker scratch for one geometry.  NULL on a bad geometry (reported) or on
// allocation failure.
dn_tensor_work *dn_tensor_work_create(const dn_tensor_opts *o);
void dn_tensor_work_free(dn_tensor_work *w);

// Scratch bytes per worker, to the dominant terms: an input to the thread cap,
// not an allocation.
size_t dn_tensor_work_bytes(const dn_tensor_opts *o);

// Denoise one patch IN PLACE, against the noise VARIANCE sigma2 (finite, >= 0).
// x holds prod(dims) values in COLUMN-MAJOR order: x[i0 + d0*(i1 + d1*(i2 + ...))].
// On return ranks[k] is the number of components kept at level k, 0 past the
// level where the recursion stopped.  Returns 0, or non-zero if the
// eigensolver failed or sigma2 was not a usable variance.
int dn_tensor_denoise(const dn_tensor_opts *o, dn_tensor_work *w, double *x,
                      double sigma2, int *ranks);

// Eigenvector fallbacks to the full solve, as dn_eig_fallbacks.
unsigned long dn_tensor_fallbacks(const dn_tensor_work *w);

// ---------------------------------------------------------------------------
// Sliding-window driver (dn_tensor_run.c)
// ---------------------------------------------------------------------------

// The three tensor modes of a multi-echo patch.
enum { DN_TMODE_SPACE = 0, DN_TMODE_TIME, DN_TMODE_ECHO };

typedef struct {
	dn_tensor_opts k;   // shrink/mean; nmodes and dims are filled in by dn_tensor_setup
	int nx, ny, nz, nt, ne;
	// Window per spatial axis.  A two-element window is win[2] = nz with the
	// patch spanning every slice; full_z records that, because center
	// assignment then writes the center COLUMN rather than one voxel.
	int win[3], stride[3], full_z;
	int order[3];         // tensor mode k is DN_TMODE_* order[k]
	int center_assign;    // write only each window's center (column)
	const float *const *echo;   // ne inputs, each nvox3d*nt, volume-major
	float *const *out;          // ne outputs, same layout
	const uint8_t *mask;        // nvox3d, or NULL: a window is denoised if it touches it
	// Noise variance: per voxel (each window uses the mean over its voxels), or
	// NULL for the single value sigma2 everywhere.
	const float *noise_var;
	double sigma2;
	float *sigma2_out;  // nvox3d mean window noise variance, or NULL
	float *rank;        // nvox3d*3, mode k at rank[k*nvox3d + v], or NULL
	float *gain;        // nvox3d SNR gain heuristic, or NULL
	int nthreads;
} dn_trun;

// Fill r->k.nmodes and r->k.dims from the geometry and order, and check that
// the patch is one the kernel can take.  Returns 0, or non-zero (reported).
int dn_tensor_setup(dn_trun *r);

// The team size dn_tensor_execute will ACTUALLY use.
int dn_tensor_threads(const dn_trun *r, int requested);

// Denoise every window position, averaging overlapping estimates.  Voxels no
// window covered keep their input values and read NaN in every map.  Output
// is byte-identical at any thread count.  Returns 0 on success.
int dn_tensor_execute(const dn_trun *r);

#ifdef __cplusplus
}
#endif

#endif // DN_TENSOR_H
