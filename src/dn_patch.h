// patch geometry and the single-patch denoising kernel.

#ifndef DN_PATCH_H
#define DN_PATCH_H

#include <stdint.h>

#include "dn_eig.h"

// Element type of the patch buffer; see the note on `pt` below.
#ifdef DN_USE_ACCELERATE
typedef double dn_pt_t;
#else
typedef float dn_pt_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Image and patch geometry.  Everything here is constant for a whole run: with
// the auto-extent rule the patch is always full size, so beta -- and therefore
// omega -- never varies from voxel to voxel and is computed exactly once.
typedef struct {
	int nx, ny, nz;      // spatial dimensions
	int nvol;            // volumes (N)
	int extent;          // patch side length k (odd)
	int m;               // k*k*k voxels per patch (M), M > N
	size_t nvox3d;       // nx*ny*nz
	double beta;         // N/M
	double omega;        // Eq. (4) coefficient
	double sigma_scale;  // 1 / sqrt(M * mu_beta); sigma = median_s * sigma_scale
} dn_geom;

// Fill geometry for the given image shape and patch extent.  Returns 0 on
// success; on failure the reason has already been reported.
int dn_geom_init(dn_geom *g, int nx, int ny, int nz, int nvol, int extent);

// The default extent: the smallest odd k with k*k*k > nvol.  Returns 0 if no
// admissible k exists within DN_MAX_EXTENT.
int dn_auto_extent(int nvol);

// Per-worker scratch.  One arena per thread; nothing here is shared.
typedef struct {
	// float on the portable path, and that is exact rather than a compromise:
	// every value in here is copied straight out of a float32 image, so widening
	// it at the point of use gives bit-for-bit the same operands a double array
	// would have held, at half the footprint -- 185 kB against 370 kB at 138x343,
	// which is what the blocked Gram loop's inner passes stream.
	//
	// DOUBLE under Accelerate, because BLAS has no mixed-precision syrk: dsyrk
	// wants the operand in the precision it accumulates in.  Paying the wider
	// gather to reach it is measured as worth it; see the Gram kernel.
	dn_pt_t *pt;     // n*m, patch transposed: pt[j*m + i] = voxel i, volume j
	double *gram;    // n*n
	double *evecs;   // n*n
	double *evals;   // n
	double *s;       // n, singular values ascending
	double *yc;      // n, the centre voxel's own vector
	double *acc;     // n, projection accumulator (kept in double, cast once)
	int32_t *offset; // m, flat spatial offsets of the patch voxels from its origin
	dn_eig *eig;
} dn_work;

dn_work *dn_work_create(const dn_geom *g);
void dn_work_free(dn_work *w);

// Denoise the voxel at spatial index (ix, iy, iz).
//
//   img       float32, VOLUME-MAJOR: img[voxel + volume*nvox3d]
//   out       N denoised values for this voxel (caller-provided, length nvol)
//   sigma     estimated noise level for this patch
//   rank      number of retained components
//
// Returns 0 on success, non-zero if the eigensolver failed to converge.
int dn_denoise_voxel(const dn_geom *g, dn_work *w, const float *img,
                     int ix, int iy, int iz,
                     float *out, float *sigma, uint16_t *rank);

#ifdef __cplusplus
}
#endif

#endif // DN_PATCH_H
