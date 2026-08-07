// shared types, diagnostics and checked allocation.

#ifndef DN_H
#define DN_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DN_NAME "svht_denoise"

// Diagnostics ALWAYS go to stderr.  One macro: an earlier dn_info was
// byte-identical to dn_err, and a name promising a distinction that does not
// exist eventually catches someone out.  Nothing in v1 writes an image to stdout, but
// keeping the discipline from the start means adding that later cannot corrupt
// the stream (niimath learned this the hard way -- see its printfx convention).
#define dn_err(...)  do { fprintf(stderr, DN_NAME ": "); fprintf(stderr, __VA_ARGS__); } while (0)

// Largest volume count the N x N dense solver will accept.  The Gram matrix and
// eigenvector block are each N*N doubles per worker, so 4096 volumes would already
// be 128 MB per thread; anything approaching this is a user error, not a workload.
#define DN_MAX_VOL 4096

// Largest patch side length.  k = 15 is 3375 voxels, far past anything useful.
#define DN_MAX_EXTENT 15

// Checked size arithmetic.  Adapted from niimath's nii_mul_size (core.c, BSD-2,
// rordenlab/niimath).  Returns 1 on overflow, 0 on success.
static inline int dn_mul_size(size_t a, size_t b, size_t *out) {
	if (a != 0 && b > SIZE_MAX / a) return 1;
	*out = a * b;
	return 0;
}

// Checked allocation.  Unlike niimath's nii_calloc these return NULL rather than
// calling exit(), because every caller here already has a cleanup path and a
// library-shaped kernel should not terminate the process.
void *dn_calloc(size_t count, size_t size);
void *dn_malloc(size_t count, size_t size);

#ifdef __cplusplus
}
#endif

#endif // DN_H
