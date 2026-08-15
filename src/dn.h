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

// Largest volume count the N x N dense solver will accept.  Several N*N double
// blocks are held per worker -- dn_run.c's worker_bytes models four -- and one of
// them alone is 128 MiB at 4096 volumes, so that many would already be half a
// gigabyte per thread; anything approaching this is a user error, not a workload.
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

// Run `fn(arg)` on `nthreads` workers and join them.  Returns the number that
// actually started; 0 means it ran inline on this thread, which is what
// nthreads == 1 asks for and also what happens if no thread could be created.
//
// Shared by dn_run.c and mrdegibbs/dg.c because a PARTIAL failure has to be
// reported: both print a thread count before starting, and continuing quietly
// with fewer workers than announced makes that number a fiction.  Having one
// copy is why degibbs cannot drift back into not saying so.
int dn_thread_run(void *(*fn)(void *), void *arg, int nthreads);

#ifdef __cplusplus
}
#endif

#endif // DN_H
