// checked allocation.

#include <stdlib.h>

#include "dn.h"

// One body: the two differ only in whether the block comes back zeroed, and
// duplicating the overflow check and both messages to express that invited them
// to drift apart.
static void *dn_alloc(size_t count, size_t size, int zero) {
	size_t bytes;
	if (dn_mul_size(count, size, &bytes)) {
		dn_err("allocation of %zu x %zu bytes overflows this build's address space\n", count, size);
		return NULL;
	}
	if (bytes == 0) bytes = 1;
	void *p = zero ? calloc(1, bytes) : malloc(bytes);
	if (!p) dn_err("out of memory (%zu bytes)\n", bytes);
	return p;
}

void *dn_calloc(size_t count, size_t size) { return dn_alloc(count, size, 1); }
void *dn_malloc(size_t count, size_t size) { return dn_alloc(count, size, 0); }
