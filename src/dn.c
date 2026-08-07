// checked allocation.

#include <stdlib.h>

#include "dn.h"

void *dn_calloc(size_t count, size_t size) {
	size_t bytes;
	if (dn_mul_size(count, size, &bytes)) {
		dn_err("allocation of %zu x %zu bytes overflows this build's address space\n", count, size);
		return NULL;
	}
	if (bytes == 0) bytes = 1;
	void *p = calloc(1, bytes);
	if (!p) dn_err("out of memory (%zu bytes)\n", bytes);
	return p;
}

void *dn_malloc(size_t count, size_t size) {
	size_t bytes;
	if (dn_mul_size(count, size, &bytes)) {
		dn_err("allocation of %zu x %zu bytes overflows this build's address space\n", count, size);
		return NULL;
	}
	if (bytes == 0) bytes = 1;
	void *p = malloc(bytes);
	if (!p) dn_err("out of memory (%zu bytes)\n", bytes);
	return p;
}
