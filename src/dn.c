// checked allocation, and the worker-pool launch both kernels share.

#include <pthread.h>
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

// One thread means no thread: it keeps the serial path usable under a debugger
// and makes it the obvious reference for the byte-identity gate.  Falling back
// to inline when NO thread could be created keeps a working run rather than
// failing one that would have succeeded slowly.
int dn_thread_run(void *(*fn)(void *), void *arg, int nthreads) {
	if (nthreads < 1) nthreads = 1;
	pthread_t *tid = (nthreads > 1)
	               ? (pthread_t *)dn_malloc((size_t)nthreads, sizeof(pthread_t)) : NULL;

	int started = 0;
	for (int i = 0; tid && i < nthreads; i++) {
		if (pthread_create(&tid[i], NULL, fn, arg) != 0) break;
		started++;
	}

	if (started == 0) {
		// The allocation failing counts here too, which is why this tests the
		// REQUEST rather than `tid`: either way the work still gets done.
		if (nthreads > 1) dn_err("could not start any worker thread; running single-threaded\n");
		fn(arg);
	} else if (started < nthreads) {
		// Partial creation used to continue silently, so a count already printed
		// by the CLI could overstate the team that really ran.  The result is
		// still correct -- the remaining workers drain the whole queue -- but the
		// reported number must not be a fiction.
		dn_err("only %d of %d worker threads could be started; continuing with %d\n",
		       started, nthreads, started);
	}

	for (int i = 0; i < started; i++) pthread_join(tid[i], NULL);
	free(tid);
	return started;
}
