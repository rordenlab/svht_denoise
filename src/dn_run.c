// whole-volume driver and worker pool.
//
// The pool is POSIX threads only, by design: the default build has no threading
// dependency to install.  Everything platform-specific is confined to this file,
// so a future Windows port replaces this one translation unit and leaves the
// kernel untouched.
//
// SCHEDULING AND DETERMINISM.  Work is handed out in chunks from a shared
// counter rather than partitioned statically up front.  That is safe for the
// the byte-identity guarantee below, for a reason worth stating explicitly: the
// order in which voxels are processed cannot affect any value, because each
// voxel reads only the immutable input and writes only its own outputs.  There
// is no reduction, no accumulator and no shared state in the arithmetic.  Static
// partitioning is the obvious alternative, but a masked run skips most
// voxels in whichever slabs fall outside the brain, so static blocks leave whole
// threads idle; dynamic chunking removes that for one mutex acquisition per
// DN_CHUNK voxels.

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#include "dn.h"
#include "dn_run.h"

// Voxels handed out per chunk.  The old fixed 2048 was ~0.9 s of work at
// 125x102, a large tail for one thread to be left holding on an irregular masked
// run.  64 is ~30 ms there and ~130 us on the smallest useful workload, so the
// mutex costs well under 0.1% either way -- which is why this is a constant
// rather than a cost model: the value only has to be small, not accurate.
#define DN_CHUNK 64

// Scratch bytes per worker, to the dominant term only.  This is a heuristic
// input to the thread cap, not an allocation, so it deliberately does NOT try to
// mirror every field in dn_work and dn_eig -- a transcription of another
// module's internals would just drift silently as those change.  The n*n and
// n*m blocks are what actually grow.
static size_t worker_bytes(const dn_geom *g) {
	const size_t n = (size_t)g->nvol, m = (size_t)g->m;
	return (n * m + 4 * n * n) * sizeof(double);
}

typedef struct {
	const dn_run *r;
	pthread_mutex_t lock;
	size_t next;        // next spatial voxel index to hand out
	size_t total;
	size_t chunk;
	int failed;         // set by any worker that could not proceed
	int eig_fail;       // count of eigensolver non-convergences
	unsigned long fallbacks;  // eigenvector fallbacks (see dn_eig_fallbacks)
} dn_pool;

int dn_default_threads(void) {
	long n = 0;
#ifdef __APPLE__
	// _SC_NPROCESSORS_ONLN is not exposed under a strict _POSIX_C_SOURCE on
	// Darwin, so ask the kernel directly.
	int ncpu = 0;
	size_t len = sizeof(ncpu);
	if (sysctlbyname("hw.logicalcpu", &ncpu, &len, NULL, 0) == 0 && ncpu > 0) n = ncpu;
#elif defined(_SC_NPROCESSORS_ONLN)
	n = sysconf(_SC_NPROCESSORS_ONLN);
#endif
	if (n < 1) return 1;
	if (n > 256) n = 256;
	return (int)n;
}

static void *worker(void *arg) {
	dn_pool *pool = (dn_pool *)arg;
	const dn_run *r = pool->r;
	const dn_geom *g = r->g;

	dn_work *w = dn_work_create(g);
	float *tmp = (float *)dn_malloc((size_t)g->nvol, sizeof(float));
	if (!w || !tmp) {
		pthread_mutex_lock(&pool->lock);
		pool->failed = 1;
		pthread_mutex_unlock(&pool->lock);
		dn_work_free(w);
		free(tmp);
		return NULL;
	}

	const int nx = g->nx, ny = g->ny;
	int local_eig_fail = 0;

	for (;;) {
		size_t lo, hi;
		pthread_mutex_lock(&pool->lock);
		if (pool->failed || pool->next >= pool->total) {
			pthread_mutex_unlock(&pool->lock);
			break;
		}
		lo = pool->next;
		hi = lo + pool->chunk;
		if (hi > pool->total) hi = pool->total;
		pool->next = hi;
		pthread_mutex_unlock(&pool->lock);

		for (size_t v = lo; v < hi; v++) {
			if (r->mask && !r->mask[v]) continue;

			// NIfTI within-volume order: x fastest, then y, then z.
			const int ix = (int)(v % (size_t)nx);
			const int iy = (int)((v / (size_t)nx) % (size_t)ny);
			const int iz = (int)(v / ((size_t)nx * ny));

			float sigma = 0.0f;
			uint16_t rank = 0;
			if (dn_denoise_voxel(g, w, r->img, ix, iy, iz, tmp, &sigma, &rank) != 0) {
				local_eig_fail++;
				continue;   // leave this voxel zero rather than write a bad basis
			}
			for (int j = 0; j < g->nvol; j++)
				r->out[v + (size_t)j * g->nvox3d] = tmp[j];
			if (r->noise) r->noise[v] = sigma;
			if (r->rank) r->rank[v] = rank;
		}
	}

	const unsigned long fb = dn_eig_fallbacks(w->eig);
	dn_work_free(w);
	free(tmp);

	if (local_eig_fail || fb) {
		pthread_mutex_lock(&pool->lock);
		pool->eig_fail += local_eig_fail;
		pool->fallbacks += fb;
		pthread_mutex_unlock(&pool->lock);
	}
	return NULL;
}

int dn_effective_threads(const dn_geom *g, size_t n_work, int requested) {
	if (!g) return 1;
	if (requested < 1) requested = 1;
	// Voxels that will actually be VISITED, not the whole grid: a one-voxel mask
	// otherwise started every core and gave each a full eigensolver arena for a
	// single eigensolve.  Thread count does not affect the output, so this is
	// free of the byte-identity promise.
	const size_t total = n_work;
	// More workers than cores never helps a CPU-bound kernel, and more workers
	// than chunks is pure overhead -- the extras allocate a full scratch arena,
	// find the queue empty and exit.
	const int hw = dn_default_threads();
	if (requested > hw) requested = hw;
	const size_t chunks = (total + DN_CHUNK - 1) / DN_CHUNK;
	if (chunks > 0 && (size_t)requested > chunks) requested = (int)chunks;
	const size_t per = worker_bytes(g);
	const size_t budget = (size_t)1 << 30;   // 1 GiB of scratch, in total
	if (per > 0 && (size_t)requested * per > budget) {
		int cap = (int)(budget / per);
		requested = (cap < 1) ? 1 : cap;
	}
	return requested;
}

int dn_run_execute(const dn_run *r) {
	if (!r || !r->g || !r->img || !r->out) return 1;
	int nthreads = dn_effective_threads(r->g, r->n_work, r->nthreads);

	dn_pool pool;
	memset(&pool, 0, sizeof(pool));
	pool.r = r;
	pool.total = r->g->nvox3d;
	pool.chunk = DN_CHUNK;

	if (pthread_mutex_init(&pool.lock, NULL) != 0) {
		dn_err("cannot create the worker mutex\n");
		return 1;
	}

	// One thread means no thread: keeps -nthreads 1 usable under a debugger and
	// makes the serial path the obvious reference for the byte-identity gate.
	if (nthreads == 1) {
		worker(&pool);
	} else {
		pthread_t *tid = (pthread_t *)dn_malloc((size_t)nthreads, sizeof(pthread_t));
		if (!tid) { pthread_mutex_destroy(&pool.lock); return 1; }
		int started = 0;
		for (int i = 0; i < nthreads; i++) {
			if (pthread_create(&tid[i], NULL, worker, &pool) != 0) break;
			started++;
		}
		if (started == 0) {
			// Fall back to running inline rather than failing the whole job.
			dn_err("could not start any worker thread; running single-threaded\n");
			worker(&pool);
		} else if (started < nthreads) {
			// Partial creation used to continue silently, so the count the CLI
			// had already printed could overstate the team that actually ran.
			// The result is still correct -- the remaining workers drain the
			// whole queue -- but the reported number must not be a fiction.
			dn_err("only %d of %d worker threads could be started; continuing with %d\n",
			       started, nthreads, started);
		}
		for (int i = 0; i < started; i++) pthread_join(tid[i], NULL);
		free(tid);
	}

	pthread_mutex_destroy(&pool.lock);

	if (pool.failed) {
		dn_err("a worker could not allocate its scratch space\n");
		return 1;
	}
	if (pool.eig_fail) {
		dn_err("the eigensolver failed to converge for %d patch%s.\n",
		       pool.eig_fail, pool.eig_fail == 1 ? "" : "es");
		dn_err("  Those voxels were left at zero rather than written from an\n");
		dn_err("  unconverged basis. This should not happen; please report it.\n");
		return 1;
	}
	if (pool.fallbacks) {
		// Not an error: the fallback IS the trusted full solve, so the values are
		// right. Reported because a non-zero count means inverse iteration is
		// struggling on this data, which is worth knowing before it gets worse.
		dn_err("note: inverse iteration fell back to the full solve for %lu patch%s\n",
		       pool.fallbacks, pool.fallbacks == 1ul ? "" : "es");
	}
	return 0;
}
