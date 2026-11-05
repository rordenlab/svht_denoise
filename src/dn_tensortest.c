// Unit checks for the tensor kernel (dn_tensor.c), on in-memory tensors.  Run by
// `make test` before test.sh.  Not part of svht_denoise.
//
// These pin properties that hold by construction and that no image-level
// comparison can isolate:
//   - a noise variance of 0 keeps every component and returns the input, with
//     and without centring (which takes the reverse pass through the means);
//   - an all-zero patch stays zero with nothing kept;
//   - given the true variance, the output is closer to a low-rank truth than
//     the input, on tall, wide and mixed shapes, under the edge cut and
//     optthresh alike;
//   - a component between the bulk edge and the optimal hard threshold is
//     kept by the one and dropped by the other, and optthresh keeps what it
//     keeps unshrunk;
//   - an unusable variance or geometry is refused rather than run.

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dn_tensor.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int pass = 0, fail = 0;

static void report(int ok, const char *what) {
	if (ok) pass++;
	else fail++;
	printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
}

// Deterministic uniform [0, 1) and standard normal, so a failure reproduces.
static double unif(uint32_t i) {
	uint32_t x = i * 2654435761u + 1013904223u;
	x ^= x >> 16; x *= 2246822519u; x ^= x >> 13; x *= 3266489917u; x ^= x >> 16;
	return ((double)x + 0.5) / 4294967296.0;
}

static double gauss(uint32_t i) {
	return sqrt(-2.0 * log(unif(2u * i))) * cos(2.0 * M_PI * unif(2u * i + 1u));
}

// A Tucker-rank (r0, r1, r2) signal of the given dims, plus noise of SD sd.
static void make(const int *d, const int *r, double sd, uint32_t seed, double *clean, double *noisy) {
	const size_t n = (size_t)d[0] * d[1] * d[2];
	for (size_t i = 0; i < n; i++) clean[i] = 0.0;
	uint32_t s = seed * 7919u;
	for (int a = 0; a < r[0]; a++)
		for (int b = 0; b < r[1]; b++)
			for (int c = 0; c < r[2]; c++) {
				const double g = 10.0 * gauss(s++) / (1.0 + a + b + c);
				for (int k = 0; k < d[2]; k++) {
					const double ck = cos(0.7 * (c + 1) * k + c);
					for (int j = 0; j < d[1]; j++) {
						const double bj = cos(0.31 * (b + 1) * j + 2.0 * b);
						for (int i = 0; i < d[0]; i++)
							clean[i + (size_t)d[0] * (j + (size_t)d[1] * k)] +=
								g * cos(0.53 * (a + 1) * i + 3.0 * a) * bj * ck;
					}
				}
			}
	for (size_t i = 0; i < n; i++) noisy[i] = clean[i] + 100.0 + sd * gauss(seed * 1000003u + (uint32_t)i);
	for (size_t i = 0; i < n; i++) clean[i] += 100.0;
}

static double rel_l2(const double *a, const double *b, size_t n) {
	double sd = 0.0, sb = 0.0;
	for (size_t i = 0; i < n; i++) {
		sd += (a[i] - b[i]) * (a[i] - b[i]);
		sb += b[i] * b[i];
	}
	return sb > 0.0 ? sqrt(sd / sb) : sqrt(sd);
}

static int run(const dn_tensor_opts *o, double *x, double sigma2, int *ranks) {
	dn_tensor_work *w = dn_tensor_work_create(o);
	if (!w) return 1;
	const int rc = dn_tensor_denoise(o, w, x, sigma2, ranks);
	dn_tensor_work_free(w);
	return rc;
}

static dn_tensor_opts opts(const int *d) {
	dn_tensor_opts o;
	memset(&o, 0, sizeof o);
	o.nmodes = 3;
	for (int k = 0; k < 3; k++) o.dims[k] = d[k];
	o.filter = DN_TFILTER_OPTSHRINK;
	return o;
}

// `dn_tensortest run d0 d1 d2 filter demean sigma2`: denoise one column-major
// float64 tensor from stdin to stdout, printing the kept ranks "r0 r1 r2" on
// stderr.  filter is optshrink, truncate or optthresh.  The kernel-level half
// of tools/tensor_parity.py.
static int run_raw(int argc, char **argv) {
	if (argc != 8) return 2;
	int d[3];
	for (int k = 0; k < 3; k++) d[k] = atoi(argv[2 + k]);
	dn_tensor_opts o = opts(d);
	if (!strcmp(argv[5], "optshrink")) o.filter = DN_TFILTER_OPTSHRINK;
	else if (!strcmp(argv[5], "truncate")) o.filter = DN_TFILTER_TRUNCATE;
	else if (!strcmp(argv[5], "optthresh")) o.filter = DN_TFILTER_OPTTHRESH;
	else return 2;
	o.subtract_mean = atoi(argv[6]);
	const double sigma2 = atof(argv[7]);
	if (dn_tensor_check(&o)) return 2;
	const size_t n = (size_t)d[0] * d[1] * d[2];
	double *x = malloc(n * sizeof(double));
	if (!x || fread(x, sizeof(double), n, stdin) != n) { free(x); return 2; }
	int r[3];
	if (run(&o, x, sigma2, r)) { free(x); return 1; }
	fprintf(stderr, "%d %d %d\n", r[0], r[1], r[2]);
	const int ok = fwrite(x, sizeof(double), n, stdout) == n;
	free(x);
	return ok ? 0 : 2;
}

int main(int argc, char **argv) {
	if (argc > 1 && !strcmp(argv[1], "run")) return run_raw(argc, argv);

	// Tall and wide in every mode: (8,60,4) has a time mode longer than the rest
	// of the tensor, so its Gram is the complementary one, and (10,2,3) has
	// later unfoldings smaller than their mode sizes.
	static const int shapes[][3] = {{3, 7, 5}, {27, 40, 3}, {8, 60, 4}, {10, 2, 3}, {64, 30, 4}};
	static const int rk[][3] = {{2, 2, 2}, {3, 2, 2}, {2, 3, 2}, {2, 1, 2}, {3, 3, 2}};
	const int nshape = (int)(sizeof shapes / sizeof shapes[0]);
	char msg[160];

	printf("tensor kernel\n");
	for (int s = 0; s < nshape; s++) {
		const int *d = shapes[s];
		const size_t n = (size_t)d[0] * d[1] * d[2];
		double *clean = malloc(n * sizeof(double)), *noisy = malloc(n * sizeof(double));
		double *x = malloc(n * sizeof(double));
		if (!clean || !noisy || !x) return 2;
		make(d, rk[s], 1.0, (uint32_t)s + 1u, clean, noisy);
		const double before = rel_l2(noisy, clean, n);
		int r[3];

		// The true variance is 1.  Denoising helps, with and without centring.
		for (int dm = 0; dm < 2; dm++) {
			dn_tensor_opts o = opts(d);
			o.subtract_mean = dm;
			memcpy(x, noisy, n * sizeof(double));
			const int rc = run(&o, x, 1.0, r);
			snprintf(msg, sizeof msg, "(%d,%d,%d)%s closer to the truth at the true variance: "
			         "%.4g -> %.4g, ranks %d,%d,%d", d[0], d[1], d[2], dm ? " demeaned," : "",
			         before, rel_l2(x, clean, n), r[0], r[1], r[2]);
			report(!rc && rel_l2(x, clean, n) < before, msg);
		}

		// sigma2 = 0 keeps everything: the output is the input.
		for (int dm = 0; dm < 2; dm++) {
			dn_tensor_opts o = opts(d);
			o.subtract_mean = dm;
			memcpy(x, noisy, n * sizeof(double));
			const int rc = run(&o, x, 0.0, r);
			const double e0 = rel_l2(x, noisy, n);
			snprintf(msg, sizeof msg, "(%d,%d,%d) sigma2 = 0%s returns the input (rel L2 %.2g)",
			         d[0], d[1], d[2], dm ? ", demeaned," : "", e0);
			report(!rc && e0 < 1e-12, msg);
		}
		free(clean); free(noisy); free(x);
	}

	// optthresh on a noisy low-rank tensor: still closer to the truth, and never
	// keeping more at level 0 than the bulk edge does, since its cut is higher.
	// Later levels see the ranks before them, so only level 0 is comparable.
	for (int s = 0; s < nshape; s++) {
		const int *d = shapes[s];
		const size_t n = (size_t)d[0] * d[1] * d[2];
		double *clean = malloc(n * sizeof(double)), *noisy = malloc(n * sizeof(double));
		double *x = malloc(n * sizeof(double));
		if (!clean || !noisy || !x) return 2;
		make(d, rk[s], 1.0, (uint32_t)s + 1u, clean, noisy);
		const double before = rel_l2(noisy, clean, n);
		int rt[3], ro[3];
		dn_tensor_opts o = opts(d);
		o.filter = DN_TFILTER_TRUNCATE;
		memcpy(x, noisy, n * sizeof(double));
		const int rc_t = run(&o, x, 1.0, rt);
		o.filter = DN_TFILTER_OPTTHRESH;
		memcpy(x, noisy, n * sizeof(double));
		const int rc_o = run(&o, x, 1.0, ro);
		snprintf(msg, sizeof msg, "(%d,%d,%d) optthresh closer to the truth: %.4g -> %.4g, "
		         "ranks %d,%d,%d (edge: %d,%d,%d)", d[0], d[1], d[2], before, rel_l2(x, clean, n),
		         ro[0], ro[1], ro[2], rt[0], rt[1], rt[2]);
		report(!rc_t && !rc_o && rel_l2(x, clean, n) < before && ro[0] <= rt[0], msg);
		free(clean); free(noisy); free(x);
	}

	// The two cuts on a rank-1 tensor of known singular value s, the same at
	// every level.  On (8,6,4) at unit variance, level 0 (8 x 24) has its edge
	// at s^2 = (sqrt 8 + sqrt 24)^2 = 59.7 and lambda_*(1/3)^2 * 24 = 81.1 for
	// optthresh; level 1 (6 x 4) cuts lower, at 19.8 and 26.5.  So s^2 = 70 is
	// kept by the edge and dropped by optthresh, and s^2 = 100 is kept by both,
	// WHOLE under optthresh and truncate and shrunk under optshrink.
	{
		const int d[3] = {8, 6, 4};
		const size_t n = 8 * 6 * 4;
		double x[8 * 6 * 4], in[8 * 6 * 4];
		static const double s2s[2] = {70.0, 100.0};
		for (int t = 0; t < 2; t++) {
			const double amp = sqrt(s2s[t] / (double)n);
			for (size_t i = 0; i < n; i++) in[i] = amp;
			int r[3];
			for (int f = DN_TFILTER_OPTSHRINK; f <= DN_TFILTER_OPTTHRESH; f++) {
				dn_tensor_opts o = opts(d);
				o.filter = f;
				memcpy(x, in, sizeof x);
				const int rc = run(&o, x, 1.0, r);
				const double e = rel_l2(x, in, n);
				const char *name = f == DN_TFILTER_OPTSHRINK ? "optshrink"
				                 : f == DN_TFILTER_TRUNCATE ? "truncate" : "optthresh";
				int ok;
				if (t == 0 && f == DN_TFILTER_OPTTHRESH) ok = r[0] == 0 && e == 1.0;
				else if (f == DN_TFILTER_OPTSHRINK) ok = r[0] == 1 && e > 1e-3 && e < 1.0;
				else ok = r[0] == 1 && e < 1e-12;
				snprintf(msg, sizeof msg, "rank-1, s^2 = %g, %s: level-0 rank %d, rel L2 to input %.3g",
				         s2s[t], name, r[0], e);
				report(!rc && ok, msg);
			}
		}
	}

	// The two shrinkers, exactly.  The rank-1 tensor above at s^2 = 70 reaches
	// the last level as a 4 x 1 unfolding (echo mode by the one kept component)
	// with r = 1, so its output is the input times the weight eta/s:
	//   gd:     eta^2 = s^2 - 2(4 + 1) + (4 - 1)^2 / s^2
	//   olesen: eta^2 = s^2 - 2(3 + 1) + (1 - 3)^2 / s^2   (R = 4-1, C = max(1-1, 1))
	{
		const int d[3] = {8, 6, 4};
		const size_t n = 8 * 6 * 4;
		double x[8 * 6 * 4], in[8 * 6 * 4];
		const double s2 = 70.0, amp = sqrt(s2 / (double)n);
		for (size_t i = 0; i < n; i++) in[i] = amp;
		const double want[2] = {sqrt(s2 - 8.0 + 4.0 / s2) / sqrt(s2), sqrt(s2 - 10.0 + 9.0 / s2) / sqrt(s2)};
		for (int obs = 0; obs < 2; obs++) {
			dn_tensor_opts o = opts(d);
			o.shrink_observed = obs;
			memcpy(x, in, sizeof x);
			int r[3];
			const int rc = run(&o, x, 1.0, r);
			double worst = 0.0;
			for (size_t i = 0; i < n; i++) worst = fmax(worst, fabs(x[i] / in[i] - want[obs]));
			snprintf(msg, sizeof msg, "rank-1 %s shrinker: weight %.12g, expected %.12g (max error %.2g)",
			         obs ? "gd" : "olesen", x[0] / in[0], want[obs], worst);
			report(!rc && worst < 1e-12, msg);
		}
	}

	// All zero: nothing kept, nothing invented, at any variance.
	{
		const int d[3] = {9, 6, 3};
		double x[9 * 6 * 3];
		memset(x, 0, sizeof x);
		dn_tensor_opts o = opts(d);
		int r[3] = {-1, -1, -1};
		const int rc = run(&o, x, 4.0, r);
		double mx = 0.0;
		for (size_t i = 0; i < sizeof x / sizeof x[0]; i++) mx = fmax(mx, fabs(x[i]));
		snprintf(msg, sizeof msg, "all-zero patch: ranks %d,%d,%d, max |out| %g", r[0], r[1], r[2], mx);
		report(!rc && !r[0] && !r[1] && !r[2] && mx == 0.0, msg);
	}

	// Refusals.
	{
		const int d[3] = {9, 6, 3};
		double x[9 * 6 * 3] = {0};
		dn_tensor_opts o = opts(d);
		int r[3];
		report(run(&o, x, -1.0, r) != 0 && run(&o, x, NAN, r) != 0, "a negative or NaN variance is refused");
		const int d1[3] = {9, 1, 3};
		dn_tensor_opts o1 = opts(d1);
		report(dn_tensor_check(&o1) != 0, "a mode of size 1 is refused");
	}

	printf("\n%d passed, %d failed\n", pass, fail);
	return fail ? 1 : 0;
}
