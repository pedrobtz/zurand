/* Streaming against a big buffer, through the C API (version 3).
 *
 * Each of K keys supplies `total` uniforms or normals, consumed by a sum of
 * squares, one key per OpenMP thread. Three ways to produce them:
 *   big      fill_normal() into a whole-stream buffer per thread, reused
 *            across keys (pages already mapped), then consume it;
 *   fresh    the same into a buffer malloc()ed per key, as R does;
 *   stream   stream_normal() through a `chunk`-value buffer per thread.
 * The values and the order of the sum are the same, so all three return
 * identical sums; only where the values live differs.
 *
 * Built by bench_stream.R with R CMD SHLIB against inst/include. */
#define R_NO_REMAP   /* Rinternals.h's `match` macro breaks LLVM's omp.h */
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <zurand.h>

/* Sum of squares in four accumulators, lane = position mod 4, so the
 * additions are four independent chains (one chain is ~1 ns per value and
 * would hide everything else) and the result does not depend on how the
 * stream is cut: chunks here are multiples of 4. */
typedef struct { double acc[4]; } sumsq;

static int consume(void *ctx, const double *x, size_t n, uint64_t start) {
    (void)start;
    sumsq *s = ctx;
    double a0 = s->acc[0], a1 = s->acc[1], a2 = s->acc[2], a3 = s->acc[3];
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        a0 += x[i] * x[i];
        a1 += x[i + 1] * x[i + 1];
        a2 += x[i + 2] * x[i + 2];
        a3 += x[i + 3] * x[i + 3];
    }
    for (; i < n; i++)
        a0 += x[i] * x[i];
    s->acc[0] = a0; s->acc[1] = a1; s->acc[2] = a2; s->acc[3] = a3;
    return 0;
}

static int thread_id(void) {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* mode 0 = big, 1 = fresh, 2 = stream; dist 0 = uniform, 1 = normal */
SEXP bs_run(SEXP keys, SEXP total_, SEXP threads_, SEXP mode_, SEXP chunk_,
            SEXP dist_) {
    const zurand_api *zr = zurand_get_api();
    int dist = Rf_asInteger(dist_);
    size_t total = (size_t)Rf_asReal(total_);
    int nt = Rf_asInteger(threads_), mode = Rf_asInteger(mode_);
    size_t chunk = (size_t)Rf_asInteger(chunk_);
    int nkey = Rf_nrows(keys);
    zurand_key *ks = (zurand_key *) R_alloc((size_t)nkey, sizeof(zurand_key));
    for (int j = 0; j < nkey; j++)
        zr->key_get(keys, j, &ks[j]);
    size_t per = mode == 2 ? chunk : total;
    double **bufs = (double **) R_alloc((size_t)nt, sizeof(double *));
    for (int t = 0; t < nt; t++) {
        bufs[t] = mode == 1 ? NULL : malloc(per * sizeof(double));
        if (mode != 1 && bufs[t] == NULL) Rf_error("out of memory");
    }
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, nkey));
    double *out = REAL(ans);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic)
#endif
    for (int j = 0; j < nkey; j++) {
        sumsq s = {{0.0, 0.0, 0.0, 0.0}};
        double *b = mode == 1 ? malloc(total * sizeof(double)) : bufs[thread_id()];
        if (mode == 2) {
            if (dist == 0)
                zr->stream_uniform(ks[j], total, 0.0, 1.0, b, chunk, consume, &s);
            else
                zr->stream_normal(ks[j], total, 0.0, 1.0, ZURAND_NORMAL_ZIGGURAT,
                                  b, chunk, consume, &s);
        } else if (b != NULL) {
            if (dist == 0)
                zr->fill_uniform(ks[j], total, 0.0, 1.0, b);
            else
                zr->fill_normal(ks[j], total, 0.0, 1.0, b);
            consume(&s, b, total, 0);
        }
        if (mode == 1) free(b);
        out[j] = (s.acc[0] + s.acc[1]) + (s.acc[2] + s.acc[3]);
    }
    for (int t = 0; t < nt; t++)
        free(bufs[t]);
    UNPROTECT(1);
    return ans;
}
