/* A downstream package using the zurand C API. The fills run inside an
 * OpenMP parallel loop over keys when OpenMP is available, because calling
 * them from worker threads is the point of the API. */
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <zurand.h>

#ifdef __has_include
#  if __has_include(<omp.h>)
#    define CLIENT_HAVE_OMP_H 1
#  endif
#else
#  define CLIENT_HAVE_OMP_H 1
#endif
#if defined(_OPENMP) && defined(CLIENT_HAVE_OMP_H)
#  define CLIENT_OPENMP 1
#endif

static void key_hex(zurand_key k, char *buf, size_t len) {
    snprintf(buf, len, "rng_key[%08x%08x%08x%08x]",
             (unsigned)(uint32_t)k.k0, (unsigned)(uint32_t)(k.k0 >> 32),
             (unsigned)(uint32_t)k.k1, (unsigned)(uint32_t)(k.k1 >> 32));
}

/* One column per key; returns an n x K matrix (a character one for bits64). */
SEXP client_fill(SEXP keys, SEXP n_, SEXP what_, SEXP a_, SEXP b_) {
    const zurand_api *zr = zurand_get_api();
    int n = Rf_asInteger(n_), what = Rf_asInteger(what_);
    double a = Rf_asReal(a_), b = Rf_asReal(b_);
    int nkey = Rf_nrows(keys);
    zurand_key *ks = (zurand_key *) R_alloc((size_t)nkey, sizeof(zurand_key));
    for (int j = 0; j < nkey; j++)
        if (zr->key_get(keys, j, &ks[j]) != ZURAND_OK)
            Rf_error("key_get failed");

    size_t total = (size_t)n * (size_t)nkey;
    SEXP ans;
    double *dout = NULL; int *iout = NULL; uint64_t *bout = NULL;
    if (what == 2) {
        ans = PROTECT(Rf_allocMatrix(INTSXP, n, nkey)); iout = INTEGER(ans);
    } else if (what == 3) {
        ans = PROTECT(Rf_allocMatrix(STRSXP, n, nkey));
        bout = (uint64_t *) R_alloc(total > 0 ? total : 1, sizeof(uint64_t));
    } else {
        ans = PROTECT(Rf_allocMatrix(REALSXP, n, nkey)); dout = REAL(ans);
    }

    int failed = 0;
#ifdef CLIENT_OPENMP
#pragma omp parallel for reduction(|:failed) schedule(static)
#endif
    for (int j = 0; j < nkey; j++) {
        size_t off = (size_t)j * (size_t)n;
        int rc;
        if (what == 0)      rc = zr->fill_uniform(ks[j], (size_t)n, a, b, dout + off);
        else if (what == 1) rc = zr->fill_normal(ks[j], (size_t)n, a, b, dout + off);
        else if (what == 2) rc = zr->fill_integer(ks[j], (size_t)n, (int)a, (int)b, iout + off);
        else if (what == 3) rc = zr->fill_bits64(ks[j], (size_t)n, bout + off);
        else                rc = zr->fill_normal_method(ks[j], (size_t)n, a, b,
                                  what == 4 ? ZURAND_NORMAL_ZIGGURAT
                                            : ZURAND_NORMAL_MCFARLAND, dout + off);
        failed |= (rc != ZURAND_OK);
    }
    if (failed)
        Rf_error("a fill returned an error code");

    if (what == 3) {
        char hex[17];
        for (size_t i = 0; i < total; i++) {
            snprintf(hex, sizeof hex, "%016llx", (unsigned long long)bout[i]);
            SET_STRING_ELT(ans, (R_xlen_t)i, Rf_mkChar(hex));
        }
    }
    UNPROTECT(1);
    return ans;
}

/* Positional fills (API version 3), one column per key, from worker
 * threads: what 0 = uniform on [a, b], 1 = normal(a, b) with `method`. */
SEXP client_fill_at(SEXP keys, SEXP n_, SEXP what_, SEXP a_, SEXP b_,
                    SEXP offset_, SEXP method_) {
    const zurand_api *zr = zurand_get_api();
    int n = Rf_asInteger(n_), what = Rf_asInteger(what_);
    int method = Rf_asInteger(method_);
    double a = Rf_asReal(a_), b = Rf_asReal(b_);
    uint64_t offset = (uint64_t)Rf_asReal(offset_);
    int nkey = Rf_nrows(keys);
    zurand_key *ks = (zurand_key *) R_alloc((size_t)nkey, sizeof(zurand_key));
    for (int j = 0; j < nkey; j++)
        if (zr->key_get(keys, j, &ks[j]) != ZURAND_OK)
            Rf_error("key_get failed");
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, nkey));
    double *out = REAL(ans);
    int failed = 0;
#ifdef CLIENT_OPENMP
#pragma omp parallel for reduction(|:failed) schedule(static)
#endif
    for (int j = 0; j < nkey; j++) {
        double *o = out + (size_t)j * (size_t)n;
        int rc = what == 0
            ? zr->fill_uniform_at(ks[j], offset, (size_t)n, a, b, o)
            : zr->fill_normal_at(ks[j], offset, (size_t)n, a, b, method, o);
        failed |= (rc != ZURAND_OK);
    }
    if (failed)
        Rf_error("a positional fill returned an error code");
    UNPROTECT(1);
    return ans;
}

/* A consumer that copies each chunk to out[start ..] and records the
 * chunk's start and length; it stops the stream after `stop_after` chunks
 * when that is positive. */
typedef struct {
    double *out;
    double *starts, *lens;
    int nchunk, max_chunks, stop_after;
    int bad;   /* a chunk arrived out of order */
    uint64_t next;
} collect_ctx;

static int collect(void *ctx_, const double *x, size_t n, uint64_t start) {
    collect_ctx *ctx = ctx_;
    if (start != ctx->next) ctx->bad = 1;
    ctx->next = start + n;
    memcpy(ctx->out + start, x, n * sizeof(double));
    if (ctx->nchunk < ctx->max_chunks) {
        ctx->starts[ctx->nchunk] = (double)start;
        ctx->lens[ctx->nchunk] = (double)n;
    }
    ctx->nchunk++;
    return ctx->stop_after > 0 && ctx->nchunk >= ctx->stop_after;
}

static int client_stream_one(const zurand_api *zr, zurand_key k, int what,
                             uint64_t total, double a, double b, int method,
                             double *buf, size_t chunk, collect_ctx *ctx) {
    return what == 0
        ? zr->stream_uniform(k, total, a, b, buf, chunk, collect, ctx)
        : zr->stream_normal(k, total, a, b, method, buf, chunk, collect, ctx);
}

/* One key streamed on the main thread: list(values, starts, lengths,
 * return code, out-of-order flag). Values past a stop stay NA. */
SEXP client_stream(SEXP key, SEXP total_, SEXP chunk_, SEXP what_, SEXP a_,
                   SEXP b_, SEXP method_, SEXP stop_after_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k;
    if (zr->key_get(key, 0, &k) != ZURAND_OK) Rf_error("key_get failed");
    int total = Rf_asInteger(total_), chunk = Rf_asInteger(chunk_);
    int max_chunks = total / chunk + 2;
    SEXP vals = PROTECT(Rf_allocVector(REALSXP, total));
    SEXP starts = PROTECT(Rf_allocVector(REALSXP, max_chunks));
    SEXP lens = PROTECT(Rf_allocVector(REALSXP, max_chunks));
    for (int i = 0; i < total; i++) REAL(vals)[i] = NA_REAL;
    double *buf = (double *) R_alloc((size_t)chunk, sizeof(double));
    collect_ctx ctx = {REAL(vals), REAL(starts), REAL(lens), 0, max_chunks,
                       Rf_asInteger(stop_after_), 0, 0};
    int rc = client_stream_one(zr, k, Rf_asInteger(what_), (uint64_t)total,
                               Rf_asReal(a_), Rf_asReal(b_),
                               Rf_asInteger(method_), buf, (size_t)chunk, &ctx);
    int nc = ctx.nchunk < max_chunks ? ctx.nchunk : max_chunks;
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, vals);
    SET_VECTOR_ELT(ans, 1, Rf_lengthgets(starts, nc));
    SET_VECTOR_ELT(ans, 2, Rf_lengthgets(lens, nc));
    SET_VECTOR_ELT(ans, 3, Rf_ScalarInteger(rc));
    SET_VECTOR_ELT(ans, 4, Rf_ScalarLogical(ctx.bad));
    UNPROTECT(4);
    return ans;
}

/* One stream per key, one key per worker thread, each with its own buffer:
 * the parallel pattern the API is for. Returns a total x K matrix. */
SEXP client_stream_par(SEXP keys, SEXP total_, SEXP chunk_, SEXP what_,
                       SEXP a_, SEXP b_, SEXP method_) {
    const zurand_api *zr = zurand_get_api();
    int total = Rf_asInteger(total_), chunk = Rf_asInteger(chunk_);
    int what = Rf_asInteger(what_), method = Rf_asInteger(method_);
    double a = Rf_asReal(a_), b = Rf_asReal(b_);
    int nkey = Rf_nrows(keys);
    zurand_key *ks = (zurand_key *) R_alloc((size_t)nkey, sizeof(zurand_key));
    for (int j = 0; j < nkey; j++)
        if (zr->key_get(keys, j, &ks[j]) != ZURAND_OK)
            Rf_error("key_get failed");
    double *bufs = (double *) R_alloc((size_t)nkey * (size_t)chunk, sizeof(double));
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, total, nkey));
    double *out = REAL(ans);
    int failed = 0;
#ifdef CLIENT_OPENMP
#pragma omp parallel for reduction(|:failed) schedule(dynamic)
#endif
    for (int j = 0; j < nkey; j++) {
        collect_ctx ctx = {out + (size_t)j * (size_t)total, NULL, NULL, 0, 0,
                           0, 0, 0};
        int rc = client_stream_one(zr, ks[j], what, (uint64_t)total, a, b,
                                   method, bufs + (size_t)j * (size_t)chunk,
                                   (size_t)chunk, &ctx);
        failed |= (rc != ZURAND_OK) | ctx.bad;
    }
    if (failed)
        Rf_error("a stream failed");
    UNPROTECT(1);
    return ans;
}

static int never(void *ctx, const double *x, size_t n, uint64_t start) {
    (void)ctx; (void)x; (void)n; (void)start;
    return 0;
}

/* format(rng_fold(key, data)) for one key and one whole number. */
SEXP client_fold(SEXP key, SEXP data_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k, out;
    if (zr->key_get(key, 0, &k) != ZURAND_OK) Rf_error("key_get failed");
    if (zr->fold_int(k, (int64_t)Rf_asReal(data_), &out) != ZURAND_OK)
        Rf_error("fold_int failed");
    char buf[48];
    key_hex(out, buf, sizeof buf);
    return Rf_mkString(buf);
}

/* Return codes for invalid calls, and the table's version. */
SEXP client_errors(SEXP key) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k, bad;
    double d[4]; int iv[4];
    if (zr->key_get(key, 0, &k) != ZURAND_OK) Rf_error("key_get failed");
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, 21));
    int *r = INTEGER(ans);
    r[0] = zr->version;
    r[1] = zr->key_get(key, 5, &bad);                         /* out of range */
    r[2] = zr->fill_uniform(k, 4, 1.0, 0.0, d);               /* min > max */
    r[3] = zr->fill_normal(k, 4, 0.0, -1.0, d);               /* sd < 0 */
    r[4] = zr->fill_uniform(k, 4, 0.0, R_PosInf, d);          /* infinite */
    r[5] = zr->fill_integer(k, 4, NA_INTEGER, 3, iv);         /* NA bound */
    bad = k; bad.engine = 7;
    r[6] = zr->fill_normal(bad, 4, 0.0, 1.0, d);              /* engine */
    bad = k; bad.stream = 2;
    r[7] = zr->fill_bits64(bad, 4, (uint64_t *) d);           /* stream */
    r[8] = zr->fill_normal(k, 0, 0.0, 1.0, d);                /* n = 0 is fine */
    r[9] = zr->fill_normal_method(k, 4, 0.0, 1.0, 2, d);      /* method */
    r[10] = zr->fill_normal_method(k, 4, 0.0, 1.0, -1, d);    /* method */
    r[11] = zr->fill_normal_method(k, 4, 0.0, -1.0, ZURAND_NORMAL_MCFARLAND, d);
    bad = k; bad.engine = 7;
    r[12] = zr->fill_normal_method(bad, 4, 0.0, 1.0, ZURAND_NORMAL_MCFARLAND, d);
    /* version 3 */
    r[13] = zr->fill_uniform_at(k, 10, 4, 1.0, 0.0, d);                  /* min > max */
    r[14] = zr->fill_normal_at(k, 10, 4, 0.0, 1.0, 2, d);                /* method */
    r[15] = zr->fill_normal_at(k, UINT64_MAX - 2, 4, 0.0, 1.0, 0, d);    /* overflow */
    r[16] = zr->stream_normal(k, 10, 0.0, 1.0, 0, d, 0, never, NULL);    /* chunk 0 */
    r[17] = zr->stream_normal(k, 10, 0.0, 1.0, 0, NULL, 4, never, NULL); /* no buf */
    r[18] = zr->stream_uniform(k, 10, 0.0, 1.0, d, 4, NULL, NULL);       /* no fn */
    r[19] = zr->stream_uniform(bad, 10, 0.0, 1.0, d, 4, never, NULL);    /* engine */
    r[20] = zr->stream_normal(k, 0, 0.0, 1.0, 1, d, 4, never, NULL);     /* total 0 */
    UNPROTECT(1);
    return ans;
}

static const R_CallMethodDef entries[] = {
    {"C_client_fill",   (DL_FUNC) &client_fill,   5},
    {"C_client_fold",   (DL_FUNC) &client_fold,   2},
    {"C_client_errors", (DL_FUNC) &client_errors, 1},
    {"C_client_fill_at", (DL_FUNC) &client_fill_at, 7},
    {"C_client_stream", (DL_FUNC) &client_stream, 8},
    {"C_client_stream_par", (DL_FUNC) &client_stream_par, 7},
    {NULL, NULL, 0}
};

void R_init_zurandclient(DllInfo *dll) {
    R_registerRoutines(dll, NULL, entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
