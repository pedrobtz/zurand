#ifdef _OPENMP
#include <omp.h>
#endif
#define R_NO_REMAP
#include <zurand.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

static uint64_t hex64(SEXP s, int i) { return strtoull(CHAR(STRING_ELT(s, i)), NULL, 0); }
static zurand_key getkey(SEXP k) {
    zurand_key z;
    z.k0 = hex64(VECTOR_ELT(k, 0), 0); z.k1 = hex64(VECTOR_ELT(k, 1), 0);
    z.engine = Rf_asInteger(VECTOR_ELT(k, 2)); z.stream = Rf_asInteger(VECTOR_ELT(k, 3));
    return z;
}
static SEXP hexstr(uint64_t v) { char b[32]; snprintf(b, sizeof b, "0x%016" PRIx64, v); return Rf_mkString(b); }
static SEXP keylist(zurand_key z) {
    const char *nm[] = {"k0", "k1", "engine", "stream", ""};
    SEXP a = PROTECT(Rf_mkNamed(VECSXP, nm));
    SET_VECTOR_ELT(a, 0, hexstr(z.k0)); SET_VECTOR_ELT(a, 1, hexstr(z.k1));
    SET_VECTOR_ELT(a, 2, Rf_ScalarInteger(z.engine)); SET_VECTOR_ELT(a, 3, Rf_ScalarInteger(z.stream));
    UNPROTECT(1); return a;
}
SEXP zc_key_get(SEXP keys, SEXP i_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key z = {0xdeadbeef, 0xdeadbeef, -7, -7};
    int rc = zr->key_get(keys, (R_xlen_t)Rf_asReal(i_), &z);
    SEXP a = PROTECT(keylist(z));
    Rf_setAttrib(a, Rf_install("rc"), Rf_ScalarInteger(rc));
    UNPROTECT(1); return a;
}
#define GUARD 8
#define SENT 0x7ff8dead0000beefULL
static SEXP result(int rc, SEXP vals, int guard_ok) {
    const char *nm[] = {"rc", "x", "guard_ok", ""};
    SEXP a = PROTECT(Rf_mkNamed(VECSXP, nm));
    SET_VECTOR_ELT(a, 0, Rf_ScalarInteger(rc)); SET_VECTOR_ELT(a, 1, vals);
    SET_VECTOR_ELT(a, 2, Rf_ScalarLogical(guard_ok));
    UNPROTECT(1); return a;
}
/* what: 0 uniform 1 normal 2 integer 3 bits64 4 normal_method. Untouched cells keep the sentinel. */
SEXP zc_fill(SEXP key, SEXP n_, SEXP what_, SEXP a_, SEXP b_, SEXP method_, SEXP imin_, SEXP imax_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key);
    size_t n = (size_t)Rf_asReal(n_); int what = Rf_asInteger(what_);
    double a = Rf_asReal(a_), b = Rf_asReal(b_);
    int rc, ok = 1;
    if (what == 2) {
        int *buf = malloc((n + 2 * GUARD) * sizeof(int));
        for (size_t i = 0; i < n + 2 * GUARD; i++) buf[i] = -777777777;
        rc = zr->fill_integer(k, n, INTEGER(imin_)[0], INTEGER(imax_)[0], buf + GUARD);
        for (int i = 0; i < GUARD; i++) if (buf[i] != -777777777 || buf[GUARD + n + i] != -777777777) ok = 0;
        SEXP v = PROTECT(Rf_allocVector(INTSXP, n));
        memcpy(INTEGER(v), buf + GUARD, n * sizeof(int)); free(buf);
        SEXP r = result(rc, v, ok); UNPROTECT(1); return r;
    }
    uint64_t *buf = malloc((n + 2 * GUARD) * 8);
    for (size_t i = 0; i < n + 2 * GUARD; i++) buf[i] = SENT;
    double *d = (double *)(buf + GUARD);
    if (what == 0) rc = zr->fill_uniform(k, n, a, b, d);
    else if (what == 1) rc = zr->fill_normal(k, n, a, b, d);
    else if (what == 3) rc = zr->fill_bits64(k, n, buf + GUARD);
    else rc = zr->fill_normal_method(k, n, a, b, Rf_asInteger(method_), d);
    for (int i = 0; i < GUARD; i++) if (buf[i] != SENT || buf[GUARD + n + i] != SENT) ok = 0;
    SEXP v;
    if (what == 3) {
        v = PROTECT(Rf_allocVector(STRSXP, n));
        char h[17];
        for (size_t i = 0; i < n; i++) { snprintf(h, sizeof h, "%016" PRIx64, buf[GUARD + i]); SET_STRING_ELT(v, i, Rf_mkChar(h)); }
    } else {
        v = PROTECT(Rf_allocVector(REALSXP, n));
        memcpy(REAL(v), d, n * 8);
    }
    free(buf);
    SEXP r = result(rc, v, ok); UNPROTECT(1); return r;
}
SEXP zc_fill_at(SEXP key, SEXP start_, SEXP n_, SEXP dist_, SEXP a_, SEXP b_, SEXP method_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key);
    uint64_t start = hex64(start_, 0);
    size_t n = (size_t)Rf_asReal(n_); int ok = 1, rc;
    uint64_t *buf = malloc((n + 2 * GUARD) * 8);
    for (size_t i = 0; i < n + 2 * GUARD; i++) buf[i] = SENT;
    double *d = (double *)(buf + GUARD);
    if (Rf_asInteger(dist_) == 0) rc = zr->fill_uniform_at(k, start, n, Rf_asReal(a_), Rf_asReal(b_), d);
    else rc = zr->fill_normal_at(k, start, n, Rf_asReal(a_), Rf_asReal(b_), Rf_asInteger(method_), d);
    for (int i = 0; i < GUARD; i++) if (buf[i] != SENT || buf[GUARD + n + i] != SENT) ok = 0;
    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(v), d, n * 8); free(buf);
    SEXP r = result(rc, v, ok); UNPROTECT(1); return r;
}
/* consecutive pieces from worker threads */
SEXP zc_par(SEXP key, SEXP start_, SEXP n_, SEXP pieces_, SEXP dist_, SEXP method_, SEXP threads_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key);
    uint64_t start = hex64(start_, 0);
    size_t n = (size_t)Rf_asReal(n_); int np = Rf_asInteger(pieces_), dist = Rf_asInteger(dist_), method = Rf_asInteger(method_);
    int nt = Rf_asInteger(threads_); (void)nt;
    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
    double *d = REAL(v); int bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic) reduction(+:bad)
#endif
    for (int p = 0; p < np; p++) {
        size_t lo = n / np * p, hi = p == np - 1 ? n : n / np * (p + 1);
        int rc = dist == 0 ? zr->fill_uniform_at(k, start + lo, hi - lo, 0.0, 1.0, d + lo)
                           : zr->fill_normal_at(k, start + lo, hi - lo, 0.0, 1.0, method, d + lo);
        bad += rc != ZURAND_OK;
    }
    SEXP r = result(bad, v, 1); UNPROTECT(1); return r;
}
typedef struct { double *keep; size_t nkeep, kept; int calls, stop_after; double *cn, *cs; int maxcalls; const double *buf; int bufbad; uint64_t expect; int seqbad; } sctx;
static int consume(void *c_, const double *x, size_t n, uint64_t start) {
    sctx *c = c_;
    if (x != c->buf) c->bufbad = 1;
    if (start != c->expect) c->seqbad = 1;
    c->expect = start + n;
    if (c->calls < c->maxcalls) { c->cn[c->calls] = (double)n; c->cs[c->calls] = (double)start; }
    c->calls++;
    for (size_t i = 0; i < n && c->kept < c->nkeep; i++) c->keep[c->kept++] = x[i];
    return c->stop_after > 0 && c->calls >= c->stop_after;
}
SEXP zc_stream(SEXP key, SEXP total_, SEXP chunk_, SEXP dist_, SEXP a_, SEXP b_, SEXP method_, SEXP stop_, SEXP keep_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key);
    uint64_t total = hex64(total_, 0);
    size_t chunk = (size_t)Rf_asReal(chunk_);
    size_t nkeep = (size_t)Rf_asReal(keep_);
    sctx c; memset(&c, 0, sizeof c);
    c.keep = malloc((nkeep + 1) * 8); c.nkeep = nkeep; c.stop_after = Rf_asInteger(stop_);
    c.maxcalls = 4096; c.cn = calloc(4096, 8); c.cs = calloc(4096, 8);
    uint64_t *buf = malloc((chunk + 2 * GUARD + 1) * 8);
    for (size_t i = 0; i < chunk + 2 * GUARD; i++) buf[i] = SENT;
    c.buf = (double *)(buf + GUARD);
    int rc = Rf_asInteger(dist_) == 0
        ? zr->stream_uniform(k, total, Rf_asReal(a_), Rf_asReal(b_), (double *)(buf + GUARD), chunk, consume, &c)
        : zr->stream_normal(k, total, Rf_asReal(a_), Rf_asReal(b_), Rf_asInteger(method_), (double *)(buf + GUARD), chunk, consume, &c);
    int ok = 1;
    for (int i = 0; i < GUARD; i++) if (buf[i] != SENT || buf[GUARD + chunk + i] != SENT) ok = 0;
    const char *nm[] = {"rc", "x", "guard_ok", "calls", "n", "start", "bufbad", "seqbad", "end", ""};
    SEXP a = PROTECT(Rf_mkNamed(VECSXP, nm));
    SET_VECTOR_ELT(a, 0, Rf_ScalarInteger(rc));
    SEXP v = Rf_allocVector(REALSXP, c.kept); SET_VECTOR_ELT(a, 1, v); memcpy(REAL(v), c.keep, c.kept * 8);
    SET_VECTOR_ELT(a, 2, Rf_ScalarLogical(ok));
    SET_VECTOR_ELT(a, 3, Rf_ScalarInteger(c.calls));
    int nc = c.calls < c.maxcalls ? c.calls : c.maxcalls;
    SEXP cn = Rf_allocVector(REALSXP, nc); SET_VECTOR_ELT(a, 4, cn); memcpy(REAL(cn), c.cn, nc * 8);
    SEXP cs = Rf_allocVector(REALSXP, nc); SET_VECTOR_ELT(a, 5, cs); memcpy(REAL(cs), c.cs, nc * 8);
    SET_VECTOR_ELT(a, 6, Rf_ScalarInteger(c.bufbad)); SET_VECTOR_ELT(a, 7, Rf_ScalarInteger(c.seqbad));
    SET_VECTOR_ELT(a, 8, hexstr(c.expect));
    free(c.keep); free(c.cn); free(c.cs); free(buf);
    UNPROTECT(1); return a;
}
SEXP zc_fold(SEXP key, SEXP data_) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key), o = {0xdeadbeef, 0xdeadbeef, -7, -7};
    int64_t data = (int64_t)strtoll(CHAR(STRING_ELT(data_, 0)), NULL, 0);
    int rc = zr->fold_int(k, data, &o);
    SEXP a = PROTECT(keylist(o));
    Rf_setAttrib(a, Rf_install("rc"), Rf_ScalarInteger(rc));
    UNPROTECT(1); return a;
}
static int never(void *c, const double *x, size_t n, uint64_t s) { (void)c; (void)x; (void)n; (void)s; return 0; }
SEXP zc_null(SEXP key) {
    const zurand_api *zr = zurand_get_api();
    zurand_key k = getkey(key);
    double d[4]; int iv[4]; uint64_t u[4];
    const char *nm[] = {"u_null_n0", "n_null_n0", "i_null_n0", "b_null_n0", "m_null_n0", "uat_null_n0", "nat_null_n0",
        "u_null_n1", "n_null_n1", "i_null_n1", "b_null_n1", "uat_null_n1", "nat_null_n1",
        "fold_null", "s_buf_null", "s_fn_null", "s_chunk0", "s_total0", "s_total0_bufnull", "s_chunk_huge", "u_n_huge", "uat_n_huge", "b_n_huge", "i_n_huge",
        "u_null_n0_const", "n_null_n0_sd0", "uat_null_n0_const", "s_const_total0", ""};
    SEXP a = PROTECT(Rf_mkNamed(INTSXP, nm)); int *r = INTEGER(a), j = 0;
    r[j++] = zr->fill_uniform(k, 0, 0, 1, NULL);
    r[j++] = zr->fill_normal(k, 0, 0, 1, NULL);
    r[j++] = zr->fill_integer(k, 0, 1, 6, NULL);
    r[j++] = zr->fill_bits64(k, 0, NULL);
    r[j++] = zr->fill_normal_method(k, 0, 0, 1, 1, NULL);
    r[j++] = zr->fill_uniform_at(k, 12345, 0, 0, 1, NULL);
    r[j++] = zr->fill_normal_at(k, 12345, 0, 0, 1, 0, NULL);
    r[j++] = zr->fill_uniform(k, 1, 0, 1, NULL);
    r[j++] = zr->fill_normal(k, 1, 0, 1, NULL);
    r[j++] = zr->fill_integer(k, 1, 1, 6, NULL);
    r[j++] = zr->fill_bits64(k, 1, NULL);
    r[j++] = zr->fill_uniform_at(k, 5, 1, 0, 1, NULL);
    r[j++] = zr->fill_normal_at(k, 5, 1, 0, 1, 0, NULL);
    r[j++] = zr->fold_int(k, 1, NULL);
    r[j++] = zr->stream_uniform(k, 10, 0, 1, NULL, 4, never, NULL);
    r[j++] = zr->stream_uniform(k, 10, 0, 1, d, 4, NULL, NULL);
    r[j++] = zr->stream_uniform(k, 10, 0, 1, d, 0, never, NULL);
    r[j++] = zr->stream_uniform(k, 0, 0, 1, d, 4, never, NULL);
    r[j++] = zr->stream_uniform(k, 0, 0, 1, NULL, 4, never, NULL);
    r[j++] = zr->stream_uniform(k, 10, 0, 1, d, (size_t)-1, never, NULL);
    r[j++] = zr->fill_uniform(k, (size_t)-1, 0, 1, d);
    r[j++] = zr->fill_uniform_at(k, 0, (size_t)-1, 0, 1, d);
    r[j++] = zr->fill_bits64(k, (size_t)-1, u);
    r[j++] = zr->fill_integer(k, (size_t)-1, 1, 6, iv);
    r[j++] = zr->fill_uniform(k, 0, 3, 3, NULL);
    r[j++] = zr->fill_normal(k, 0, 3, 0, NULL);
    r[j++] = zr->fill_uniform_at(k, 7, 0, 3, 3, NULL);
    r[j++] = zr->stream_uniform(k, 0, 3, 3, d, 4, never, NULL);
    UNPROTECT(1); return a;
}
static const R_CallMethodDef entries[] = {
    {"C_zc_key_get", (DL_FUNC)&zc_key_get, 2}, {"C_zc_fill", (DL_FUNC)&zc_fill, 8},
    {"C_zc_fill_at", (DL_FUNC)&zc_fill_at, 7}, {"C_zc_stream", (DL_FUNC)&zc_stream, 9},
    {"C_zc_fold", (DL_FUNC)&zc_fold, 2}, {"C_zc_null", (DL_FUNC)&zc_null, 1},
    {"C_zc_par", (DL_FUNC)&zc_par, 7}, {NULL, NULL, 0}};
void R_init_zc3(DllInfo *dll) { R_registerRoutines(dll, NULL, entries, NULL, NULL); R_useDynamicSymbols(dll, FALSE); }
