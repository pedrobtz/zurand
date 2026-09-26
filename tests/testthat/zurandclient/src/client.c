/* A downstream package using the zurand C API. The fills run inside an
 * OpenMP parallel loop over keys when OpenMP is available, because calling
 * them from worker threads is the point of the API. */
#include <stdio.h>
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
        else                rc = zr->fill_bits64(ks[j], (size_t)n, bout + off);
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
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, 9));
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
    UNPROTECT(1);
    return ans;
}

static const R_CallMethodDef entries[] = {
    {"C_client_fill",   (DL_FUNC) &client_fill,   5},
    {"C_client_fold",   (DL_FUNC) &client_fold,   2},
    {"C_client_errors", (DL_FUNC) &client_errors, 1},
    {NULL, NULL, 0}
};

void R_init_zurandclient(DllInfo *dll) {
    R_registerRoutines(dll, NULL, entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
