/* A zurand stream tested through zurngtest's C callable. The fill folds the
 * key by a running block counter, so successive blocks are successive
 * substreams of one stream, and writes each 64-bit word as two 32-bit
 * words, low half first. */
#include <stdint.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <zurand.h>
#include <zurngtest.h>

#define BLOCK 2048 /* 64-bit words per folded substream */

typedef struct {
    const zurand_api *zr;
    zurand_key key;
    int64_t block;
} stream;

static void fill(void *state, uint32_t *out, size_t n) {
    stream *s = state;
    uint64_t buf[BLOCK];
    while (n > 0) {
        zurand_key k;
        if (s->zr->fold_int(s->key, s->block++, &k) != ZURAND_OK ||
            s->zr->fill_bits64(k, BLOCK, buf) != ZURAND_OK)
            Rf_error("zurand failed to fill a block");
        size_t words = n < 2 * BLOCK ? n : 2 * BLOCK;
        for (size_t i = 0; i < words; i++)
            out[i] = (uint32_t)(buf[i / 2] >> (32 * (i % 2)));
        out += words;
        n -= words;
    }
}

SEXP pract_zurand_c(SEXP keys, SEXP options) {
    stream s;
    s.zr = zurand_get_api();
    s.block = 0;
    if (s.zr->key_get(keys, 0, &s.key) != ZURAND_OK)
        Rf_error("`keys` must hold at least one key");
    return zurngtest_test_cfun(fill, &s, options);
}

static const R_CallMethodDef calls[] = {
    {"pract_zurand_c", (DL_FUNC) &pract_zurand_c, 2},
    {NULL, NULL, 0}
};

void R_init_zurandpract(DllInfo *dll) {
    R_registerRoutines(dll, NULL, calls, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
