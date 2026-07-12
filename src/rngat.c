#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>

#include "Random123/philox.h"

/*
 * rngat: random access ("at") random numbers built on Philox4x64-10.
 *
 * Every value is a pure function of (key, index, domain). One Philox
 * evaluation yields a 256-bit block (four uint64 words), and we pack
 * four consecutive logical positions into it:
 *
 *   ctr    = {index >> 2, domain, purpose, 0}   (block = index / 4)
 *   block  = philox4x64(ctr, key)               (128-bit key, 256-bit out)
 *   value  = block.v[index & 3]                 (word = index % 4)
 *
 * So value(index) is still O(1) random access for any index, but a
 * contiguous run of indices costs ~1 Philox call per 4 values instead
 * of 1 per value. Because Philox output words are statistically
 * independent, using all four words as consecutive stream values is
 * the standard Random123 usage and passes the same test batteries.
 *
 * The `purpose` counter word gives domain separation between ordinary
 * draws (PURPOSE_DRAW) and key derivation via fold_in (PURPOSE_FOLD),
 * so a key derived with fold_in(key, i) shares no bits with any value
 * drawn at index i (draw and fold blocks never coincide).
 *
 * Keys are 16-byte raw vectors (class "rngat_key") holding the two
 * 64-bit key words in little-endian order, so keys are portable
 * across architectures.
 */

#define RNGAT_PURPOSE_DRAW ((uint64_t)0)
#define RNGAT_PURPOSE_FOLD ((uint64_t)1)

#define RNGAT_KEY_BYTES 16

/* ---- key (de)serialization, explicitly little-endian ---- */

static uint64_t load_le64(const Rbyte *p) {
    uint64_t x = 0;
    for (int i = 7; i >= 0; i--) x = (x << 8) | (uint64_t)p[i];
    return x;
}

static void store_le64(Rbyte *p, uint64_t x) {
    for (int i = 0; i < 8; i++) { p[i] = (Rbyte)(x & 0xff); x >>= 8; }
}

static philox4x64_key_t key_from_sexp(SEXP key) {
    if (TYPEOF(key) != RAWSXP || Rf_xlength(key) != RNGAT_KEY_BYTES ||
        !Rf_inherits(key, "rngat_key"))
        Rf_error("`key` must be an <rngat_key>; create one with rng_key()");
    philox4x64_key_t k;
    const Rbyte *p = RAW(key);
    k.v[0] = load_le64(p);
    k.v[1] = load_le64(p + 8);
    return k;
}

static SEXP key_to_sexp(uint64_t k0, uint64_t k1) {
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, RNGAT_KEY_BYTES));
    store_le64(RAW(ans), k0);
    store_le64(RAW(ans) + 8, k1);
    Rf_classgets(ans, Rf_mkString("rngat_key"));
    UNPROTECT(1);
    return ans;
}

/* ---- scalar conversions ---- */

static uint64_t u64_from_double(double x, const char *what) {
    if (!R_FINITE(x) || x != trunc(x))
        Rf_error("`%s` must be a whole number, not %g", what, x);
    if (x < -9223372036854775808.0 || x >= 9223372036854775808.0)
        Rf_error("`%s` is out of the representable integer range", what);
    return (uint64_t)(int64_t)x;
}

static void check_numeric(SEXP x, const char *what) {
    if (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)
        Rf_error("`%s` must be an integer or double vector", what);
}

/*
 * A numeric vector opened for fast element access: the integer/double
 * type is resolved and the data pointer fetched once, so the hot loop
 * pays neither a per-element TYPEOF branch nor an accessor call.
 */
typedef struct {
    int is_int;
    const int *ip;
    const double *dp;
    R_xlen_t n;
} num_vec;

static num_vec num_vec_open(SEXP x, const char *what) {
    check_numeric(x, what);
    num_vec v;
    v.n = Rf_xlength(x);
    v.is_int = (TYPEOF(x) == INTSXP);
    v.ip = v.is_int ? INTEGER(x) : NULL;
    v.dp = v.is_int ? NULL : REAL(x);
    return v;
}

static uint64_t num_vec_get(num_vec v, R_xlen_t i, const char *what) {
    if (v.is_int) {
        int e = v.ip[i];
        if (e == NA_INTEGER) Rf_error("`%s` must not contain missing values", what);
        return (uint64_t)(int64_t)e;
    }
    return u64_from_double(v.dp[i], what);
}

static uint64_t u64_scalar(SEXP x, const char *what) {
    num_vec v = num_vec_open(x, what);
    if (v.n != 1)
        Rf_error("`%s` must be a single value", what);
    return num_vec_get(v, 0, what);
}

/* ---- core generator ---- */

/*
 * Force-inlined with the same idiom philox.h uses for its own rounds:
 * a plain `static inline` is only a hint, and clang declines it here
 * (philox4x64's forced inlining makes this function large), leaving a
 * call per block whose 32-byte struct return goes through memory on
 * arm64 -- measured at ~60% overhead on the bulk path. Inlined into the
 * fill loops, the constant purpose and zero fourth word also let the
 * compiler fold part of the first Philox round.
 */
R123_STATIC_INLINE R123_FORCE_INLINE(philox4x64_ctr_t rngat_block(
    philox4x64_key_t key, uint64_t index, uint64_t domain, uint64_t purpose));
R123_STATIC_INLINE philox4x64_ctr_t rngat_block(philox4x64_key_t key,
                                                uint64_t index,
                                                uint64_t domain,
                                                uint64_t purpose) {
    philox4x64_ctr_t ctr = {{index, domain, purpose, 0}};
    return philox4x64(ctr, key);
}

/* 53 random bits -> double strictly inside (0, 1) */
static double u01_open(uint64_t bits) {
    return ((double)(bits >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/* ---- .Call entry points ---- */

SEXP C_rng_key(SEXP seed) {
    uint64_t s = u64_scalar(seed, "seed");
    return key_to_sexp(s, 0);
}

SEXP C_fold_in(SEXP key, SEXP identity) {
    philox4x64_key_t k = key_from_sexp(key);
    uint64_t id = u64_scalar(identity, "identity");
    philox4x64_ctr_t out = rngat_block(k, id, 0, RNGAT_PURPOSE_FOLD);
    return key_to_sexp(out.v[0], out.v[1]);
}

typedef enum { DRAW_BITS, DRAW_UNIF, DRAW_NORM } draw_kind;

/*
 * Gather loop for arbitrary index/domain vectors, stamped once per draw
 * kind (like fill_seq_*) so no per-value dispatch survives. The block
 * cache makes contiguous stretches cost ~1 Philox call per 4 values;
 * correctness never depends on it -- value(index) is a pure function of
 * (index & 3) and the block at (index >> 2, domain) -- it only skips
 * recomputation.
 */
#define RNGAT_DEFINE_FILL_AT(SUFFIX, EMIT)                                   \
static void fill_at_##SUFFIX(double *out, R_xlen_t n, philox4x64_key_t k,   \
                             num_vec iv, num_vec dv,                         \
                             uint64_t idx0, uint64_t dom0) {                 \
    int idx_scalar = (iv.n == 1), dom_scalar = (dv.n == 1);                  \
    philox4x64_ctr_t block = {{0, 0, 0, 0}};                                 \
    uint64_t cur_bidx = 0, cur_dom = 0;                                      \
    int have = 0;                                                            \
    for (R_xlen_t i = 0; i < n; i++) {                                       \
        uint64_t idx = idx_scalar ? idx0 : num_vec_get(iv, i, "index");      \
        uint64_t dom = dom_scalar ? dom0 : num_vec_get(dv, i, "domain");     \
        uint64_t bidx = idx >> 2;                                            \
        unsigned w = (unsigned)(idx & 3u);                                   \
        if (!have || bidx != cur_bidx || dom != cur_dom) {                   \
            block = rngat_block(k, bidx, dom, RNGAT_PURPOSE_DRAW);           \
            cur_bidx = bidx; cur_dom = dom; have = 1;                        \
        }                                                                    \
        out[i] = EMIT(block.v[w]);                                           \
    }                                                                        \
}

#define RNGAT_EMIT_BITS(bits) ((double)(uint32_t)(bits))
#define RNGAT_EMIT_UNIF(bits) u01_open(bits)
#define RNGAT_EMIT_NORM(bits) qnorm(u01_open(bits), 0.0, 1.0, TRUE, FALSE)

RNGAT_DEFINE_FILL_AT(bits, RNGAT_EMIT_BITS)
RNGAT_DEFINE_FILL_AT(unif, RNGAT_EMIT_UNIF)
RNGAT_DEFINE_FILL_AT(norm, RNGAT_EMIT_NORM)

static SEXP draw_at(SEXP key, SEXP index, SEXP domain, draw_kind kind) {
    philox4x64_key_t k = key_from_sexp(key);
    num_vec iv = num_vec_open(index, "index");
    num_vec dv = num_vec_open(domain, "domain");

    R_xlen_t ni = iv.n, nd = dv.n;
    R_xlen_t n = (ni == 0 || nd == 0) ? 0 : (ni > nd ? ni : nd);
    if (n > 0 && ((ni != 1 && ni != n) || (nd != 1 && nd != n)))
        Rf_error("`index` (length %lld) and `domain` (length %lld) must have "
                 "the same length, or one of them length 1",
                 (long long)ni, (long long)nd);

    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    double *out = REAL(ans);

    /* Scalar operands are read and validated once, outside the loop. */
    uint64_t idx0 = (ni == 1) && n ? num_vec_get(iv, 0, "index") : 0;
    uint64_t dom0 = (nd == 1) && n ? num_vec_get(dv, 0, "domain") : 0;

    switch (kind) {
    case DRAW_BITS: fill_at_bits(out, n, k, iv, dv, idx0, dom0); break;
    case DRAW_UNIF: fill_at_unif(out, n, k, iv, dv, idx0, dom0); break;
    case DRAW_NORM: fill_at_norm(out, n, k, iv, dv, idx0, dom0); break;
    }

    UNPROTECT(1);
    return ans;
}

/*
 * Fill out[0..n) with values at indices start, start+1, ..., start+n-1,
 * emitting four values per Philox call. Stamped once per draw kind with
 * the conversion baked in: relying on the inliner to specialize a
 * `draw_kind` parameter left a shared function with per-value dispatch
 * (verified in the disassembly), which cost ~25% on the bulk path.
 */
#define RNGAT_DEFINE_FILL_SEQ(SUFFIX, EMIT)                                  \
static void fill_seq_##SUFFIX(double *out, R_xlen_t n, philox4x64_key_t k,  \
                              uint64_t start, uint64_t dom) {                \
    uint64_t bidx = start >> 2;                                              \
    R_xlen_t i = 0;                                                          \
    /* head: to the next block boundary (index wraps mod 2^64 by design) */ \
    if ((start & 3u) != 0) {                                                 \
        philox4x64_ctr_t b = rngat_block(k, bidx, dom, RNGAT_PURPOSE_DRAW);  \
        for (unsigned w = start & 3u; w < 4 && i < n; w++)                   \
            out[i++] = EMIT(b.v[w]);                                         \
        bidx++;                                                              \
    }                                                                        \
    /* middle: whole blocks, four values per Philox call */                 \
    for (; n - i >= 4; i += 4, bidx++) {                                     \
        philox4x64_ctr_t b = rngat_block(k, bidx, dom, RNGAT_PURPOSE_DRAW);  \
        out[i]     = EMIT(b.v[0]);                                           \
        out[i + 1] = EMIT(b.v[1]);                                           \
        out[i + 2] = EMIT(b.v[2]);                                           \
        out[i + 3] = EMIT(b.v[3]);                                           \
    }                                                                        \
    /* tail: at most three leftover values from one final block */           \
    if (i < n) {                                                             \
        philox4x64_ctr_t b = rngat_block(k, bidx, dom, RNGAT_PURPOSE_DRAW);  \
        for (unsigned w = 0; i < n; w++, i++) out[i] = EMIT(b.v[w]);         \
    }                                                                        \
}

RNGAT_DEFINE_FILL_SEQ(bits, RNGAT_EMIT_BITS)
RNGAT_DEFINE_FILL_SEQ(unif, RNGAT_EMIT_UNIF)
RNGAT_DEFINE_FILL_SEQ(norm, RNGAT_EMIT_NORM)

/*
 * Sequential fast path: bit-identical to draw_at() over the index vector
 * start + 0:(n-1), but with no index vector to allocate, read, or validate.
 */
static SEXP draw_seq(SEXP key, SEXP start, SEXP len, SEXP domain, draw_kind kind) {
    philox4x64_key_t k = key_from_sexp(key);
    uint64_t s = u64_scalar(start, "start");
    uint64_t dom = u64_scalar(domain, "domain");

    double lend = Rf_asReal(len);
    if (ISNAN(lend) || lend != trunc(lend) || lend < 0)
        Rf_error("`len` must be a single non-negative whole number");
    if (lend > (double)R_XLEN_T_MAX)
        Rf_error("`len` is too large for an R vector");
    R_xlen_t n = (R_xlen_t)lend;

    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    double *out = REAL(ans);

    switch (kind) {
    case DRAW_BITS: fill_seq_bits(out, n, k, s, dom); break;
    case DRAW_UNIF: fill_seq_unif(out, n, k, s, dom); break;
    case DRAW_NORM: fill_seq_norm(out, n, k, s, dom); break;
    }

    UNPROTECT(1);
    return ans;
}

SEXP C_bits_at(SEXP key, SEXP index, SEXP domain) {
    return draw_at(key, index, domain, DRAW_BITS);
}

SEXP C_runif_at(SEXP key, SEXP index, SEXP domain) {
    return draw_at(key, index, domain, DRAW_UNIF);
}

SEXP C_rnorm_at(SEXP key, SEXP index, SEXP domain) {
    return draw_at(key, index, domain, DRAW_NORM);
}

SEXP C_bits_seq(SEXP key, SEXP start, SEXP len, SEXP domain) {
    return draw_seq(key, start, len, domain, DRAW_BITS);
}

SEXP C_runif_seq(SEXP key, SEXP start, SEXP len, SEXP domain) {
    return draw_seq(key, start, len, domain, DRAW_UNIF);
}

SEXP C_rnorm_seq(SEXP key, SEXP start, SEXP len, SEXP domain) {
    return draw_seq(key, start, len, domain, DRAW_NORM);
}

/* ---- registration ---- */

static const R_CallMethodDef CallEntries[] = {
    {"C_rng_key",  (DL_FUNC) &C_rng_key,  1},
    {"C_fold_in",  (DL_FUNC) &C_fold_in,  2},
    {"C_bits_at",   (DL_FUNC) &C_bits_at,   3},
    {"C_runif_at",  (DL_FUNC) &C_runif_at,  3},
    {"C_rnorm_at",  (DL_FUNC) &C_rnorm_at,  3},
    {"C_bits_seq",  (DL_FUNC) &C_bits_seq,  4},
    {"C_runif_seq", (DL_FUNC) &C_runif_seq, 4},
    {"C_rnorm_seq", (DL_FUNC) &C_rnorm_seq, 4},
    {NULL, NULL, 0}
};

void R_init_rngat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
