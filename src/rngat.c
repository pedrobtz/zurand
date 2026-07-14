#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Random123/philox.h"

/* NumPy's ziggurat tables (BSD 3-clause, see src/numpyzig/LICENSE). The
 * header is vendored verbatim and also carries float/exponential tables
 * this package does not use. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-const-variable"
#endif
#include "numpyzig/ziggurat_constants.h"
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
/* Fixed-point brackets of the wedge test, generated from the NumPy
 * tables by tools/generate-zig-bounds.R: most wedge decisions need one
 * wide multiply instead of exp(). */
#include "zigbounds.h"

/* Width of the deferred band on the chord side of the wedge shortcut, in
 * the same 52-bit fixed point as rngat_zig_gap. Must dominate the ~7-unit
 * chord crossing that ki rounding causes at layer edges (measured in
 * tools/generate-zig-bounds.R) plus the fallback's own double-rounding
 * (~tens of units); 4096 leaves two orders of magnitude of headroom at a
 * hit rate of ~2^-40 per wedge draw. */
#define RNGAT_ZIG_GUARD ((uint64_t)4096)

/*
 * rngat: stateless random numbers built on Philox4x64-10.
 *
 * A public key vector is an opaque integer matrix with one key per row and
 * four 32-bit words per row:
 *
 *   {k0 low, k0 high, k1 low, k1 high}
 *
 * R integers are signed and have NA as a distinguished bit pattern, so the
 * words are copied by bit pattern rather than interpreted as R values. Users
 * should treat keys as opaque S3 objects and use format()/print() for display.
 *
 * Every sampler is a pure function of (key, n, args). Internally the key is
 * the Philox key and the output position 0:(n - 1) is the counter. Separate
 * purpose values domain-separate folding and distribution draws.
 */

#define RNGAT_KEY_WORDS 4
#define RNGAT_ENGINE "philox4x64"
#define RNGAT_MAX_EXACT_INT 9007199254740991.0

#define RNGAT_PURPOSE_BITS    ((uint64_t)0)
#define RNGAT_PURPOSE_UNIFORM ((uint64_t)1)
#define RNGAT_PURPOSE_NORMAL  ((uint64_t)2)
#define RNGAT_PURPOSE_INTEGER ((uint64_t)3)
#define RNGAT_PURPOSE_FOLD    ((uint64_t)4)

#define RNGAT_OMP_MIN_VALUES ((R_xlen_t)32768)

/* Package-local OpenMP thread cap set by rng_threads(); 0 means "no cap"
 * (use omp_get_max_threads()). Applied per pragma via num_threads(), so it
 * governs only rngat's own fills, not other OpenMP code in the process. */
static int rngat_thread_cap = 0;

static int rngat_threads(void) {
#ifdef _OPENMP
    int max = omp_get_max_threads();
    return (rngat_thread_cap > 0 && rngat_thread_cap < max) ? rngat_thread_cap
                                                            : max;
#else
    return 1;
#endif
}

/* ---- key representation ---- */

static SEXP engine_symbol(void) {
    return Rf_install("engine");
}

static void check_engine(SEXP engine) {
    if (TYPEOF(engine) != STRSXP || Rf_xlength(engine) != 1 ||
        strcmp(CHAR(STRING_ELT(engine, 0)), RNGAT_ENGINE) != 0)
        Rf_error("`engine` must be \"%s\"", RNGAT_ENGINE);
}

/* `words` is INTEGER(key), hoisted by the caller so that no R API is
 * touched when keys are re-read inside OpenMP worker threads. */
static uint32_t load_word_at(const int *words, R_xlen_t nkey, R_xlen_t i, int word) {
    uint32_t out;
    memcpy(&out, words + i + nkey * (R_xlen_t)word, sizeof out);
    return out;
}

static void store_word_at(SEXP x, R_xlen_t nkey, R_xlen_t i, int word, uint32_t value) {
    memcpy(INTEGER(x) + i + nkey * (R_xlen_t)word, &value, sizeof value);
}

static R_xlen_t key_count(SEXP key) {
    if (TYPEOF(key) != INTSXP ||
        !Rf_inherits(key, "rng_key"))
        Rf_error("`key` must be an <rng_key>; create one with rng_key()");

    SEXP dim = Rf_getAttrib(key, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || Rf_xlength(dim) != 2 ||
        INTEGER(dim)[0] < 0 || INTEGER(dim)[1] != RNGAT_KEY_WORDS)
        Rf_error("`key` has an invalid internal shape");
    if (Rf_xlength(key) != (R_xlen_t)INTEGER(dim)[0] * RNGAT_KEY_WORDS)
        Rf_error("`key` has an invalid internal length");

    SEXP engine = Rf_getAttrib(key, engine_symbol());
    check_engine(engine);
    return (R_xlen_t)INTEGER(dim)[0];
}

static SEXP key_engine(SEXP key) {
    (void) key_count(key);
    return Rf_getAttrib(key, engine_symbol());
}

static philox4x64_key_t key_from_words(const int *words, R_xlen_t nkey, R_xlen_t i) {
    philox4x64_key_t k;
    uint64_t w0 = (uint64_t)load_word_at(words, nkey, i, 0);
    uint64_t w1 = (uint64_t)load_word_at(words, nkey, i, 1);
    uint64_t w2 = (uint64_t)load_word_at(words, nkey, i, 2);
    uint64_t w3 = (uint64_t)load_word_at(words, nkey, i, 3);
    k.v[0] = w0 | (w1 << 32);
    k.v[1] = w2 | (w3 << 32);
    return k;
}

static SEXP alloc_key_vector(R_xlen_t nkey, SEXP engine) {
    check_engine(engine);
    if (nkey > R_XLEN_T_MAX / RNGAT_KEY_WORDS)
        Rf_error("too many keys requested");
    if (nkey > INT_MAX)
        Rf_error("too many keys for a matrix-shaped key vector");

    SEXP ans = PROTECT(Rf_allocVector(INTSXP, nkey * RNGAT_KEY_WORDS));
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = (int)nkey;
    INTEGER(dim)[1] = RNGAT_KEY_WORDS;
    Rf_setAttrib(ans, R_DimSymbol, dim);
    Rf_setAttrib(ans, engine_symbol(), engine);
    SEXP cls = PROTECT(Rf_mkString("rng_key"));
    Rf_classgets(ans, cls);
    UNPROTECT(3);
    return ans;
}

static void set_key_words(SEXP ans, R_xlen_t nkey, R_xlen_t i, uint64_t k0, uint64_t k1) {
    store_word_at(ans, nkey, i, 0, (uint32_t)(k0 & 0xffffffffu));
    store_word_at(ans, nkey, i, 1, (uint32_t)(k0 >> 32));
    store_word_at(ans, nkey, i, 2, (uint32_t)(k1 & 0xffffffffu));
    store_word_at(ans, nkey, i, 3, (uint32_t)(k1 >> 32));
}

SEXP C_rng_key_format(SEXP key) {
    R_xlen_t nkey = key_count(key);
    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(Rf_allocVector(STRSXP, nkey));
    char buf[42];
    for (R_xlen_t i = 0; i < nkey; i++) {
        snprintf(buf, sizeof buf, "rng_key[%08" PRIx32 "%08" PRIx32 "%08" PRIx32 "%08" PRIx32 "]",
                 load_word_at(kw, nkey, i, 0), load_word_at(kw, nkey, i, 1),
                 load_word_at(kw, nkey, i, 2), load_word_at(kw, nkey, i, 3));
        SET_STRING_ELT(ans, i, Rf_mkChar(buf));
    }
    UNPROTECT(1);
    return ans;
}

/* ---- scalar validation ---- */

static double numeric_scalar(SEXP x, const char *what) {
    if (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)
        Rf_error("`%s` must be an integer or double scalar", what);
    if (Rf_xlength(x) != 1)
        Rf_error("`%s` must be a single value", what);

    double out;
    if (TYPEOF(x) == INTSXP) {
        int value = INTEGER(x)[0];
        if (value == NA_INTEGER)
            Rf_error("`%s` must not be missing", what);
        out = (double)value;
    } else {
        out = REAL(x)[0];
        if (!R_FINITE(out))
            Rf_error("`%s` must be finite", what);
    }
    return out;
}

static uint64_t u64_seed(SEXP x) {
    double value = numeric_scalar(x, "seed");
    if (value != trunc(value))
        Rf_error("`seed` must be a whole number");
    if (value < 0 || value > RNGAT_MAX_EXACT_INT)
        Rf_error("`seed` must be between 0 and 2^53 - 1");
    return (uint64_t)value;
}

static R_xlen_t length_scalar(SEXP x, const char *what) {
    double value = numeric_scalar(x, what);
    if (value != trunc(value) || value < 0)
        Rf_error("`%s` must be a single non-negative whole number", what);
    if (value > (double)R_XLEN_T_MAX)
        Rf_error("`%s` is too large for an R vector", what);
    return (R_xlen_t)value;
}

static int bits_scalar(SEXP x) {
    double value = numeric_scalar(x, "bits");
    if (value != 32.0 && value != 64.0)
        Rf_error("`bits` must be either 32 or 64");
    return (int)value;
}

static double finite_scalar(SEXP x, const char *what) {
    return numeric_scalar(x, what);
}

static int integer_bound(SEXP x, const char *what) {
    double value = numeric_scalar(x, what);
    if (value != trunc(value))
        Rf_error("`%s` must be a whole number", what);
    if (value < (double)INT_MIN + 1.0 || value > (double)INT_MAX)
        Rf_error("`%s` must be in R's non-missing integer range", what);
    return (int)value;
}

static int64_t exact_i64_from_double(double value, const char *what) {
    if (!R_FINITE(value) || value != trunc(value))
        Rf_error("`%s` must contain only finite whole numbers", what);
    if (value < -RNGAT_MAX_EXACT_INT || value > RNGAT_MAX_EXACT_INT)
        Rf_error("`%s` values must be exactly representable integers", what);
    return (int64_t)value;
}

/* ---- core generator ---- */

static uint64_t splitmix64_next(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

R123_STATIC_INLINE R123_FORCE_INLINE(philox4x64_ctr_t rngat_block(
    philox4x64_key_t key, uint64_t index, uint64_t domain, uint64_t purpose));
R123_STATIC_INLINE philox4x64_ctr_t rngat_block(philox4x64_key_t key,
                                                uint64_t index,
                                                uint64_t domain,
                                                uint64_t purpose) {
    philox4x64_ctr_t ctr = {{index, domain, purpose, 0}};
    return philox4x64(ctr, key);
}

static uint64_t rngat_word(philox4x64_key_t key, uint64_t index,
                           uint64_t domain, uint64_t purpose) {
    philox4x64_ctr_t block = rngat_block(key, index >> 2, domain, purpose);
    return block.v[index & 3u];
}

/* Top 52 random bits, centered in their bucket -> exactly (m + 0.5) *
 * 2^-52, strictly inside (0, 1): stuff the bits into the mantissa of a
 * double in [1, 2), then shift the interval; the subtraction is exact.
 * Integer ops plus one FP subtract, so fill loops vectorize. (The former
 * ((bits >> 11) + 0.5) * 2^-53 form rounded for bits >= 2^52 and could
 * even yield exactly 1.0, breaking the open interval.) */
static double u01_open(uint64_t bits) {
    uint64_t stuffed = (bits >> 12) | UINT64_C(0x3ff0000000000000);
    double d;
    memcpy(&d, &stuffed, sizeof d);
    return d - (1.0 - 0x1.0p-53);
}

/* 53-bit uniform in [0, 1), numpy's next_double mapping. */
#define RNGAT_U64_TO_DOUBLE(u) (((u) >> 11) * 0x1.0p-53)

/* Retry-word stream for the slow path: group g >= 1 is one Philox block
 * at counter {index, g, purpose}, consumed word by word, so all four
 * words of each retry block are used and draw `index` stays a pure
 * function of (key, index). Distinct from the fast-path counters, whose
 * second word is always 0. (The former one-word-per-attempt scheme paid
 * a full Philox call per retry word and discarded the other three;
 * batching measured ~3% faster on bulk rnorm overall.) */
typedef struct {
    philox4x64_ctr_t blk;
    uint64_t group;
    unsigned w;
} zig_stream;

R123_STATIC_INLINE uint64_t zig_next(zig_stream *s, philox4x64_key_t key,
                                     uint64_t index) {
    if (s->w == 4) {
        s->blk = rngat_block(key, index, ++s->group, RNGAT_PURPOSE_NORMAL);
        s->w = 0;
    }
    return s->blk.v[s->w++];
}

/* Cold continuation of zig_normal_at once the one-word fast path has
 * rejected: wedge acceptance and the layer-0 tail, redrawing words from
 * the zig_stream retry blocks. Most wedge decisions resolve in fixed
 * point against the Dnorm bracket; only the narrow ambiguous band pays
 * exp(). Kept out of line so the fast path inlines into the sampler
 * loops. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, cold))
#endif
static double zig_normal_slow(philox4x64_key_t key, uint64_t index, uint64_t r) {
    zig_stream s = {.group = 0, .w = 4};
    int idx = (int)(r & 0xff);
    int sign = (int)((r >> 8) & 0x1);
    uint64_t rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
    double x = (double)rabs * wi_double[idx];
    if (sign)
        x = -x;
    /* the caller established rabs >= ki_double[idx] for the entry word */

    for (;;) {
        uint64_t Y = zig_next(&s, key, index);

        if (idx == 0) {
            /* layer-0 tail; the first ordinate reuses Y */
            double yy = -log1p(-RNGAT_U64_TO_DOUBLE(Y));
            for (;;) {
                double xx = -ziggurat_nor_inv_r *
                    log1p(-RNGAT_U64_TO_DOUBLE(zig_next(&s, key, index)));
                if (yy + yy > xx * xx)
                    return sign ? -(ziggurat_nor_r + xx)
                                : ziggurat_nor_r + xx;
                yy = -log1p(-RNGAT_U64_TO_DOUBLE(zig_next(&s, key, index)));
            }
        }

        /* wedge: compare Y * (width of layer idx) against the position in
         * the layer, with rngat_zig_gap bracketing the exp curve around its
         * chord; f is concave below the inflection layer (x < 1, curve
         * above the chord) and convex above it. RNGAT_ZIG_GUARD widens the
         * ambiguous band on the chord side: ki rounding lets the true curve
         * cross the chord by a few fixed-point units, and the fallback's
         * own double rounding lives there too, so the hairline band defers
         * to the exp() test instead of deciding. With the guard, shortcut
         * decisions never contradict the fallback. */
        uint64_t L = (UINT64_C(1) << 52) - ki_double[idx];
        uint64_t R = (UINT64_C(1) << 52) - rabs;
        uint64_t YL;
        (void)mulhilo64(Y, L, &YL);
        int accept, reject;
        if (idx > RNGAT_ZIG_INFLECTION) {
            reject = YL > R + RNGAT_ZIG_GUARD;
            accept = !reject && YL + rngat_zig_gap[idx] < R;
        } else if (idx < RNGAT_ZIG_INFLECTION) {
            accept = YL + RNGAT_ZIG_GUARD < R;
            reject = !accept && YL > R + rngat_zig_gap[idx];
        } else {
            reject = YL > R + rngat_zig_gap_hi52;
            accept = !reject && YL + rngat_zig_gap[idx] < R;
        }
        if (accept)
            return x;
        if (!reject) {
            double u = RNGAT_U64_TO_DOUBLE(Y);
            if ((fi_double[idx - 1] - fi_double[idx]) * u + fi_double[idx] <
                exp(-0.5 * x * x))
                return x;
        }

        r = zig_next(&s, key, index);
        idx = (int)(r & 0xff);
        sign = (int)((r >> 8) & 0x1);
        rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
        x = (double)rabs * wi_double[idx];
        if (sign)
            x = -x;
        if (rabs < ki_double[idx])
            return x;
    }
}

/* One standard normal deviate for output position `index`, following
 * numpy's random_standard_normal: 8 bits of layer index, 1 sign bit and a
 * 52-bit magnitude from a single word decide ~99% of draws with one table
 * compare. `r` is the attempt-0 word from the shared Philox block. */
R123_STATIC_INLINE double zig_normal_at(philox4x64_key_t key, uint64_t index,
                                        uint64_t r) {
    int idx = (int)(r & 0xff);
    uint64_t rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
    /* Negate before the int->double convert: (-rabs) * wi and
     * -(rabs * wi) are bit-identical, and the signed form costs one
     * cneg instead of a second multiply feeding a select. Computed
     * unconditionally so the convert/multiply chain starts before the
     * acceptance compare resolves. */
    int64_t s = (r >> 8) & 0x1 ? -(int64_t)rabs : (int64_t)rabs;
    double x = (double)s * wi_double[idx];
    if (R123_BUILTIN_EXPECT(rabs < ki_double[idx], 1))
        return x;
    return zig_normal_slow(key, index, r);
}

/* Each Philox block yields four outputs and depends only on its counter,
 * so the block loops below parallelize with bit-identical results.
 * `threads` gates the inner parallel region; callers pass 0 when they
 * already parallelize over key columns. */
/* Fills standard normals; C_rng_normal applies mean/sd in a separate
 * vectorizable pass so the hot loop stays load/compare/multiply only.
 *
 * Works in two passes over a small stack chunk: a tight Philox-only loop
 * (independent iterations the compiler can pipeline across the 10-round
 * dependency chain), then the ziggurat transform over the buffered words.
 * Same words, same decisions, same output as a fused loop. (Transforming
 * in place over the output slice instead — randompack's layout — measured
 * ~12% slower here: it turns the output stream's write-once pattern into
 * write-read-write.) */
#define RNGAT_CHUNK_BLOCKS 128 /* 512 words, 4 KiB per thread */

static void fill_normal_column(double *out, R_xlen_t n,
                               philox4x64_key_t key, int threads) {
    R_xlen_t nblock = n >> 2;
    R_xlen_t nchunk = (nblock + RNGAT_CHUNK_BLOCKS - 1) / RNGAT_CHUNK_BLOCKS;
#ifdef _OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, nblock, nchunk, key) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        uint64_t buf[RNGAT_CHUNK_BLOCKS * 4];
        R_xlen_t b0 = c * RNGAT_CHUNK_BLOCKS;
        int nb = (int)(nblock - b0 < RNGAT_CHUNK_BLOCKS ? nblock - b0
                                                        : RNGAT_CHUNK_BLOCKS);
        for (int j = 0; j < nb; j++) {
            philox4x64_ctr_t block = rngat_block(key, (uint64_t)(b0 + j), 0,
                                                 RNGAT_PURPOSE_NORMAL);
            memcpy(buf + 4 * j, block.v, sizeof block.v);
        }

        double *o = out + (b0 << 2);
        uint64_t base = (uint64_t)(b0 << 2);
        for (int j = 0; j < nb * 4; j++)
            o[j] = zig_normal_at(key, base + (uint64_t)j, buf[j]);
    }

    R_xlen_t i = nblock << 2;
    if (i < n) {
        philox4x64_ctr_t block = rngat_block(key, (uint64_t)nblock, 0,
                                             RNGAT_PURPOSE_NORMAL);
        for (unsigned w = 0; i < n; w++, i++)
            out[i] = zig_normal_at(key, (uint64_t)i, block.v[w]);
    }
}

/* Fills uniforms on (0, 1); C_rng_uniform applies min/span in a separate
 * vectorizable pass so the hot loop stays Philox plus the bit trick.
 * `threads` is the thread count for the inner parallel region; callers
 * pass 0 to keep the column serial (e.g. when parallelizing over key
 * columns). */
static void fill_uniform_column(double *out, R_xlen_t n,
                                philox4x64_key_t key, int threads) {
    R_xlen_t nblock = n >> 2;
#ifdef _OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, nblock, key) schedule(static)
#endif
    for (R_xlen_t b = 0; b < nblock; b++) {
        philox4x64_ctr_t block = rngat_block(key, (uint64_t)b, 0,
                                             RNGAT_PURPOSE_UNIFORM);
        double *o = out + (b << 2);
        o[0] = u01_open(block.v[0]);
        o[1] = u01_open(block.v[1]);
        o[2] = u01_open(block.v[2]);
        o[3] = u01_open(block.v[3]);
    }

    R_xlen_t i = nblock << 2;
    if (i < n) {
        philox4x64_ctr_t block = rngat_block(key, (uint64_t)nblock, 0,
                                             RNGAT_PURPOSE_UNIFORM);
        for (unsigned w = 0; i < n; w++, i++)
            out[i] = u01_open(block.v[w]);
    }
}

/* ---- stable fold hashing ---- */

static void hash_byte(uint64_t *h, unsigned char byte) {
    *h ^= (uint64_t)byte;
    *h *= UINT64_C(1099511628211);
}

static void hash_bytes(uint64_t *h, const unsigned char *bytes, size_t n) {
    for (size_t i = 0; i < n; i++) hash_byte(h, bytes[i]);
}

static void hash_u64(uint64_t *h, uint64_t x) {
    for (int i = 0; i < 8; i++) {
        hash_byte(h, (unsigned char)(x & 0xffu));
        x >>= 8;
    }
}

static uint64_t hash_data(SEXP data) {
    uint64_t h = UINT64_C(1469598103934665603);
    R_xlen_t n = Rf_xlength(data);

    switch (TYPEOF(data)) {
    case STRSXP:
        hash_byte(&h, 1);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++) {
            SEXP s = STRING_ELT(data, i);
            if (s == NA_STRING)
                Rf_error("`data` must not contain missing values");
            const char *utf8 = Rf_translateCharUTF8(s);
            size_t len = strlen(utf8);
            hash_u64(&h, (uint64_t)len);
            hash_bytes(&h, (const unsigned char *)utf8, len);
        }
        break;

    /* Integer and double data share one tag and one width so that folding
     * by 1L and by 1 derives the same key; doubles are already restricted
     * to exact whole numbers. */
    case INTSXP:
        hash_byte(&h, 2);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++) {
            int value = INTEGER(data)[i];
            if (value == NA_INTEGER)
                Rf_error("`data` must not contain missing values");
            hash_u64(&h, (uint64_t)(int64_t)value);
        }
        break;

    case REALSXP:
        hash_byte(&h, 2);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++) {
            int64_t value = exact_i64_from_double(REAL(data)[i], "data");
            hash_u64(&h, (uint64_t)value);
        }
        break;

    case LGLSXP:
        hash_byte(&h, 4);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++) {
            int value = LOGICAL(data)[i];
            if (value == NA_LOGICAL)
                Rf_error("`data` must not contain missing values");
            hash_byte(&h, (unsigned char)(value ? 1 : 0));
        }
        break;

    case RAWSXP:
        hash_byte(&h, 5);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++)
            hash_byte(&h, (unsigned char)RAW(data)[i]);
        break;

    default:
        Rf_error("`data` must be a character, integer, double, logical or raw vector");
    }

    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
}

/* ---- .Call entry points: keys ---- */

static R_xlen_t checked_product(R_xlen_t a, R_xlen_t b, const char *what) {
    if (a != 0 && b > R_XLEN_T_MAX / a)
        Rf_error("requested `%s` output is too large", what);
    return a * b;
}

static void set_sample_dim(SEXP ans, R_xlen_t n, R_xlen_t nkey) {
    if (nkey == 1)
        return;
    if (n > INT_MAX || nkey > INT_MAX)
        Rf_error("multi-key output dimensions are too large");
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = (int)n;
    INTEGER(dim)[1] = (int)nkey;
    Rf_setAttrib(ans, R_DimSymbol, dim);
    UNPROTECT(1);
}

SEXP C_rng_key(SEXP seed, SEXP n_, SEXP engine) {
    check_engine(engine);
    uint64_t state = u64_seed(seed);
    R_xlen_t n = length_scalar(n_, "n");

    SEXP ans = PROTECT(alloc_key_vector(n, engine));
    for (R_xlen_t i = 0; i < n; i++) {
        uint64_t k0 = splitmix64_next(&state);
        uint64_t k1 = splitmix64_next(&state);
        set_key_words(ans, n, i, k0, k1);
    }
    UNPROTECT(1);
    return ans;
}

static uint32_t r_random_u32(void) {
    return (uint32_t)(unif_rand() * 4294967296.0);
}

SEXP C_rng_key_from_r(SEXP n_, SEXP engine) {
    check_engine(engine);
    R_xlen_t n = length_scalar(n_, "n");

    SEXP ans = PROTECT(alloc_key_vector(n, engine));
    GetRNGstate();
    for (R_xlen_t i = 0; i < n; i++) {
        uint64_t k0 = ((uint64_t)r_random_u32() << 32) | (uint64_t)r_random_u32();
        uint64_t k1 = ((uint64_t)r_random_u32() << 32) | (uint64_t)r_random_u32();
        set_key_words(ans, n, i, k0, k1);
    }
    PutRNGstate();
    UNPROTECT(1);
    return ans;
}

SEXP C_rng_fold(SEXP key, SEXP data) {
    R_xlen_t nkey = key_count(key);
    uint64_t h = hash_data(data);
    uint64_t domain = h ^ UINT64_C(0x9e3779b97f4a7c15);
    SEXP engine = key_engine(key);

    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(alloc_key_vector(nkey, engine));
    for (R_xlen_t i = 0; i < nkey; i++) {
        philox4x64_key_t k = key_from_words(kw, nkey, i);
        philox4x64_ctr_t out = rngat_block(k, h, domain, RNGAT_PURPOSE_FOLD);
        set_key_words(ans, nkey, i, out.v[0], out.v[1]);
    }
    UNPROTECT(1);
    return ans;
}

/* NULL queries the effective maximum. A value sets the cap (0 removes it)
 * and returns the previous *raw* cap, 0 when none was set, so that
 * rng_threads(old) restores the uncapped state rather than pinning
 * whatever the OpenMP maximum happened to be at save time. */
SEXP C_rng_threads(SEXP threads_) {
    if (Rf_isNull(threads_))
        return Rf_ScalarInteger(rngat_threads());

    double value = numeric_scalar(threads_, "threads");
    if (value != trunc(value) || value < 0)
        Rf_error("`threads` must be a single whole number that is at least 1, "
                 "or 0 to remove the cap");
    int prev = rngat_thread_cap;
    rngat_thread_cap = value > (double)INT_MAX ? INT_MAX : (int)value;
    return Rf_ScalarInteger(prev);
}

/* ---- .Call entry points: samplers ---- */

SEXP C_rng_uniform(SEXP key, SEXP n_, SEXP min_, SEXP max_) {
    R_xlen_t nkey = key_count(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    double min = finite_scalar(min_, "min");
    double max = finite_scalar(max_, "max");
    if (min > max)
        Rf_error("`min` must be less than or equal to `max`");

    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
    set_sample_dim(ans, n, nkey);
    double *out = REAL(ans);
    double span = max - min;
    int nt = rngat_threads();
    if (span == 0.0) {
#ifdef _OPENMP
#pragma omp parallel for if(nt > 1 && total >= RNGAT_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, min) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = min;
        UNPROTECT(1);
        return ans;
    }
    int par_cols = nt > 1 && nkey > 1 && total >= RNGAT_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= RNGAT_OMP_MIN_VALUES ? nt : 0;
#ifdef _OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, par_rows) schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_from_words(kw, nkey, col);
        fill_uniform_column(out + n * col, n, k, par_rows);
    }
    if (min != 0.0 || span != 1.0) {
#ifdef _OPENMP
#pragma omp parallel for if(nt > 1 && total >= RNGAT_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, min, span) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = min + span * out[i];
    }
    UNPROTECT(1);
    return ans;
}

SEXP C_rng_normal(SEXP key, SEXP n_, SEXP mean_, SEXP sd_) {
    R_xlen_t nkey = key_count(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    double mean = finite_scalar(mean_, "mean");
    double sd = finite_scalar(sd_, "sd");
    if (sd < 0)
        Rf_error("`sd` must be non-negative");

    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
    set_sample_dim(ans, n, nkey);
    double *out = REAL(ans);
    int nt = rngat_threads();
    if (sd == 0.0) {
#ifdef _OPENMP
#pragma omp parallel for if(nt > 1 && total >= RNGAT_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, mean) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = mean;
        UNPROTECT(1);
        return ans;
    }
    int par_cols = nt > 1 && nkey > 1 && total >= RNGAT_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= RNGAT_OMP_MIN_VALUES ? nt : 0;
#ifdef _OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, par_rows) schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_from_words(kw, nkey, col);
        fill_normal_column(out + n * col, n, k, par_rows);
    }
    if (mean != 0.0 || sd != 1.0) {
#ifdef _OPENMP
#pragma omp parallel for if(nt > 1 && total >= RNGAT_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, mean, sd) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = mean + sd * out[i];
    }
    UNPROTECT(1);
    return ans;
}

/* Lemire bounded sampling: accept unless the low half of x * range falls
 * below the bias threshold, which the caller computes once per range. */
static int lemire_accept(uint32_t x, uint32_t range, uint32_t threshold,
                         uint32_t *offset) {
    uint64_t product = (uint64_t)x * (uint64_t)range;
    if ((uint32_t)product < threshold)
        return 0;
    *offset = (uint32_t)(product >> 32);
    return 1;
}

/* Rejection continuation for one index; attempt 0 was consumed from the
 * shared block by fill_integer_column. */
static uint32_t bounded_u32_retry(philox4x64_key_t k, uint64_t index,
                                  uint32_t range, uint32_t threshold) {
    for (uint64_t attempt = 1; ; attempt++) {
        uint32_t x = (uint32_t)rngat_word(k, index, attempt, RNGAT_PURPOSE_INTEGER);
        uint32_t offset;
        if (lemire_accept(x, range, threshold, &offset))
            return offset;
    }
}

static void fill_integer_column(int *out, R_xlen_t n, philox4x64_key_t key,
                                int min, uint32_t range, uint32_t threshold,
                                int threads) {
    R_xlen_t nblock = n >> 2;
#ifdef _OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, nblock, key, min, range, threshold) schedule(static)
#endif
    for (R_xlen_t b = 0; b < nblock; b++) {
        philox4x64_ctr_t block = rngat_block(key, (uint64_t)b, 0,
                                             RNGAT_PURPOSE_INTEGER);
        int *o = out + (b << 2);
        for (unsigned w = 0; w < 4; w++) {
            uint32_t offset;
            if (!lemire_accept((uint32_t)block.v[w], range, threshold, &offset))
                offset = bounded_u32_retry(key, (uint64_t)((b << 2) + w),
                                           range, threshold);
            o[w] = min + (int)offset;
        }
    }

    R_xlen_t i = nblock << 2;
    if (i < n) {
        philox4x64_ctr_t block = rngat_block(key, (uint64_t)nblock, 0,
                                             RNGAT_PURPOSE_INTEGER);
        for (unsigned w = 0; i < n; w++, i++) {
            uint32_t offset;
            if (!lemire_accept((uint32_t)block.v[w], range, threshold, &offset))
                offset = bounded_u32_retry(key, (uint64_t)i, range, threshold);
            out[i] = min + (int)offset;
        }
    }
}

SEXP C_rng_integer(SEXP key, SEXP n_, SEXP min_, SEXP max_) {
    R_xlen_t nkey = key_count(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    int min = integer_bound(min_, "min");
    int max = integer_bound(max_, "max");
    if (min > max)
        Rf_error("`min` must be less than or equal to `max`");

    uint32_t range = (uint32_t)((int64_t)max - (int64_t)min + 1);
    uint32_t threshold = (uint32_t)((UINT64_C(0x100000000) - range) % range);
    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, total));
    set_sample_dim(ans, n, nkey);
    int *out = INTEGER(ans);
    int nt = rngat_threads();
    int par_cols = nt > 1 && nkey > 1 && total >= RNGAT_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= RNGAT_OMP_MIN_VALUES ? nt : 0;
#ifdef _OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, min, range, threshold, par_rows) schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_from_words(kw, nkey, col);
        fill_integer_column(out + n * col, n, k, min, range, threshold, par_rows);
    }
    UNPROTECT(1);
    return ans;
}

SEXP C_rng_bits(SEXP key, SEXP n_, SEXP bits_) {
    R_xlen_t nkey = key_count(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    int bits = bits_scalar(bits_);
    const int *kw = INTEGER(key);

    if (bits == 32) {
        SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
        set_sample_dim(ans, n, nkey);
        double *out = REAL(ans);
        for (R_xlen_t col = 0; col < nkey; col++) {
            philox4x64_key_t k = key_from_words(kw, nkey, col);
            double *col_out = out + n * col;
            for (R_xlen_t i = 0; i < n; i += 4) {
                philox4x64_ctr_t block = rngat_block(k, (uint64_t)(i >> 2), 0,
                                                     RNGAT_PURPOSE_BITS);
                R_xlen_t stop = n - i < 4 ? n - i : 4;
                for (R_xlen_t w = 0; w < stop; w++)
                    col_out[i + w] = (double)(uint32_t)block.v[w];
            }
        }
        UNPROTECT(1);
        return ans;
    }

    SEXP ans = PROTECT(Rf_allocVector(STRSXP, total));
    set_sample_dim(ans, n, nkey);
    char buf[17];
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_from_words(kw, nkey, col);
        for (R_xlen_t i = 0; i < n; i += 4) {
            philox4x64_ctr_t block = rngat_block(k, (uint64_t)(i >> 2), 0,
                                                 RNGAT_PURPOSE_BITS);
            R_xlen_t stop = n - i < 4 ? n - i : 4;
            for (R_xlen_t w = 0; w < stop; w++) {
                snprintf(buf, sizeof buf, "%016" PRIx64, block.v[w]);
                SET_STRING_ELT(ans, i + w + n * col, Rf_mkChar(buf));
            }
        }
    }
    UNPROTECT(1);
    return ans;
}

/* ---- registration ---- */

static const R_CallMethodDef CallEntries[] = {
    {"C_rng_key",        (DL_FUNC) &C_rng_key,        3},
    {"C_rng_key_from_r", (DL_FUNC) &C_rng_key_from_r, 2},
    {"C_rng_key_format", (DL_FUNC) &C_rng_key_format, 1},
    {"C_rng_fold",       (DL_FUNC) &C_rng_fold,       2},
    {"C_rng_uniform",    (DL_FUNC) &C_rng_uniform,    4},
    {"C_rng_normal",     (DL_FUNC) &C_rng_normal,     4},
    {"C_rng_integer",    (DL_FUNC) &C_rng_integer,    4},
    {"C_rng_bits",       (DL_FUNC) &C_rng_bits,       3},
    {"C_rng_threads",    (DL_FUNC) &C_rng_threads,    1},
    {NULL, NULL, 0}
};

void R_init_rngat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
