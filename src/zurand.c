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
/* WRE's idiom is a bare `#ifdef _OPENMP` around <omp.h>, which assumes that a
 * compiler defining _OPENMP also ships the header. The clang sanitizer
 * containers disprove that: R's Makeconf there puts -fopenmp in
 * SHLIB_OPENMP_CFLAGS, so _OPENMP is defined, while libomp-dev is absent and
 * the include is a fatal error. Gate every OpenMP use on one macro instead,
 * so that build is simply serial rather than broken. */
#ifdef __has_include
#  if __has_include(<omp.h>)
#    define ZURAND_HAVE_OMP_H 1
#  endif
#else
/* No __has_include: assume the header is present, the previous behaviour. */
#  define ZURAND_HAVE_OMP_H 1
#endif

#if defined(_OPENMP) && defined(ZURAND_HAVE_OMP_H)
#  include <omp.h>
#  define ZURAND_OPENMP 1
#endif

#include "Random123/philox.h"
#include "Random123/threefry.h"

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

/* No fused multiply-add in this translation unit.
 *
 * The scaling passes are `min + span * x` and `mean + sd * x`. A compiler
 * with FMA available -- every arm64 target, and x86_64 built with -mfma --
 * is free to contract each into a single fused instruction that rounds once
 * instead of twice. The result differs from the unfused form by an ulp.
 *
 * That silently breaks the package's central promise. CI caught it: macOS
 * arm64 disagreed with x86_64 on exactly the two golden expectations that
 * pass a non-default min/max or mean/sd, and on no others. Reproducibility
 * is the product here, so the ulp is not an acceptable trade for one fused
 * instruction in a pass that is memory-bound anyway.
 *
 * clang honours this pragma. GCC has only implemented it recently, and
 * cannot contract on an x86_64 target without FMA enabled in any case; a
 * GCC build for an FMA-capable target should also pass -ffp-contract=off,
 * which is what the planned configure script will probe for. */
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 14)
#pragma STDC FP_CONTRACT OFF
#endif

/* Width of the deferred band on the chord side of the wedge shortcut, in
 * the same 52-bit fixed point as zurand_zig_gap. Must dominate the ~7-unit
 * chord crossing that ki rounding causes at layer edges (measured in
 * tools/generate-zig-bounds.R) plus the fallback's own double-rounding
 * (~tens of units); 4096 leaves two orders of magnitude of headroom at a
 * hit rate of ~2^-40 per wedge draw. */
#define ZURAND_ZIG_GUARD ((uint64_t)4096)

/*
 * zurand: stateless random numbers built on Philox4x64-10.
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

#define ZURAND_KEY_WORDS 4
/* Engine names, as stored on every key's `engine` attribute. Adding one is
 * additive by design: an existing key keeps its name and therefore its
 * stream forever, and a caller opts into a new engine explicitly. */
#define ZURAND_ENGINE_PHILOX   "philox4x64"
#define ZURAND_ENGINE_THREEFRY "threefry4x64"
#define ZURAND_ENGINE_XOSHIRO  "xoshiro256pp"
#define ZURAND_ENGINE ZURAND_ENGINE_PHILOX  /* the default */

/* Both engines produce the same 256-bit block type; only their key types
 * differ. The sampler core is written against this. */
typedef philox4x64_ctr_t zurand_ctr_t;

typedef enum { ZURAND_ENG_PHILOX = 0, ZURAND_ENG_THREEFRY = 1,
               ZURAND_ENG_XOSHIRO = 2 } zurand_engine_t;
#define ZURAND_MAX_EXACT_INT 9007199254740991.0

#define ZURAND_PURPOSE_BITS    ((uint64_t)0)
#define ZURAND_PURPOSE_UNIFORM ((uint64_t)1)
#define ZURAND_PURPOSE_NORMAL  ((uint64_t)2)
#define ZURAND_PURPOSE_INTEGER ((uint64_t)3)
#define ZURAND_PURPOSE_FOLD    ((uint64_t)4)

#define ZURAND_OMP_MIN_VALUES ((R_xlen_t)32768)

/* Package-local OpenMP thread cap set by rng_threads(); 0 means "no cap"
 * (use omp_get_max_threads()). Applied per pragma via num_threads(), so it
 * governs only zurand's own fills, not other OpenMP code in the process. */
static int zurand_thread_cap = 0;

static int zurand_threads(void) {
#ifdef ZURAND_OPENMP
    int max = omp_get_max_threads();
    return (zurand_thread_cap > 0 && zurand_thread_cap < max) ? zurand_thread_cap
                                                            : max;
#else
    return 1;
#endif
}

/* ---- key representation ---- */

static SEXP engine_symbol(void) {
    return Rf_install("engine");
}

static int engine_code(SEXP engine) {
    if (TYPEOF(engine) == STRSXP && Rf_xlength(engine) == 1) {
        const char *s = CHAR(STRING_ELT(engine, 0));
        if (strcmp(s, ZURAND_ENGINE_PHILOX) == 0)   return ZURAND_ENG_PHILOX;
        if (strcmp(s, ZURAND_ENGINE_THREEFRY) == 0) return ZURAND_ENG_THREEFRY;
        if (strcmp(s, ZURAND_ENGINE_XOSHIRO) == 0)  return ZURAND_ENG_XOSHIRO;
    }
    return -1;
}

static void check_engine(SEXP engine) {
    if (engine_code(engine) < 0)
        Rf_error("`engine` must be \"%s\", \"%s\" or \"%s\"",
                 ZURAND_ENGINE_PHILOX, ZURAND_ENGINE_THREEFRY,
                 ZURAND_ENGINE_XOSHIRO);
}

/* The engine a key was created with; every sampler dispatches on this. */
static zurand_engine_t key_engine_code(SEXP key) {
    SEXP e = Rf_getAttrib(key, engine_symbol());
    int c = engine_code(e);
    if (c < 0)
        Rf_error("<rng_key> has an unrecognised engine");
    return (zurand_engine_t)c;
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
        INTEGER(dim)[0] < 0 || INTEGER(dim)[1] != ZURAND_KEY_WORDS)
        Rf_error("`key` has an invalid internal shape");
    if (Rf_xlength(key) != (R_xlen_t)INTEGER(dim)[0] * ZURAND_KEY_WORDS)
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

/* Threefry takes a 256-bit key; a zurand key stores 128. The upper half is
 * derived from the lower by xor with the same two constants Threefry itself
 * uses as its Weyl increments, so all four key words depend on the stored
 * key rather than two of them being fixed for every user.
 *
 * This does not invent entropy: the key space stays 128 bits, which is
 * exactly what rng_key(seed) produces for philox too. It only avoids handing
 * Threefry's schedule two constant words. */
static threefry4x64_key_t key_from_words_threefry(const int *words,
                                                  R_xlen_t nkey, R_xlen_t i) {
    philox4x64_key_t p = key_from_words(words, nkey, i);
    threefry4x64_key_t k;
    k.v[0] = p.v[0];
    k.v[1] = p.v[1];
    k.v[2] = p.v[0] ^ UINT64_C(0x9e3779b97f4a7c15);
    k.v[3] = p.v[1] ^ UINT64_C(0xbb67ae8584caa73b);
    return k;
}

static SEXP alloc_key_vector(R_xlen_t nkey, SEXP engine) {
    check_engine(engine);
    if (nkey > R_XLEN_T_MAX / ZURAND_KEY_WORDS)
        Rf_error("too many keys requested");
    if (nkey > INT_MAX)
        Rf_error("too many keys for a matrix-shaped key vector");

    SEXP ans = PROTECT(Rf_allocVector(INTSXP, nkey * ZURAND_KEY_WORDS));
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = (int)nkey;
    INTEGER(dim)[1] = ZURAND_KEY_WORDS;
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
    if (value < 0 || value > ZURAND_MAX_EXACT_INT)
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
    if (value < -ZURAND_MAX_EXACT_INT || value > ZURAND_MAX_EXACT_INT)
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
#define ZURAND_U64_TO_DOUBLE(u) (((u) >> 11) * 0x1.0p-53)

static int lemire_accept(uint32_t x, uint32_t range, uint32_t threshold,
                         uint32_t *offset) {
    uint64_t product = (uint64_t)x * (uint64_t)range;
    if ((uint32_t)product < threshold)
        return 0;
    *offset = (uint32_t)(product >> 32);
    return 1;
}

/* ---- engine-parameterised sampler cores ----
 *
 * Both engines are Random123 generators, so both give counter-based random
 * access: a value is a pure function of (key, counter) and any index can be
 * reached directly without iterating. That property is zurand's whole point
 * and is preserved whichever engine a key names.
 *
 * Round counts are the lowest that Random123's authors report passing
 * BigCrush, not the library defaults, which carry an extra safety margin:
 * philox4x64 defaults to 10 and threefry4x64 to 20. philox stays at 10
 * because it is the shipped default stream and must never change; threefry
 * is new, so it takes 13. Both round counts have published KAT vectors,
 * which tests/testthat/test-kat.R checks against.
 */
#define ZURAND_CHUNK_BLOCKS 128 /* 512 words, 4 KiB per thread */
#define ZURAND_CHUNK_WORDS (ZURAND_CHUNK_BLOCKS * 4)

#define ZE_SUFFIX philox
#define ZE_KEY_T  philox4x64_key_t
#define ZE_GEN(c, k) philox4x64_R(10, (c), (k))
#include "zurand_engine.h"

#define ZE_SUFFIX threefry
#define ZE_KEY_T  threefry4x64_key_t
#define ZE_GEN(c, k) threefry4x64_R(13, (c), (k))
#include "zurand_engine.h"

/* ---- xoshiro256++, seeded per chunk by Philox ----
 *
 * The counter-based engines recompute a value from its counter, which
 * costs 2.5-3 ns per word however the fill is arranged. This one pays that
 * once per 512-word chunk to derive a 256-bit xoshiro state, then runs a
 * recurrence of about five operations per word. Measured in this package's
 * own fill it is roughly twice philox on uniform and 1.7x on normal.
 *
 * What it keeps, which is everything the API actually promises:
 *
 *   - chunk c depends only on (key, purpose, c), so the same key gives the
 *     same values, draw i does not depend on n, and a threaded fill is
 *     bit-identical to a serial one because threads split on chunks.
 *   - rng_fold(), key vectors and substreams are untouched: Philox still
 *     derives all of them, and the ziggurat and rejection retry paths
 *     still draw from Philox at (index, attempt), so a retried draw is
 *     still a pure function of (key, index).
 *
 * What it gives up is O(1) access to an arbitrary index -- reaching word w
 * inside a chunk costs up to 511 recurrence steps. No exported function
 * offers indexed access, so nothing user-visible changes; the philox
 * engine remains for anyone who wants that property held open.
 */
R123_STATIC_INLINE uint64_t zurand_rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static void chunk_words_xoshiro(philox4x64_key_t key, uint64_t c,
                                uint64_t purpose, uint64_t *buf, int nwords) {
    philox4x64_ctr_t ctr = {{c, 0, purpose, 0}};
    philox4x64_ctr_t st = philox4x64_R(10, ctr, key);
    uint64_t s0 = st.v[0], s1 = st.v[1], s2 = st.v[2], s3 = st.v[3];
    /* xoshiro cannot leave the all-zero state. Philox reaching it has
     * probability 2^-256, but the branch is once per chunk and free. */
    if ((s0 | s1 | s2 | s3) == 0) s0 = 1;
    for (int j = 0; j < nwords; j++) {
        uint64_t r = zurand_rotl64(s0 + s3, 23) + s0;   /* the "++" scrambler */
        uint64_t t = s1 << 17;
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t;
        s3 = zurand_rotl64(s3, 45);
        buf[j] = r;
    }
}

/* Philox for the key type and for the retry paths; xoshiro for bulk words. */
#define ZE_SUFFIX xoshiro
#define ZE_KEY_T  philox4x64_key_t
#define ZE_GEN(c, k) philox4x64_R(10, (c), (k))
#define ZE_CUSTOM_CHUNK 1
#include "zurand_engine.h"


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
    zurand_engine_t eng = key_engine_code(key);
    uint64_t h = hash_data(data);
    uint64_t domain = h ^ UINT64_C(0x9e3779b97f4a7c15);
    SEXP engine = key_engine(key);

    const int *kw = INTEGER(key);
    SEXP ans = PROTECT(alloc_key_vector(nkey, engine));
    for (R_xlen_t i = 0; i < nkey; i++) {
        zurand_ctr_t out;
        if (eng == ZURAND_ENG_PHILOX)
            out = zurand_block_philox(key_from_words(kw, nkey, i), h, domain,
                                      ZURAND_PURPOSE_FOLD);
        else
            out = zurand_block_threefry(key_from_words_threefry(kw, nkey, i), h,
                                        domain, ZURAND_PURPOSE_FOLD);
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
        return Rf_ScalarInteger(zurand_threads());

    double value = numeric_scalar(threads_, "threads");
    if (value != trunc(value) || value < 0)
        Rf_error("`threads` must be a single whole number that is at least 1, "
                 "or 0 to remove the cap");
    int prev = zurand_thread_cap;
    zurand_thread_cap = value > (double)INT_MAX ? INT_MAX : (int)value;
    return Rf_ScalarInteger(prev);
}

/* ---- .Call entry points: samplers ---- */

SEXP C_rng_uniform(SEXP key, SEXP n_, SEXP min_, SEXP max_) {
    R_xlen_t nkey = key_count(key);
    zurand_engine_t eng = key_engine_code(key);
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
    int nt = zurand_threads();
    if (span == 0.0) {
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, min) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = min;
        UNPROTECT(1);
        return ans;
    }
    int par_cols = nt > 1 && nkey > 1 && total >= ZURAND_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= ZURAND_OMP_MIN_VALUES ? nt : 0;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, par_rows, eng) schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        if (eng == ZURAND_ENG_THREEFRY)
            fill_uniform_column_threefry(out + n * col, n,
                                         key_from_words_threefry(kw, nkey, col),
                                         par_rows);
        else if (eng == ZURAND_ENG_XOSHIRO)
            fill_uniform_column_xoshiro(out + n * col, n,
                                        key_from_words(kw, nkey, col), par_rows);
        else
            fill_uniform_column_philox(out + n * col, n,
                                       key_from_words(kw, nkey, col), par_rows);
    }
    if (min != 0.0 || span != 1.0) {
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
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
    zurand_engine_t eng = key_engine_code(key);
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
    int nt = zurand_threads();
    if (sd == 0.0) {
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, mean) schedule(static)
#endif
        for (R_xlen_t i = 0; i < total; i++)
            out[i] = mean;
        UNPROTECT(1);
        return ans;
    }
    int par_cols = nt > 1 && nkey > 1 && total >= ZURAND_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= ZURAND_OMP_MIN_VALUES ? nt : 0;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, par_rows, eng) schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        if (eng == ZURAND_ENG_THREEFRY)
            fill_normal_column_threefry(out + n * col, n,
                                        key_from_words_threefry(kw, nkey, col),
                                        par_rows);
        else if (eng == ZURAND_ENG_XOSHIRO)
            fill_normal_column_xoshiro(out + n * col, n,
                                       key_from_words(kw, nkey, col), par_rows);
        else
            fill_normal_column_philox(out + n * col, n,
                                      key_from_words(kw, nkey, col), par_rows);
    }
    if (mean != 0.0 || sd != 1.0) {
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
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

SEXP C_rng_integer(SEXP key, SEXP n_, SEXP min_, SEXP max_) {
    R_xlen_t nkey = key_count(key);
    zurand_engine_t eng = key_engine_code(key);
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
    int nt = zurand_threads();
    int par_cols = nt > 1 && nkey > 1 && total >= ZURAND_OMP_MIN_VALUES;
    int par_rows = !par_cols && n >= ZURAND_OMP_MIN_VALUES ? nt : 0;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(par_cols) num_threads(nt) default(none) \
    shared(kw, out, n, nkey, min, range, threshold, par_rows, eng) \
    schedule(static)
#endif
    for (R_xlen_t col = 0; col < nkey; col++) {
        if (eng == ZURAND_ENG_THREEFRY)
            fill_integer_column_threefry(out + n * col, n,
                                         key_from_words_threefry(kw, nkey, col),
                                         min, range, threshold, par_rows);
        else if (eng == ZURAND_ENG_XOSHIRO)
            fill_integer_column_xoshiro(out + n * col, n,
                                        key_from_words(kw, nkey, col),
                                        min, range, threshold, par_rows);
        else
            fill_integer_column_philox(out + n * col, n,
                                       key_from_words(kw, nkey, col),
                                       min, range, threshold, par_rows);
    }
    UNPROTECT(1);
    return ans;
}

/* rng_bits() exposes the raw stream, so it must follow whichever engine
 * the key names. xoshiro produces words a chunk at a time, so one block of
 * four is taken from the front of its chunk. */
static void bits_block(zurand_engine_t eng, philox4x64_key_t kp,
                       threefry4x64_key_t kt, uint64_t b, zurand_ctr_t *out) {
    if (eng == ZURAND_ENG_THREEFRY) {
        *out = zurand_block_threefry(kt, b, 0, ZURAND_PURPOSE_BITS);
    } else if (eng == ZURAND_ENG_XOSHIRO) {
        uint64_t buf[ZURAND_CHUNK_WORDS];
        uint64_t c = b / ZURAND_CHUNK_BLOCKS;
        int off = (int)((b % ZURAND_CHUNK_BLOCKS) * 4);
        chunk_words_xoshiro(kp, c, ZURAND_PURPOSE_BITS, buf, off + 4);
        memcpy(out->v, buf + off, sizeof out->v);
    } else {
        *out = zurand_block_philox(kp, b, 0, ZURAND_PURPOSE_BITS);
    }
}

SEXP C_rng_bits(SEXP key, SEXP n_, SEXP bits_) {
    R_xlen_t nkey = key_count(key);
    zurand_engine_t eng = key_engine_code(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    int bits = bits_scalar(bits_);
    const int *kw = INTEGER(key);

    if (bits == 32) {
        SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
        set_sample_dim(ans, n, nkey);
        double *out = REAL(ans);
        for (R_xlen_t col = 0; col < nkey; col++) {
            philox4x64_key_t kp = key_from_words(kw, nkey, col);
            threefry4x64_key_t kt = key_from_words_threefry(kw, nkey, col);
            double *col_out = out + n * col;
            for (R_xlen_t i = 0; i < n; i += 4) {
                zurand_ctr_t block;
                    bits_block(eng, kp, kt, (uint64_t)(i >> 2), &block);
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
        philox4x64_key_t kp = key_from_words(kw, nkey, col);
        threefry4x64_key_t kt = key_from_words_threefry(kw, nkey, col);
        for (R_xlen_t i = 0; i < n; i += 4) {
            zurand_ctr_t block;
                bits_block(eng, kp, kt, (uint64_t)(i >> 2), &block);
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

void R_init_zurand(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
