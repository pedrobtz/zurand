#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>
#include "zurand.h"   /* the public C API, from inst/include */
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
 * The pragma alone is not enough. GCC 13 ignores it, GCC contracts by
 * default on every arm64 target, and both compilers ignore it under
 * -ffp-contract=fast, which a user's ~/.R/Makevars may set and which R
 * places after the package's own flags. A compiler flag of our own could
 * be overridden the same way, and R CMD check reports any -f flag in
 * src/Makevars as non-portable. So every product that feeds an addition
 * goes through zurand_rounded(), which no compiler can see through, and
 * the pragma stays as a second layer.
 *
 * Checking it: `clang -O2 -mfma -ffp-contract=fast` (x86_64, i386, and
 * --target=arm64-apple-macos) must emit no vfmadd/vfmsub/fmadd/fmla at all;
 * that is how the call sites were found. -ffast-math is outside the
 * contract: it licenses reassociation as well as contraction, and still
 * fuses inside the ziggurat's exp() argument. */
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 14)
#pragma STDC FP_CONTRACT OFF
#endif

/* Returns x unchanged, after the optimiser has lost track of how it was
 * computed: the product is rounded to double here and cannot be fused into
 * the addition that follows. The empty asm emits no instruction; it only
 * pins x in a vector register. The fallback's volatile store also rounds
 * away x87 excess precision on 32-bit x86. */
static inline double zurand_rounded(double x) {
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__SSE2_MATH__))
    __asm__("" : "+x"(x));
#elif defined(__GNUC__) && defined(__aarch64__)
    __asm__("" : "+w"(x));
#else
    volatile double v = x;
    x = v;
#endif
    return x;
}

/* The same barrier on two doubles at once, so the whole-array scaling
 * passes stay two-wide SIMD (SSE2 / NEON). With the scalar barrier the
 * compiler can no longer vectorise them. */
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__SSE2_MATH__) || \
                          defined(__aarch64__))
#  define ZURAND_HAVE_V2D 1
typedef double zurand_v2d __attribute__((vector_size(16)));
static inline zurand_v2d zurand_rounded2(zurand_v2d x) {
#  if defined(__aarch64__)
    __asm__("" : "+w"(x));
#  else
    __asm__("" : "+x"(x));
#  endif
    return x;
}
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

/* Both engines produce the same 256-bit block type; only their key types
 * differ. The sampler core is written against this. */
typedef philox4x64_ctr_t zurand_ctr_t;

typedef enum { ZURAND_ENG_PHILOX = ZURAND_ENGINE_PHILOX4X64,
               ZURAND_ENG_THREEFRY = ZURAND_ENGINE_THREEFRY4X64,
               ZURAND_ENG_XOSHIRO = ZURAND_ENGINE_XOSHIRO256PP } zurand_engine_t;
#define ZURAND_MAX_EXACT_INT 9007199254740991.0

#define ZURAND_PURPOSE_BITS    ((uint64_t)0)
#define ZURAND_PURPOSE_UNIFORM ((uint64_t)1)
#define ZURAND_PURPOSE_NORMAL  ((uint64_t)2)
#define ZURAND_PURPOSE_INTEGER ((uint64_t)3)
#define ZURAND_PURPOSE_FOLD    ((uint64_t)4)

/* Counter word 3: which engine a Philox evaluation is made on behalf of.
 * The xoshiro engine is keyed by Philox keys and draws its sub-chunk
 * seeds and its retry words from Philox; without its own tag those were
 * the philox engine's own blocks, so for any key the xoshiro stream was a
 * function of the philox stream (xoshiro's first word for rng_key(42) was
 * rotl(s0 + s3, 23) + s0 of philox's first block). Threefry is a
 * different permutation and needs no tag. See dev/design.md, section 3.2. */
#define ZURAND_TAG_PHILOX      ((uint64_t)0)
#define ZURAND_TAG_XOSHIRO     ((uint64_t)1)

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

/* Stream-format version, recorded on every key next to its engine. It is
 * what lets a stream be fixed after release without breaking saved keys:
 * a bug found after v0.1.0 becomes stream 2 for new keys, and stream-1 keys
 * keep producing stream-1 values. Only stream 1 exists. A key without the
 * attribute predates it and is stream 1. See dev/design.md, section 3.1. */
#define ZURAND_STREAM_CURRENT 1

static SEXP stream_symbol(void) {
    return Rf_install("stream");
}

static int key_stream(SEXP key) {
    SEXP s = Rf_getAttrib(key, stream_symbol());
    if (s == R_NilValue)
        return 1;
    int v = NA_INTEGER;
    if (TYPEOF(s) == INTSXP && Rf_xlength(s) == 1)
        v = INTEGER(s)[0];
    else if (TYPEOF(s) == REALSXP && Rf_xlength(s) == 1 &&
             REAL(s)[0] == (double)(int)REAL(s)[0])
        v = (int)REAL(s)[0];
    if (v == NA_INTEGER || v < 1)
        Rf_error("<rng_key> has an invalid `stream` attribute");
    if (v > ZURAND_STREAM_CURRENT)
        Rf_error("<rng_key> uses stream version %d, but this version of "
                 "zurand implements streams up to %d; update zurand",
                 v, ZURAND_STREAM_CURRENT);
    return v;
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
    (void) key_stream(key);
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
static threefry4x64_key_t threefry_key_from_philox(philox4x64_key_t p) {
    threefry4x64_key_t k;
    k.v[0] = p.v[0];
    k.v[1] = p.v[1];
    k.v[2] = p.v[0] ^ UINT64_C(0x9e3779b97f4a7c15);
    k.v[3] = p.v[1] ^ UINT64_C(0xbb67ae8584caa73b);
    return k;
}

static threefry4x64_key_t key_from_words_threefry(const int *words,
                                                  R_xlen_t nkey, R_xlen_t i) {
    return threefry_key_from_philox(key_from_words(words, nkey, i));
}

static SEXP alloc_key_vector(R_xlen_t nkey, SEXP engine, int stream) {
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
    SEXP st = PROTECT(Rf_ScalarInteger(stream));
    Rf_setAttrib(ans, stream_symbol(), st);
    SEXP cls = PROTECT(Rf_mkString("rng_key"));
    Rf_classgets(ans, cls);
    UNPROTECT(4);
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
/* The counter-based chunk_words copies whole blocks, so a request for
 * nwords may write up to 3 words past it. Every chunk buffer is sized
 * ZURAND_CHUNK_WORDS + ZURAND_CHUNK_SLACK for that reason. */
#define ZURAND_CHUNK_SLACK 3
/* An engine may take a wider chunk than the default: xoshiro groups four
 * 512-word sub-chunks so four can run in parallel SIMD lanes. Every chunk
 * buffer is sized for the widest, so one constant governs stack use. */
#define ZURAND_MAX_CHUNK_WORDS (ZURAND_CHUNK_WORDS * 10)

/* exp() and log1p() for the ziggurat's slow path: fdlibm, so the result
 * does not depend on the platform's libm. See the file for why. */
#include "zurand_fdlibm.h"

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

/* Sub-chunk: the unit a single xoshiro state covers, and the unit the
 * stream is defined in. A fill asks for four of these at a time so four
 * can run in parallel lanes, but sub-chunk s is always seeded at
 * {s, 0, purpose, 0}, so global word w comes from sub-chunk w / 512 at
 * offset w % 512 whatever the grouping. Grouping is a fill decision; it
 * is not visible in the output. */
#define ZURAND_XOSHIRO_SUB ZURAND_CHUNK_WORDS
/* Sub-chunks per fill chunk. The grouping is not part of the stream -- each
 * sub-chunk is seeded from its own index -- so it only sets how much work
 * a vector path gets at once: two AVX2 calls of 4 lanes, or one NEON call
 * of 8 (four 2-lane registers, which 32 vector registers hold without
 * spilling). Was 4 until the NEON path; output did not change. */
#if defined(__aarch64__) && defined(__ARM_NEON)
/* NEON: eight sub-chunks in four registers plus two on the scalar ALUs,
 * which are otherwise idle while the vector pipes work (see below). */
#  define ZURAND_XOSHIRO_LANES 10
#else
#  define ZURAND_XOSHIRO_LANES 8
#endif
#define ZURAND_AVX2_LANES 4

R123_STATIC_INLINE void xoshiro_seed(philox4x64_key_t key, uint64_t sub,
                                     uint64_t purpose, uint64_t *st) {
    philox4x64_ctr_t ctr = {{sub, 0, purpose, ZURAND_TAG_XOSHIRO}};
    philox4x64_ctr_t r = philox4x64_R(10, ctr, key);
    st[0] = r.v[0]; st[1] = r.v[1]; st[2] = r.v[2]; st[3] = r.v[3];
    /* xoshiro cannot leave the all-zero state; Philox reaching it has
     * probability 2^-256, but the branch is once per sub-chunk. */
    if ((st[0] | st[1] | st[2] | st[3]) == 0) st[0] = 1;
}

/* One sub-chunk, scalar. The reference definition of the stream. */
static void xoshiro_sub_scalar(philox4x64_key_t key, uint64_t sub,
                               uint64_t purpose, uint64_t *buf, int nwords) {
    uint64_t s[4];
    xoshiro_seed(key, sub, purpose, s);
    for (int j = 0; j < nwords; j++) {
        uint64_t r = zurand_rotl64(s[0] + s[3], 23) + s[0];  /* "++" scrambler */
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t;
        s[3] = zurand_rotl64(s[3], 45);
        buf[j] = r;
    }
}

/* ---- optional AVX2 path ----
 *
 * Four *sub-chunks* in the four lanes of one register, each running the
 * scalar recurrence unchanged. Vectorising across sub-chunks rather than
 * within one is what makes this a pure speedup: the output is the scalar
 * output, bit for bit, so no stream depends on whether it ran.
 *
 * Interleaving lanes within a sub-chunk was measured first and rejected --
 * four lanes are sixteen live state words on sixteen registers, and it
 * cost ~20% on baseline builds while losing even under AVX2. See
 * dev/simd/ and Phase 5 of dev/roadmap.md.
 *
 * Built with a target attribute rather than a separate -mavx2 object, so
 * src/Makevars stays portable and R CMD check sees no unusual flags. */
#if defined(__x86_64__) || defined(__i386__)
#  if defined(__GNUC__) || defined(__clang__)
#    define ZURAND_X86_DISPATCH 1
#    include <immintrin.h>
#  endif
#endif

#ifdef ZURAND_X86_DISPATCH
__attribute__((target("avx2")))
static void xoshiro_group_avx2(philox4x64_key_t key, uint64_t sub0,
                               uint64_t purpose, uint64_t *buf) {
    uint64_t st[ZURAND_AVX2_LANES][4];
    for (int l = 0; l < ZURAND_AVX2_LANES; l++)
        xoshiro_seed(key, sub0 + (uint64_t)l, purpose, st[l]);

    __m256i s0 = _mm256_set_epi64x((long long)st[3][0], (long long)st[2][0],
                                   (long long)st[1][0], (long long)st[0][0]);
    __m256i s1 = _mm256_set_epi64x((long long)st[3][1], (long long)st[2][1],
                                   (long long)st[1][1], (long long)st[0][1]);
    __m256i s2 = _mm256_set_epi64x((long long)st[3][2], (long long)st[2][2],
                                   (long long)st[1][2], (long long)st[0][2]);
    __m256i s3 = _mm256_set_epi64x((long long)st[3][3], (long long)st[2][3],
                                   (long long)st[1][3], (long long)st[0][3]);
#define ZR_VROTL(x, k) _mm256_or_si256(_mm256_slli_epi64((x), (k)), \
                                       _mm256_srli_epi64((x), 64 - (k)))
    /* Four steps are buffered and transposed so each sub-chunk receives
     * contiguous stores; storing per step would scatter across four
     * 4 KiB-apart destinations. */
    for (int j = 0; j < ZURAND_XOSHIRO_SUB; j += 4) {
        __m256i v[4];
        for (int k = 0; k < 4; k++) {
            v[k] = _mm256_add_epi64(
                ZR_VROTL(_mm256_add_epi64(s0, s3), 23), s0);
            __m256i t = _mm256_slli_epi64(s1, 17);
            s2 = _mm256_xor_si256(s2, s0); s3 = _mm256_xor_si256(s3, s1);
            s1 = _mm256_xor_si256(s1, s2); s0 = _mm256_xor_si256(s0, s3);
            s2 = _mm256_xor_si256(s2, t);  s3 = ZR_VROTL(s3, 45);
        }
        __m256i t0 = _mm256_unpacklo_epi64(v[0], v[1]);
        __m256i t1 = _mm256_unpackhi_epi64(v[0], v[1]);
        __m256i t2 = _mm256_unpacklo_epi64(v[2], v[3]);
        __m256i t3 = _mm256_unpackhi_epi64(v[2], v[3]);
        _mm256_storeu_si256((__m256i *)(buf + 0 * ZURAND_XOSHIRO_SUB + j),
                            _mm256_permute2x128_si256(t0, t2, 0x20));
        _mm256_storeu_si256((__m256i *)(buf + 1 * ZURAND_XOSHIRO_SUB + j),
                            _mm256_permute2x128_si256(t1, t3, 0x20));
        _mm256_storeu_si256((__m256i *)(buf + 2 * ZURAND_XOSHIRO_SUB + j),
                            _mm256_permute2x128_si256(t0, t2, 0x31));
        _mm256_storeu_si256((__m256i *)(buf + 3 * ZURAND_XOSHIRO_SUB + j),
                            _mm256_permute2x128_si256(t1, t3, 0x31));
    }
#undef ZR_VROTL
}

/* The same four sub-chunks, converted to uniforms in registers and stored
 * straight into the output: no word buffer, and no second pass reading it
 * back. The conversion is u01_open() lane by lane -- the same shift, OR and
 * subtraction -- so the values are identical to the scalar path.
 *
 * Measured (2026-09-26, uniform, ratio to dqrng in the same bench::mark()):
 * +17% at n = 1e7 on an i5-8500B; on an EPYC 7763 900 M/s one thread and
 * 1995 M/s on four, from 569 and 1370. Non-temporal stores on top were
 * tried and dropped: slower at the R level on both machines, one thread
 * and four, because at these sizes the cost is page faults on freshly
 * allocated vectors, which streaming stores do not avoid. */
__attribute__((target("avx2")))
static void xoshiro_group_avx2_u01(philox4x64_key_t key, uint64_t sub0,
                                   uint64_t purpose, double *out) {
    uint64_t st[ZURAND_AVX2_LANES][4];
    for (int l = 0; l < ZURAND_AVX2_LANES; l++)
        xoshiro_seed(key, sub0 + (uint64_t)l, purpose, st[l]);

    __m256i s0 = _mm256_set_epi64x((long long)st[3][0], (long long)st[2][0],
                                   (long long)st[1][0], (long long)st[0][0]);
    __m256i s1 = _mm256_set_epi64x((long long)st[3][1], (long long)st[2][1],
                                   (long long)st[1][1], (long long)st[0][1]);
    __m256i s2 = _mm256_set_epi64x((long long)st[3][2], (long long)st[2][2],
                                   (long long)st[1][2], (long long)st[0][2]);
    __m256i s3 = _mm256_set_epi64x((long long)st[3][3], (long long)st[2][3],
                                   (long long)st[1][3], (long long)st[0][3]);
    const __m256i exponent = _mm256_set1_epi64x(0x3ff0000000000000LL);
    const __m256d bias = _mm256_set1_pd(1.0 - 0x1.0p-53);
#define ZR_VROTL(x, k) _mm256_or_si256(_mm256_slli_epi64((x), (k)), \
                                       _mm256_srli_epi64((x), 64 - (k)))
#define ZR_U01(w) _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256( \
                      _mm256_srli_epi64((w), 12), exponent)), bias)
#define ZR_STORE(p, v) _mm256_storeu_pd((p), (v))
    for (int j = 0; j < ZURAND_XOSHIRO_SUB; j += 4) {
        __m256i v[4];
        for (int k = 0; k < 4; k++) {
            v[k] = _mm256_add_epi64(
                ZR_VROTL(_mm256_add_epi64(s0, s3), 23), s0);
            __m256i t = _mm256_slli_epi64(s1, 17);
            s2 = _mm256_xor_si256(s2, s0); s3 = _mm256_xor_si256(s3, s1);
            s1 = _mm256_xor_si256(s1, s2); s0 = _mm256_xor_si256(s0, s3);
            s2 = _mm256_xor_si256(s2, t);  s3 = ZR_VROTL(s3, 45);
        }
        __m256i t0 = _mm256_unpacklo_epi64(v[0], v[1]);
        __m256i t1 = _mm256_unpackhi_epi64(v[0], v[1]);
        __m256i t2 = _mm256_unpacklo_epi64(v[2], v[3]);
        __m256i t3 = _mm256_unpackhi_epi64(v[2], v[3]);
        ZR_STORE(out + 0 * ZURAND_XOSHIRO_SUB + j,
                 ZR_U01(_mm256_permute2x128_si256(t0, t2, 0x20)));
        ZR_STORE(out + 1 * ZURAND_XOSHIRO_SUB + j,
                 ZR_U01(_mm256_permute2x128_si256(t1, t3, 0x20)));
        ZR_STORE(out + 2 * ZURAND_XOSHIRO_SUB + j,
                 ZR_U01(_mm256_permute2x128_si256(t0, t2, 0x31)));
        ZR_STORE(out + 3 * ZURAND_XOSHIRO_SUB + j,
                 ZR_U01(_mm256_permute2x128_si256(t1, t3, 0x31)));
    }
#undef ZR_STORE
#undef ZR_U01
#undef ZR_VROTL
}
#endif

/* What the CPU offers (-1 until probed) and whether the user has turned it
 * off. Separated so that rng_simd(TRUE) restores detection rather than
 * asserting a capability the machine may not have. */
/* ---- NEON (every aarch64 CPU has it; no runtime check) ----
 *
 * Eight sub-chunks per call, two per 128-bit register, four registers --
 * randompack's width, which on arm64 does not spill. After two steps a
 * register holds [x0 y0], [x1 y1] for sub-chunks (p, p+1); zip1/zip2 turn
 * that into [x0 x1] and [y0 y1], contiguous in each sub-chunk.
 *
 * Two further sub-chunks run on the scalar integer ALUs in the same loop:
 * an Apple core has six of them, idle while four vector pipes do the
 * work, so they add throughput almost for free.
 *
 * Rotations are SHL, USHR, ORR. SHL + SRI is one instruction shorter and
 * about 30% slower on Apple M1 (dev/simd/m1_compare.c): SRI reads its
 * destination, which lengthens the loop-carried chain of the recurrence.
 *
 * Measured on GitHub's macos-26-arm64 runner, dev/simd/neon.c: uniform
 * 2.06-2.27x the scalar path (1300-1500 M/s), and bit-identical to it. On
 * the ubuntu-24.04-arm runner (Neoverse) it is level with scalar: that
 * machine is store-bound, so the gain is Apple Silicon's. */
#if defined(__aarch64__) && defined(__ARM_NEON)
#  define ZURAND_NEON 1
#  include <arm_neon.h>

#define ZR_NROTL(x, k) vorrq_u64(vshlq_n_u64((x), (k)), vshrq_n_u64((x), 64 - (k)))
#define ZR_NSTEP(s0, s1, s2, s3, r) do {                                   \
        (r) = vaddq_u64(ZR_NROTL(vaddq_u64((s0), (s3)), 23), (s0));         \
        uint64x2_t t_ = vshlq_n_u64((s1), 17);                             \
        (s2) = veorq_u64((s2), (s0)); (s3) = veorq_u64((s3), (s1));       \
        (s1) = veorq_u64((s1), (s2)); (s0) = veorq_u64((s0), (s3));       \
        (s2) = veorq_u64((s2), t_);   (s3) = ZR_NROTL((s3), 45);           \
    } while (0)
#define ZR_NSEED(p, sub) do {                                              \
        uint64_t a_[4], b_[4];                                             \
        xoshiro_seed(key, (sub), purpose, a_);                             \
        xoshiro_seed(key, (sub) + 1, purpose, b_);                         \
        uint64_t l_[2];                                                    \
        l_[0] = a_[0]; l_[1] = b_[0]; p##0 = vld1q_u64(l_);                \
        l_[0] = a_[1]; l_[1] = b_[1]; p##1 = vld1q_u64(l_);                \
        l_[0] = a_[2]; l_[1] = b_[2]; p##2 = vld1q_u64(l_);                \
        l_[0] = a_[3]; l_[1] = b_[3]; p##3 = vld1q_u64(l_);                \
    } while (0)

#define ZR_SSTEP(s, r) do {                                                \
        (r) = zurand_rotl64((s)[0] + (s)[3], 23) + (s)[0];                  \
        uint64_t t_ = (s)[1] << 17;                                        \
        (s)[2] ^= (s)[0]; (s)[3] ^= (s)[1]; (s)[1] ^= (s)[2];              \
        (s)[0] ^= (s)[3]; (s)[2] ^= t_; (s)[3] = zurand_rotl64((s)[3], 45); \
    } while (0)

static inline float64x2_t zr_neon_u01(uint64x2_t w) {
    uint64x2_t m = vorrq_u64(vshrq_n_u64(w, 12),
                             vdupq_n_u64(UINT64_C(0x3ff0000000000000)));
    return vsubq_f64(vreinterpretq_f64_u64(m), vdupq_n_f64(1.0 - 0x1.0p-53));
}

/* Words for 10 sub-chunks into buf (sub-chunk l at buf + l * SUB). */
static void xoshiro_group_neon(philox4x64_key_t key, uint64_t sub0,
                               uint64_t purpose, uint64_t *buf) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    ZR_NSEED(a, sub0);     ZR_NSEED(b, sub0 + 2);
    ZR_NSEED(c, sub0 + 4); ZR_NSEED(d, sub0 + 6);
    uint64_t e[4], f[4];
    xoshiro_seed(key, sub0 + 8, purpose, e);
    xoshiro_seed(key, sub0 + 9, purpose, f);
#define ZR_NPAIR(p, l) do { uint64x2_t r0_, r1_;                             \
        ZR_NSTEP(p##0, p##1, p##2, p##3, r0_);                             \
        ZR_NSTEP(p##0, p##1, p##2, p##3, r1_);                             \
        vst1q_u64(buf + (l) * ZURAND_XOSHIRO_SUB + j, vzip1q_u64(r0_, r1_)); \
        vst1q_u64(buf + ((l) + 1) * ZURAND_XOSHIRO_SUB + j, vzip2q_u64(r0_, r1_)); \
    } while (0)
    uint64_t *be = buf + 8 * ZURAND_XOSHIRO_SUB, *bf = buf + 9 * ZURAND_XOSHIRO_SUB;
    for (int j = 0; j < ZURAND_XOSHIRO_SUB; j += 2) {
        ZR_SSTEP(e, be[j]); ZR_SSTEP(f, bf[j]);
        ZR_NPAIR(a, 0); ZR_NPAIR(b, 2);
        ZR_SSTEP(e, be[j + 1]); ZR_SSTEP(f, bf[j + 1]);
        ZR_NPAIR(c, 4); ZR_NPAIR(d, 6);
    }
#undef ZR_NPAIR
}

/* The same, converted to uniforms in registers and stored into `out`. */
static void xoshiro_group_neon_u01(philox4x64_key_t key, uint64_t sub0,
                                   uint64_t purpose, double *out) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    ZR_NSEED(a, sub0);     ZR_NSEED(b, sub0 + 2);
    ZR_NSEED(c, sub0 + 4); ZR_NSEED(d, sub0 + 6);
    uint64_t e[4], f[4];
    xoshiro_seed(key, sub0 + 8, purpose, e);
    xoshiro_seed(key, sub0 + 9, purpose, f);
    double *oe = out + 8 * ZURAND_XOSHIRO_SUB, *of = out + 9 * ZURAND_XOSHIRO_SUB;
#define ZR_NPAIR(p, l) do { uint64x2_t r0_, r1_;                             \
        ZR_NSTEP(p##0, p##1, p##2, p##3, r0_);                             \
        ZR_NSTEP(p##0, p##1, p##2, p##3, r1_);                             \
        vst1q_f64(out + (l) * ZURAND_XOSHIRO_SUB + j,                      \
                  zr_neon_u01(vzip1q_u64(r0_, r1_)));                      \
        vst1q_f64(out + ((l) + 1) * ZURAND_XOSHIRO_SUB + j,                \
                  zr_neon_u01(vzip2q_u64(r0_, r1_)));                      \
    } while (0)
    for (int j = 0; j < ZURAND_XOSHIRO_SUB; j += 2) {
        uint64_t e0, e1, f0, f1;
        ZR_SSTEP(e, e0); ZR_SSTEP(f, f0);
        ZR_NPAIR(a, 0); ZR_NPAIR(b, 2);
        ZR_SSTEP(e, e1); ZR_SSTEP(f, f1);
        ZR_NPAIR(c, 4); ZR_NPAIR(d, 6);
        oe[j] = u01_open(e0); oe[j + 1] = u01_open(e1);
        of[j] = u01_open(f0); of[j + 1] = u01_open(f1);
    }
#undef ZR_NPAIR
}
#undef ZR_SSTEP
#undef ZR_NSEED
#undef ZR_NSTEP
#undef ZR_NROTL
#endif /* NEON */

/* DC ZVA on arm64: zero whole cache blocks without reading them first.
 *
 * zurand's stream puts each 512-value sub-chunk in its own contiguous
 * 4 KiB block, so a NEON group writes ten blocks at once. That pattern
 * does not get the "whole line is being overwritten, skip the read" path a
 * single sequential stream gets on Apple cores, so every output line was
 * first read from memory. Zeroing the group's region with DC ZVA first
 * gives each store a line it already owns. Measured on the macos-latest
 * (Apple M1) runner, dev/simd/m1_compare.c: 23-25% less time per value on
 * large fills, reused or freshly allocated; ~20% more in cache, where the
 * lines never left anyway -- hence the size threshold. Apple only: on
 * the ubuntu-24.04-arm runner (Neoverse N2) it measured level or slower,
 * so other arm64 cores keep ordinary stores.
 *
 * Only blocks lying entirely inside the region are zeroed. R aligns vector
 * data to 16 bytes, not to the block, and a partial block at either end
 * holds bytes of the neighbouring chunk or of the vector's header; those
 * edges keep ordinary stores. Whole blocks inside one chunk are never
 * touched by another thread. The block size is read once at load (0 if
 * DC ZVA is prohibited), so worker threads only read it. */
static size_t zurand_zva_bytes = 0;

#ifndef ZURAND_ZVA_MIN_VALUES
#  define ZURAND_ZVA_MIN_VALUES ((R_xlen_t)1 << 17)   /* 1 MiB of doubles */
#endif

static void zurand_zva_init(void) {
#if defined(__APPLE__) && defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    uint64_t v;
    __asm__ volatile("mrs %0, dczid_el0" : "=r"(v));
    zurand_zva_bytes = (v & 16) ? 0 : (size_t)4 << (v & 15);
#endif
}

static inline void zurand_zva_region(double *p, size_t bytes) {
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    size_t b = zurand_zva_bytes;
    if (!b)
        return;
    uintptr_t lo = ((uintptr_t)p + b - 1) & ~(uintptr_t)(b - 1);
    uintptr_t hi = ((uintptr_t)p + bytes) & ~(uintptr_t)(b - 1);
    for (uintptr_t c = lo; c < hi; c += b)
        __asm__ volatile("dc zva, %0" : : "r"(c) : "memory");
#else
    (void)p; (void)bytes;
#endif
}

static int zurand_avx2_available = -1;
static int zurand_simd_off = 0;

static int zurand_use_neon(void) {
#ifdef ZURAND_NEON
    return !zurand_simd_off;
#else
    return 0;
#endif
}

static int zurand_use_avx2(void) {
    if (zurand_avx2_available < 0) {
#ifdef ZURAND_X86_DISPATCH
        zurand_avx2_available = __builtin_cpu_supports("avx2") ? 1 : 0;
#else
        zurand_avx2_available = 0;
#endif
    }
    return zurand_avx2_available && !zurand_simd_off;
}

/* Querying returns the active path; setting FALSE forces the scalar one.
 * The point of exposing this is that "both paths emit identical bits" is
 * the whole basis for dispatching at all, and a claim nobody can check is
 * not worth making -- tests/testthat/test-simd.R checks it. */
SEXP C_rng_simd(SEXP enable) {
    const char *prev = zurand_use_avx2() ? "avx2"
                     : zurand_use_neon() ? "neon" : "none";
    if (enable != R_NilValue) {
        int e = Rf_asLogical(enable);
        if (e == NA_LOGICAL)
            Rf_error("`enable` must be TRUE or FALSE");
        zurand_simd_off = !e;
    }
    return Rf_mkString(prev);
}

static void chunk_words_xoshiro(philox4x64_key_t key, uint64_t c,
                                uint64_t purpose, uint64_t *buf, int nwords) {
    uint64_t sub0 = c * ZURAND_XOSHIRO_LANES;
    if (nwords == ZURAND_XOSHIRO_SUB * ZURAND_XOSHIRO_LANES) {
#ifdef ZURAND_X86_DISPATCH
        if (zurand_use_avx2()) {
            xoshiro_group_avx2(key, sub0, purpose, buf);
            xoshiro_group_avx2(key, sub0 + ZURAND_AVX2_LANES, purpose,
                               buf + ZURAND_AVX2_LANES * ZURAND_XOSHIRO_SUB);
            return;
        }
#endif
#ifdef ZURAND_NEON
        if (zurand_use_neon()) {
            xoshiro_group_neon(key, sub0, purpose, buf);
            return;
        }
#endif
    }
    /* Scalar: the same sub-chunks, one at a time. Also the path taken for
     * the final short chunk of a fill, where the group is not full. */
    for (int l = 0; l < ZURAND_XOSHIRO_LANES; l++) {
        int have = nwords - l * ZURAND_XOSHIRO_SUB;
        if (have <= 0) break;
        xoshiro_sub_scalar(key, sub0 + (uint64_t)l, purpose,
                           buf + l * ZURAND_XOSHIRO_SUB,
                           have < ZURAND_XOSHIRO_SUB ? have : ZURAND_XOSHIRO_SUB);
    }
}

/* One full xoshiro chunk of uniforms, written straight into `o` when the
 * AVX2 path is on; returns 0 to fall back to words-then-convert. */
static int uniform_chunk_xoshiro_fast(philox4x64_key_t key, uint64_t c,
                                      double *o, int m, int large) {
    if (m != ZURAND_XOSHIRO_SUB * ZURAND_XOSHIRO_LANES)
        return 0;
    uint64_t sub0 = c * ZURAND_XOSHIRO_LANES;
#ifdef ZURAND_X86_DISPATCH
    if (zurand_use_avx2()) {
        xoshiro_group_avx2_u01(key, sub0, ZURAND_PURPOSE_UNIFORM, o);
        xoshiro_group_avx2_u01(key, sub0 + ZURAND_AVX2_LANES,
                               ZURAND_PURPOSE_UNIFORM,
                               o + ZURAND_AVX2_LANES * ZURAND_XOSHIRO_SUB);
        return 1;
    }
#endif
#ifdef ZURAND_NEON
    if (zurand_use_neon()) {
        if (large)
            zurand_zva_region(o, (size_t)m * sizeof(double));
        xoshiro_group_neon_u01(key, sub0, ZURAND_PURPOSE_UNIFORM, o);
        return 1;
    }
#endif
    (void)key; (void)o; (void)sub0; (void)large;
    return 0;
}

/* Philox for the key type and for the retry paths; xoshiro for bulk words. */
#define ZE_UNIFORM_FAST(key, c, o, m, large) \
    uniform_chunk_xoshiro_fast((key), (c), (o), (m), (large))
#define ZE_SUFFIX xoshiro
#define ZE_KEY_T  philox4x64_key_t
#define ZE_CHUNK_WORDS (ZURAND_CHUNK_WORDS * ZURAND_XOSHIRO_LANES)
#define ZE_GEN(c, k) philox4x64_R(10, (c), (k))
#define ZE_TAG ZURAND_TAG_XOSHIRO
#define ZE_CUSTOM_CHUNK 1
#include "zurand_engine.h"

/* ---- per-key fills, shared by the .Call entry points and the C API ----
 *
 * One key, one output buffer, no R API, so they are safe on any thread.
 * par_rows > 1 splits one long fill over OpenMP threads, which the R
 * samplers use; the C API passes 0 and leaves threading to its caller. */

static void fill_uniform_key(zurand_engine_t eng, philox4x64_key_t kp,
                             double *out, R_xlen_t n, int par_rows) {
    if (eng == ZURAND_ENG_THREEFRY)
        fill_uniform_column_threefry(out, n, threefry_key_from_philox(kp), par_rows);
    else if (eng == ZURAND_ENG_XOSHIRO)
        fill_uniform_column_xoshiro(out, n, kp, par_rows);
    else
        fill_uniform_column_philox(out, n, kp, par_rows);
}

static void fill_normal_key(zurand_engine_t eng, philox4x64_key_t kp,
                            double *out, R_xlen_t n, int par_rows) {
    if (eng == ZURAND_ENG_THREEFRY)
        fill_normal_column_threefry(out, n, threefry_key_from_philox(kp), par_rows);
    else if (eng == ZURAND_ENG_XOSHIRO)
        fill_normal_column_xoshiro(out, n, kp, par_rows);
    else
        fill_normal_column_philox(out, n, kp, par_rows);
}

static void fill_integer_key(zurand_engine_t eng, philox4x64_key_t kp,
                             int *out, R_xlen_t n, int min, uint32_t range,
                             uint32_t threshold, int par_rows) {
    if (eng == ZURAND_ENG_THREEFRY)
        fill_integer_column_threefry(out, n, threefry_key_from_philox(kp),
                                     min, range, threshold, par_rows);
    else if (eng == ZURAND_ENG_XOSHIRO)
        fill_integer_column_xoshiro(out, n, kp, min, range, threshold, par_rows);
    else
        fill_integer_column_philox(out, n, kp, min, range, threshold, par_rows);
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

/* Final avalanche (MurmurHash3's fmix64), shared by rng_fold() and the C
 * API's fold_int() so the two cannot derive different keys. */
static uint64_t hash_finish(uint64_t h) {
    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
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

    return hash_finish(h);
}

/* out[i] = a + s * out[i], with the product rounded before the add on
 * every compiler (see zurand_rounded). Shared by the min/max and mean/sd
 * passes; elementwise, so the threaded result equals the serial one. */
static void affine_pass(double *out, R_xlen_t total, double a, double s,
                        int nt) {
    (void)nt;  /* the thread count, used only by the OpenMP pragmas */
#ifdef ZURAND_HAVE_V2D
    const R_xlen_t npair = total / 2;
    const zurand_v2d av = {a, a}, sv = {s, s};
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, npair, av, sv) schedule(static)
#endif
    for (R_xlen_t p = 0; p < npair; p++) {
        zurand_v2d x;
        memcpy(&x, out + 2 * p, sizeof x);
        x = av + zurand_rounded2(sv * x);
        memcpy(out + 2 * p, &x, sizeof x);
    }
    if (total & 1)
        out[total - 1] = a + zurand_rounded(s * out[total - 1]);
#else
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(nt > 1 && total >= ZURAND_OMP_MIN_VALUES) \
    num_threads(nt) default(none) shared(out, total, a, s) schedule(static)
#endif
    for (R_xlen_t i = 0; i < total; i++)
        out[i] = a + zurand_rounded(s * out[i]);
#endif
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

    SEXP ans = PROTECT(alloc_key_vector(n, engine, ZURAND_STREAM_CURRENT));
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

    SEXP ans = PROTECT(alloc_key_vector(n, engine, ZURAND_STREAM_CURRENT));
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
    SEXP ans = PROTECT(alloc_key_vector(nkey, engine, key_stream(key)));
    for (R_xlen_t i = 0; i < nkey; i++) {
        zurand_ctr_t out;
        /* xoshiro keys are Philox keys -- Philox derives their chunk seeds
         * and their retry words -- so they fold with Philox too. The earlier
         * two-way branch sent them down the threefry path with a
         * threefry-expanded key, which contradicted the documentation and
         * tied xoshiro's substreams to another engine's key schedule. */
        if (eng == ZURAND_ENG_THREEFRY)
            out = zurand_block_threefry(key_from_words_threefry(kw, nkey, i), h,
                                        domain, ZURAND_PURPOSE_FOLD);
        else
            out = zurand_block_philox(key_from_words(kw, nkey, i), h, domain,
                                      ZURAND_PURPOSE_FOLD);
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
    for (R_xlen_t col = 0; col < nkey; col++)
        fill_uniform_key(eng, key_from_words(kw, nkey, col), out + n * col, n,
                         par_rows);
    if (min != 0.0 || span != 1.0)
        affine_pass(out, total, min, span, nt);
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
    for (R_xlen_t col = 0; col < nkey; col++)
        fill_normal_key(eng, key_from_words(kw, nkey, col), out + n * col, n,
                        par_rows);
    if (mean != 0.0 || sd != 1.0)
        affine_pass(out, total, mean, sd, nt);
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
    for (R_xlen_t col = 0; col < nkey; col++)
        fill_integer_key(eng, key_from_words(kw, nkey, col), out + n * col, n,
                         min, range, threshold, par_rows);
    UNPROTECT(1);
    return ans;
}

/* Words per chunk for an engine, mirroring ZE_CHUNK_WORDS inside the
 * engine header. Code outside the header -- rng_bits() -- has to agree
 * with it: passing a 512-word chunk index to an engine that groups four
 * of them silently reads the wrong sub-chunks. */
static R_xlen_t engine_chunk_words_raw(zurand_engine_t eng) {
    return eng == ZURAND_ENG_XOSHIRO
        ? (R_xlen_t)ZURAND_CHUNK_WORDS * ZURAND_XOSHIRO_LANES
        : (R_xlen_t)ZURAND_CHUNK_WORDS;
}

static R_xlen_t engine_chunk_words(zurand_engine_t eng) {
    R_xlen_t w = engine_chunk_words_raw(eng);
    /* Every caller writes this many words into a buffer sized by
     * ZURAND_MAX_CHUNK_WORDS. Leaving one at the narrow width while this
     * returned the wide one overflowed the stack by 1533 words, which the
     * test suite ran straight past and only R CMD check caught, as an
     * abort with no output. */
    if (w > ZURAND_MAX_CHUNK_WORDS)
        Rf_error("internal: chunk width %lld exceeds buffer capacity %d",
                 (long long)w, ZURAND_MAX_CHUNK_WORDS);
    return w;
}

/* rng_bits() exposes the raw stream and must follow whichever engine the
 * key names, a chunk at a time. The previous version fetched one block per
 * call, which for xoshiro meant regenerating from the chunk start every
 * four words -- 9x slower than the other engines, and the statistical
 * audit's `bits` sampler goes through exactly this path. */
static void bits_chunk(zurand_engine_t eng, philox4x64_key_t kp,
                       threefry4x64_key_t kt, uint64_t c, uint64_t *buf,
                       int nwords) {
    if (eng == ZURAND_ENG_THREEFRY)
        chunk_words_threefry(kt, c, ZURAND_PURPOSE_BITS, buf, nwords);
    else if (eng == ZURAND_ENG_XOSHIRO)
        chunk_words_xoshiro(kp, c, ZURAND_PURPOSE_BITS, buf, nwords);
    else
        chunk_words_philox(kp, c, ZURAND_PURPOSE_BITS, buf, nwords);
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
            R_xlen_t cw = engine_chunk_words(eng);
            R_xlen_t nchunk = (n + cw - 1) / cw;
            for (R_xlen_t c = 0; c < nchunk; c++) {
                uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
                R_xlen_t w0 = c * cw;
                int m = (int)(n - w0 < cw ? n - w0 : cw);
                bits_chunk(eng, kp, kt, (uint64_t)c, buf, m);
                for (int j = 0; j < m; j++)
                    col_out[w0 + j] = (double)(uint32_t)buf[j];
            }
        }
        UNPROTECT(1);
        return ans;
    }

    SEXP ans = PROTECT(Rf_allocVector(STRSXP, total));
    set_sample_dim(ans, n, nkey);
    char hex[17];
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t kp = key_from_words(kw, nkey, col);
        threefry4x64_key_t kt = key_from_words_threefry(kw, nkey, col);
        R_xlen_t cw = engine_chunk_words(eng);
            R_xlen_t nchunk = (n + cw - 1) / cw;
        for (R_xlen_t c = 0; c < nchunk; c++) {
            uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
            R_xlen_t w0 = c * cw;
            int m = (int)(n - w0 < cw ? n - w0 : cw);
            bits_chunk(eng, kp, kt, (uint64_t)c, buf, m);
            for (int j = 0; j < m; j++) {
                snprintf(hex, sizeof hex, "%016" PRIx64, buf[j]);
                SET_STRING_ELT(ans, w0 + j + n * col, Rf_mkChar(hex));
            }
        }
    }
    UNPROTECT(1);
    return ans;
}

/* ---- C API (inst/include/zurand.h) ----
 *
 * Everything below except api_key_get() runs without the R API, so the
 * fills can be called from a caller's worker threads. Arguments are checked
 * as the R samplers check them, but reported as return codes: an R error
 * raised on a worker thread would take the process down. */

/* The API's engine codes are the internal ones; zurand_engine_t is defined
 * from the header's constants, so they cannot drift apart. */

static int api_key_ok(zurand_key k) {
    return (k.engine == ZURAND_ENG_PHILOX || k.engine == ZURAND_ENG_THREEFRY ||
            k.engine == ZURAND_ENG_XOSHIRO) &&
           k.stream >= 1 && k.stream <= ZURAND_STREAM_CURRENT;
}

static philox4x64_key_t api_philox_key(zurand_key k) {
    philox4x64_key_t p;
    p.v[0] = k.k0;
    p.v[1] = k.k1;
    return p;
}

static int api_key_get(SEXP keys, R_xlen_t i, zurand_key *out) {
    R_xlen_t nkey = key_count(keys);          /* R errors if not a key vector */
    if (i < 0 || i >= nkey)
        return ZURAND_EINVAL;
    philox4x64_key_t p = key_from_words(INTEGER(keys), nkey, i);
    out->k0 = p.v[0];
    out->k1 = p.v[1];
    out->engine = (int)key_engine_code(keys);
    out->stream = key_stream(keys);
    return ZURAND_OK;
}

/* hash_data() for a length-one integer or double vector holding `data`. */
static int api_fold_int(zurand_key k, int64_t data, zurand_key *out) {
    if (!api_key_ok(k))
        return ZURAND_EKEY;
    uint64_t h = UINT64_C(1469598103934665603);
    hash_byte(&h, 2);
    hash_u64(&h, 1);
    hash_u64(&h, (uint64_t)data);
    h = hash_finish(h);
    uint64_t domain = h ^ UINT64_C(0x9e3779b97f4a7c15);
    philox4x64_key_t p = api_philox_key(k);
    zurand_ctr_t r = k.engine == ZURAND_ENG_THREEFRY
        ? zurand_block_threefry(threefry_key_from_philox(p), h, domain,
                                ZURAND_PURPOSE_FOLD)
        : zurand_block_philox(p, h, domain, ZURAND_PURPOSE_FOLD);
    out->k0 = r.v[0];
    out->k1 = r.v[1];
    out->engine = k.engine;
    out->stream = k.stream;
    return ZURAND_OK;
}

static int api_fill_uniform(zurand_key k, size_t n, double min, double max,
                            double *out) {
    if (!api_key_ok(k))
        return ZURAND_EKEY;
    if (!isfinite(min) || !isfinite(max) || min > max ||
        n > (size_t)R_XLEN_T_MAX)
        return ZURAND_EINVAL;
    R_xlen_t nn = (R_xlen_t)n;
    double span = max - min;
    if (span == 0.0) {
        for (R_xlen_t i = 0; i < nn; i++)
            out[i] = min;
        return ZURAND_OK;
    }
    fill_uniform_key((zurand_engine_t)k.engine, api_philox_key(k), out, nn, 0);
    if (min != 0.0 || span != 1.0)
        affine_pass(out, nn, min, span, 1);
    return ZURAND_OK;
}

static int api_fill_normal(zurand_key k, size_t n, double mean, double sd,
                           double *out) {
    if (!api_key_ok(k))
        return ZURAND_EKEY;
    if (!isfinite(mean) || !isfinite(sd) || sd < 0 || n > (size_t)R_XLEN_T_MAX)
        return ZURAND_EINVAL;
    R_xlen_t nn = (R_xlen_t)n;
    if (sd == 0.0) {
        for (R_xlen_t i = 0; i < nn; i++)
            out[i] = mean;
        return ZURAND_OK;
    }
    fill_normal_key((zurand_engine_t)k.engine, api_philox_key(k), out, nn, 0);
    if (mean != 0.0 || sd != 1.0)
        affine_pass(out, nn, mean, sd, 1);
    return ZURAND_OK;
}

static int api_fill_integer(zurand_key k, size_t n, int min, int max,
                            int *out) {
    if (!api_key_ok(k))
        return ZURAND_EKEY;
    if (min == NA_INTEGER || max == NA_INTEGER || min > max ||
        n > (size_t)R_XLEN_T_MAX)
        return ZURAND_EINVAL;
    uint32_t range = (uint32_t)((int64_t)max - (int64_t)min + 1);
    uint32_t threshold = (uint32_t)((UINT64_C(0x100000000) - range) % range);
    fill_integer_key((zurand_engine_t)k.engine, api_philox_key(k), out,
                     (R_xlen_t)n, min, range, threshold, 0);
    return ZURAND_OK;
}

static int api_fill_bits64(zurand_key k, size_t n, uint64_t *out) {
    if (!api_key_ok(k))
        return ZURAND_EKEY;
    if (n > (size_t)R_XLEN_T_MAX)
        return ZURAND_EINVAL;
    zurand_engine_t eng = (zurand_engine_t)k.engine;
    philox4x64_key_t kp = api_philox_key(k);
    threefry4x64_key_t kt = threefry_key_from_philox(kp);
    R_xlen_t nn = (R_xlen_t)n;
    R_xlen_t cw = engine_chunk_words_raw(eng);
    R_xlen_t nchunk = (nn + cw - 1) / cw;
    for (R_xlen_t c = 0; c < nchunk; c++) {
        uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
        R_xlen_t w0 = c * cw;
        int m = (int)(nn - w0 < cw ? nn - w0 : cw);
        bits_chunk(eng, kp, kt, (uint64_t)c, buf, m);
        memcpy(out + w0, buf, (size_t)m * sizeof(uint64_t));
    }
    return ZURAND_OK;
}

static const zurand_api zurand_api_table = {
    ZURAND_API_VERSION,
    sizeof(zurand_api),
    api_key_get,
    api_fold_int,
    api_fill_uniform,
    api_fill_normal,
    api_fill_integer,
    api_fill_bits64
};

static const zurand_api *zurand_api_get(void) {
    return &zurand_api_table;
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
    {"C_rng_simd",       (DL_FUNC) &C_rng_simd,       1},
    {NULL, NULL, 0}
};

void R_init_zurand(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_RegisterCCallable("zurand", "zurand_api", (DL_FUNC) &zurand_api_get);
    /* Settle the SIMD dispatch now, on the main thread, so C API calls from
     * worker threads only ever read it. */
    (void) zurand_use_avx2();
    zurand_zva_init();
}
