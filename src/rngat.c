#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "Random123/philox.h"

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

/* ---- key representation ---- */

static SEXP engine_symbol(void) {
    return Rf_install("engine");
}

static void check_engine(SEXP engine) {
    if (TYPEOF(engine) != STRSXP || Rf_xlength(engine) != 1 ||
        strcmp(CHAR(STRING_ELT(engine, 0)), RNGAT_ENGINE) != 0)
        Rf_error("`engine` must be \"%s\"", RNGAT_ENGINE);
}

static uint32_t load_word_at(SEXP x, R_xlen_t nkey, R_xlen_t i, int word) {
    uint32_t out;
    memcpy(&out, INTEGER(x) + i + nkey * (R_xlen_t)word, sizeof out);
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

static philox4x64_key_t key_at(SEXP key, R_xlen_t i) {
    R_xlen_t nkey = key_count(key);
    if (i < 0 || i >= nkey)
        Rf_error("internal key index out of range");

    philox4x64_key_t k;
    uint64_t w0 = (uint64_t)load_word_at(key, nkey, i, 0);
    uint64_t w1 = (uint64_t)load_word_at(key, nkey, i, 1);
    uint64_t w2 = (uint64_t)load_word_at(key, nkey, i, 2);
    uint64_t w3 = (uint64_t)load_word_at(key, nkey, i, 3);
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
    Rf_classgets(ans, Rf_mkString("rng_key"));
    UNPROTECT(2);
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
    SEXP ans = PROTECT(Rf_allocVector(STRSXP, nkey));
    char buf[42];
    for (R_xlen_t i = 0; i < nkey; i++) {
        snprintf(buf, sizeof buf, "rng_key[%08" PRIx32 "%08" PRIx32 "%08" PRIx32 "%08" PRIx32 "]",
                 load_word_at(key, nkey, i, 0), load_word_at(key, nkey, i, 1),
                 load_word_at(key, nkey, i, 2), load_word_at(key, nkey, i, 3));
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

/* 53 random bits -> double strictly inside (0, 1) */
static double u01_open(uint64_t bits) {
    return ((double)(bits >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/* ---- stable fold hashing ---- */

static void hash_byte(uint64_t *h, unsigned char byte) {
    *h ^= (uint64_t)byte;
    *h *= UINT64_C(1099511628211);
}

static void hash_bytes(uint64_t *h, const unsigned char *bytes, size_t n) {
    for (size_t i = 0; i < n; i++) hash_byte(h, bytes[i]);
}

static void hash_u32(uint64_t *h, uint32_t x) {
    for (int i = 0; i < 4; i++) {
        hash_byte(h, (unsigned char)(x & 0xffu));
        x >>= 8;
    }
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

    case INTSXP:
        hash_byte(&h, 2);
        hash_u64(&h, (uint64_t)n);
        for (R_xlen_t i = 0; i < n; i++) {
            int value = INTEGER(data)[i];
            if (value == NA_INTEGER)
                Rf_error("`data` must not contain missing values");
            hash_u32(&h, (uint32_t)value);
        }
        break;

    case REALSXP:
        hash_byte(&h, 3);
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

    SEXP ans = PROTECT(alloc_key_vector(nkey, engine));
    for (R_xlen_t i = 0; i < nkey; i++) {
        philox4x64_key_t k = key_at(key, i);
        philox4x64_ctr_t out = rngat_block(k, h, domain, RNGAT_PURPOSE_FOLD);
        set_key_words(ans, nkey, i, out.v[0], out.v[1]);
    }
    UNPROTECT(1);
    return ans;
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

    SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
    set_sample_dim(ans, n, nkey);
    double *out = REAL(ans);
    double span = max - min;
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_at(key, col);
        for (R_xlen_t i = 0; i < n; i++) {
            uint64_t bits = rngat_word(k, (uint64_t)i, 0, RNGAT_PURPOSE_UNIFORM);
            out[i + n * col] = min + span * u01_open(bits);
        }
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

    SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
    set_sample_dim(ans, n, nkey);
    double *out = REAL(ans);
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_at(key, col);
        for (R_xlen_t i = 0; i < n; i++) {
            uint64_t bits = rngat_word(k, (uint64_t)i, 0, RNGAT_PURPOSE_NORMAL);
            out[i + n * col] = mean + sd * qnorm(u01_open(bits), 0.0, 1.0, TRUE, FALSE);
        }
    }
    UNPROTECT(1);
    return ans;
}

static uint32_t bounded_u32(philox4x64_key_t k, uint64_t index, uint32_t range) {
    uint32_t threshold = (uint32_t)((UINT64_C(0x100000000) - range) % range);
    for (uint64_t attempt = 0; ; attempt++) {
        uint32_t x = (uint32_t)rngat_word(k, index, attempt, RNGAT_PURPOSE_INTEGER);
        uint64_t product = (uint64_t)x * (uint64_t)range;
        uint32_t low = (uint32_t)product;
        if (low >= threshold)
            return (uint32_t)(product >> 32);
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
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, total));
    set_sample_dim(ans, n, nkey);
    int *out = INTEGER(ans);
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_at(key, col);
        for (R_xlen_t i = 0; i < n; i++) {
            uint32_t offset = bounded_u32(k, (uint64_t)i, range);
            out[i + n * col] = min + (int)offset;
        }
    }
    UNPROTECT(1);
    return ans;
}

SEXP C_rng_bits(SEXP key, SEXP n_, SEXP bits_) {
    R_xlen_t nkey = key_count(key);
    R_xlen_t n = length_scalar(n_, "n");
    R_xlen_t total = checked_product(n, nkey, "sample");
    int bits = bits_scalar(bits_);

    if (bits == 32) {
        SEXP ans = PROTECT(Rf_allocVector(REALSXP, total));
        set_sample_dim(ans, n, nkey);
        double *out = REAL(ans);
        for (R_xlen_t col = 0; col < nkey; col++) {
            philox4x64_key_t k = key_at(key, col);
            for (R_xlen_t i = 0; i < n; i++) {
                uint64_t word = rngat_word(k, (uint64_t)i, 0, RNGAT_PURPOSE_BITS);
                out[i + n * col] = (double)(uint32_t)word;
            }
        }
        UNPROTECT(1);
        return ans;
    }

    SEXP ans = PROTECT(Rf_allocVector(STRSXP, total));
    set_sample_dim(ans, n, nkey);
    char buf[17];
    for (R_xlen_t col = 0; col < nkey; col++) {
        philox4x64_key_t k = key_at(key, col);
        for (R_xlen_t i = 0; i < n; i++) {
            uint64_t word = rngat_word(k, (uint64_t)i, 0, RNGAT_PURPOSE_BITS);
            snprintf(buf, sizeof buf, "%016" PRIx64, word);
            SET_STRING_ELT(ans, i + n * col, Rf_mkChar(buf));
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
    {NULL, NULL, 0}
};

void R_init_rngat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
