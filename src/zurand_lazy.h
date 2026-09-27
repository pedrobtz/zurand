/* Lazy random vectors: an ALTREP double vector whose elements are computed
 * when read. Included once from zurand.c, after the fills it uses.
 *
 * Every value is a pure function of (key, position), so element i of a lazy
 * vector is simply position i of rng_uniform() or rng_normal() for its key,
 * computed on demand. Reading it all gives identically the materialised
 * vector.
 *
 * State. data1 is list(key, meta, st): the rng_key (for serialisation and
 * printing), meta = c(n, dist, a, s), and st, a raw vector holding a
 * lazy_state -- the decoded key, the scaling, and a one-chunk cache.
 * data2 is the materialised vector once something has asked for the data
 * pointer, R_NilValue until then.
 *
 * Reads. The cache holds one engine chunk (512 values; 4096 or 5120 for
 * xoshiro256pp), filled by the column fills' own per-chunk code. A single
 * read that misses fills the chunk only when reads walk forward (the last
 * one plus one, or a chunk's first position); a scattered read computes the
 * value alone -- one Philox block, or for xoshiro256pp a reseed and up to
 * 512 steps -- because filling a chunk for it measured 50x slower. Regions
 * (sum(), x[a:b]) take whole chunks straight into R's buffer, partial ones
 * through the cache, and long ones through the parallel range fill; x[idx]
 * reads runs of consecutive indices as regions.
 *
 * Materialisation. Anything that asks for the data pointer gets the whole
 * vector, filled once with the ordinary parallel fill and kept in data2;
 * from then on every read is served from it.
 *
 * Scaling is applied per element with the formula affine_pass() uses, so
 * scaled values are identical however they are reached. */

#include <R_ext/Altrep.h>

/* Regions this long go to fill_range_key(), threads and all. */
#define ZURAND_LAZY_BIG_REGION ((R_xlen_t)1 << 15)
/* Runs of consecutive indices this long in x[idx] are read as a region. */
#define ZURAND_LAZY_RUN 16

typedef struct {
    philox4x64_key_t kp;
    threefry4x64_key_t kt;
    int engine;
    int dist;
    double a, s;               /* x = a + s * u, or a when s == 0 */
    R_xlen_t n;
    R_xlen_t cw;               /* the engine's chunk: 512, or 4096/5120 */
    R_xlen_t cchunk;           /* chunk in the cache, -1 for none */
    R_xlen_t clen;
    R_xlen_t last;             /* the last position read one at a time */
    int written;               /* the data pointer went out writable */
    double cache[ZURAND_MAX_CHUNK_WORDS];
} lazy_state;

static R_altrep_class_t zurand_lazy_class;
/* How a single read that misses the cache is served, for measurement:
 * 0 computes the value alone, 2 fills the chunk around it, 1 (default)
 * fills the chunk only when reads move forward one at a time. */
static int zurand_lazy_cache_mode = 1;

static lazy_state *lazy_st(SEXP x) {
    return (lazy_state *) RAW(VECTOR_ELT(R_altrep_data1(x), 2));
}

static uint64_t lazy_purpose(int dist) {
    return dist == ZURAND_DIST_NORMAL ? ZURAND_PURPOSE_NORMAL
         : dist == ZURAND_DIST_NORMAL_MCF ? ZURAND_PURPOSE_NORMAL_MCF
         : ZURAND_PURPOSE_UNIFORM;
}

static double lazy_scale1(const lazy_state *st, double u) {
    return st->s == 0.0 ? st->a
         : (st->a != 0.0 || st->s != 1.0) ? st->a + zurand_rounded(st->s * u)
         : u;
}

static R_xlen_t lazy_chunk_len(const lazy_state *st, R_xlen_t c) {
    R_xlen_t left = st->n - c * st->cw;
    return left < st->cw ? left : st->cw;
}

/* Engine chunk c, scaled, into out: the column fills' own per-chunk code,
 * vector paths included. */
static void lazy_chunk(const lazy_state *st, R_xlen_t c, double *out,
                       R_xlen_t m) {
    if (st->s == 0.0) {
        for (R_xlen_t j = 0; j < m; j++)
            out[j] = st->a;
        return;
    }
    dist_chunk_key((zurand_engine_t)st->engine, st->kp, st->kt, (uint64_t)c,
                   out, (int)m, st->dist);
    if (st->a != 0.0 || st->s != 1.0)
        affine_pass(out, m, st->a, st->s, 1);
}

static void lazy_fill_cache(lazy_state *st, R_xlen_t c) {
    st->clen = lazy_chunk_len(st, c);
    lazy_chunk(st, c, st->cache, st->clen);
    st->cchunk = c;
}

/* Position i alone: one Philox block for the counter engines; a reseed and
 * i % 512 + 1 steps for xoshiro256pp. */
static double lazy_value(const lazy_state *st, R_xlen_t i) {
    if (st->s == 0.0)
        return st->a;
    uint64_t p = (uint64_t)i, purpose = lazy_purpose(st->dist);
    double u;
    if (st->engine == ZURAND_ENG_XOSHIRO) {
        uint64_t words[ZURAND_XOSHIRO_SUB];
        int j = (int)(p % ZURAND_XOSHIRO_SUB);
        xoshiro_sub_scalar(st->kp, p / ZURAND_XOSHIRO_SUB, purpose, words, j + 1);
        u = st->dist == ZURAND_DIST_NORMAL ? zig_normal_at_xoshiro(st->kp, p, words[j])
          : st->dist == ZURAND_DIST_NORMAL_MCF ? mcf_normal_at_xoshiro(st->kp, p, words[j])
          : u01_open(words[j]);
    } else if (st->engine == ZURAND_ENG_THREEFRY) {
        uint64_t w = zurand_word_threefry(st->kt, p, 0, purpose);
        u = st->dist == ZURAND_DIST_NORMAL ? zig_normal_at_threefry(st->kt, p, w)
          : st->dist == ZURAND_DIST_NORMAL_MCF ? mcf_normal_at_threefry(st->kt, p, w)
          : u01_open(w);
    } else {
        uint64_t w = zurand_word_philox(st->kp, p, 0, purpose);
        u = st->dist == ZURAND_DIST_NORMAL ? zig_normal_at_philox(st->kp, p, w)
          : st->dist == ZURAND_DIST_NORMAL_MCF ? mcf_normal_at_philox(st->kp, p, w)
          : u01_open(w);
    }
    return lazy_scale1(st, u);
}

/* One read: from the cache, or -- on a miss -- by filling the chunk when
 * the reads are walking forward, or by computing the value alone when
 * they are scattered, where a chunk would be mostly wasted. */
static double lazy_elt(lazy_state *st, R_xlen_t i) {
    R_xlen_t c = i / st->cw;
    if (c == st->cchunk) {
        st->last = i;
        return st->cache[i - c * st->cw];
    }
    int fill = zurand_lazy_cache_mode == 2 ||
        (zurand_lazy_cache_mode == 1 && (i == st->last + 1 || i == c * st->cw));
    st->last = i;
    if (!fill)
        return lazy_value(st, i);
    lazy_fill_cache(st, c);
    return st->cache[i - c * st->cw];
}

/* Positions i .. i + n - 1, scaled, into buf: whole chunks straight in,
 * partial ones through the cache, long regions to the parallel fill. */
static void lazy_region(lazy_state *st, R_xlen_t i, R_xlen_t n, double *buf) {
    if (n >= ZURAND_LAZY_BIG_REGION && st->s != 0.0) {
        int nt = zurand_threads();
        fill_range_key((zurand_engine_t)st->engine, st->kp, buf, (uint64_t)i, n,
                       st->dist, n >= ZURAND_OMP_MIN_VALUES ? nt : 0);
        if (st->a != 0.0 || st->s != 1.0)
            affine_pass(buf, n, st->a, st->s, nt);
        return;
    }
    R_xlen_t done = 0;
    while (done < n) {
        R_xlen_t pos = i + done;
        R_xlen_t c = pos / st->cw;
        R_xlen_t off = pos - c * st->cw;
        R_xlen_t m = lazy_chunk_len(st, c);
        if (off == 0 && n - done >= m && c != st->cchunk) {
            lazy_chunk(st, c, buf + done, m);
            done += m;
            continue;
        }
        if (c != st->cchunk)
            lazy_fill_cache(st, c);
        R_xlen_t take = m - off < n - done ? m - off : n - done;
        memcpy(buf + done, st->cache + off, (size_t)take * sizeof(double));
        done += take;
    }
    st->last = i + n - 1;
}

/* ---- construction ---- */

static SEXP lazy_make(SEXP key, SEXP meta) {
    const double *m = REAL(meta);
    SEXP st_raw = PROTECT(Rf_allocVector(RAWSXP, sizeof(lazy_state)));
    lazy_state *st = (lazy_state *) RAW(st_raw);
    memset(st, 0, sizeof *st);
    st->kp = key_from_words(INTEGER(key), key_count(key), 0);
    st->kt = threefry_key_from_philox(st->kp);
    st->engine = (int)key_engine_code(key);
    st->n = (R_xlen_t)m[0];
    st->dist = (int)m[1];
    st->a = m[2];
    st->s = m[3];
    st->cw = engine_chunk_words_raw((zurand_engine_t)st->engine);
    st->cchunk = -1;
    st->last = -2;
    SEXP data1 = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(data1, 0, key);
    SET_VECTOR_ELT(data1, 1, meta);
    SET_VECTOR_ELT(data1, 2, st_raw);
    SEXP ans = R_new_altrep(zurand_lazy_class, data1, R_NilValue);
    UNPROTECT(2);
    return ans;
}

/* rng_lazy_uniform() / rng_lazy_normal(): arguments as the samplers take
 * them, dist from ZURAND_DIST_*. */
SEXP C_rng_lazy(SEXP key, SEXP n_, SEXP dist_, SEXP a_, SEXP b_) {
    if (key_count(key) != 1)
        Rf_error("`key` must be a single key");
    (void) key_engine_code(key);
    key_stream(key);
    R_xlen_t n = length_scalar(n_, "n");
    int dist = Rf_asInteger(dist_);
    double a, s;
    if (dist == ZURAND_DIST_UNIFORM) {
        double min = finite_scalar(a_, "min"), max = finite_scalar(b_, "max");
        if (min > max)
            Rf_error("`min` must be less than or equal to `max`");
        a = min;
        s = max - min;
    } else {
        a = finite_scalar(a_, "mean");
        s = finite_scalar(b_, "sd");
        if (s < 0)
            Rf_error("`sd` must be non-negative");
    }
    SEXP meta = PROTECT(Rf_allocVector(REALSXP, 4));
    REAL(meta)[0] = (double)n;
    REAL(meta)[1] = (double)dist;
    REAL(meta)[2] = a;
    REAL(meta)[3] = s;
    SEXP ans = lazy_make(key, meta);
    UNPROTECT(1);
    return ans;
}

/* ---- materialisation ---- */

static SEXP lazy_materialise(SEXP x) {
    SEXP v = R_altrep_data2(x);
    if (v != R_NilValue)
        return v;
    lazy_state *st = lazy_st(x);
    R_xlen_t n = st->n;
    v = PROTECT(Rf_allocVector(REALSXP, n));
    double *out = REAL(v);
    zurand_advise_output(out, (size_t)n * sizeof(double));
    int nt = zurand_threads();
    if (st->s == 0.0) {
        for (R_xlen_t i = 0; i < n; i++)
            out[i] = st->a;
    } else {
        int par_rows = n >= ZURAND_OMP_MIN_VALUES ? nt : 0;
        if (st->dist == ZURAND_DIST_UNIFORM)
            fill_uniform_key((zurand_engine_t)st->engine, st->kp, out, n, par_rows);
        else
            fill_normal_key((zurand_engine_t)st->engine, st->kp, out, n, par_rows,
                            st->dist == ZURAND_DIST_NORMAL_MCF
                                ? ZURAND_NORMAL_MCFARLAND : ZURAND_NORMAL_ZIGGURAT);
        if (st->a != 0.0 || st->s != 1.0)
            affine_pass(out, n, st->a, st->s, nt);
    }
    R_set_altrep_data2(x, v);
    UNPROTECT(1);
    return v;
}

SEXP C_rng_lazy_materialised(SEXP x) {
    if (!ALTREP(x) || !R_altrep_inherits(x, zurand_lazy_class))
        return Rf_ScalarLogical(NA_LOGICAL);
    return Rf_ScalarLogical(R_altrep_data2(x) != R_NilValue);
}

SEXP C_rng_lazy_cache(SEXP mode) {
    int prev = zurand_lazy_cache_mode;
    if (mode != R_NilValue)
        zurand_lazy_cache_mode = Rf_asInteger(mode);
    return Rf_ScalarInteger(prev);
}

/* ---- ALTREP methods ---- */

static R_xlen_t lazy_Length(SEXP x) {
    return lazy_st(x)->n;
}

static double lazy_Elt(SEXP x, R_xlen_t i) {
    SEXP v = R_altrep_data2(x);
    if (v != R_NilValue)
        return REAL(v)[i];
    return lazy_elt(lazy_st(x), i);
}

static R_xlen_t lazy_Get_region(SEXP x, R_xlen_t i, R_xlen_t n, double *buf) {
    SEXP v = R_altrep_data2(x);
    lazy_state *st = lazy_st(x);
    if (i >= st->n)
        return 0;
    if (n > st->n - i)
        n = st->n - i;
    if (v != R_NilValue)
        memcpy(buf, REAL(v) + i, (size_t)n * sizeof(double));
    else
        lazy_region(st, i, n, buf);
    return n;
}

/* x[idx]: runs of consecutive indices are read as regions, the rest one by
 * one. Indices arrive 1-based, as R's ExtractSubset() takes them; NA and
 * out-of-range ones give NA. */
static R_xlen_t lazy_index(SEXP indx, R_xlen_t k, R_xlen_t n) {
    if (TYPEOF(indx) == INTSXP) {
        int v = INTEGER_ELT(indx, k);
        return v == NA_INTEGER || v < 1 || v > n ? -1 : (R_xlen_t)v - 1;
    }
    double d = REAL_ELT(indx, k);
    if (!R_FINITE(d))
        return -1;
    R_xlen_t v = (R_xlen_t)d - 1;
    return v < 0 || v >= n ? -1 : v;
}

static SEXP lazy_Extract_subset(SEXP x, SEXP indx, SEXP call) {
    (void)call;
    if (R_altrep_data2(x) != R_NilValue ||
        (TYPEOF(indx) != INTSXP && TYPEOF(indx) != REALSXP))
        return NULL;
    lazy_state *st = lazy_st(x);
    R_xlen_t m = XLENGTH(indx);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, m));
    double *out = REAL(ans);
    R_xlen_t k = 0;
    while (k < m) {
        R_xlen_t i = lazy_index(indx, k, st->n);
        if (i < 0) {
            out[k++] = NA_REAL;
            continue;
        }
        R_xlen_t r = 1;
        while (k + r < m && r < st->n - i && lazy_index(indx, k + r, st->n) == i + r)
            r++;
        if (r >= ZURAND_LAZY_RUN)
            lazy_region(st, i, r, out + k);
        else
            for (R_xlen_t j = 0; j < r; j++)
                out[k + j] = lazy_elt(st, i + j);
        k += r;
    }
    UNPROTECT(1);
    return ans;
}

static void *lazy_Dataptr(SEXP x, Rboolean writable) {
    SEXP v = lazy_materialise(x);
    if (writable)
        lazy_st(x)->written = 1;   /* the values may now differ from the key's */
    return REAL(v);
}

static const void *lazy_Dataptr_or_null(SEXP x) {
    SEXP v = R_altrep_data2(x);
    return v == R_NilValue ? NULL : (const void *) REAL(v);
}

static int lazy_No_NA(SEXP x) {
    (void)x;
    return 1;
}

/* The key and arguments, so a lazy vector saves as a few hundred bytes --
 * unless its data pointer has gone out writable: then it may have been
 * modified, and NULL makes R save the values as a plain vector. */
static SEXP lazy_Serialized_state(SEXP x) {
    if (lazy_st(x)->written)
        return NULL;
    SEXP d1 = R_altrep_data1(x);
    SEXP s = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(s, 0, VECTOR_ELT(d1, 0));
    SET_VECTOR_ELT(s, 1, VECTOR_ELT(d1, 1));
    UNPROTECT(1);
    return s;
}

static SEXP lazy_Unserialize(SEXP cls, SEXP state) {
    (void)cls;
    return lazy_make(VECTOR_ELT(state, 0), VECTOR_ELT(state, 1));
}

static SEXP lazy_Duplicate(SEXP x, Rboolean deep) {
    (void)deep;
    SEXP v = R_altrep_data2(x);
    if (v != R_NilValue)
        return Rf_duplicate(v);           /* a plain vector from here on */
    SEXP d1 = R_altrep_data1(x);
    return lazy_make(VECTOR_ELT(d1, 0), VECTOR_ELT(d1, 1));
}

static Rboolean lazy_Inspect(SEXP x, int pre, int deep, int pvec,
                             void (*inspect_subtree)(SEXP, int, int, int)) {
    (void)pre; (void)deep; (void)pvec; (void)inspect_subtree;
    const lazy_state *st = lazy_st(x);
    static const char *dists[] = {"uniform", "normal", "normal (mcfarland)"};
    static const char *engines[] = {"philox4x64", "threefry4x64", "xoshiro256pp"};
    Rprintf(" zurand lazy %s, n = %.0f, engine %s, %s\n", dists[st->dist],
            (double)st->n, engines[st->engine],
            R_altrep_data2(x) == R_NilValue ? "not materialised" : "materialised");
    return TRUE;
}

static void zurand_lazy_init(DllInfo *dll) {
    R_altrep_class_t c = R_make_altreal_class("zurand_lazy", "zurand", dll);
    R_set_altrep_Length_method(c, lazy_Length);
    R_set_altrep_Inspect_method(c, lazy_Inspect);
    R_set_altrep_Serialized_state_method(c, lazy_Serialized_state);
    R_set_altrep_Unserialize_method(c, lazy_Unserialize);
    R_set_altrep_Duplicate_method(c, lazy_Duplicate);
    R_set_altvec_Dataptr_method(c, lazy_Dataptr);
    R_set_altvec_Dataptr_or_null_method(c, lazy_Dataptr_or_null);
    R_set_altreal_Elt_method(c, lazy_Elt);
    R_set_altreal_Get_region_method(c, lazy_Get_region);
    R_set_altvec_Extract_subset_method(c, lazy_Extract_subset);
    R_set_altreal_No_NA_method(c, lazy_No_NA);
    zurand_lazy_class = c;
}
