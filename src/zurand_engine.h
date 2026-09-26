/* Engine-parameterised sampler core -- included once per engine.
 *
 * Every Random123 generator zurand offers has a 256-bit counter and a
 * 256-bit output, so between engines only the key type and the generation
 * call differ. This file is included once per engine from zurand.c with
 * ZE_SUFFIX, ZE_KEY_T and ZE_GEN defined, yielding one monomorphic copy of
 * the sampler core each.
 *
 * A runtime branch inside zurand_block() would have been far simpler, and
 * would also have destroyed the forced inlining the hot loops depend on --
 * see the inlining note in CLAUDE.md. Generating a copy per engine keeps
 * every inner loop exactly as tight as the single-engine version was.
 *
 * Engine-independent helpers (u01_open, ZURAND_U64_TO_DOUBLE,
 * lemire_accept, ZURAND_CHUNK_BLOCKS) deliberately live in zurand.c so
 * they are defined once; this file uses them and must be included after.
 */
#ifndef ZE_SUFFIX
#  error "include zurand_engine.h from zurand.c only, with ZE_* defined"
#endif

/* Words per chunk for this engine. The default is ZURAND_CHUNK_WORDS; an
 * engine defines ZE_CHUNK_WORDS before including to take a wider one. The
 * value only changes how much the fill asks for at a time -- for the
 * counter-based engines word w is block w/4 regardless, so their output
 * does not depend on it. */
#ifndef ZE_CHUNK_WORDS
#  define ZE_CHUNK_WORDS ZURAND_CHUNK_WORDS
#endif

#define ZE_CAT2(a, b) a##_##b
#define ZE_CAT(a, b) ZE_CAT2(a, b)
#define ZE_N(name) ZE_CAT(name, ZE_SUFFIX)

R123_STATIC_INLINE R123_FORCE_INLINE(zurand_ctr_t ZE_N(zurand_block)(
    ZE_KEY_T key, uint64_t index, uint64_t domain, uint64_t purpose));
R123_STATIC_INLINE zurand_ctr_t ZE_N(zurand_block)(ZE_KEY_T key,
                                                uint64_t index,
                                                uint64_t domain,
                                                uint64_t purpose) {
    zurand_ctr_t ctr = {{index, domain, purpose, 0}};
    return ZE_GEN(ctr, key);
}

static uint64_t ZE_N(zurand_word)(ZE_KEY_T key, uint64_t index,
                           uint64_t domain, uint64_t purpose) {
    zurand_ctr_t block = ZE_N(zurand_block)(key, index >> 2, domain, purpose);
    return block.v[index & 3u];
}



/* Retry-word stream for the slow path: group g >= 1 is one Philox block
 * at counter {index, g, purpose}, consumed word by word, so all four
 * words of each retry block are used and draw `index` stays a pure
 * function of (key, index). Distinct from the fast-path counters, whose
 * second word is always 0. (The former one-word-per-attempt scheme paid
 * a full Philox call per retry word and discarded the other three;
 * batching measured ~3% faster on bulk rnorm overall.) */
typedef struct {
    zurand_ctr_t blk;
    uint64_t group;
    unsigned w;
} ZE_N(zig_stream);

R123_STATIC_INLINE uint64_t ZE_N(zig_next)(ZE_N(zig_stream) *s, ZE_KEY_T key,
                                     uint64_t index) {
    if (s->w == 4) {
        s->blk = ZE_N(zurand_block)(key, index, ++s->group, ZURAND_PURPOSE_NORMAL);
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
static double ZE_N(zig_normal_slow)(ZE_KEY_T key, uint64_t index, uint64_t r) {
    ZE_N(zig_stream) s = {.group = 0, .w = 4};
    int idx = (int)(r & 0xff);
    int sign = (int)((r >> 8) & 0x1);
    uint64_t rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
    double x = (double)rabs * wi_double[idx];
    if (sign)
        x = -x;
    /* the caller established rabs >= ki_double[idx] for the entry word */

    for (;;) {
        uint64_t Y = ZE_N(zig_next)(&s, key, index);

        if (idx == 0) {
            /* layer-0 tail; the first ordinate reuses Y */
            double yy = -log1p(-ZURAND_U64_TO_DOUBLE(Y));
            for (;;) {
                double xx = zurand_rounded(-ziggurat_nor_inv_r *
                    log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(zig_next)(&s, key, index))));
                if (yy + yy > xx * xx)
                    return sign ? -(ziggurat_nor_r + xx)
                                : ziggurat_nor_r + xx;
                yy = -log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(zig_next)(&s, key, index)));
            }
        }

        /* wedge: compare Y * (width of layer idx) against the position in
         * the layer, with zurand_zig_gap bracketing the exp curve around its
         * chord; f is concave below the inflection layer (x < 1, curve
         * above the chord) and convex above it. ZURAND_ZIG_GUARD widens the
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
        if (idx > ZURAND_ZIG_INFLECTION) {
            reject = YL > R + ZURAND_ZIG_GUARD;
            accept = !reject && YL + zurand_zig_gap[idx] < R;
        } else if (idx < ZURAND_ZIG_INFLECTION) {
            accept = YL + ZURAND_ZIG_GUARD < R;
            reject = !accept && YL > R + zurand_zig_gap[idx];
        } else {
            reject = YL > R + zurand_zig_gap_hi52;
            accept = !reject && YL + zurand_zig_gap[idx] < R;
        }
        if (accept)
            return x;
        if (!reject) {
            double u = ZURAND_U64_TO_DOUBLE(Y);
            double rise = zurand_rounded((fi_double[idx - 1] - fi_double[idx]) * u);
            if (rise + fi_double[idx] < exp(-0.5 * x * x))
                return x;
        }

        r = ZE_N(zig_next)(&s, key, index);
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
R123_STATIC_INLINE double ZE_N(zig_normal_at)(ZE_KEY_T key, uint64_t index,
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
    return ZE_N(zig_normal_slow)(key, index, r);
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

/* Fill `nwords` words for chunk `c`. This is the one place an engine
 * decides how bits are produced in bulk, and it is a *chunk* rather than a
 * block because that is the granularity an engine can amortise over: the
 * xoshiro engine pays one Philox call here to seed a chunk and then runs a
 * cheap recurrence, which a per-block hook could not express.
 *
 * `nwords` rather than always the full chunk, so that rng_uniform(key, 1)
 * still costs one block instead of a whole chunk.
 *
 * Word w of a column is chunk w / ZURAND_CHUNK_WORDS at offset
 * w % ZURAND_CHUNK_WORDS, which for the counter-based engines is exactly
 * block w / 4 word w % 4 -- the mapping they had before this indirection
 * existed, so their output is unchanged. */
#ifndef ZE_CUSTOM_CHUNK
static void ZE_N(chunk_words)(ZE_KEY_T key, uint64_t c, uint64_t purpose,
                              uint64_t *buf, int nwords) {
    /* Always copy the whole 32-byte block. A variable-length memcpy here --
     * trimming the last block to nwords -- cost philox 18% (306 -> 249 M/s):
     * the size is not a compile-time constant, so every block in the hot
     * loop became a library call instead of two register moves. The
     * caller's buffer carries ZURAND_CHUNK_SLACK words for the over-copy. */
    int nb = (nwords + 3) >> 2;
    for (int j = 0; j < nb; j++) {
        zurand_ctr_t block = ZE_N(zurand_block)(
            key, c * (ZE_CHUNK_WORDS / 4) + (uint64_t)j, 0, purpose);
        memcpy(buf + 4 * j, block.v, sizeof block.v);
    }
}
#endif

static void ZE_N(fill_normal_column)(double *out, R_xlen_t n,
                               ZE_KEY_T key, int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(chunk_words)(key, (uint64_t)c, ZURAND_PURPOSE_NORMAL, buf, m);

        double *o = out + w0;
        for (int j = 0; j < m; j++)
            o[j] = ZE_N(zig_normal_at)(key, (uint64_t)(w0 + j), buf[j]);
    }
}

/* Fills uniforms on (0, 1); C_rng_uniform applies min/span in a separate
 * vectorizable pass so the hot loop stays Philox plus the bit trick.
 *
 * Two passes over a stack chunk, for the same reason fill_normal_column
 * uses them: the fused form makes each iteration's stores wait on that
 * iteration's 10-round Philox chain, and splitting generation from
 * transformation lets the generator loop pipeline across independent
 * blocks. Same counters, same words, same order, same output.
 * `threads` is the thread count for the inner parallel region; callers
 * pass 0 to keep the column serial (e.g. when parallelizing over key
 * columns). */
static void ZE_N(fill_uniform_column)(double *out, R_xlen_t n,
                                ZE_KEY_T key, int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(chunk_words)(key, (uint64_t)c, ZURAND_PURPOSE_UNIFORM, buf, m);

        double *o = out + w0;
        for (int j = 0; j < m; j++)
            o[j] = u01_open(buf[j]);
    }
}

/* Rejection continuation for one index; attempt 0 was consumed from the
 * shared block by fill_integer_column. */
static uint32_t ZE_N(bounded_u32_retry)(ZE_KEY_T k, uint64_t index,
                                  uint32_t range, uint32_t threshold) {
    for (uint64_t attempt = 1; ; attempt++) {
        uint32_t x = (uint32_t)ZE_N(zurand_word)(k, index, attempt, ZURAND_PURPOSE_INTEGER);
        uint32_t offset;
        if (lemire_accept(x, range, threshold, &offset))
            return offset;
    }
}

static void ZE_N(fill_integer_column)(int *out, R_xlen_t n, ZE_KEY_T key,
                                int min, uint32_t range, uint32_t threshold,
                                int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key, min, range, threshold) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(chunk_words)(key, (uint64_t)c, ZURAND_PURPOSE_INTEGER, buf, m);

        int *o = out + w0;
        for (int j = 0; j < m; j++) {
            uint32_t offset;
            if (!lemire_accept((uint32_t)buf[j], range, threshold, &offset))
                offset = ZE_N(bounded_u32_retry)(key, (uint64_t)(w0 + j),
                                                 range, threshold);
            o[j] = min + (int)offset;
        }
    }
}

#undef ZE_N
#undef ZE_CAT
#undef ZE_CAT2
#undef ZE_SUFFIX
#undef ZE_KEY_T
#undef ZE_GEN
#undef ZE_CUSTOM_CHUNK
#undef ZE_CHUNK_WORDS
