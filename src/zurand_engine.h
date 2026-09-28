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

/* Counter word 3 for this engine's Philox/Threefry calls (see
 * ZURAND_TAG_* in zurand.c). */
#ifndef ZE_TAG
#  define ZE_TAG ZURAND_TAG_PHILOX
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
    zurand_ctr_t ctr = {{index, domain, purpose, ZE_TAG}};
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
 * loops. `noinline` only, not `cold`: GCC compiles a cold function for
 * size, and that cost 2-5% of the whole normal fill on Zen 3 (three
 * runners, 10 of 10 rounds each). */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
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
            double yy = -zurand_log1p(-ZURAND_U64_TO_DOUBLE(Y));
            for (;;) {
                double xx = zurand_rounded(-ziggurat_nor_inv_r *
                    zurand_log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(zig_next)(&s, key, index))));
                if (yy + yy > xx * xx)
                    return sign ? -(ziggurat_nor_r + xx)
                                : ziggurat_nor_r + xx;
                yy = -zurand_log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(zig_next)(&s, key, index)));
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
        int shortcut = zig_wedge_shortcut(idx, rabs, Y);
        if (shortcut > 0)
            return x;
        if (shortcut == 0) {
            double u = ZURAND_U64_TO_DOUBLE(Y);
            double rise = zurand_rounded((fi_double[idx - 1] - fi_double[idx]) * u);
            if (rise + fi_double[idx] < zurand_exp(-0.5 * x * x))
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
/* Entry to the slow path from the signed-table fast path below: settles
 * the one case the new accept test sends here that the old one did not,
 * rabs == 0 in a layer with ki > 0, which the old code accepted as +0.0.
 * `noinline`, not `cold`: see zig_normal_slow. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
static double ZE_N(zig_normal_reject)(ZE_KEY_T key, uint64_t index, uint64_t r) {
    if (((r >> 9) & UINT64_C(0x000fffffffffffff)) == 0 && ki_double[r & 0xff] != 0)
        return 0.0;
    return ZE_N(zig_normal_slow)(key, index, r);
}

/* The fast path: one table multiply and one compare per draw. Bit-identical
 * to the original (double)(sign ? -rabs : rabs) * wi[idx] accepted when
 * rabs < ki[idx]; see zurand_wis / zurand_kim1 in zurand.c. */
R123_STATIC_INLINE double ZE_N(zig_normal_at)(ZE_KEY_T key, uint64_t index,
                                        uint64_t r) {
    uint64_t rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
    double x = (double)rabs * zurand_wis[r & 0x1ff];
    if (R123_BUILTIN_EXPECT(rabs - 1 < zurand_kim1[r & 0xff], 1))
        return x;
    return ZE_N(zig_normal_reject)(key, index, r);
}

/* ---- McFarland's modified ziggurat: rng_normal(method = "mcfarland") ----
 *
 * C. D. McFarland (2016), "A modified ziggurat algorithm for generating
 * exponentially and normally distributed pseudorandom numbers", J. Stat.
 * Comput. Simul. 86(7), 1281-1294; https://github.com/cd-mcfarland/fast_prng
 * (MIT). Tables and bounds: mcfarland_normal.h.
 *
 * Fast path, 253/256 of draws: the word's low byte picks layer i; below
 * ZURAND_MCF_LAYERS the draw is X_i * (the word as a signed integer). The
 * accept test needs no table load and the sign comes with the word; the
 * low byte is squashed by the conversion to double. Measured 21-39% less
 * time per draw than the NumPy-table fast path on eight CPUs
 * (dev/simd/mcfarland.c).
 *
 * Edge, 3/256: alias sampling picks the tail (j = 0) or overhang j, then
 * rejection inside it. Every word after the first comes from retry blocks
 * at (index, attempt >= 1) under this method's own purpose value, so draw
 * i is a pure function of (key, i), as for the NumPy-table sampler. The
 * sign is the first word's top bit, and its other 63 bits are the first
 * overhang position, as in McFarland's C and Boyce's Fortran (BiF-lib):
 * the low byte, known to be >= 253 here, is squashed by the conversion to
 * double. That fits 88% of edge draws in one retry block instead of 71%
 * (2,477 fewer Philox calls per 1e6 draws). Further uniforms u in
 * [0, 2^63) are a word's top 63 bits. */
R123_STATIC_INLINE uint64_t ZE_N(mcf_next)(ZE_N(zig_stream) *s, ZE_KEY_T key,
                                     uint64_t index) {
    if (s->w == 4) {
        s->blk = ZE_N(zurand_block)(key, index, ++s->group, ZURAND_PURPOSE_NORMAL_MCF);
        s->w = 0;
    }
    return s->blk.v[s->w++];
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))   /* not cold: see zig_normal_slow */
#endif
static double ZE_N(mcf_normal_edge)(ZE_KEY_T key, uint64_t index, uint64_t w) {
    ZE_N(zig_stream) s = {.group = 0, .w = 4};
    const double *X = zurand_mcf_x, *Y = zurand_mcf_y;
    double sign = (w >> 63) ? -1.0 : 1.0;

    uint64_t a = ZE_N(mcf_next)(&s, key, index);
    unsigned j = (unsigned)(a & 0xff);
    if ((a >> 8) >= zurand_mcf_alias_t[j])
        j = zurand_mcf_alias_j[j];

    double x;
    if (j > ZURAND_MCF_J_INFLECTION) {
        /* x < 1, the curve above the chord: below the chord (u2 >= u1) is
         * inside; far above it is outside; the band between asks exp(). */
        uint64_t u1 = w & UINT64_C(0x7fffffffffffffff);
        for (;;) {
            uint64_t u2 = ZE_N(mcf_next)(&s, key, index) >> 1;
            x = zurand_mcf_interp(X, j, u1);
            int shortcut = mcf_edge_shortcut(j, u1, u2);
            if (shortcut > 0 || (shortcut == 0 &&
                zurand_mcf_interp(Y, j, u2) < zurand_exp(-0.5 * x * x)))
                break;
            u1 = ZE_N(mcf_next)(&s, key, index) >> 1;
        }
    } else if (j == 0) {
        /* tail beyond X_0: Marsaglia's method, as the NumPy-table sampler */
        for (;;) {
            double xx = zurand_rounded(-zurand_mcf_inv_x0 *
                zurand_log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(mcf_next)(&s, key, index))));
            double yy = -zurand_log1p(-ZURAND_U64_TO_DOUBLE(ZE_N(mcf_next)(&s, key, index)));
            if (yy + yy > xx * xx) {
                x = zurand_mcf_x0 + xx;
                break;
            }
        }
    } else if (j < ZURAND_MCF_J_INFLECTION) {
        /* x > 1, the curve below the chord: reflect into the lower
         * triangle; well below the chord is inside. */
        uint64_t u1 = w & UINT64_C(0x7fffffffffffffff);
        for (;;) {
            uint64_t u2 = ZE_N(mcf_next)(&s, key, index) >> 1;
            if (u2 < u1) {
                uint64_t t = u1; u1 = u2; u2 = t;
            }
            x = zurand_mcf_interp(X, j, u1);
            if (mcf_edge_shortcut(j, u1, u2) > 0 ||
                zurand_mcf_interp(Y, j, u2) < zurand_exp(-0.5 * x * x))
                break;
            u1 = ZE_N(mcf_next)(&s, key, index) >> 1;
        }
    } else {
        /* the overhang straddling x = 1: the whole rectangle */
        uint64_t u1 = w & UINT64_C(0x7fffffffffffffff);
        for (;;) {
            uint64_t u2 = ZE_N(mcf_next)(&s, key, index) >> 1;
            x = zurand_mcf_interp(X, j, u1);
            if (zurand_mcf_interp(Y, j, u2) < zurand_exp(-0.5 * x * x))
                break;
            u1 = ZE_N(mcf_next)(&s, key, index) >> 1;
        }
    }
    return sign * x;
}

R123_STATIC_INLINE double ZE_N(mcf_normal_at)(ZE_KEY_T key, uint64_t index,
                                        uint64_t w) {
    unsigned i = (unsigned)(w & 0xff);
    if (R123_BUILTIN_EXPECT(i < ZURAND_MCF_LAYERS, 1))
        return zurand_mcf_x[i] * (double)(int64_t)w;
    return ZE_N(mcf_normal_edge)(key, index, w);
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

/* One chunk of a sampler: positions c * ZE_CHUNK_WORDS onward, m of them
 * (m < ZE_CHUNK_WORDS only for the last chunk of a fill), into o[0 .. m-1].
 * Every fill -- a whole column, a range from an offset, a stream -- is a
 * loop over these, so all of them give the same value at each position.
 * Forced inline so the column loops compile exactly as they did when the
 * body was written in them. */
R123_STATIC_INLINE R123_FORCE_INLINE(void ZE_N(normal_chunk)(
    ZE_KEY_T key, uint64_t c, double *o, int m));
R123_STATIC_INLINE void ZE_N(normal_chunk)(ZE_KEY_T key, uint64_t c,
                                           double *o, int m) {
    uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
    uint64_t w0 = c * (uint64_t)ZE_CHUNK_WORDS;
    ZE_N(chunk_words)(key, c, ZURAND_PURPOSE_NORMAL, buf, m);

    /* Unrolled by four: measured ~20% faster on Apple M1, level on x86. */
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j]     = ZE_N(zig_normal_at)(key, w0 + (uint64_t)j,       buf[j]);
        o[j + 1] = ZE_N(zig_normal_at)(key, w0 + (uint64_t)j + 1,   buf[j + 1]);
        o[j + 2] = ZE_N(zig_normal_at)(key, w0 + (uint64_t)j + 2,   buf[j + 2]);
        o[j + 3] = ZE_N(zig_normal_at)(key, w0 + (uint64_t)j + 3,   buf[j + 3]);
    }
    for (; j < m; j++)
        o[j] = ZE_N(zig_normal_at)(key, w0 + (uint64_t)j, buf[j]);
}

static void ZE_N(fill_normal_column)(double *out, R_xlen_t n,
                               ZE_KEY_T key, int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(normal_chunk)(key, (uint64_t)c, out + w0, m);
    }
}

/* The same chunk and fill for McFarland's method, under its own purpose
 * value. */
R123_STATIC_INLINE R123_FORCE_INLINE(void ZE_N(normal_mcf_chunk)(
    ZE_KEY_T key, uint64_t c, double *o, int m));
R123_STATIC_INLINE void ZE_N(normal_mcf_chunk)(ZE_KEY_T key, uint64_t c,
                                               double *o, int m) {
    uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
    uint64_t w0 = c * (uint64_t)ZE_CHUNK_WORDS;
    ZE_N(chunk_words)(key, c, ZURAND_PURPOSE_NORMAL_MCF, buf, m);

    /* Unrolled by four: 12-20% faster on arm64, level on x86. */
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j]     = ZE_N(mcf_normal_at)(key, w0 + (uint64_t)j,     buf[j]);
        o[j + 1] = ZE_N(mcf_normal_at)(key, w0 + (uint64_t)j + 1, buf[j + 1]);
        o[j + 2] = ZE_N(mcf_normal_at)(key, w0 + (uint64_t)j + 2, buf[j + 2]);
        o[j + 3] = ZE_N(mcf_normal_at)(key, w0 + (uint64_t)j + 3, buf[j + 3]);
    }
    for (; j < m; j++)
        o[j] = ZE_N(mcf_normal_at)(key, w0 + (uint64_t)j, buf[j]);
}

static void ZE_N(fill_normal_mcf_column)(double *out, R_xlen_t n,
                                   ZE_KEY_T key, int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(normal_mcf_chunk)(key, (uint64_t)c, out + w0, m);
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
R123_STATIC_INLINE R123_FORCE_INLINE(void ZE_N(uniform_chunk)(
    ZE_KEY_T key, uint64_t c, double *o, int m, int large));
R123_STATIC_INLINE void ZE_N(uniform_chunk)(ZE_KEY_T key, uint64_t c,
                                            double *o, int m, int large) {
#ifdef ZE_UNIFORM_FAST
    /* An engine may write a whole chunk of uniforms itself. `large` is
     * decided on the whole fill, not on one thread's share, and only
     * chooses how the stores are made. */
    if (ZE_UNIFORM_FAST(key, c, o, m, large))
        return;
#else
    (void)large;
#endif
    uint64_t buf[ZURAND_MAX_CHUNK_WORDS + ZURAND_CHUNK_SLACK];
    ZE_N(chunk_words)(key, c, ZURAND_PURPOSE_UNIFORM, buf, m);

    for (int j = 0; j < m; j++)
        o[j] = u01_open(buf[j]);
}

static void ZE_N(fill_uniform_column)(double *out, R_xlen_t n,
                                ZE_KEY_T key, int threads) {
    R_xlen_t nchunk = (n + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZE_UNIFORM_FAST
    int large = n >= ZURAND_ZVA_MIN_VALUES;
#else
    int large = 0;
#endif
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(out, n, nchunk, key, large) schedule(static)
#endif
    for (R_xlen_t c = 0; c < nchunk; c++) {
        R_xlen_t w0 = c * ZE_CHUNK_WORDS;
        int m = (int)(n - w0 < ZE_CHUNK_WORDS ? n - w0 : ZE_CHUNK_WORDS);
        ZE_N(uniform_chunk)(key, (uint64_t)c, out + w0, m, large);
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
            /* offset can exceed INT_MAX when max - min >= 2^31; the sum
             * itself lies in [min, max], so widen, add, then narrow. */
            o[j] = (int)((int64_t)min + (int64_t)offset);
        }
    }
}

/* One chunk of any double sampler (ZURAND_DIST_*). */
static void ZE_N(dist_chunk)(ZE_KEY_T key, uint64_t c, double *o, int m,
                             int dist, int large) {
    if (dist == ZURAND_DIST_NORMAL)
        ZE_N(normal_chunk)(key, c, o, m);
    else if (dist == ZURAND_DIST_NORMAL_MCF)
        ZE_N(normal_mcf_chunk)(key, c, o, m);
    else
        ZE_N(uniform_chunk)(key, c, o, m, large);
}

/* Positions start .. start + n - 1 of a double sampler into out[0 .. n-1]:
 * the column fill from an offset. Chunks are generated from their first
 * position, so a range that starts inside a chunk generates that chunk's
 * head into a stack buffer and copies the part it needs; every later chunk
 * is written in place, and the last may be short, as in a column fill.
 * `stage` holds ZE_CHUNK_WORDS doubles. */
static void ZE_N(fill_range)(double *out, uint64_t start, R_xlen_t n,
                             ZE_KEY_T key, int dist, int threads,
                             double *stage) {
    if (n <= 0)
        return;
    uint64_t c0 = start / (uint64_t)ZE_CHUNK_WORDS;
    int head = (int)(start % (uint64_t)ZE_CHUNK_WORDS);
    R_xlen_t done = 0;
#ifdef ZE_UNIFORM_FAST
    int large = n >= ZURAND_ZVA_MIN_VALUES;
#else
    int large = 0;
#endif
    if (head) {
        R_xlen_t room = ZE_CHUNK_WORDS - head;
        done = n < room ? n : room;
        ZE_N(dist_chunk)(key, c0, stage, head + (int)done, dist, 0);
        memcpy(out, stage + head, (size_t)done * sizeof(double));
        c0++;
    }
    double *rest = out + done;
    R_xlen_t nrest = n - done;
    R_xlen_t nchunk = (nrest + ZE_CHUNK_WORDS - 1) / ZE_CHUNK_WORDS;
#ifdef ZURAND_OPENMP
#pragma omp parallel for if(threads > 1) \
    num_threads(threads > 0 ? threads : 1) default(none) \
    shared(rest, nrest, nchunk, key, dist, large, c0) schedule(static)
#endif
    for (R_xlen_t k = 0; k < nchunk; k++) {
        R_xlen_t w0 = k * ZE_CHUNK_WORDS;
        int m = (int)(nrest - w0 < ZE_CHUNK_WORDS ? nrest - w0 : ZE_CHUNK_WORDS);
        ZE_N(dist_chunk)(key, c0 + (uint64_t)k, rest + w0, m, dist, large);
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
#undef ZE_TAG
#undef ZE_UNIFORM_FAST
