/* Vectorised Box-Muller in the style of VectorizedRNG.jl, made
 * deterministic: every lane runs the same sequence of IEEE operations
 * (+ - * / sqrt, all correctly rounded) as the scalar reference, with no
 * fused multiply-add anywhere, so SIMD and scalar agree bit for bit and so
 * do all platforms. Words in cache, normals out; compared with McFarland's
 * fast path and the shipped NumPy-table ziggurat.
 *
 * Pair (w1, w2) -> (z1, z2) = R * (sin t, cos t), words 2k and 2k+1 giving
 * values 2k and 2k+1:
 *   angle   a = top 52 bits of w1 as [0, 1); sin(pi a / 2) and
 *           cos(pi a / 2) = sin(pi (1 - a) / 2) by one odd polynomial
 *           (VectorizedRNG.jl's approx_sin8 coefficients, MIT); the two
 *           low bits of w1 pick the quadrant signs. No range reduction.
 *   radius  u = ((double)w2 + 0.5) * 2^-64 in (0, 1], conversion
 *           correctly rounded via exact 32-bit halves; R = sqrt(-2 ln u)
 *           with fdlibm's log, both of its final forms computed and one
 *           selected, so it is branch-free.
 *
 *   BS  scalar reference
 *   BV  vector: AVX2 4 lanes / AVX-512 8 lanes (x86, run-time checked),
 *       NEON 2 lanes (arm64); no FMA, identical to BS
 *   BF  the same vector code with fused multiply-add in the polynomials
 *       and the log, as VectorizedRNG.jl does: what giving up
 *       cross-platform identity would buy
 *   M1  McFarland fast path, unrolled x4 (for comparison)
 *   Z6  shipped NumPy-table ziggurat fast path (for comparison)
 *
 *   cc -O2 -ffp-contract=off -I../src -o /tmp/bm ../dev/simd/boxmuller.c -lm && /tmp/bm
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "numpyzig/ziggurat_constants.h"
#if defined(__aarch64__)
#  include <arm_neon.h>
#elif defined(__x86_64__)
#  include <immintrin.h>
#endif

#define NW 4096
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }

/* no-FMA barrier, as zurand_rounded() */
static inline double RB(double x) {
#if defined(__x86_64__)
    __asm__("" : "+x"(x));
#elif defined(__aarch64__)
    __asm__("" : "+w"(x));
#endif
    return x;
}
static inline double d_of(uint64_t b) { double d; memcpy(&d, &b, 8); return d; }
static inline uint64_t b_of(double d) { uint64_t b; memcpy(&b, &d, 8); return b; }

/* sin(pi x / 2) / x as a degree-8 polynomial in x^2 on [0, 1], fitted in
 * mpmath for relative error (max 2.6e-19; see the roadmap entry) */
static const double C[9] = {
    0x1.921fb54442d18p+0, -0x1.4abbce625be52p-1, 0x1.466bc6775aa6ep-4,
    -0x1.32d2cce627536p-8, 0x1.5078348518a8ep-13, -0x1.e3074ddce95eep-19,
    0x1.e8f35e8e4077cp-25, -0x1.6f7974c32ffe1p-31, 0x1.9d1b3e747fae1p-38 };

/* fdlibm e_log.c constants */
static const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10,
    Lg1 = 6.666666666666735130e-01, Lg2 = 3.999999999940941908e-01, Lg3 = 2.857142874366239149e-01,
    Lg4 = 2.222219843214978396e-01, Lg5 = 1.818357216161805012e-01, Lg6 = 1.531383769920937332e-01,
    Lg7 = 1.479819860511658591e-01;

static inline double sinq(double x) {   /* sin(pi x / 2), x in [0, 1] */
    double t = x * x, p = C[8];
    for (int k = 7; k >= 0; k--) p = RB(p * t) + C[k];
    return p * x;
}
/* ln(u) for u in [2^-65, 1], fdlibm's algorithm without the small-f shortcut */
static inline double logu(double u) {
    uint64_t b = b_of(u);
    int32_t hx = (int32_t)(b >> 32);
    int32_t k = (hx >> 20) - 1023;
    hx &= 0x000fffff;
    int32_t i = (hx + 0x95f64) & 0x100000;
    double x = d_of(((uint64_t)(uint32_t)(hx | (i ^ 0x3ff00000)) << 32) | (b & 0xffffffffu));
    k += i >> 20;
    double f = x - 1.0;
    double s = f / (2.0 + f), dk = (double)k, z = s * s, w = z * z;
    double t1 = RB(w * (Lg2 + RB(w * (Lg4 + RB(w * Lg6)))));
    double t2 = RB(z * (Lg1 + RB(w * (Lg3 + RB(w * (Lg5 + RB(w * Lg7)))))));
    int32_t sel = (hx - 0x6147a) | (0x6b851 - hx);
    double R = t2 + t1;
    double hfsq = RB(RB(0.5 * f) * f);
    double A = RB(dk * ln2_hi) - ((hfsq - (RB(s * (hfsq + R)) + RB(dk * ln2_lo))) - f);
    double B = RB(dk * ln2_hi) - ((RB(s * (f - R)) - RB(dk * ln2_lo)) - f);
    return sel > 0 ? A : B;
}
static inline double u64_to_d(uint64_t w) {   /* correctly rounded, via exact halves */
    return RB((double)(uint32_t)(w >> 32) * 4294967296.0) + (double)(uint32_t)w;
}
static inline void bm_pair(uint64_t w1, uint64_t w2, double *z1, double *z2) {
    double a = d_of((w1 >> 12) | 0x3ff0000000000000ULL) - 1.0;   /* [0, 1) exact */
    double b = 1.0 - a;                                            /* (0, 1] exact */
    double s = sinq(a), c = sinq(b);
    s = (w1 & 1) ? -s : s;
    c = (w1 & 2) ? -c : c;
    double u = (u64_to_d(w2) + 0.5) * 0x1p-64;
    double r = sqrt(-2.0 * logu(u));
    *z1 = r * s; *z2 = r * c;
}
__attribute__((noinline)) static void bs(const uint64_t *buf, double *o, int m) {
    for (int j = 0; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}

/* ---------------- vector paths ---------------- */
#if defined(__x86_64__)
#define V_BODY(VT, VI, P, I, W)                                                          \
    const VT one = P##set1_pd(1.0), two = P##set1_pd(2.0), half = P##set1_pd(0.5);       \
    VT a = P##sub_pd(P##castsi##W##_pd(I##or_si##W(I##srli_epi64(w1, 12),                \
                     I##set1_epi64x(0x3ff0000000000000LL))), one);                       \
    VT b = P##sub_pd(one, a);                                                            \
    VT ta = P##mul_pd(a, a), tb = P##mul_pd(b, b);                                       \
    VT ps = P##set1_pd(C[8]), pc = P##set1_pd(C[8]);                                     \
    for (int k = 7; k >= 0; k--) {                                                       \
        ps = MA(ps, ta, P##set1_pd(C[k]));                                               \
        pc = MA(pc, tb, P##set1_pd(C[k]));                                               \
    }                                                                                    \
    VT s = P##mul_pd(ps, a), c = P##mul_pd(pc, b);                                       \
    VT sgn = P##castsi##W##_pd(I##slli_epi64(w1, 63));                                   \
    s = P##xor_pd(s, sgn);                                                               \
    c = P##xor_pd(c, P##castsi##W##_pd(I##slli_epi64(I##srli_epi64(w1, 1), 63)));        \
    /* u64 -> double: hi * 2^32 + lo, each half exact via magic numbers */               \
    VT hi = P##sub_pd(P##castsi##W##_pd(I##or_si##W(I##srli_epi64(w2, 32),               \
                      I##set1_epi64x(0x4530000000000000LL))), P##set1_pd(0x1p84));       \
    VT lo = P##sub_pd(P##castsi##W##_pd(I##or_si##W(I##and_si##W(w2,                     \
                      I##set1_epi64x(0xffffffffLL)), I##set1_epi64x(0x4330000000000000LL))), \
                      P##set1_pd(0x1p52));                                               \
    VT u = P##mul_pd(P##add_pd(P##add_pd(hi, lo), half), P##set1_pd(0x1p-64));           \
    /* log */                                                                            \
    VI ub = P##castpd_si##W(u);                                                          \
    VI hx = I##srli_epi64(ub, 32);                                                       \
    VI k = I##sub_epi64(I##srli_epi64(hx, 20), I##set1_epi64x(1023));                    \
    hx = I##and_si##W(hx, I##set1_epi64x(0x000fffff));                                   \
    VI i = I##and_si##W(I##add_epi64(hx, I##set1_epi64x(0x95f64)), I##set1_epi64x(0x100000)); \
    VI newhi = I##or_si##W(hx, I##xor_si##W(i, I##set1_epi64x(0x3ff00000)));             \
    VT x = P##castsi##W##_pd(I##or_si##W(I##slli_epi64(newhi, 32),                       \
                      I##and_si##W(ub, I##set1_epi64x(0xffffffffLL))));                  \
    k = I##add_epi64(k, I##srli_epi64(i, 20));                                           \
    VT f = P##sub_pd(x, one);                                                            \
    VT sl = P##div_pd(f, P##add_pd(two, f));                                             \
    VT dk = P##sub_pd(P##castsi##W##_pd(I##add_epi64(k,                                  \
                      I##set1_epi64x(0x4330000000100000LL))), P##set1_pd(0x1p52 + 0x1p20)); \
    VT z = P##mul_pd(sl, sl), ww = P##mul_pd(z, z);                                      \
    VT t1 = P##mul_pd(ww, MA(ww, MA(ww, P##set1_pd(Lg6), P##set1_pd(Lg4)),              \
                              P##set1_pd(Lg2)));                                         \
    VT t2 = P##mul_pd(z, MA(ww, MA(ww, MA(ww, P##set1_pd(Lg7), P##set1_pd(Lg5)),         \
                              P##set1_pd(Lg3)), P##set1_pd(Lg1)));                       \
    VT R = P##add_pd(t2, t1);                                                            \
    VT hfsq = P##mul_pd(P##mul_pd(half, f), f);                                          \
    VT dh = P##mul_pd(dk, P##set1_pd(ln2_hi)), dl = P##mul_pd(dk, P##set1_pd(ln2_lo));   \
    VT A = P##sub_pd(dh, P##sub_pd(P##sub_pd(hfsq, MA(sl, P##add_pd(hfsq, R), dl)), f));  \
    VT B = P##sub_pd(dh, P##sub_pd(P##sub_pd(P##mul_pd(sl, P##sub_pd(f, R)), dl), f));   \
    VI selv = I##or_si##W(I##sub_epi32(hx, I##set1_epi64x(0x6147a)),                     \
                          I##sub_epi32(I##set1_epi64x(0x6b851), hx));

#define MA(a, b, c) _mm256_add_pd(_mm256_mul_pd(a, b), c)
__attribute__((target("avx2"), noinline))
static void bv_avx2(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 8 <= m; j += 8) {
        __m256i v0 = _mm256_loadu_si256((const __m256i *)(buf + j));
        __m256i v1 = _mm256_loadu_si256((const __m256i *)(buf + j + 4));
        __m256i w1 = _mm256_permute4x64_epi64(_mm256_unpacklo_epi64(v0, v1), 0xD8);
        __m256i w2 = _mm256_permute4x64_epi64(_mm256_unpackhi_epi64(v0, v1), 0xD8);
        V_BODY(__m256d, __m256i, _mm256_, _mm256_, 256)
        /* sel > 0 as a signed 32-bit value in the low half of each lane */
        __m256i pos = _mm256_cmpgt_epi32(_mm256_and_si256(selv, _mm256_set1_epi64x(0xffffffffLL)),
                                         _mm256_setzero_si256());
        pos = _mm256_shuffle_epi32(pos, 0xA0);   /* broadcast each lane's low dword */
        __m256d ln = _mm256_blendv_pd(B, A, _mm256_castsi256_pd(pos));
        __m256d r = _mm256_sqrt_pd(_mm256_mul_pd(_mm256_set1_pd(-2.0), ln));
        __m256d z1 = _mm256_mul_pd(r, s), z2 = _mm256_mul_pd(r, c);
        __m256d e0 = _mm256_unpacklo_pd(z1, z2), e1 = _mm256_unpackhi_pd(z1, z2);
        _mm256_storeu_pd(o + j, _mm256_permute2f128_pd(e0, e1, 0x20));
        _mm256_storeu_pd(o + j + 4, _mm256_permute2f128_pd(e0, e1, 0x31));
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef MA
#define MA(a, b, c) _mm512_add_pd(_mm512_mul_pd(a, b), c)
__attribute__((target("avx512f,avx512dq"), noinline))
static void bv_avx512(const uint64_t *buf, double *o, int m) {
    int j = 0;
    const __m512i ev = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0), od = _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);
    const __m512i lo_i = _mm512_set_epi64(11, 3, 10, 2, 9, 1, 8, 0), hi_i = _mm512_set_epi64(15, 7, 14, 6, 13, 5, 12, 4);
    for (; j + 16 <= m; j += 16) {
        __m512i v0 = _mm512_loadu_si512((const void *)(buf + j));
        __m512i v1 = _mm512_loadu_si512((const void *)(buf + j + 8));
        __m512i w1 = _mm512_permutex2var_epi64(v0, ev, v1);
        __m512i w2 = _mm512_permutex2var_epi64(v0, od, v1);
#define _mm512_set1_epi64x _mm512_set1_epi64
        V_BODY(__m512d, __m512i, _mm512_, _mm512_, 512)
#undef _mm512_set1_epi64x
        __mmask8 pos = _mm512_cmpgt_epi64_mask(_mm512_srai_epi64(_mm512_slli_epi64(selv, 32), 32),
                                               _mm512_setzero_si512());
        __m512d ln = _mm512_mask_blend_pd(pos, B, A);
        __m512d r = _mm512_sqrt_pd(_mm512_mul_pd(_mm512_set1_pd(-2.0), ln));
        __m512d z1 = _mm512_mul_pd(r, s), z2 = _mm512_mul_pd(r, c);
        _mm512_storeu_pd(o + j, _mm512_permutex2var_pd(z1, lo_i, z2));
        _mm512_storeu_pd(o + j + 8, _mm512_permutex2var_pd(z1, hi_i, z2));
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef MA

#define MA(a, b, c) _mm256_fmadd_pd(a, b, c)
__attribute__((target("avx2,fma"), noinline))
static void bf_avx2(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 8 <= m; j += 8) {
        __m256i v0 = _mm256_loadu_si256((const __m256i *)(buf + j));
        __m256i v1 = _mm256_loadu_si256((const __m256i *)(buf + j + 4));
        __m256i w1 = _mm256_permute4x64_epi64(_mm256_unpacklo_epi64(v0, v1), 0xD8);
        __m256i w2 = _mm256_permute4x64_epi64(_mm256_unpackhi_epi64(v0, v1), 0xD8);
        V_BODY(__m256d, __m256i, _mm256_, _mm256_, 256)
        /* sel > 0 as a signed 32-bit value in the low half of each lane */
        __m256i pos = _mm256_cmpgt_epi32(_mm256_and_si256(selv, _mm256_set1_epi64x(0xffffffffLL)),
                                         _mm256_setzero_si256());
        pos = _mm256_shuffle_epi32(pos, 0xA0);   /* broadcast each lane's low dword */
        __m256d ln = _mm256_blendv_pd(B, A, _mm256_castsi256_pd(pos));
        __m256d r = _mm256_sqrt_pd(_mm256_mul_pd(_mm256_set1_pd(-2.0), ln));
        __m256d z1 = _mm256_mul_pd(r, s), z2 = _mm256_mul_pd(r, c);
        __m256d e0 = _mm256_unpacklo_pd(z1, z2), e1 = _mm256_unpackhi_pd(z1, z2);
        _mm256_storeu_pd(o + j, _mm256_permute2f128_pd(e0, e1, 0x20));
        _mm256_storeu_pd(o + j + 4, _mm256_permute2f128_pd(e0, e1, 0x31));
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef MA
#define MA(a, b, c) _mm512_fmadd_pd(a, b, c)
__attribute__((target("avx512f,avx512dq"), noinline))
static void bf_avx512(const uint64_t *buf, double *o, int m) {
    int j = 0;
    const __m512i ev = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0), od = _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);
    const __m512i lo_i = _mm512_set_epi64(11, 3, 10, 2, 9, 1, 8, 0), hi_i = _mm512_set_epi64(15, 7, 14, 6, 13, 5, 12, 4);
    for (; j + 16 <= m; j += 16) {
        __m512i v0 = _mm512_loadu_si512((const void *)(buf + j));
        __m512i v1 = _mm512_loadu_si512((const void *)(buf + j + 8));
        __m512i w1 = _mm512_permutex2var_epi64(v0, ev, v1);
        __m512i w2 = _mm512_permutex2var_epi64(v0, od, v1);
#define _mm512_set1_epi64x _mm512_set1_epi64
        V_BODY(__m512d, __m512i, _mm512_, _mm512_, 512)
#undef _mm512_set1_epi64x
        __mmask8 pos = _mm512_cmpgt_epi64_mask(_mm512_srai_epi64(_mm512_slli_epi64(selv, 32), 32),
                                               _mm512_setzero_si512());
        __m512d ln = _mm512_mask_blend_pd(pos, B, A);
        __m512d r = _mm512_sqrt_pd(_mm512_mul_pd(_mm512_set1_pd(-2.0), ln));
        __m512d z1 = _mm512_mul_pd(r, s), z2 = _mm512_mul_pd(r, c);
        _mm512_storeu_pd(o + j, _mm512_permutex2var_pd(z1, lo_i, z2));
        _mm512_storeu_pd(o + j + 8, _mm512_permutex2var_pd(z1, hi_i, z2));
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef MA

static int has512(void) { __builtin_cpu_init(); return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"); }
static int has2(void) { __builtin_cpu_init(); return __builtin_cpu_supports("avx2"); }
static int hasfma(void) { __builtin_cpu_init(); return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"); }
#elif defined(__aarch64__)
#define NMA(a, b, c) vaddq_f64(vmulq_f64(a, b), c)
__attribute__((noinline)) static void bv_neon(const uint64_t *buf, double *o, int m) {
    int j = 0;
    const float64x2_t one = vdupq_n_f64(1.0), two = vdupq_n_f64(2.0), half = vdupq_n_f64(0.5);
    for (; j + 4 <= m; j += 4) {
        uint64x2x2_t d = vld2q_u64(buf + j);      /* de-interleaves: val[0] evens, val[1] odds */
        uint64x2_t w1 = d.val[0], w2 = d.val[1];
        float64x2_t a = vsubq_f64(vreinterpretq_f64_u64(vorrq_u64(vshrq_n_u64(w1, 12),
                                  vdupq_n_u64(0x3ff0000000000000ULL))), one);
        float64x2_t b = vsubq_f64(one, a);
        float64x2_t ta = vmulq_f64(a, a), tb = vmulq_f64(b, b);
        float64x2_t ps = vdupq_n_f64(C[8]), pc = ps;
        for (int k = 7; k >= 0; k--) {
            ps = NMA(ps, ta, vdupq_n_f64(C[k]));
            pc = NMA(pc, tb, vdupq_n_f64(C[k]));
        }
        float64x2_t s = vmulq_f64(ps, a), c = vmulq_f64(pc, b);
        s = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(s), vshlq_n_u64(w1, 63)));
        c = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(c), vshlq_n_u64(vshrq_n_u64(w1, 1), 63)));
        float64x2_t hi = vmulq_f64(vcvtq_f64_u64(vshrq_n_u64(w2, 32)), vdupq_n_f64(4294967296.0));
        float64x2_t lo = vcvtq_f64_u64(vandq_u64(w2, vdupq_n_u64(0xffffffffu)));
        float64x2_t u = vmulq_f64(vaddq_f64(vaddq_f64(hi, lo), half), vdupq_n_f64(0x1p-64));
        uint64x2_t ub = vreinterpretq_u64_f64(u);
        int64x2_t hx = vreinterpretq_s64_u64(vshrq_n_u64(ub, 32));
        int64x2_t k = vsubq_s64(vshrq_n_s64(hx, 20), vdupq_n_s64(1023));
        hx = vandq_s64(hx, vdupq_n_s64(0x000fffff));
        int64x2_t i = vandq_s64(vaddq_s64(hx, vdupq_n_s64(0x95f64)), vdupq_n_s64(0x100000));
        int64x2_t newhi = vorrq_s64(hx, veorq_s64(i, vdupq_n_s64(0x3ff00000)));
        float64x2_t x = vreinterpretq_f64_u64(vorrq_u64(vshlq_n_u64(vreinterpretq_u64_s64(newhi), 32),
                                               vandq_u64(ub, vdupq_n_u64(0xffffffffu))));
        k = vaddq_s64(k, vshrq_n_s64(i, 20));
        float64x2_t f = vsubq_f64(x, one);
        float64x2_t sl = vdivq_f64(f, vaddq_f64(two, f));
        float64x2_t dk = vcvtq_f64_s64(k);
        float64x2_t z = vmulq_f64(sl, sl), ww = vmulq_f64(z, z);
        float64x2_t t1 = vmulq_f64(ww, NMA(ww, NMA(ww, vdupq_n_f64(Lg6), vdupq_n_f64(Lg4)), vdupq_n_f64(Lg2)));
        float64x2_t t2 = vmulq_f64(z, NMA(ww, NMA(ww, NMA(ww, vdupq_n_f64(Lg7), vdupq_n_f64(Lg5)),
                                           vdupq_n_f64(Lg3)), vdupq_n_f64(Lg1)));
        float64x2_t R = vaddq_f64(t2, t1);
        float64x2_t hfsq = vmulq_f64(vmulq_f64(half, f), f);
        float64x2_t dh = vmulq_f64(dk, vdupq_n_f64(ln2_hi)), dl = vmulq_f64(dk, vdupq_n_f64(ln2_lo));
        float64x2_t A = vsubq_f64(dh, vsubq_f64(vsubq_f64(hfsq, NMA(sl, vaddq_f64(hfsq, R), dl)), f));
        float64x2_t B = vsubq_f64(dh, vsubq_f64(vsubq_f64(vmulq_f64(sl, vsubq_f64(f, R)), dl), f));
        int64x2_t sel = vorrq_s64(vsubq_s64(hx, vdupq_n_s64(0x6147a)), vsubq_s64(vdupq_n_s64(0x6b851), hx));
        uint64x2_t pos = vcgtq_s64(vshrq_n_s64(vshlq_n_s64(sel, 32), 32), vdupq_n_s64(0));
        float64x2_t ln = vbslq_f64(pos, A, B);
        float64x2_t r = vsqrtq_f64(vmulq_f64(vdupq_n_f64(-2.0), ln));
        float64x2x2_t out = {{vmulq_f64(r, s), vmulq_f64(r, c)}};
        vst2q_f64(o + j, out);                    /* re-interleaves */
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef NMA

#define NMA(a, b, c) vfmaq_f64(c, a, b)
__attribute__((noinline)) static void bf_neon(const uint64_t *buf, double *o, int m) {
    int j = 0;
    const float64x2_t one = vdupq_n_f64(1.0), two = vdupq_n_f64(2.0), half = vdupq_n_f64(0.5);
    for (; j + 4 <= m; j += 4) {
        uint64x2x2_t d = vld2q_u64(buf + j);      /* de-interleaves: val[0] evens, val[1] odds */
        uint64x2_t w1 = d.val[0], w2 = d.val[1];
        float64x2_t a = vsubq_f64(vreinterpretq_f64_u64(vorrq_u64(vshrq_n_u64(w1, 12),
                                  vdupq_n_u64(0x3ff0000000000000ULL))), one);
        float64x2_t b = vsubq_f64(one, a);
        float64x2_t ta = vmulq_f64(a, a), tb = vmulq_f64(b, b);
        float64x2_t ps = vdupq_n_f64(C[8]), pc = ps;
        for (int k = 7; k >= 0; k--) {
            ps = NMA(ps, ta, vdupq_n_f64(C[k]));
            pc = NMA(pc, tb, vdupq_n_f64(C[k]));
        }
        float64x2_t s = vmulq_f64(ps, a), c = vmulq_f64(pc, b);
        s = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(s), vshlq_n_u64(w1, 63)));
        c = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(c), vshlq_n_u64(vshrq_n_u64(w1, 1), 63)));
        float64x2_t hi = vmulq_f64(vcvtq_f64_u64(vshrq_n_u64(w2, 32)), vdupq_n_f64(4294967296.0));
        float64x2_t lo = vcvtq_f64_u64(vandq_u64(w2, vdupq_n_u64(0xffffffffu)));
        float64x2_t u = vmulq_f64(vaddq_f64(vaddq_f64(hi, lo), half), vdupq_n_f64(0x1p-64));
        uint64x2_t ub = vreinterpretq_u64_f64(u);
        int64x2_t hx = vreinterpretq_s64_u64(vshrq_n_u64(ub, 32));
        int64x2_t k = vsubq_s64(vshrq_n_s64(hx, 20), vdupq_n_s64(1023));
        hx = vandq_s64(hx, vdupq_n_s64(0x000fffff));
        int64x2_t i = vandq_s64(vaddq_s64(hx, vdupq_n_s64(0x95f64)), vdupq_n_s64(0x100000));
        int64x2_t newhi = vorrq_s64(hx, veorq_s64(i, vdupq_n_s64(0x3ff00000)));
        float64x2_t x = vreinterpretq_f64_u64(vorrq_u64(vshlq_n_u64(vreinterpretq_u64_s64(newhi), 32),
                                               vandq_u64(ub, vdupq_n_u64(0xffffffffu))));
        k = vaddq_s64(k, vshrq_n_s64(i, 20));
        float64x2_t f = vsubq_f64(x, one);
        float64x2_t sl = vdivq_f64(f, vaddq_f64(two, f));
        float64x2_t dk = vcvtq_f64_s64(k);
        float64x2_t z = vmulq_f64(sl, sl), ww = vmulq_f64(z, z);
        float64x2_t t1 = vmulq_f64(ww, NMA(ww, NMA(ww, vdupq_n_f64(Lg6), vdupq_n_f64(Lg4)), vdupq_n_f64(Lg2)));
        float64x2_t t2 = vmulq_f64(z, NMA(ww, NMA(ww, NMA(ww, vdupq_n_f64(Lg7), vdupq_n_f64(Lg5)),
                                           vdupq_n_f64(Lg3)), vdupq_n_f64(Lg1)));
        float64x2_t R = vaddq_f64(t2, t1);
        float64x2_t hfsq = vmulq_f64(vmulq_f64(half, f), f);
        float64x2_t dh = vmulq_f64(dk, vdupq_n_f64(ln2_hi)), dl = vmulq_f64(dk, vdupq_n_f64(ln2_lo));
        float64x2_t A = vsubq_f64(dh, vsubq_f64(vsubq_f64(hfsq, NMA(sl, vaddq_f64(hfsq, R), dl)), f));
        float64x2_t B = vsubq_f64(dh, vsubq_f64(vsubq_f64(vmulq_f64(sl, vsubq_f64(f, R)), dl), f));
        int64x2_t sel = vorrq_s64(vsubq_s64(hx, vdupq_n_s64(0x6147a)), vsubq_s64(vdupq_n_s64(0x6b851), hx));
        uint64x2_t pos = vcgtq_s64(vshrq_n_s64(vshlq_n_s64(sel, 32), 32), vdupq_n_s64(0));
        float64x2_t ln = vbslq_f64(pos, A, B);
        float64x2_t r = vsqrtq_f64(vmulq_f64(vdupq_n_f64(-2.0), ln));
        float64x2x2_t out = {{vmulq_f64(r, s), vmulq_f64(r, c)}};
        vst2q_f64(o + j, out);                    /* re-interleaves */
    }
    for (; j + 2 <= m; j += 2) bm_pair(buf[j], buf[j + 1], o + j, o + j + 1);
}
#undef NMA
#endif

/* ---------------- comparisons ---------------- */
static const double *MX;   /* McFarland layer lengths, built at start */
static double mx[254];
__attribute__((noinline)) static double slow(uint64_t index, uint64_t r) {
    double x = (double)((r >> 9) & 0xfffffffffffffULL) * wi_double[r & 0xff];
    for (int i = 0; i < 10; i++) x = x * 0.999 + (double)(index & 7) * 1e-9;
    return (r >> 8) & 1 ? -x : x;
}
static inline double fast_m(uint64_t index, uint64_t r) {
    unsigned i = (unsigned)(r & 0xff);
    if (__builtin_expect(i < 253, 1)) return MX[i] * (double)(int64_t)r;
    return slow(index, r);
}
__attribute__((noinline)) static void m1(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j] = fast_m(j, buf[j]); o[j + 1] = fast_m(j + 1, buf[j + 1]);
        o[j + 2] = fast_m(j + 2, buf[j + 2]); o[j + 3] = fast_m(j + 3, buf[j + 3]);
    }
}
static double wis[512]; static uint64_t kim1[256];
static inline double fast_z(uint64_t index, uint64_t r) {
    uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
    double x = (double)rabs * wis[r & 0x1ff];
    if (__builtin_expect(rabs - 1 < kim1[r & 0xff], 1)) return x;
    return slow(index, r);
}
__attribute__((noinline)) static void z6(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j] = fast_z(j, buf[j]); o[j + 1] = fast_z(j + 1, buf[j + 1]);
        o[j + 2] = fast_z(j + 2, buf[j + 2]); o[j + 3] = fast_z(j + 3, buf[j + 3]);
    }
}

typedef void (*fn)(const uint64_t *, double *, int);

int main(int argc, char **argv) {
    int reps = argc > 1 ? atoi(argv[1]) : 10000;
    for (int i = 0; i < 256; i++) { wis[i] = wi_double[i]; wis[256 + i] = -wi_double[i]; kim1[i] = ki_double[i] ? ki_double[i] - 1 : 0; }
    for (int i = 0; i < 253; i++) mx[i] = wi_double[i] * 0x1p-11;   /* timing stand-in with the right shape */
    MX = mx;
    uint64_t *buf = malloc(NW * 8); double *o = malloc(NW * 8), *ref = malloc(NW * 8);
    uint64_t s = 0x9e3779b97f4a7c15ULL;
    for (int j = 0; j < NW; j++) {
        uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        buf[j] = z ^ (z >> 31);
    }
    /* accuracy of the pieces against libm (statistics only; libm is not the reference) */
    double me = 0, ml = 0;
    for (int k = 0; k <= 100000; k++) {
        double x = k / 100000.0;
        double e = fabs(sinq(x) - sin(M_PI / 2 * x)); if (e > me) me = e;
        double u = (k + 0.5) / 100001.0, l = fabs(logu(u) - log(u)) / fabs(log(u)); if (l > ml) ml = l;
    }
    printf("max |sinq - sin| = %.2e;  max rel |logu - log| = %.2e\n", me, ml);
    bs(buf, ref, NW);
    double sm = 0, sq = 0; for (int j = 0; j < NW; j++) { sm += ref[j]; sq += ref[j] * ref[j]; }
    printf("sanity on %d values: mean %+.4f  var %.4f\n", NW, sm / NW, sq / NW - (sm / NW) * (sm / NW));
    struct { const char *name; fn f; int ok; } v[8]; int nv = 0;
    v[nv++] = (typeof(v[0])){"Z6 shipped ziggurat fast path", z6, 1};
    v[nv++] = (typeof(v[0])){"M1 McFarland fast path", m1, 1};
    v[nv++] = (typeof(v[0])){"BS Box-Muller scalar", bs, 1};
#if defined(__x86_64__)
    v[nv++] = (typeof(v[0])){"BV Box-Muller AVX2, 4 lanes", bv_avx2, has2()};
    v[nv++] = (typeof(v[0])){"BV Box-Muller AVX-512, 8 lanes", bv_avx512, has512()};
    v[nv++] = (typeof(v[0])){"BF Box-Muller AVX2+FMA, 4 lanes", bf_avx2, hasfma()};
    v[nv++] = (typeof(v[0])){"BF Box-Muller AVX-512+FMA, 8 lanes", bf_avx512, has512()};
#elif defined(__aarch64__)
    v[nv++] = (typeof(v[0])){"BV Box-Muller NEON, 2 lanes", bv_neon, 1};
    v[nv++] = (typeof(v[0])){"BF Box-Muller NEON+FMA, 2 lanes", bf_neon, 1};
#endif
    for (int i = 3; i < nv; i++) if (v[i].ok) {
        v[i].f(buf, o, NW);
        int bad = 0; double md = 0;
        for (int j = 0; j < NW; j++) { bad += b_of(o[j]) != b_of(ref[j]); double d = fabs(o[j] - ref[j]); if (d > md) md = d; }
        printf("%-36s identical to scalar: %s (%d of %d differ, max |diff| %.1e)\n", v[i].name, bad ? "NO" : "yes", bad, NW, md);
    }
    for (int round = 1; round <= 3; round++) {
        printf("-- round %d (ns per normal, words in cache)\n", round);
        for (int i = 0; i < nv; i++) {
            if (!v[i].ok) { printf("  %-36s skipped\n", v[i].name); continue; }
            v[i].f(buf, o, NW);
            double t = now();
            for (int r = 0; r < reps; r++) v[i].f(buf, o, NW);
            printf("  %-36s %.3f\n", v[i].name, (now() - t) / reps / NW * 1e9);
        }
    }
    return 0;
}
