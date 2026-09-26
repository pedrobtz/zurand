/* AVX-512 for the xoshiro256pp generator: 8 sub-chunks in 8 lanes, the
 * native 64-bit rotate (vprolq) instead of shift-shift-or, and 32 vector
 * registers. One 8-lane group is exactly one chunk (8 x 512 words), so it
 * slots into the shipped chunk structure with no stream change: words are
 * asserted equal to the AVX2 path's before timing. randompack ships an
 * AVX-512 path; zurand does not, yet.
 *
 *   U2  uniform, AVX2 4 lanes, fused u01 (shipped)
 *   U5  uniform, AVX-512 8 lanes, fused u01
 *   A2  normal, AVX2 words -> buffer -> ziggurat transform (shipped)
 *   A5  normal, AVX-512 words -> buffer -> the same transform
 *
 * Each function carries its own target attribute and the file is compiled
 * without -m flags, as the package is; the AVX-512 variants run only if
 * the CPU has avx512f and avx512dq.
 *
 *   cc -O2 -DR123_USE_MULHILO64_C99=1 -I../src -o /tmp/a512 ../dev/simd/avx512.c && /tmp/a512 10000000
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>
#include "Random123/philox.h"
#include "numpyzig/ziggurat_constants.h"

#define SUB 512
#define CHUNK (8 * SUB)

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }
static philox4x64_key_t PK = {{0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL}};

static double wis[512];
static uint64_t kim1[256];
static void tables_init(void) {
    for (int i = 0; i < 256; i++) {
        wis[i] = wi_double[i]; wis[256 + i] = -wi_double[i];
        kim1[i] = ki_double[i] ? ki_double[i] - 1 : 0;
    }
}
__attribute__((noinline)) static double slow(uint64_t index, uint64_t r) {
    double x = (double)((r >> 9) & 0xfffffffffffffULL) * wi_double[r & 0xff];
    for (int i = 0; i < 10; i++) x = x * 0.999 + (double)(index & 7) * 1e-9;
    return (r >> 8) & 1 ? -x : x;
}
__attribute__((noinline, cold)) static double reject(uint64_t index, uint64_t r) {
    if (((r >> 9) & 0xfffffffffffffULL) == 0 && ki_double[r & 0xff] != 0) return 0.0;
    return slow(index, r);
}
static inline double fast(uint64_t index, uint64_t r) {
    uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
    double x = (double)rabs * wis[r & 0x1ff];
    if (__builtin_expect(rabs - 1 < kim1[r & 0xff], 1)) return x;
    return reject(index, r);
}
static void transform(const uint64_t *buf, double *o, int m, uint64_t w0) {
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j]     = fast(w0 + j,     buf[j]);
        o[j + 1] = fast(w0 + j + 1, buf[j + 1]);
        o[j + 2] = fast(w0 + j + 2, buf[j + 2]);
        o[j + 3] = fast(w0 + j + 3, buf[j + 3]);
    }
    for (; j < m; j++) o[j] = fast(w0 + j, buf[j]);
}
static inline void seed(uint64_t sub, uint64_t *s) {
    philox4x64_ctr_t c = {{sub, 0, 2, 1}};
    philox4x64_ctr_t r = philox4x64(c, PK);
    s[0] = r.v[0]; s[1] = r.v[1]; s[2] = r.v[2]; s[3] = r.v[3];
    if (!(s[0] | s[1] | s[2] | s[3])) s[0] = 1;
}

/* ---------------- AVX2, 4 lanes (shipped) ---------------- */
#define VROTL2(x, k) _mm256_or_si256(_mm256_slli_epi64((x), (k)), _mm256_srli_epi64((x), 64 - (k)))
#define LOAD2(sub0)                                                          \
    uint64_t st[4][4];                                                       \
    for (int l = 0; l < 4; l++) seed((sub0) + l, st[l]);                     \
    __m256i s0 = _mm256_set_epi64x(st[3][0], st[2][0], st[1][0], st[0][0]);  \
    __m256i s1 = _mm256_set_epi64x(st[3][1], st[2][1], st[1][1], st[0][1]);  \
    __m256i s2 = _mm256_set_epi64x(st[3][2], st[2][2], st[1][2], st[0][2]);  \
    __m256i s3 = _mm256_set_epi64x(st[3][3], st[2][3], st[1][3], st[0][3])
#define STEP2(v)                                                             \
    for (int k = 0; k < 4; k++) {                                            \
        v[k] = _mm256_add_epi64(VROTL2(_mm256_add_epi64(s0, s3), 23), s0);   \
        __m256i t = _mm256_slli_epi64(s1, 17);                               \
        s2 = _mm256_xor_si256(s2, s0); s3 = _mm256_xor_si256(s3, s1);        \
        s1 = _mm256_xor_si256(s1, s2); s0 = _mm256_xor_si256(s0, s3);        \
        s2 = _mm256_xor_si256(s2, t);  s3 = VROTL2(s3, 45);                  \
    }
#define TR2(v, r0, r1, r2, r3)                                               \
    __m256i t0 = _mm256_unpacklo_epi64(v[0], v[1]);                          \
    __m256i t1 = _mm256_unpackhi_epi64(v[0], v[1]);                          \
    __m256i t2 = _mm256_unpacklo_epi64(v[2], v[3]);                          \
    __m256i t3 = _mm256_unpackhi_epi64(v[2], v[3]);                          \
    __m256i r0 = _mm256_permute2x128_si256(t0, t2, 0x20);                    \
    __m256i r1 = _mm256_permute2x128_si256(t1, t3, 0x20);                    \
    __m256i r2 = _mm256_permute2x128_si256(t0, t2, 0x31);                    \
    __m256i r3 = _mm256_permute2x128_si256(t1, t3, 0x31)

__attribute__((target("avx2")))
static void words_avx2(uint64_t sub0, uint64_t *buf) {
    LOAD2(sub0);
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP2(v); TR2(v, r0, r1, r2, r3);
        _mm256_storeu_si256((__m256i *)(buf + 0 * SUB + j), r0);
        _mm256_storeu_si256((__m256i *)(buf + 1 * SUB + j), r1);
        _mm256_storeu_si256((__m256i *)(buf + 2 * SUB + j), r2);
        _mm256_storeu_si256((__m256i *)(buf + 3 * SUB + j), r3);
    }
}
__attribute__((target("avx2")))
static void u01_avx2(uint64_t sub0, double *out) {
    LOAD2(sub0);
    const __m256i exponent = _mm256_set1_epi64x(0x3ff0000000000000LL);
    const __m256d bias = _mm256_set1_pd(1.0 - 0x1.0p-53);
#define U01(w) _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_srli_epi64((w), 12), exponent)), bias)
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP2(v); TR2(v, r0, r1, r2, r3);
        _mm256_storeu_pd(out + 0 * SUB + j, U01(r0));
        _mm256_storeu_pd(out + 1 * SUB + j, U01(r1));
        _mm256_storeu_pd(out + 2 * SUB + j, U01(r2));
        _mm256_storeu_pd(out + 3 * SUB + j, U01(r3));
    }
#undef U01
}

/* ---------------- AVX-512, 8 lanes ---------------- */
#define LOAD5(sub0)                                                          \
    uint64_t st[8][4];                                                       \
    for (int l = 0; l < 8; l++) seed((sub0) + l, st[l]);                     \
    __m512i s0 = _mm512_set_epi64(st[7][0], st[6][0], st[5][0], st[4][0],    \
                                  st[3][0], st[2][0], st[1][0], st[0][0]);   \
    __m512i s1 = _mm512_set_epi64(st[7][1], st[6][1], st[5][1], st[4][1],    \
                                  st[3][1], st[2][1], st[1][1], st[0][1]);   \
    __m512i s2 = _mm512_set_epi64(st[7][2], st[6][2], st[5][2], st[4][2],    \
                                  st[3][2], st[2][2], st[1][2], st[0][2]);   \
    __m512i s3 = _mm512_set_epi64(st[7][3], st[6][3], st[5][3], st[4][3],    \
                                  st[3][3], st[2][3], st[1][3], st[0][3])
#define STEP5(v)                                                             \
    for (int k = 0; k < 8; k++) {                                            \
        v[k] = _mm512_add_epi64(_mm512_rol_epi64(_mm512_add_epi64(s0, s3), 23), s0); \
        __m512i t = _mm512_slli_epi64(s1, 17);                               \
        s2 = _mm512_xor_si512(s2, s0); s3 = _mm512_xor_si512(s3, s1);        \
        s1 = _mm512_xor_si512(s1, s2); s0 = _mm512_xor_si512(s0, s3);        \
        s2 = _mm512_xor_si512(s2, t);  s3 = _mm512_rol_epi64(s3, 45);        \
    }
/* 8x8 transpose of 64-bit lanes: rows = steps in, rows = sub-chunks out.
 * unpack pairs, then two rounds of 128-bit lane shuffles. */
#define TR5(v, r)                                                            \
    __m512i t[8], u[8];                                                      \
    for (int i = 0; i < 4; i++) {                                            \
        t[2 * i]     = _mm512_unpacklo_epi64(v[2 * i], v[2 * i + 1]);        \
        t[2 * i + 1] = _mm512_unpackhi_epi64(v[2 * i], v[2 * i + 1]);        \
    }                                                                        \
    u[0] = _mm512_shuffle_i64x2(t[0], t[2], 0x88);                           \
    u[1] = _mm512_shuffle_i64x2(t[0], t[2], 0xDD);                           \
    u[2] = _mm512_shuffle_i64x2(t[1], t[3], 0x88);                           \
    u[3] = _mm512_shuffle_i64x2(t[1], t[3], 0xDD);                           \
    u[4] = _mm512_shuffle_i64x2(t[4], t[6], 0x88);                           \
    u[5] = _mm512_shuffle_i64x2(t[4], t[6], 0xDD);                           \
    u[6] = _mm512_shuffle_i64x2(t[5], t[7], 0x88);                           \
    u[7] = _mm512_shuffle_i64x2(t[5], t[7], 0xDD);                           \
    r[0] = _mm512_shuffle_i64x2(u[0], u[4], 0x88);                           \
    r[4] = _mm512_shuffle_i64x2(u[0], u[4], 0xDD);                           \
    r[2] = _mm512_shuffle_i64x2(u[1], u[5], 0x88);                           \
    r[6] = _mm512_shuffle_i64x2(u[1], u[5], 0xDD);                           \
    r[1] = _mm512_shuffle_i64x2(u[2], u[6], 0x88);                           \
    r[5] = _mm512_shuffle_i64x2(u[2], u[6], 0xDD);                           \
    r[3] = _mm512_shuffle_i64x2(u[3], u[7], 0x88);                           \
    r[7] = _mm512_shuffle_i64x2(u[3], u[7], 0xDD)

__attribute__((target("avx512f,avx512dq")))
static void words_avx512(uint64_t sub0, uint64_t *buf) {
    LOAD5(sub0);
    for (int j = 0; j < SUB; j += 8) {
        __m512i v[8], r[8]; STEP5(v); TR5(v, r);
        for (int l = 0; l < 8; l++)
            _mm512_storeu_si512((void *)(buf + l * SUB + j), r[l]);
    }
}
__attribute__((target("avx512f,avx512dq")))
static void u01_avx512(uint64_t sub0, double *out) {
    LOAD5(sub0);
    const __m512i exponent = _mm512_set1_epi64(0x3ff0000000000000LL);
    const __m512d bias = _mm512_set1_pd(1.0 - 0x1.0p-53);
    for (int j = 0; j < SUB; j += 8) {
        __m512i v[8], r[8]; STEP5(v); TR5(v, r);
        for (int l = 0; l < 8; l++)
            _mm512_storeu_pd(out + l * SUB + j,
                _mm512_sub_pd(_mm512_castsi512_pd(_mm512_or_si512(
                    _mm512_srli_epi64(r[l], 12), exponent)), bias));
    }
}

typedef void (*fill_fn)(double *, long);
static void fill_U2(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        u01_avx2(c * 8, out + c * CHUNK);
        u01_avx2(c * 8 + 4, out + c * CHUNK + 4 * SUB);
    }
}
static void fill_U5(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) u01_avx512(c * 8, out + c * CHUNK);
}
static void fill_A2(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        uint64_t buf[CHUNK];
        words_avx2(c * 8, buf); words_avx2(c * 8 + 4, buf + 4 * SUB);
        transform(buf, out + c * CHUNK, CHUNK, (uint64_t)c * CHUNK);
    }
}
static void fill_A5(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        uint64_t buf[CHUNK];
        words_avx512(c * 8, buf);
        transform(buf, out + c * CHUNK, CHUNK, (uint64_t)c * CHUNK);
    }
}

int main(int argc, char **argv) {
    long n = argc > 1 ? atol(argv[1]) : 10000000L; n -= n % CHUNK;
    int reps = argc > 2 ? atoi(argv[2]) : 3;
    __builtin_cpu_init();
    if (!(__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))) {
        printf("no AVX-512 on this CPU\n"); return 0;
    }
    tables_init();
    { uint64_t a[CHUNK], b[CHUNK];
      words_avx2(800, a); words_avx2(804, a + 4 * SUB); words_avx512(800, b);
      printf("AVX-512 words == AVX2 words, chunk 100: %s\n", memcmp(a, b, sizeof a) ? "NO" : "yes"); }
    double *ref = malloc((size_t)n * 8), *out = malloc((size_t)n * 8);
    fill_U2(ref, n); fill_U5(out, n);
    printf("U5 == U2: %s\n", memcmp(ref, out, (size_t)n * 8) ? "NO" : "yes");
    fill_A2(ref, n); fill_A5(out, n);
    printf("A5 == A2: %s\n", memcmp(ref, out, (size_t)n * 8) ? "NO" : "yes");
    struct { const char *name; fill_fn f; } v[] = {
        {"U2 uniform AVX2 (shipped)", fill_U2}, {"U5 uniform AVX-512", fill_U5},
        {"A2 normal AVX2 words + transform", fill_A2}, {"A5 normal AVX-512 words + transform", fill_A5}};
    /* round robin over the variants so no variant owns a warm or cold
     * phase of the process; median over rounds of the mean of `reps` */
    enum { ROUNDS = 7 };
    double tr[4][ROUNDS], tf[4][ROUNDS], acc = 0;
    for (int i = 0; i < 4; i++) v[i].f(out, n);
    for (int rd = 0; rd < ROUNDS; rd++)
        for (int i = 0; i < 4; i++) {
            double t = now(); for (int r = 0; r < reps; r++) v[i].f(out, n); tr[i][rd] = (now() - t) / reps / n * 1e9;
            t = now();
            for (int r = 0; r < reps; r++) { double *p = malloc((size_t)n * 8); v[i].f(p, n); acc += p[n - 1]; free(p); }
            tf[i][rd] = (now() - t) / reps / n * 1e9;
        }
    printf("n = %ld, %d rounds x %d reps; ns per value, median [min-max]\n%-40s %22s %22s\n", n, ROUNDS, reps, "", "reused", "fresh malloc");
    for (int i = 0; i < 4; i++) {
        double mr = 0, mf = 0, lr = 1e9, hr = 0, lf = 1e9, hf = 0;
        for (int rd = 0; rd < ROUNDS; rd++) {   /* insertion sort copies for the medians */
            if (tr[i][rd] < lr) lr = tr[i][rd]; if (tr[i][rd] > hr) hr = tr[i][rd];
            if (tf[i][rd] < lf) lf = tf[i][rd]; if (tf[i][rd] > hf) hf = tf[i][rd];
        }
        double cr[ROUNDS], cf[ROUNDS];
        memcpy(cr, tr[i], sizeof cr); memcpy(cf, tf[i], sizeof cf);
        for (int x = 1; x < ROUNDS; x++) { double k = cr[x]; int y = x - 1; while (y >= 0 && cr[y] > k) { cr[y + 1] = cr[y]; y--; } cr[y + 1] = k; }
        for (int x = 1; x < ROUNDS; x++) { double k = cf[x]; int y = x - 1; while (y >= 0 && cf[y] > k) { cf[y + 1] = cf[y]; y--; } cf[y + 1] = k; }
        mr = cr[ROUNDS / 2]; mf = cf[ROUNDS / 2];
        printf("%-40s %7.3f [%.3f-%.3f] %7.3f [%.3f-%.3f]\n", v[i].name, mr, lr, hr, mf, lf, hf);
    }
    printf("(%g)\n", acc);
    free(ref); free(out);
    return 0;
}
