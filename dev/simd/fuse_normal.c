/* Where the normal fill's time goes on x86, and whether fusing the AVX2
 * generator with the ziggurat transform buys anything. Same words, same
 * bits, same fast path (the shipped N6 structure); the rare slow path is a
 * stand-in shared by every variant, as in zig_compare.c.
 *
 *   A  shipped structure: 8 sub-chunks (two AVX2 groups) -> 32 KiB word
 *      buffer -> transform pass -> out
 *   B  fused: transform straight from the transposed vectors, no buffer;
 *      rejects recorded and fixed after the group (the AVX2 code stays a
 *      leaf, which MinGW needs)
 *   D  strip-mined: the buffer shrunk to 2 KiB (16 steps x 4 lanes) so
 *      buffer and tables both stay in L1; plain calls, no leaf constraint
 *   U  uniform fill (the shipped fused u01 loop), to place the numbers
 *
 * Each variant runs into a reused buffer (generation only) and into a
 * fresh malloc per repetition (what R does: page faults included).
 *
 *   clang -O2 -mavx2 -DR123_USE_MULHILO64_C99=1 -I../src \
 *         -o /tmp/fuse ../dev/simd/fuse_normal.c && /tmp/fuse 10000000
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
#define LANES 4
#define CHUNK (2 * LANES * SUB)          /* 4096 words, as shipped */

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
/* N6: the shipped transform */
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
#define VROTL(x, k) _mm256_or_si256(_mm256_slli_epi64((x), (k)), _mm256_srli_epi64((x), 64 - (k)))
#define LOAD_STATE(sub0)                                                     \
    uint64_t st[LANES][4];                                                   \
    for (int l = 0; l < LANES; l++) seed((sub0) + l, st[l]);                 \
    __m256i s0 = _mm256_set_epi64x(st[3][0], st[2][0], st[1][0], st[0][0]);  \
    __m256i s1 = _mm256_set_epi64x(st[3][1], st[2][1], st[1][1], st[0][1]);  \
    __m256i s2 = _mm256_set_epi64x(st[3][2], st[2][2], st[1][2], st[0][2]);  \
    __m256i s3 = _mm256_set_epi64x(st[3][3], st[2][3], st[1][3], st[0][3])
#define STEP4(v)                                                             \
    for (int k = 0; k < 4; k++) {                                            \
        v[k] = _mm256_add_epi64(VROTL(_mm256_add_epi64(s0, s3), 23), s0);    \
        __m256i t = _mm256_slli_epi64(s1, 17);                               \
        s2 = _mm256_xor_si256(s2, s0); s3 = _mm256_xor_si256(s3, s1);        \
        s1 = _mm256_xor_si256(s1, s2); s0 = _mm256_xor_si256(s0, s3);        \
        s2 = _mm256_xor_si256(s2, t);  s3 = VROTL(s3, 45);                   \
    }
#define TRANSPOSE(v, r0, r1, r2, r3)                                         \
    __m256i t0 = _mm256_unpacklo_epi64(v[0], v[1]);                          \
    __m256i t1 = _mm256_unpackhi_epi64(v[0], v[1]);                          \
    __m256i t2 = _mm256_unpacklo_epi64(v[2], v[3]);                          \
    __m256i t3 = _mm256_unpackhi_epi64(v[2], v[3]);                          \
    __m256i r0 = _mm256_permute2x128_si256(t0, t2, 0x20);                    \
    __m256i r1 = _mm256_permute2x128_si256(t1, t3, 0x20);                    \
    __m256i r2 = _mm256_permute2x128_si256(t0, t2, 0x31);                    \
    __m256i r3 = _mm256_permute2x128_si256(t1, t3, 0x31)

/* words of 4 sub-chunks -> buf, sub-chunk l at buf + l*SUB (shipped) */
static void group_words(uint64_t sub0, uint64_t *buf) {
    LOAD_STATE(sub0);
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP4(v); TRANSPOSE(v, r0, r1, r2, r3);
        _mm256_storeu_si256((__m256i *)(buf + 0 * SUB + j), r0);
        _mm256_storeu_si256((__m256i *)(buf + 1 * SUB + j), r1);
        _mm256_storeu_si256((__m256i *)(buf + 2 * SUB + j), r2);
        _mm256_storeu_si256((__m256i *)(buf + 3 * SUB + j), r3);
    }
}
/* shipped fused u01 */
static void group_u01(uint64_t sub0, double *out) {
    LOAD_STATE(sub0);
    const __m256i exponent = _mm256_set1_epi64x(0x3ff0000000000000LL);
    const __m256d bias = _mm256_set1_pd(1.0 - 0x1.0p-53);
#define U01(w) _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_srli_epi64((w), 12), exponent)), bias)
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP4(v); TRANSPOSE(v, r0, r1, r2, r3);
        _mm256_storeu_pd(out + 0 * SUB + j, U01(r0));
        _mm256_storeu_pd(out + 1 * SUB + j, U01(r1));
        _mm256_storeu_pd(out + 2 * SUB + j, U01(r2));
        _mm256_storeu_pd(out + 3 * SUB + j, U01(r3));
    }
#undef U01
}
/* B: fused normal. Leaf: rejects go to a list, fixed by the caller. */
static int group_normal_fused(uint64_t sub0, double *out, uint16_t *rej, uint64_t *rejw) {
    LOAD_STATE(sub0);
    int nrej = 0;
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP4(v); TRANSPOSE(v, r0, r1, r2, r3);
        __m256i rows[4] = {r0, r1, r2, r3};
        for (int l = 0; l < 4; l++) {
            uint64_t w[4];
            _mm256_storeu_si256((__m256i *)w, rows[l]);   /* 4 words, stays in registers/L1 */
            double *o = out + l * SUB + j;
            for (int k = 0; k < 4; k++) {
                uint64_t r = w[k];
                uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
                o[k] = (double)rabs * wis[r & 0x1ff];
                if (__builtin_expect(rabs - 1 >= kim1[r & 0xff], 0)) {
                    rej[nrej] = (uint16_t)(l * SUB + j + k); rejw[nrej] = r; nrej++;
                }
            }
        }
    }
    return nrej;
}
/* D: strip-mined buffer, 16 steps per strip */
#define STRIP 16
static void group_normal_strip(uint64_t sub0, double *out, uint64_t w0) {
    LOAD_STATE(sub0);
    uint64_t buf[LANES][STRIP * 4];
    for (int j = 0; j < SUB; j += STRIP * 4) {
        for (int s = 0; s < STRIP; s++) {
            __m256i v[4]; STEP4(v); TRANSPOSE(v, r0, r1, r2, r3);
            _mm256_storeu_si256((__m256i *)(buf[0] + 4 * s), r0);
            _mm256_storeu_si256((__m256i *)(buf[1] + 4 * s), r1);
            _mm256_storeu_si256((__m256i *)(buf[2] + 4 * s), r2);
            _mm256_storeu_si256((__m256i *)(buf[3] + 4 * s), r3);
        }
        for (int l = 0; l < 4; l++)
            transform(buf[l], out + l * SUB + j, STRIP * 4, w0 + l * SUB + j);
    }
}

/* uniform with non-temporal 128-bit stores (R's data is 16-byte aligned) */
static void group_u01_nt(uint64_t sub0, double *out) {
    LOAD_STATE(sub0);
    const __m256i exponent = _mm256_set1_epi64x(0x3ff0000000000000LL);
    const __m256d bias = _mm256_set1_pd(1.0 - 0x1.0p-53);
#define U01(w) _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_srli_epi64((w), 12), exponent)), bias)
#define NT(p, v) do { __m256d v_ = (v); _mm_stream_pd((p), _mm256_castpd256_pd128(v_)); \
                      _mm_stream_pd((p) + 2, _mm256_extractf128_pd(v_, 1)); } while (0)
    for (int j = 0; j < SUB; j += 4) {
        __m256i v[4]; STEP4(v); TRANSPOSE(v, r0, r1, r2, r3);
        NT(out + 0 * SUB + j, U01(r0));
        NT(out + 1 * SUB + j, U01(r1));
        NT(out + 2 * SUB + j, U01(r2));
        NT(out + 3 * SUB + j, U01(r3));
    }
    _mm_sfence();
#undef NT
#undef U01
}

typedef void (*fill_fn)(double *, long);
static void fill_A(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        uint64_t buf[CHUNK];
        group_words(c * 2 * LANES, buf);
        group_words(c * 2 * LANES + LANES, buf + LANES * SUB);
        transform(buf, out + c * CHUNK, CHUNK, (uint64_t)c * CHUNK);
    }
}
static void fill_B(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        for (int g = 0; g < 2; g++) {
            uint16_t rej[LANES * SUB]; uint64_t rejw[LANES * SUB];
            double *o = out + c * CHUNK + g * LANES * SUB;
            int nr = group_normal_fused(c * 2 * LANES + g * LANES, o, rej, rejw);
            for (int k = 0; k < nr; k++)
                o[rej[k]] = reject((uint64_t)(o - out) + rej[k], rejw[k]);
        }
    }
}
static void fill_D(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++)
        for (int g = 0; g < 2; g++) {
            double *o = out + c * CHUNK + g * LANES * SUB;
            group_normal_strip(c * 2 * LANES + g * LANES, o, (uint64_t)(o - out));
        }
}
static void fill_U(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        group_u01(c * 2 * LANES, out + c * CHUNK);
        group_u01(c * 2 * LANES + LANES, out + c * CHUNK + LANES * SUB);
    }
}
static void fill_UNT(double *out, long n) {
    for (long c = 0; c < n / CHUNK; c++) {
        group_u01_nt(c * 2 * LANES, out + c * CHUNK);
        group_u01_nt(c * 2 * LANES + LANES, out + c * CHUNK + LANES * SUB);
    }
}
static void fill_W(double *out, long n) {   /* words only, to the buffer; touches out once per chunk */
    for (long c = 0; c < n / CHUNK; c++) {
        uint64_t buf[CHUNK];
        group_words(c * 2 * LANES, buf);
        group_words(c * 2 * LANES + LANES, buf + LANES * SUB);
        out[c * CHUNK] = (double)buf[CHUNK - 1];
    }
}

int main(int argc, char **argv) {
    long n = argc > 1 ? atol(argv[1]) : 10000000L; n -= n % CHUNK;
    int reps = argc > 2 ? atoi(argv[2]) : 10;
    tables_init();
    double *ref = malloc((size_t)n * 8), *out = malloc((size_t)n * 8);
    fill_A(ref, n);
    fill_B(out, n); printf("B == A: %s\n", memcmp(ref, out, (size_t)n * 8) ? "NO" : "yes");
    fill_D(out, n); printf("D == A: %s\n", memcmp(ref, out, (size_t)n * 8) ? "NO" : "yes");
    fill_U(ref, n); fill_UNT(out, n); printf("UN == U: %s\n", memcmp(ref, out, (size_t)n * 8) ? "NO" : "yes");
    struct { const char *name; fill_fn f; } v[] = {
        {"W  words only (generator)", fill_W}, {"U  uniform, fused (shipped)", fill_U},
        {"UN uniform, fused, non-temporal 128-bit", fill_UNT},
        {"A  normal, buffered (shipped)", fill_A}, {"B  normal, fused, rejects listed", fill_B},
        {"D  normal, 2 KiB strips", fill_D}};
    printf("n = %ld, %d reps; ns per value\n%-36s %10s %12s\n", n, reps, "", "reused", "fresh malloc");
    for (int i = 0; i < 6; i++) {
        double acc = 0;
        v[i].f(out, n);
        double t = now(); for (int r = 0; r < reps; r++) v[i].f(out, n); double e1 = (now() - t) / reps;
        t = now();
        for (int r = 0; r < reps; r++) { double *p = malloc((size_t)n * 8); v[i].f(p, n); acc += p[n - 1]; free(p); }
        double e2 = (now() - t) / reps;
        printf("%-36s %10.3f %12.3f   (%g)\n", v[i].name, e1 / n * 1e9, e2 / n * 1e9, acc);
    }
    free(ref); free(out);
    return 0;
}
