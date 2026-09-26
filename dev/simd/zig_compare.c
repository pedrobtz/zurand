/* The ziggurat's fast path in isolation: 64-bit words already in cache,
 * turned into normals. Every variant must give the same bits as N0, the
 * structure zurand ships; the rare slow path is a stand-in of comparable
 * cost, shared by all variants.
 *
 *   N0  one draw per iteration, branch to the slow path inline (shipped)
 *   N1  branch-free: every draw computed, a reject flag recorded, rejects
 *       fixed afterwards; restrict pointers
 *   N2  N1 unrolled by four, independent chains
 *   N3  vector fast path: NEON 2 lanes (lane loads) / AVX2 4 lanes (plain
 *       loads, no gathers), rejects fixed afterwards
 *   N4  randompack's structure: words in the output array, converted in
 *       place, backwards
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "numpyzig/ziggurat_constants.h"
#if defined(__aarch64__)
#  include <arm_neon.h>
#elif defined(__AVX2__)
#  include <immintrin.h>
#endif

#define NW 4096
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }

__attribute__((noinline)) static double slow(uint64_t index, uint64_t r) {
  /* stand-in for the wedge/tail: ~40 dependent operations */
  double x = (double)((r >> 9) & 0xfffffffffffffULL) * wi_double[r & 0xff];
  for (int i = 0; i < 10; i++) x = x * 0.999 + (double)(index & 7) * 1e-9;
  return (r >> 8) & 1 ? -x : x;
}

static inline double fast_or_slow(uint64_t index, uint64_t r) {
  int idx = (int)(r & 0xff);
  uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
  int64_t s = (r >> 8) & 1 ? -(int64_t)rabs : (int64_t)rabs;
  double x = (double)s * wi_double[idx];
  if (__builtin_expect(rabs < ki_double[idx], 1)) return x;
  return slow(index, r);
}

__attribute__((noinline)) static void n0(const uint64_t *buf, double *o, int m) {
  for (int j = 0; j < m; j++) o[j] = fast_or_slow((uint64_t)j, buf[j]);
}

__attribute__((noinline)) static void n1(const uint64_t *restrict buf, double *restrict o, int m) {
  unsigned char rej[NW];
  int any = 0;
  for (int j = 0; j < m; j++) {
    uint64_t r = buf[j]; int idx = (int)(r & 0xff);
    uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
    int64_t s = (r >> 8) & 1 ? -(int64_t)rabs : (int64_t)rabs;
    o[j] = (double)s * wi_double[idx];
    unsigned char b = rabs >= ki_double[idx];
    rej[j] = b; any |= b;
  }
  if (any) for (int j = 0; j < m; j++) if (rej[j]) o[j] = slow((uint64_t)j, buf[j]);
}

__attribute__((noinline)) static void n2(const uint64_t *restrict buf, double *restrict o, int m) {
  unsigned char rej[NW];
  int any = 0;
  for (int j = 0; j < m; j += 4) {
#define ONE(k) do { uint64_t r = buf[j + k]; int idx = (int)(r & 0xff);               \
      uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;                             \
      int64_t s = (r >> 8) & 1 ? -(int64_t)rabs : (int64_t)rabs;                   \
      o[j + k] = (double)s * wi_double[idx];                                        \
      unsigned char b = rabs >= ki_double[idx]; rej[j + k] = b; any |= b; } while (0)
    ONE(0); ONE(1); ONE(2); ONE(3);
#undef ONE
  }
  if (any) for (int j = 0; j < m; j++) if (rej[j]) o[j] = slow((uint64_t)j, buf[j]);
}

#if defined(__aarch64__)
__attribute__((noinline)) static void n3(const uint64_t *restrict buf, double *restrict o, int m) {
  unsigned char rej[NW / 2];
  int any = 0;
  const uint64x2_t mask52 = vdupq_n_u64(0x000fffffffffffffULL), exp52 = vdupq_n_u64(0x4330000000000000ULL);
  const float64x2_t two52 = vdupq_n_f64(0x1p52);
  for (int j = 0; j < m; j += 2) {
    uint64x2_t w = vld1q_u64(buf + j);
    uint64x2_t rabs = vandq_u64(vshrq_n_u64(w, 9), mask52);
    uint64x2_t sign = vshlq_n_u64(vshrq_n_u64(w, 8), 63);
    sign = vbicq_u64(sign, vceqzq_u64(rabs));                    /* +0.0 at rabs == 0 */
    float64x2_t d = vsubq_f64(vreinterpretq_f64_u64(vorrq_u64(rabs, exp52)), two52);
    int i0 = (int)(buf[j] & 0xff), i1 = (int)(buf[j + 1] & 0xff);
    float64x2_t wi = vld1q_lane_f64(&wi_double[i1], vld1q_dup_f64(&wi_double[i0]), 1);
    uint64x2_t ki = vld1q_lane_u64(&ki_double[i1], vld1q_dup_u64(&ki_double[i0]), 1);
    float64x2_t x = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(vmulq_f64(d, wi)), sign));
    vst1q_f64(o + j, x);
    uint64x2_t bad = vcgeq_u64(rabs, ki);
    unsigned char b = (unsigned char)((vgetq_lane_u64(bad, 0) & 1) | (vgetq_lane_u64(bad, 1) & 2));
    rej[j >> 1] = b; any |= b;
  }
  if (any) for (int j = 0; j < m; j++) if (rej[j >> 1] & (1 << (j & 1))) o[j] = slow((uint64_t)j, buf[j]);
}
#elif defined(__AVX2__)
__attribute__((noinline)) static void n3(const uint64_t *restrict buf, double *restrict o, int m) {
  unsigned char rej[NW / 4];
  int any = 0;
  const __m256i mask52 = _mm256_set1_epi64x(0x000fffffffffffffLL), exp52 = _mm256_set1_epi64x(0x4330000000000000LL);
  const __m256d two52 = _mm256_set1_pd(0x1p52);
  for (int j = 0; j < m; j += 4) {
    __m256i w = _mm256_loadu_si256((const __m256i *)(buf + j));
    __m256i rabs = _mm256_and_si256(_mm256_srli_epi64(w, 9), mask52);
    __m256i sign = _mm256_slli_epi64(_mm256_srli_epi64(w, 8), 63);
    sign = _mm256_andnot_si256(_mm256_cmpeq_epi64(rabs, _mm256_setzero_si256()), sign);
    __m256d d = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(rabs, exp52)), two52);
    uint64_t i0 = buf[j] & 0xff, i1 = buf[j + 1] & 0xff, i2 = buf[j + 2] & 0xff, i3 = buf[j + 3] & 0xff;
    __m256d wi = _mm256_set_pd(wi_double[i3], wi_double[i2], wi_double[i1], wi_double[i0]);
    __m256i ki = _mm256_set_epi64x((long long)ki_double[i3], (long long)ki_double[i2], (long long)ki_double[i1], (long long)ki_double[i0]);
    _mm256_storeu_pd(o + j, _mm256_xor_pd(_mm256_mul_pd(d, wi), _mm256_castsi256_pd(sign)));
    unsigned char b = (unsigned char)(~_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(ki, rabs))) & 0xf);
    rej[j >> 2] = b; any |= b;
  }
  if (any) for (int j = 0; j < m; j++) if (rej[j >> 2] & (1 << (j & 3))) o[j] = slow((uint64_t)j, buf[j]);
}
#else
static void n3(const uint64_t *buf, double *o, int m) { n2(buf, o, m); }
#endif

/* N4: randompack's structure, in place over the output array, backwards */
__attribute__((noinline)) static void n4(uint64_t *wo, int m) {
  for (int j = m; j-- > 0;) {
    uint64_t r; memcpy(&r, wo + j, 8);
    double x = fast_or_slow((uint64_t)j, r);
    memcpy(wo + j, &x, 8);
  }
}

int main(void) {
  static uint64_t words[NW], scratch[NW];
  static double o0[NW], o[NW];
  uint64_t z = 0x9e3779b97f4a7c15ULL;
  for (int i = 0; i < NW; i++) { z += 0x9e3779b97f4a7c15ULL; uint64_t t = z;
    t = (t ^ (t >> 30)) * 0xbf58476d1ce4e5b9ULL; t = (t ^ (t >> 27)) * 0x94d049bb133111ebULL; words[i] = t ^ (t >> 31); }
  n0(words, o0, NW);
  int rej = 0; for (int i = 0; i < NW; i++) rej += (words[i] >> 9 & 0xfffffffffffffULL) >= ki_double[words[i] & 0xff];
  printf("rejected in fast path: %d of %d\n", rej, NW);
  void (*fs[3])(const uint64_t *, double *, int) = {n1, n2, n3};
  const char *nm[5] = {"N0 shipped: branch per draw", "N1 branch-free, restrict, fix-up",
                       "N2 N1 unrolled x4", "N3 vector fast path", "N4 randompack: in place, backwards"};
  for (int k = 0; k < 3; k++) { memset(o, 0, sizeof o); fs[k](words, o, NW);
    printf("%-38s identical to N0: %s\n", nm[k + 1], memcmp(o, o0, sizeof o) ? "NO" : "yes"); }
  memcpy(scratch, words, sizeof words); n4(scratch, NW);
  printf("%-38s identical to N0: %s\n", nm[4], memcmp(scratch, o0, sizeof o0) ? "NO" : "yes");
  const long REPS = 20000;
  for (int round = 0; round < 3; round++) {
    printf("-- round %d  (ns per draw, words in cache)\n", round + 1);
    double best[5] = {1e9, 1e9, 1e9, 1e9, 1e9};
    for (int t = 0; t < 5; t++) {
      double s = now(); for (long r = 0; r < REPS; r++) n0(words, o, NW); s = now() - s; if (s < best[0]) best[0] = s;
      for (int k = 0; k < 3; k++) { s = now(); for (long r = 0; r < REPS; r++) fs[k](words, o, NW); s = now() - s; if (s < best[k + 1]) best[k + 1] = s; }
      s = now(); for (long r = 0; r < REPS; r++) { memcpy(scratch, words, sizeof words); n4(scratch, NW); } s = now() - s; if (s < best[4]) best[4] = s;
    }
    for (int k = 0; k < 5; k++) printf("  %-38s %.3f%s\n", nm[k], best[k] / (REPS * (double)NW) * 1e9, k == 4 ? "  (includes copying the words in)" : "");
  }
  return 0;
}
