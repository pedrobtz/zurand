/* NEON uniform fill for xoshiro256pp, in zurand's layout: sub-chunk s (512
 * values) is one xoshiro256++ stream, stored contiguously; the fill is a
 * sequence of groups of 4 sub-chunks. Every variant must match the scalar
 * reference bit for bit, and says so.
 *
 *   scalar   today's arm64 path: one lane, words to a buffer, then convert
 *   neon4    2 register sets x 2 lanes = one group of 4 sub-chunks, fused
 *   neon8    4 register sets x 2 lanes = two groups at once, fused
 *
 * Rotations use SHL + SRI (shift right and insert): 2 instructions, not 3.
 * Seeds are a cheap stand-in for Philox, one per 512 values in all
 * variants, so they cost the same everywhere. */
#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SUB 512
#define N (1L << 24)                      /* 16.8M values, a multiple of 8*SUB */
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }
static inline uint64_t mix(uint64_t z) { z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31); }
static void seed(uint64_t s, uint64_t st[4]) { for (int i = 0; i < 4; i++) st[i] = mix(s * 4 + i); }
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline double u01(uint64_t b) { uint64_t s = (b >> 12) | 0x3ff0000000000000ULL;
  double d; memcpy(&d, &s, 8); return d - (1.0 - 0x1p-53); }

__attribute__((noinline)) static void fill_scalar(double *out) {
  uint64_t buf[4 * SUB];
  for (long g = 0; g < N / (4 * SUB); g++) {
    for (int l = 0; l < 4; l++) {
      uint64_t s[4]; seed(g * 4 + l, s);
      for (int j = 0; j < SUB; j++) {
        uint64_t r = rotl(s[0] + s[3], 23) + s[0], t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        buf[l * SUB + j] = r;
      }
    }
    double *o = out + g * 4 * SUB;
    for (int j = 0; j < 4 * SUB; j++) o[j] = u01(buf[j]);
  }
}

#define VROTL(x, k) vsriq_n_u64(vshlq_n_u64((x), (k)), (x), 64 - (k))
#define VSTEP(s0, s1, s2, s3, r) do { \
  r = vaddq_u64(VROTL(vaddq_u64(s0, s3), 23), s0); uint64x2_t t = vshlq_n_u64(s1, 17); \
  s2 = veorq_u64(s2, s0); s3 = veorq_u64(s3, s1); s1 = veorq_u64(s1, s2); s0 = veorq_u64(s0, s3); \
  s2 = veorq_u64(s2, t); s3 = VROTL(s3, 45); } while (0)
static inline float64x2_t vu01(uint64x2_t w) {
  uint64x2_t m = vorrq_u64(vshrq_n_u64(w, 12), vdupq_n_u64(0x3ff0000000000000ULL));
  return vsubq_f64(vreinterpretq_f64_u64(m), vdupq_n_f64(1.0 - 0x1p-53)); }
/* Registers hold lanes (sub-chunk p, sub-chunk p+1). Two steps, then zip:
 * [x0 y0],[x1 y1] -> [x0 x1] for sub-chunk p, [y0 y1] for p+1. */
static inline void load2(uint64_t sub, uint64x2_t *s0, uint64x2_t *s1, uint64x2_t *s2, uint64x2_t *s3) {
  uint64_t a[4], b[4]; seed(sub, a); seed(sub + 1, b);
  *s0 = (uint64x2_t){a[0], b[0]}; *s1 = (uint64x2_t){a[1], b[1]};
  *s2 = (uint64x2_t){a[2], b[2]}; *s3 = (uint64x2_t){a[3], b[3]}; }
#define PAIR(s0, s1, s2, s3, dst0, dst1, j) do { uint64x2_t r0, r1; \
  VSTEP(s0, s1, s2, s3, r0); VSTEP(s0, s1, s2, s3, r1); \
  vst1q_f64((dst0) + (j), vu01(vzip1q_u64(r0, r1))); \
  vst1q_f64((dst1) + (j), vu01(vzip2q_u64(r0, r1))); } while (0)

__attribute__((noinline)) static void fill_neon4(double *out) {
  for (long g = 0; g < N / (4 * SUB); g++) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3;
    load2(g * 4, &a0, &a1, &a2, &a3); load2(g * 4 + 2, &b0, &b1, &b2, &b3);
    double *o = out + g * 4 * SUB;
    for (int j = 0; j < SUB; j += 2) {
      PAIR(a0, a1, a2, a3, o + 0 * SUB, o + 1 * SUB, j);
      PAIR(b0, b1, b2, b3, o + 2 * SUB, o + 3 * SUB, j);
    }
  }
}

__attribute__((noinline)) static void fill_neon8(double *out) {
  for (long g = 0; g < N / (8 * SUB); g++) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    load2(g * 8, &a0, &a1, &a2, &a3); load2(g * 8 + 2, &b0, &b1, &b2, &b3);
    load2(g * 8 + 4, &c0, &c1, &c2, &c3); load2(g * 8 + 6, &d0, &d1, &d2, &d3);
    double *o = out + g * 8 * SUB;
    for (int j = 0; j < SUB; j += 2) {
      PAIR(a0, a1, a2, a3, o + 0 * SUB, o + 1 * SUB, j);
      PAIR(b0, b1, b2, b3, o + 2 * SUB, o + 3 * SUB, j);
      PAIR(c0, c1, c2, c3, o + 4 * SUB, o + 5 * SUB, j);
      PAIR(d0, d1, d2, d3, o + 6 * SUB, o + 7 * SUB, j);
    }
  }
}

int main(void) {
  double *ref = malloc(N * 8), *o = malloc(N * 8);
  void (*f[3])(double *) = {fill_scalar, fill_neon4, fill_neon8};
  const char *nm[3] = {"scalar (today)", "neon4 fused", "neon8 fused"};
  fill_scalar(ref);
  for (int k = 1; k < 3; k++) { memset(o, 0, N * 8); f[k](o);
    printf("%-16s identical to scalar: %s\n", nm[k], memcmp(o, ref, N * 8) ? "NO" : "yes"); }
  for (int mode = 0; mode < 2; mode++) {
    printf(mode ? "-- fresh allocation per call\n" : "-- reused buffer\n");
    double best[3] = {1e9, 1e9, 1e9};
    for (int rep = 0; rep < 15; rep++) for (int k = 0; k < 3; k++) {
      double *d = mode ? malloc(N * 8) : o; double t = now(); f[k](d); t = now() - t;
      if (mode) free(d); if (t < best[k]) best[k] = t; }
    for (int k = 0; k < 3; k++)
      printf("  %-16s %6.0f M/s  %.2fx scalar\n", nm[k], N / best[k] / 1e6, best[0] / best[k]);
  }
  return 0;
}
