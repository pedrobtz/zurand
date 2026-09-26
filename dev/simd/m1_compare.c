/* Where does the time go in zurand's and randompack's NEON uniform fills?
 *
 * Both run xoshiro256++ in 8 lanes (four 2-lane registers). They differ in
 * four ways, and each variant below changes one of them:
 *   seeding    zurand reseeds each 512-value sub-chunk from Philox4x64-10;
 *              randompack seeds once and keeps a running state
 *   rotations  zurand: SHL + SRI (2 instructions); randompack: SHL, SHR, ORR
 *   layout     zurand zips lane pairs so each sub-chunk is contiguous;
 *              randompack stores lanes interleaved, as they come
 *   passes     zurand converts in registers and stores once; randompack
 *              fills a 144-word buffer (state loaded and stored per refill),
 *              then converts it in a second pass
 * The randompack variants reimplement the loops of randompack 0.1.10
 * (src/x256pp_neon_scalar.inc, src/rand_dble.inc; BSD 3-clause, Kristjan
 * Jonasson) for measurement only; the conversion keeps its FMA form.
 *
 * Reports ns per value in cache (4096 values, reused) and in memory
 * (16.8M values, reused, so page faults are excluded). */
#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Random123/philox.h"

#define SUB 512
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }
static const philox4x64_key_t KEY = {{0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL}};

/* ---- seeding ---- */
static inline void seed_philox(uint64_t sub, uint64_t st[4]) {
  philox4x64_ctr_t c = {{sub, 0, 1, 1}};
  philox4x64_ctr_t r = philox4x64_R(10, c, KEY);
  st[0] = r.v[0]; st[1] = r.v[1]; st[2] = r.v[2]; st[3] = r.v[3];
  if (!(st[0] | st[1] | st[2] | st[3])) st[0] = 1;
}
static inline void seed_cheap(uint64_t sub, uint64_t st[4]) {
  st[0] = sub ^ 0x9e3779b97f4a7c15ULL; st[1] = sub * 3 + 1;
  st[2] = ~sub; st[3] = sub + 0x2545F4914F6CDD1DULL;
}

/* ---- rotations ---- */
#define ROT_SRI(x, k) vsriq_n_u64(vshlq_n_u64((x), (k)), (x), 64 - (k))
#define ROT_ORR(x, k) vorrq_u64(vshlq_n_u64((x), (k)), vshrq_n_u64((x), 64 - (k)))
#define STEP(ROT, s0, s1, s2, s3, r) do {                                  \
    (r) = vaddq_u64(ROT(vaddq_u64((s0), (s3)), 23), (s0));                 \
    uint64x2_t t_ = vshlq_n_u64((s1), 17);                                 \
    (s2) = veorq_u64((s2), (s0)); (s3) = veorq_u64((s3), (s1));           \
    (s1) = veorq_u64((s1), (s2)); (s0) = veorq_u64((s0), (s3));           \
    (s2) = veorq_u64((s2), t_);   (s3) = ROT((s3), 45); } while (0)

static inline float64x2_t u01_open(uint64x2_t w) {   /* zurand: (m + 0.5) 2^-52 */
  uint64x2_t m = vorrq_u64(vshrq_n_u64(w, 12), vdupq_n_u64(0x3ff0000000000000ULL));
  return vsubq_f64(vreinterpretq_f64_u64(m), vdupq_n_f64(1.0 - 0x1p-53)); }

/* ---- zurand: 8 sub-chunks per call, fused ---- */
#define LOADPAIR(SEED, p, sub) do { uint64_t a_[4], b_[4], l_[2];          \
    SEED((sub), a_); SEED((sub) + 1, b_);                                  \
    l_[0] = a_[0]; l_[1] = b_[0]; p##0 = vld1q_u64(l_);                    \
    l_[0] = a_[1]; l_[1] = b_[1]; p##1 = vld1q_u64(l_);                    \
    l_[0] = a_[2]; l_[1] = b_[2]; p##2 = vld1q_u64(l_);                    \
    l_[0] = a_[3]; l_[1] = b_[3]; p##3 = vld1q_u64(l_); } while (0)
#define ZURAND_FILL(NAME, SEED, ROT)                                        \
static void NAME(double *out, long n) {                                     \
  for (long g = 0; g < n / (8 * SUB); g++) {                                \
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3; \
    uint64_t s0 = (uint64_t)g * 8;                                          \
    LOADPAIR(SEED, a, s0); LOADPAIR(SEED, b, s0 + 2);                       \
    LOADPAIR(SEED, c, s0 + 4); LOADPAIR(SEED, d, s0 + 6);                   \
    double *o = out + g * 8 * SUB;                                          \
    for (int j = 0; j < SUB; j += 2) {                                      \
      uint64x2_t r0, r1;                                                    \
      STEP(ROT, a0, a1, a2, a3, r0); STEP(ROT, a0, a1, a2, a3, r1);         \
      vst1q_f64(o + 0 * SUB + j, u01_open(vzip1q_u64(r0, r1)));             \
      vst1q_f64(o + 1 * SUB + j, u01_open(vzip2q_u64(r0, r1)));             \
      STEP(ROT, b0, b1, b2, b3, r0); STEP(ROT, b0, b1, b2, b3, r1);         \
      vst1q_f64(o + 2 * SUB + j, u01_open(vzip1q_u64(r0, r1)));             \
      vst1q_f64(o + 3 * SUB + j, u01_open(vzip2q_u64(r0, r1)));             \
      STEP(ROT, c0, c1, c2, c3, r0); STEP(ROT, c0, c1, c2, c3, r1);         \
      vst1q_f64(o + 4 * SUB + j, u01_open(vzip1q_u64(r0, r1)));             \
      vst1q_f64(o + 5 * SUB + j, u01_open(vzip2q_u64(r0, r1)));             \
      STEP(ROT, d0, d1, d2, d3, r0); STEP(ROT, d0, d1, d2, d3, r1);         \
      vst1q_f64(o + 6 * SUB + j, u01_open(vzip1q_u64(r0, r1)));             \
      vst1q_f64(o + 7 * SUB + j, u01_open(vzip2q_u64(r0, r1)));             \
    }                                                                       \
  }                                                                         \
}
ZURAND_FILL(z_real,        seed_philox, ROT_SRI)   /* zurand as shipped */
ZURAND_FILL(z_cheapseed,   seed_cheap,  ROT_SRI)   /* minus Philox */
ZURAND_FILL(z_orr,         seed_philox, ROT_ORR)   /* with randompack's rotations */

/* ---- tricks, all producing zurand's exact stream ---- */
/* T1: conversion in one SRI: insert w >> 12 under the exponent of 1.0 */
static inline float64x2_t u01_sri(uint64x2_t w) {
  uint64x2_t m = vsriq_n_u64(vdupq_n_u64(0x3ff0000000000000ULL), w, 12);
  return vsubq_f64(vreinterpretq_f64_u64(m), vdupq_n_f64(1.0 - 0x1p-53)); }
#define ZURAND_FILL_T1(NAME, SEED)                                          \
static void NAME(double *out, long n) {                                     \
  for (long g = 0; g < n / (8 * SUB); g++) {                                \
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3; \
    uint64_t s0 = (uint64_t)g * 8;                                          \
    LOADPAIR(SEED, a, s0); LOADPAIR(SEED, b, s0 + 2);                       \
    LOADPAIR(SEED, c, s0 + 4); LOADPAIR(SEED, d, s0 + 6);                   \
    double *o = out + g * 8 * SUB;                                          \
    for (int j = 0; j < SUB; j += 2) {                                      \
      uint64x2_t r0, r1;                                                    \
      STEP(ROT_SRI, a0, a1, a2, a3, r0); STEP(ROT_SRI, a0, a1, a2, a3, r1); \
      vst1q_f64(o + 0 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));              \
      vst1q_f64(o + 1 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));              \
      STEP(ROT_SRI, b0, b1, b2, b3, r0); STEP(ROT_SRI, b0, b1, b2, b3, r1); \
      vst1q_f64(o + 2 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));              \
      vst1q_f64(o + 3 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));              \
      STEP(ROT_SRI, c0, c1, c2, c3, r0); STEP(ROT_SRI, c0, c1, c2, c3, r1); \
      vst1q_f64(o + 4 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));              \
      vst1q_f64(o + 5 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));              \
      STEP(ROT_SRI, d0, d1, d2, d3, r0); STEP(ROT_SRI, d0, d1, d2, d3, r1); \
      vst1q_f64(o + 6 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));              \
      vst1q_f64(o + 7 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));              \
    }                                                                       \
  }                                                                         \
}
ZURAND_FILL_T1(z_t1, seed_philox)

/* T2: T1 plus store pairs -- four steps per register, then one 32-byte
 * STP (vst1q_f64_x2) per sub-chunk instead of two 16-byte stores. */
#define QUAD(p, l) do { uint64x2_t r0, r1, r2, r3;                          \
    STEP(ROT_SRI, p##0, p##1, p##2, p##3, r0); STEP(ROT_SRI, p##0, p##1, p##2, p##3, r1); \
    STEP(ROT_SRI, p##0, p##1, p##2, p##3, r2); STEP(ROT_SRI, p##0, p##1, p##2, p##3, r3); \
    float64x2x2_t x_, y_;                                                   \
    x_.val[0] = u01_sri(vzip1q_u64(r0, r1)); x_.val[1] = u01_sri(vzip1q_u64(r2, r3)); \
    y_.val[0] = u01_sri(vzip2q_u64(r0, r1)); y_.val[1] = u01_sri(vzip2q_u64(r2, r3)); \
    vst1q_f64_x2(o + (l) * SUB + j, x_); vst1q_f64_x2(o + ((l) + 1) * SUB + j, y_); } while (0)
static void z_t2(double *out, long n) {
  for (long g = 0; g < n / (8 * SUB); g++) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    uint64_t s0 = (uint64_t)g * 8;
    LOADPAIR(seed_philox, a, s0); LOADPAIR(seed_philox, b, s0 + 2);
    LOADPAIR(seed_philox, c, s0 + 4); LOADPAIR(seed_philox, d, s0 + 6);
    double *o = out + g * 8 * SUB;
    for (int j = 0; j < SUB; j += 4) { QUAD(a, 0); QUAD(b, 2); QUAD(c, 4); QUAD(d, 6); }
  }
}

/* T3: T1 plus two scalar sub-chunks in the same loop, so the integer ALUs
 * work alongside the NEON pipes. Groups of 10 sub-chunks (the grouping is
 * not part of the stream). */
static inline uint64_t rotl_(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline double u01s(uint64_t r) { uint64_t b = (r >> 12) | 0x3ff0000000000000ULL;
  double d; memcpy(&d, &b, 8); return d - (1.0 - 0x1p-53); }
#define SSTEP(s, r) do { r = rotl_(s[0] + s[3], 23) + s[0]; uint64_t t_ = s[1] << 17; \
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t_; s[3] = rotl_(s[3], 45); } while (0)
static void z_t3(double *out, long n) {
  long groups = n / (10 * SUB);
  for (long g = 0; g < groups; g++) {
    uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    uint64_t s0 = (uint64_t)g * 10;
    LOADPAIR(seed_philox, a, s0); LOADPAIR(seed_philox, b, s0 + 2);
    LOADPAIR(seed_philox, c, s0 + 4); LOADPAIR(seed_philox, d, s0 + 6);
    uint64_t e[4], f[4]; seed_philox(s0 + 8, e); seed_philox(s0 + 9, f);
    double *o = out + g * 10 * SUB;
    for (int j = 0; j < SUB; j += 2) {
      uint64x2_t r0, r1; uint64_t e0, e1, f0, f1;
      SSTEP(e, e0); SSTEP(f, f0);
      STEP(ROT_SRI, a0, a1, a2, a3, r0); STEP(ROT_SRI, a0, a1, a2, a3, r1);
      vst1q_f64(o + 0 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));
      vst1q_f64(o + 1 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));
      STEP(ROT_SRI, b0, b1, b2, b3, r0); STEP(ROT_SRI, b0, b1, b2, b3, r1);
      vst1q_f64(o + 2 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));
      vst1q_f64(o + 3 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));
      SSTEP(e, e1); SSTEP(f, f1);
      STEP(ROT_SRI, c0, c1, c2, c3, r0); STEP(ROT_SRI, c0, c1, c2, c3, r1);
      vst1q_f64(o + 4 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));
      vst1q_f64(o + 5 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));
      STEP(ROT_SRI, d0, d1, d2, d3, r0); STEP(ROT_SRI, d0, d1, d2, d3, r1);
      vst1q_f64(o + 6 * SUB + j, u01_sri(vzip1q_u64(r0, r1)));
      vst1q_f64(o + 7 * SUB + j, u01_sri(vzip2q_u64(r0, r1)));
      o[8 * SUB + j] = u01s(e0); o[8 * SUB + j + 1] = u01s(e1);
      o[9 * SUB + j] = u01s(f0); o[9 * SUB + j + 1] = u01s(f1);
    }
  }
}

/* ---- round 2: combine what worked ---- */
/* V1: ORR rotations (measured faster than SRI on M1: SRI reads its
 * destination, lengthening the loop-carried chain). */
static inline float64x2_t u01o(uint64x2_t w) { return u01_open(w); }
#define PAIRV(ROT, p, dst, l, j) do { uint64x2_t r0, r1;                     \
    STEP(ROT, p##0, p##1, p##2, p##3, r0); STEP(ROT, p##0, p##1, p##2, p##3, r1); \
    vst1q_f64((dst) + (l) * STRIDE + (j), u01o(vzip1q_u64(r0, r1)));          \
    vst1q_f64((dst) + ((l) + 1) * STRIDE + (j), u01o(vzip2q_u64(r0, r1))); } while (0)
#define SEEDS8(s0) uint64x2_t a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3; \
    LOADPAIR(seed_philox, a, (s0)); LOADPAIR(seed_philox, b, (s0) + 2);    \
    LOADPAIR(seed_philox, c, (s0) + 4); LOADPAIR(seed_philox, d, (s0) + 6)

/* V3: ORR, staged: 8 sub-chunks into a cache-resident buffer, then copied
 * out sequentially (one write stream instead of eight 4 KiB apart).
 * STRIDE = SUB is the naive layout; SUB + 8 pads away the 4 KiB aliasing. */
#define STAGED(NAME, STR)                                                   \
static void NAME(double *out, long n) {                                     \
  enum { STRIDE = STR };                                                    \
  static double stage[8 * STR] __attribute__((aligned(64)));                \
  for (long g = 0; g < n / (8 * SUB); g++) {                                \
    SEEDS8((uint64_t)g * 8);                                                \
    for (int j = 0; j < SUB; j += 2) {                                      \
      PAIRV(ROT_ORR, a, stage, 0, j); PAIRV(ROT_ORR, b, stage, 2, j);       \
      PAIRV(ROT_ORR, c, stage, 4, j); PAIRV(ROT_ORR, d, stage, 6, j);       \
    }                                                                       \
    double *o = out + g * 8 * SUB;                                          \
    for (int l = 0; l < 8; l++) memcpy(o + l * SUB, stage + l * STRIDE, SUB * sizeof(double)); \
  }                                                                         \
}
STAGED(z_v3_naive, SUB)
STAGED(z_v3_pad, SUB + 8)

/* V2: ORR + two scalar sub-chunks alongside (groups of 10), direct stores */
static void z_v2(double *out, long n) {
  enum { STRIDE = SUB };
  for (long g = 0; g < n / (10 * SUB); g++) {
    uint64_t s0 = (uint64_t)g * 10;
    SEEDS8(s0);
    uint64_t e[4], f[4]; seed_philox(s0 + 8, e); seed_philox(s0 + 9, f);
    double *o = out + g * 10 * SUB;
    for (int j = 0; j < SUB; j += 2) {
      uint64_t e0, e1, f0, f1;
      SSTEP(e, e0); SSTEP(f, f0);
      PAIRV(ROT_ORR, a, o, 0, j); PAIRV(ROT_ORR, b, o, 2, j);
      SSTEP(e, e1); SSTEP(f, f1);
      PAIRV(ROT_ORR, c, o, 4, j); PAIRV(ROT_ORR, d, o, 6, j);
      o[8 * SUB + j] = u01s(e0); o[8 * SUB + j + 1] = u01s(e1);
      o[9 * SUB + j] = u01s(f0); o[9 * SUB + j + 1] = u01s(f1);
    }
  }
}

/* V4: ORR + scalar sub-chunks + padded staging */
static void z_v4(double *out, long n) {
  enum { STRIDE = SUB + 8 };
  static double stage[10 * (SUB + 8)] __attribute__((aligned(64)));
  for (long g = 0; g < n / (10 * SUB); g++) {
    uint64_t s0 = (uint64_t)g * 10;
    SEEDS8(s0);
    uint64_t e[4], f[4]; seed_philox(s0 + 8, e); seed_philox(s0 + 9, f);
    for (int j = 0; j < SUB; j += 2) {
      uint64_t e0, e1, f0, f1;
      SSTEP(e, e0); SSTEP(f, f0);
      PAIRV(ROT_ORR, a, stage, 0, j); PAIRV(ROT_ORR, b, stage, 2, j);
      SSTEP(e, e1); SSTEP(f, f1);
      PAIRV(ROT_ORR, c, stage, 4, j); PAIRV(ROT_ORR, d, stage, 6, j);
      stage[8 * STRIDE + j] = u01s(e0); stage[8 * STRIDE + j + 1] = u01s(e1);
      stage[9 * STRIDE + j] = u01s(f0); stage[9 * STRIDE + j + 1] = u01s(f1);
    }
    double *o = out + g * 10 * SUB;
    for (int l = 0; l < 10; l++) memcpy(o + l * SUB, stage + l * STRIDE, SUB * sizeof(double));
  }
}

/* zurand stream reference: scalar, sub-chunk by sub-chunk */
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static void z_scalar(double *out, long n) {
  for (long s = 0; s < n / SUB; s++) {
    uint64_t st[4]; seed_philox((uint64_t)s, st);
    for (int j = 0; j < SUB; j++) {
      uint64_t r = rotl(st[0] + st[3], 23) + st[0], t = st[1] << 17;
      st[2] ^= st[0]; st[3] ^= st[1]; st[1] ^= st[2]; st[0] ^= st[3]; st[2] ^= t;
      st[3] = rotl(st[3], 45);
      uint64_t b = (r >> 12) | 0x3ff0000000000000ULL; double d; memcpy(&d, &b, 8);
      out[s * SUB + j] = d - (1.0 - 0x1p-53);
    }
  }
}

/* ---- randompack style: 8 lanes interleaved, state in memory ---- */
#define BUFSIZE 144
typedef struct { uint64_t s0[8], s1[8], s2[8], s3[8]; } xo256;
static void rp_init(xo256 *st) { for (int l = 0; l < 8; l++) { uint64_t s[4];
  seed_philox(1000 + l, s); st->s0[l] = s[0]; st->s1[l] = s[1]; st->s2[l] = s[2]; st->s3[l] = s[3]; } }
#define RP_FILL(NAME, ROT)                                                  \
static void NAME(uint64_t *buf, size_t len, xo256 *st) {                    \
  uint64x2_t a0 = vld1q_u64(&st->s0[0]), a1 = vld1q_u64(&st->s1[0]),       \
             a2 = vld1q_u64(&st->s2[0]), a3 = vld1q_u64(&st->s3[0]);       \
  uint64x2_t b0 = vld1q_u64(&st->s0[2]), b1 = vld1q_u64(&st->s1[2]),       \
             b2 = vld1q_u64(&st->s2[2]), b3 = vld1q_u64(&st->s3[2]);       \
  uint64x2_t c0 = vld1q_u64(&st->s0[4]), c1 = vld1q_u64(&st->s1[4]),       \
             c2 = vld1q_u64(&st->s2[4]), c3 = vld1q_u64(&st->s3[4]);       \
  uint64x2_t d0 = vld1q_u64(&st->s0[6]), d1 = vld1q_u64(&st->s1[6]),       \
             d2 = vld1q_u64(&st->s2[6]), d3 = vld1q_u64(&st->s3[6]);       \
  for (size_t i = 0; i < len; i += 8) {                                     \
    uint64x2_t ra, rb, rc, rd;                                              \
    STEP(ROT, a0, a1, a2, a3, ra); STEP(ROT, b0, b1, b2, b3, rb);           \
    STEP(ROT, c0, c1, c2, c3, rc); STEP(ROT, d0, d1, d2, d3, rd);           \
    vst1q_u64(buf + i, ra); vst1q_u64(buf + i + 2, rb);                     \
    vst1q_u64(buf + i + 4, rc); vst1q_u64(buf + i + 6, rd);                 \
  }                                                                         \
  vst1q_u64(&st->s0[0], a0); vst1q_u64(&st->s1[0], a1); vst1q_u64(&st->s2[0], a2); vst1q_u64(&st->s3[0], a3); \
  vst1q_u64(&st->s0[2], b0); vst1q_u64(&st->s1[2], b1); vst1q_u64(&st->s2[2], b2); vst1q_u64(&st->s3[2], b3); \
  vst1q_u64(&st->s0[4], c0); vst1q_u64(&st->s1[4], c1); vst1q_u64(&st->s2[4], c2); vst1q_u64(&st->s3[4], c3); \
  vst1q_u64(&st->s0[6], d0); vst1q_u64(&st->s1[6], d1); vst1q_u64(&st->s2[6], d2); vst1q_u64(&st->s3[6], d3); \
}
RP_FILL(rp_gen_orr, ROT_ORR)
RP_FILL(rp_gen_sri, ROT_SRI)
/* rand_unif_4x52 for [0, 1): shift = a - b = -1, scale b = 1, as FMA */
static void rp_convert(double *x, size_t len, const uint64_t *u) {
  uint64x2_t expo = vdupq_n_u64(0x3ff0000000000000ULL);
  float64x2_t s = vdupq_n_f64(1.0), c = vdupq_n_f64(-1.0);
  for (size_t i = 0; i < len; i += 4) {
    uint64x2_t r0 = vorrq_u64(vshrq_n_u64(vld1q_u64(u + i), 12), expo);
    uint64x2_t r1 = vorrq_u64(vshrq_n_u64(vld1q_u64(u + i + 2), 12), expo);
    vst1q_f64(x + i, vfmaq_f64(c, vreinterpretq_f64_u64(r0), s));
    vst1q_f64(x + i + 2, vfmaq_f64(c, vreinterpretq_f64_u64(r1), s));
  }
}
#define RP_FULL(NAME, GEN)                                                  \
static void NAME(double *out, long n) {                                     \
  static xo256 st; static int init = 0; if (!init) { rp_init(&st); init = 1; } \
  uint64_t buf[BUFSIZE];                                                    \
  long i = 0;                                                               \
  while (i + BUFSIZE <= n) { GEN(buf, BUFSIZE, &st); rp_convert(out + i, BUFSIZE, buf); i += BUFSIZE; } \
  if (i < n) { GEN(buf, BUFSIZE, &st); rp_convert(out + i, (size_t)(n - i) & ~3UL, buf); } \
}
RP_FULL(rp_real, rp_gen_orr)       /* randompack as shipped (NEON) */
RP_FULL(rp_sri,  rp_gen_sri)       /* with SRI rotations */
/* randompack's generator with the conversion fused: no buffer pass */
static void rp_fused(double *out, long n) {
  static xo256 st; static int init = 0; if (!init) { rp_init(&st); init = 1; }
  uint64x2_t a0 = vld1q_u64(&st.s0[0]), a1 = vld1q_u64(&st.s1[0]), a2 = vld1q_u64(&st.s2[0]), a3 = vld1q_u64(&st.s3[0]);
  uint64x2_t b0 = vld1q_u64(&st.s0[2]), b1 = vld1q_u64(&st.s1[2]), b2 = vld1q_u64(&st.s2[2]), b3 = vld1q_u64(&st.s3[2]);
  uint64x2_t c0 = vld1q_u64(&st.s0[4]), c1 = vld1q_u64(&st.s1[4]), c2 = vld1q_u64(&st.s2[4]), c3 = vld1q_u64(&st.s3[4]);
  uint64x2_t d0 = vld1q_u64(&st.s0[6]), d1 = vld1q_u64(&st.s1[6]), d2 = vld1q_u64(&st.s2[6]), d3 = vld1q_u64(&st.s3[6]);
  for (long i = 0; i + 8 <= n; i += 8) {
    uint64x2_t ra, rb, rc, rd;
    STEP(ROT_ORR, a0, a1, a2, a3, ra); STEP(ROT_ORR, b0, b1, b2, b3, rb);
    STEP(ROT_ORR, c0, c1, c2, c3, rc); STEP(ROT_ORR, d0, d1, d2, d3, rd);
    vst1q_f64(out + i, u01_open(ra)); vst1q_f64(out + i + 2, u01_open(rb));
    vst1q_f64(out + i + 4, u01_open(rc)); vst1q_f64(out + i + 6, u01_open(rd));
  }
  vst1q_u64(&st.s0[0], a0);   /* keep the state live */
}

typedef void (*fill_fn)(double *, long);
static double best_ns(fill_fn f, double *buf, long n, long total) {
  long reps = total / n; double best = 1e30;
  for (int trial = 0; trial < 7; trial++) {
    double t = now(); for (long r = 0; r < reps; r++) f(buf, n); t = now() - t;
    if (t < best) best = t; }
  return best / (double)(reps * n) * 1e9;
}

int main(void) {
  const long SMALL = 40 * SUB, LARGE = 40L * SUB * 820;   /* 20480 values, 16.8M */
  double *big = malloc(LARGE * 8), *ref = malloc(LARGE * 8);
  const long CHK = 40 * SUB;
  z_scalar(ref, CHK);
  fill_fn chk[6] = {z_real, z_orr, z_v2, z_v3_naive, z_v3_pad, z_v4}; const char *cn[6] = {"shipped", "V1", "V2", "V3n", "V3p", "V4"};
  for (int k = 0; k < 6; k++) { memset(big, 0, CHK * 8); chk[k](big, CHK);
    printf("zurand %-8s == scalar stream: %s\n", cn[k], memcmp(big, ref, CHK * 8) ? "NO" : "yes"); }
  struct { const char *name; fill_fn f; } v[] = {
    {"zurand as shipped (Philox seeds, SRI, zip, fused)", z_real},
    {"  V1 ORR rotations",                              z_orr},
    {"  V2 ORR + 2 scalar sub-chunks",                  z_v2},
    {"  V3 ORR + staged copy, naive layout",            z_v3_naive},
    {"  V3 ORR + staged copy, padded layout",           z_v3_pad},
    {"  V4 ORR + scalar + padded staging",              z_v4},
    {"randompack as shipped (ORR, buffer + FMA pass)",    rp_real},
    {"  fused conversion, no buffer",                    rp_fused},
  };
  int nv = sizeof v / sizeof v[0];
  for (int round = 0; round < 3; round++) {
    printf("-- round %d        ns/value in cache (20480)  ns/value in memory (16.8M)\n", round + 1);
    for (int k = 0; k < nv; k++) {
      double s = best_ns(v[k].f, big, SMALL, 1L << 26);
      double l = best_ns(v[k].f, big, LARGE, LARGE * 4);
      printf("  %-52s %6.3f            %6.3f\n", v[k].name, s, l);
    }
  }
  return 0;
}
