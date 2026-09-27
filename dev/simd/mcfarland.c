/* McFarland's modified ziggurat (fast_prng, MIT) against the shipped
 * NumPy-table ziggurat, fast paths only: 64-bit words already in cache,
 * turned into normals. The rare edge is a stand-in of the same cost for
 * both, called at each method's own rate (NumPy tables ~1.49%, McFarland
 * 3/256 = 1.17%).
 *
 *   Z6  shipped: signed wi table, rabs - 1 < ki - 1, unrolled by 4
 *   M0  McFarland: i = w & 0xff; i < 253 -> X[i] * (double)(int64)w
 *   M1  M0 unrolled by 4
 *   M2  branch-free: every draw computed, rejects (i >= 253) listed and
 *       fixed afterwards -- possible because the accept test needs no load
 *   M3  M2 vectorised where int64 -> double exists in hardware: NEON
 *       (scvtf, 2 lanes) or AVX-512DQ (vcvtqq2pd, 8 lanes); lane loads for
 *       X[i], no gathers
 *
 *   cc -O2 -I../src -o /tmp/mcf ../dev/simd/mcfarland.c && /tmp/mcf
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

/* McFarland's layer lengths, scaled by 2^-63 (normal.h, 12 digits; the
 * timing does not depend on their last digits) */
static const double MX[254] = { 3.94216628254e-19, 3.72049450041e-19, 3.58270244806e-19, 3.48074762365e-19, 3.39901771719e-19, 3.33037783603e-19, 3.27094388176e-19, 3.21835771325e-19, 3.17107585418e-19, 3.1280307407e-19, 3.08845206558e-19, 3.05176506241e-19, 3.01752902926e-19, 2.98539834407e-19, 2.95509674628e-19, 2.92639979885e-19, 2.899122587e-19, 2.87311087802e-19, 2.84823463271e-19, 2.82438315352e-19, 2.80146139647e-19, 2.77938712618e-19, 2.75808869214e-19, 2.73750326983e-19, 2.71757545434e-19, 2.69825612475e-19, 2.67950151888e-19, 2.66127247304e-19, 2.6435337928e-19, 2.6262537282e-19, 2.60940353352e-19, 2.59295709543e-19, 2.57689061732e-19, 2.56118234977e-19, 2.54581235934e-19, 2.53076232924e-19, 2.51601538678e-19, 2.50155595336e-19, 2.48736961354e-19, 2.47344300031e-19, 2.45976369429e-19, 2.44632013479e-19, 2.43310154111e-19, 2.42009784271e-19, 2.40729961704e-19, 2.39469803409e-19, 2.38228480673e-19, 2.37005214619e-19, 2.35799272207e-19, 2.34609962621e-19, 2.33436634011e-19, 2.32278670547e-19, 2.31135489743e-19, 2.30006540027e-19, 2.28891298528e-19, 2.27789269059e-19, 2.26699980275e-19, 2.25622983985e-19, 2.24557853607e-19, 2.23504182749e-19, 2.22461583905e-19, 2.21429687253e-19, 2.20408139549e-19, 2.19396603103e-19, 2.18394754837e-19, 2.17402285409e-19, 2.164188984e-19, 2.15444309566e-19, 2.14478246135e-19, 2.13520446164e-19, 2.12570657924e-19, 2.11628639347e-19, 2.10694157491e-19, 2.09766988055e-19, 2.08846914916e-19, 2.079337297e-19, 2.0702723138e-19, 2.06127225897e-19, 2.05233525809e-19, 2.04345949953e-19, 2.03464323137e-19, 2.02588475842e-19, 2.01718243948e-19, 2.00853468469e-19, 1.99993995309e-19, 1.9913967503e-19, 1.9829036263e-19, 1.97445917335e-19, 1.96606202405e-19, 1.95771084943e-19, 1.94940435722e-19, 1.9411412902e-19, 1.93292042452e-19, 1.92474056827e-19, 1.91660056003e-19, 1.90849926746e-19, 1.90043558606e-19, 1.89240843788e-19, 1.88441677035e-19, 1.87645955517e-19, 1.86853578721e-19, 1.8606444835e-19, 1.85278468221e-19, 1.84495544175e-19, 1.83715583984e-19, 1.82938497262e-19, 1.82164195388e-19, 1.81392591419e-19, 1.80623600019e-19, 1.7985713738e-19, 1.79093121154e-19, 1.78331470384e-19, 1.77572105435e-19, 1.76814947933e-19, 1.76059920701e-19, 1.753069477e-19, 1.74555953971e-19, 1.73806865576e-19, 1.73059609547e-19, 1.72314113829e-19, 1.71570307233e-19, 1.70828119379e-19, 1.7008748065e-19, 1.69348322146e-19, 1.68610575631e-19, 1.67874173493e-19, 1.67139048692e-19, 1.66405134721e-19, 1.6567236556e-19, 1.64940675631e-19, 1.64209999755e-19, 1.63480273116e-19, 1.62751431209e-19, 1.62023409806e-19, 1.61296144913e-19, 1.60569572726e-19, 1.59843629593e-19, 1.59118251972e-19, 1.58393376391e-19, 1.57668939404e-19, 1.56944877552e-19, 1.56221127324e-19, 1.55497625108e-19, 1.54774307158e-19, 1.54051109542e-19, 1.53327968107e-19, 1.52604818431e-19, 1.51881595777e-19, 1.51158235054e-19, 1.50434670764e-19, 1.49710836959e-19, 1.48986667191e-19, 1.48262094465e-19, 1.47537051186e-19, 1.46811469107e-19, 1.46085279278e-19, 1.4535841199e-19, 1.44630796717e-19, 1.43902362058e-19, 1.43173035676e-19, 1.42442744238e-19, 1.41711413344e-19, 1.40978967466e-19, 1.40245329873e-19, 1.39510422558e-19, 1.38774166165e-19, 1.38036479905e-19, 1.37297281475e-19, 1.36556486972e-19, 1.35814010798e-19, 1.35069765568e-19, 1.34323662007e-19, 1.33575608847e-19, 1.32825512715e-19, 1.32073278015e-19, 1.31318806805e-19, 1.30561998669e-19, 1.29802750579e-19, 1.29040956749e-19, 1.28276508483e-19, 1.2750929401e-19, 1.26739198313e-19, 1.25966102948e-19, 1.25189885844e-19, 1.24410421101e-19, 1.23627578765e-19, 1.22841224598e-19, 1.2205121982e-19, 1.21257420848e-19, 1.20459679002e-19, 1.19657840201e-19, 1.18851744634e-19, 1.18041226403e-19, 1.17226113142e-19, 1.16406225609e-19, 1.15581377245e-19, 1.14751373693e-19, 1.13916012285e-19, 1.13075081485e-19, 1.12228360281e-19, 1.11375617531e-19, 1.10516611251e-19, 1.09651087832e-19, 1.08778781199e-19, 1.07899411881e-19, 1.07012685997e-19, 1.06118294148e-19, 1.05215910191e-19, 1.043051899e-19, 1.0338576948e-19, 1.02457263929e-19, 1.01519265222e-19, 1.00571340295e-19, 9.96130287997e-20, 9.86438405995e-20, 9.76632529648e-20, 9.66707074276e-20, 9.56656062409e-20, 9.46473083804e-20, 9.36151250173e-20, 9.25683143709e-20, 9.15060758376e-20, 9.04275432677e-20, 8.93317772338e-20, 8.82177561023e-20, 8.70843656749e-20, 8.59303871096e-20, 8.47544827642e-20, 8.35551795085e-20, 8.23308489336e-20, 8.10796837291e-20, 7.97996692841e-20, 7.84885492861e-20, 7.71437837009e-20, 7.57624969795e-20, 7.43414135785e-20, 7.28767768074e-20, 7.13642454435e-20, 6.97987602408e-20, 6.81743689448e-20, 6.64839929862e-20, 6.47191103452e-20, 6.28693148131e-20, 6.09216875483e-20, 5.88598735756e-20, 5.66626751161e-20, 5.43018136309e-20, 5.17381717445e-20, 4.89150317224e-20, 4.57447418908e-20, 4.20788025686e-20, 3.76259867224e-20, 3.16285898059e-20, 0.0 };
#define M_IMAX 253

/* ---- shared stand-in for the rare edge ---- */
__attribute__((noinline)) static double slow(uint64_t index, uint64_t r) {
    double x = (double)((r >> 9) & 0xfffffffffffffULL) * wi_double[r & 0xff];
    for (int i = 0; i < 10; i++) x = x * 0.999 + (double)(index & 7) * 1e-9;
    return (r >> 8) & 1 ? -x : x;
}

/* ---- Z6: shipped ---- */
static double wis[512];
static uint64_t kim1[256];
static void tables_init(void) {
    for (int i = 0; i < 256; i++) {
        wis[i] = wi_double[i]; wis[256 + i] = -wi_double[i];
        kim1[i] = ki_double[i] ? ki_double[i] - 1 : 0;
    }
}
__attribute__((noinline)) static double reject_z(uint64_t index, uint64_t r) {
    if (((r >> 9) & 0xfffffffffffffULL) == 0 && ki_double[r & 0xff] != 0) return 0.0;
    return slow(index, r);
}
static inline double fast_z(uint64_t index, uint64_t r) {
    uint64_t rabs = (r >> 9) & 0x000fffffffffffffULL;
    double x = (double)rabs * wis[r & 0x1ff];
    if (__builtin_expect(rabs - 1 < kim1[r & 0xff], 1)) return x;
    return reject_z(index, r);
}
__attribute__((noinline)) static void z6(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j] = fast_z(j, buf[j]); o[j + 1] = fast_z(j + 1, buf[j + 1]);
        o[j + 2] = fast_z(j + 2, buf[j + 2]); o[j + 3] = fast_z(j + 3, buf[j + 3]);
    }
    for (; j < m; j++) o[j] = fast_z(j, buf[j]);
}

/* ---- McFarland ---- */
static inline double fast_m(uint64_t index, uint64_t r) {
    unsigned i = (unsigned)(r & 0xff);
    if (__builtin_expect(i < M_IMAX, 1)) return MX[i] * (double)(int64_t)r;
    return slow(index, r);
}
__attribute__((noinline)) static void m0(const uint64_t *buf, double *o, int m) {
    for (int j = 0; j < m; j++) o[j] = fast_m(j, buf[j]);
}
__attribute__((noinline)) static void m1(const uint64_t *buf, double *o, int m) {
    int j = 0;
    for (; j + 4 <= m; j += 4) {
        o[j] = fast_m(j, buf[j]); o[j + 1] = fast_m(j + 1, buf[j + 1]);
        o[j + 2] = fast_m(j + 2, buf[j + 2]); o[j + 3] = fast_m(j + 3, buf[j + 3]);
    }
    for (; j < m; j++) o[j] = fast_m(j, buf[j]);
}
static void fix(const uint64_t *buf, double *o, const uint16_t *rej, int nr) {
    for (int k = 0; k < nr; k++) o[rej[k]] = slow(rej[k], buf[rej[k]]);
}
__attribute__((noinline)) static void m2(const uint64_t *restrict buf, double *restrict o, int m) {
    uint16_t rej[NW]; int nr = 0;
    for (int j = 0; j < m; j++) {
        uint64_t r = buf[j]; unsigned i = (unsigned)(r & 0xff);
        o[j] = MX[i < M_IMAX ? i : 0] * (double)(int64_t)r;
        rej[nr] = (uint16_t)j; nr += i >= M_IMAX;
    }
    fix(buf, o, rej, nr);
}
#if defined(__aarch64__)
#  define M3_NAME "M3 NEON 2 lanes, scvtf"
__attribute__((noinline)) static void m3(const uint64_t *restrict buf, double *restrict o, int m) {
    uint16_t rej[NW]; int nr = 0;
    int j = 0;
    for (; j + 2 <= m; j += 2) {
        int64x2_t w = vld1q_s64((const int64_t *)buf + j);
        unsigned i0 = (unsigned)(buf[j] & 0xff), i1 = (unsigned)(buf[j + 1] & 0xff);
        float64x2_t x = vsetq_lane_f64(MX[i1], vdupq_n_f64(MX[i0]), 1);
        vst1q_f64(o + j, vmulq_f64(x, vcvtq_f64_s64(w)));
        rej[nr] = (uint16_t)j; nr += i0 >= M_IMAX;
        rej[nr] = (uint16_t)(j + 1); nr += i1 >= M_IMAX;
    }
    for (; j < m; j++) { unsigned i = buf[j] & 0xff; o[j] = MX[i] * (double)(int64_t)buf[j]; rej[nr] = j; nr += i >= M_IMAX; }
    fix(buf, o, rej, nr);
}
static int m3_ok(void) { return 1; }
#elif defined(__x86_64__)
#  define M3_NAME "M3 AVX-512DQ 8 lanes, vcvtqq2pd"
__attribute__((target("avx512f,avx512dq"), noinline))
static void m3(const uint64_t *restrict buf, double *restrict o, int m) {
    uint16_t rej[NW]; int nr = 0;
    int j = 0;
    for (; j + 8 <= m; j += 8) {
        __m512i w = _mm512_loadu_si512((const void *)(buf + j));
        double t[8];
        for (int k = 0; k < 8; k++) {
            unsigned i = (unsigned)(buf[j + k] & 0xff);
            t[k] = MX[i];
            rej[nr] = (uint16_t)(j + k); nr += i >= M_IMAX;
        }
        _mm512_storeu_pd(o + j, _mm512_mul_pd(_mm512_loadu_pd(t), _mm512_cvtepi64_pd(w)));
    }
    for (; j < m; j++) { unsigned i = buf[j] & 0xff; o[j] = MX[i] * (double)(int64_t)buf[j]; rej[nr] = j; nr += i >= M_IMAX; }
    fix(buf, o, rej, nr);
}
static int m3_ok(void) { __builtin_cpu_init(); return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"); }
#else
#  define M3_NAME "M3 (none on this arch)"
static void m3(const uint64_t *buf, double *o, int m) { m2(buf, o, m); }
static int m3_ok(void) { return 1; }
#endif

typedef void (*fn)(const uint64_t *, double *, int);

int main(int argc, char **argv) {
    int reps = argc > 1 ? atoi(argv[1]) : 20000;
    tables_init();
    uint64_t *buf = malloc(NW * 8); double *o = malloc(NW * 8), *ref = malloc(NW * 8);
    uint64_t s = 0x9e3779b97f4a7c15ULL;
    for (int j = 0; j < NW; j++) {
        uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        buf[j] = z ^ (z >> 31);
    }
    int rz = 0, rm = 0;
    for (int j = 0; j < NW; j++) {
        uint64_t r = buf[j], rabs = (r >> 9) & 0xfffffffffffffULL;
        rz += rabs - 1 >= kim1[r & 0xff]; rm += (r & 0xff) >= M_IMAX;
    }
    printf("edge rate on these %d words: NumPy tables %.2f%%, McFarland %.2f%%\n", NW, 100.0 * rz / NW, 100.0 * rm / NW);
    m0(buf, ref, NW);
    struct { const char *name; fn f; int ok; } v[] = {
        {"Z6 shipped (NumPy tables, signed wi)", z6, 1},
        {"M0 McFarland, branch per draw", m0, 1}, {"M1 M0 unrolled x4", m1, 1},
        {"M2 branch-free, rejects listed", m2, 1}, {M3_NAME, m3, m3_ok()}};
    for (int i = 2; i < 5; i++) if (v[i].ok) {
        v[i].f(buf, o, NW);
        printf("%-40s identical to M0: %s\n", v[i].name, memcmp(o, ref, NW * 8) ? "NO" : "yes");
    }
    for (int round = 1; round <= 3; round++) {
        printf("-- round %d (ns per draw, words in cache)\n", round);
        for (int i = 0; i < 5; i++) {
            if (!v[i].ok) { printf("  %-40s skipped\n", v[i].name); continue; }
            v[i].f(buf, o, NW);
            double t = now();
            for (int r = 0; r < reps; r++) v[i].f(buf, o, NW);
            printf("  %-40s %.3f\n", v[i].name, (now() - t) / reps / NW * 1e9);
        }
    }
    return 0;
}
