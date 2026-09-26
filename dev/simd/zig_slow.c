/* The ziggurat slow path, taken apart: how often a rejected draw reaches
 * exp() after the fixed-point wedge brackets, and what zurand_exp costs
 * against the platform libm. Same tables, same bracket logic and the same
 * fdlibm port as src/zurand.c; the words come from splitmix64 because only
 * their distribution matters here.
 *
 *   cc -O2 -I../src -o /tmp/zs ../dev/simd/zig_slow.c -lm && /tmp/zs
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "numpyzig/ziggurat_constants.h"
#include "zigbounds.h"

static inline double zurand_rounded(double x) {
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__SSE2_MATH__))
    __asm__("" : "+x"(x));
#elif defined(__GNUC__) && defined(__aarch64__)
    __asm__("" : "+w"(x));
#else
    volatile double v = x; x = v;
#endif
    return x;
}
#include "zurand_fdlibm.h"
#define ZURAND_ZIG_GUARD ((uint64_t)4096)

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }
static uint64_t sm = 0x9e3779b97f4a7c15ULL;
static inline uint64_t next(void) {
    uint64_t z = (sm += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint64_t mulhi(uint64_t a, uint64_t b) { return (uint64_t)(((__uint128_t)a * b) >> 64); }

int main(void) {
    /* 1. where slow-path entries go */
    long draws = 20000000, entries = 0, tail = 0, layer1 = 0, amb = 0, acc = 0, rej = 0;
    long exp_calls = 0, loops = 0;
    for (long i = 0; i < draws; i++) {
        uint64_t r = next();
        int idx = (int)(r & 0xff);
        uint64_t rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
        if (rabs < ki_double[idx]) continue;
        entries++;
        if (idx == 1) layer1++;
        /* follow the real loop until acceptance, counting exp() calls */
        for (;;) {
            loops++;
            uint64_t Y = next();
            if (idx == 0) { tail++; break; }   /* tail: two log1p per attempt, rare */
            uint64_t L = (UINT64_C(1) << 52) - ki_double[idx];
            uint64_t R = (UINT64_C(1) << 52) - rabs;
            uint64_t YL = mulhi(Y, L);
            int a, j;
            if (idx > ZURAND_ZIG_INFLECTION) {
                j = YL > R + ZURAND_ZIG_GUARD; a = !j && YL + zurand_zig_gap[idx] < R;
            } else if (idx < ZURAND_ZIG_INFLECTION) {
                a = YL + ZURAND_ZIG_GUARD < R; j = !a && YL > R + zurand_zig_gap[idx];
            } else {
                j = YL > R + zurand_zig_gap_hi52; a = !j && YL + zurand_zig_gap[idx] < R;
            }
            if (a) { acc++; break; }
            if (!j) {
                amb++; exp_calls++;
                double x = (double)rabs * wi_double[idx];
                double u = (double)(Y >> 11) * 0x1.0p-53;
                double rise = zurand_rounded((fi_double[idx - 1] - fi_double[idx]) * u);
                if (rise + fi_double[idx] < zurand_exp(-0.5 * x * x)) break;
            } else rej++;
            r = next(); idx = (int)(r & 0xff);
            rabs = (r >> 9) & UINT64_C(0x000fffffffffffff);
            if (rabs < ki_double[idx]) break;
        }
    }
    printf("slow-path entries: %.3f%% of draws; of these layer 1: %.1f%%, tail: %.2f%%\n",
           100.0 * entries / draws, 100.0 * layer1 / entries, 100.0 * tail / entries);
    printf("wedge tests: %ld; resolved by bracket: accept %.1f%%, reject %.1f%%; ambiguous -> exp(): %.1f%%\n",
           acc + rej + amb, 100.0 * acc / (acc + rej + amb), 100.0 * rej / (acc + rej + amb), 100.0 * amb / (acc + rej + amb));
    printf("exp() calls per slow-path entry: %.3f; per draw: %.5f; loop iterations per entry: %.3f\n",
           (double)exp_calls / entries, (double)exp_calls / draws, (double)loops / entries);

    /* 2. what exp costs: independent calls (throughput) and a chain (latency) */
    long n = 5000000;
    double *args = malloc(n * sizeof *args);
    for (long i = 0; i < n; i++) args[i] = -6.5 * (double)(next() >> 11) * 0x1.0p-53;
    double s = 0, t;
    t = now(); for (long i = 0; i < n; i++) s += zurand_exp(args[i]); double tz = (now() - t) / n * 1e9;
    t = now(); for (long i = 0; i < n; i++) s += exp(args[i]);        double tl = (now() - t) / n * 1e9;
    double c = 0;
    t = now(); for (long i = 0; i < n; i++) c = zurand_exp(args[i] + c * 1e-300); double tzc = (now() - t) / n * 1e9;
    t = now(); for (long i = 0; i < n; i++) c = exp(args[i] + c * 1e-300);        double tlc = (now() - t) / n * 1e9;
    printf("exp ns/call   throughput: zurand_exp %.1f, libm %.1f   latency chain: zurand_exp %.1f, libm %.1f   (%g %g)\n",
           tz, tl, tzc, tlc, s, c);
    t = now(); for (long i = 0; i < n; i++) s += zurand_log1p(-args[i] / 7.0); double tp = (now() - t) / n * 1e9;
    t = now(); for (long i = 0; i < n; i++) s += log1p(-args[i] / 7.0);        double tq = (now() - t) / n * 1e9;
    printf("log1p ns/call throughput: zurand_log1p %.1f, libm %.1f  (%g)\n", tp, tq, s);
    return 0;
}
