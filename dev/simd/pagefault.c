/* What a fresh 80 MB vector costs on this kernel, and whether asking for
 * huge pages or pre-faulting the region before the fill changes it. R
 * allocates large vectors with malloc(), which glibc serves with a fresh
 * mmap() above its threshold, so every call to a sampler pays these page
 * faults; so does every competitor. Linux only for the madvise variants;
 * elsewhere only the baselines run.
 *
 *   cc -O2 -o /tmp/pf ../dev/simd/pagefault.c && /tmp/pf 10000000
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef __linux__
#  include <sys/mman.h>
#endif
#ifndef MADV_POPULATE_WRITE
#  define MADV_POPULATE_WRITE 23
#endif

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
                          return t.tv_sec + t.tv_nsec * 1e-9; }

/* the fill: a cheap sequential store stream, like a fast generator */
static void fill(double *p, long n) {
    uint64_t s = 0x9e3779b97f4a7c15ULL;
    for (long i = 0; i < n; i++) { s += 0x9e3779b97f4a7c15ULL; p[i] = (double)(s >> 11) * 0x1.0p-53; }
}

static void advise(double *p, long n, int hugepage, int populate, int nohuge) {
#ifdef __linux__
    uintptr_t lo = ((uintptr_t)p + 4095) & ~(uintptr_t)4095;
    uintptr_t hi = ((uintptr_t)p + (size_t)n * 8) & ~(uintptr_t)4095;
    if (hi <= lo) return;
    if (nohuge) madvise((void *)lo, hi - lo, MADV_NOHUGEPAGE);
    if (hugepage) madvise((void *)lo, hi - lo, MADV_HUGEPAGE);
    if (populate) madvise((void *)lo, hi - lo, MADV_POPULATE_WRITE);
#else
    (void)p; (void)n; (void)hugepage; (void)populate; (void)nohuge;
#endif
}

int main(int argc, char **argv) {
    long n = argc > 1 ? atol(argv[1]) : 10000000L;
    int reps = 8;
#ifdef __linux__
    FILE *f = fopen("/sys/kernel/mm/transparent_hugepage/enabled", "r");
    char line[128] = "?";
    if (f) { if (!fgets(line, sizeof line, f)) line[0] = 0; fclose(f); }
    printf("THP enabled: %s", line);
#endif
    struct { const char *name; int reuse, hp, pop, memset_only, nohp; } v[] = {
        {"reused buffer (no faults)", 1, 0, 0, 0},
        {"fresh malloc, memset 0 only", 0, 0, 0, 1},
        {"fresh malloc, fill", 0, 0, 0, 0},
        {"fresh malloc, MADV_HUGEPAGE, fill", 0, 1, 0, 0},
        {"fresh malloc, MADV_POPULATE_WRITE, fill", 0, 0, 1, 0},
        {"fresh malloc, both, fill", 0, 1, 1, 0},
        {"fresh malloc, MADV_NOHUGEPAGE (4 KiB pages), fill", 0, 0, 0, 0, 1},
    };
    double *keep = malloc((size_t)n * 8); fill(keep, n);
    double acc = 0;
    printf("n = %ld doubles (%.0f MB), %d reps; ns per value\n", n, n * 8 / 1e6, reps);
    for (int i = 0; i < 7; i++) {
        double t = now();
        for (int r = 0; r < reps; r++) {
            double *p = v[i].reuse ? keep : malloc((size_t)n * 8);
            if (v[i].hp || v[i].pop || v[i].nohp) advise(p, n, v[i].hp, v[i].pop, v[i].nohp);
            if (v[i].memset_only) memset(p, 0, (size_t)n * 8); else fill(p, n);
            acc += p[n - 1];
            if (!v[i].reuse) free(p);
        }
        printf("%-50s %7.3f\n", v[i].name, (now() - t) / reps / n * 1e9);
    }
    printf("(%g)\n", acc);
    return 0;
}
