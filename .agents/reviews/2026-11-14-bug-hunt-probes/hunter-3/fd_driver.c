#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
static inline double zurand_rounded(double x) {
#if defined(__x86_64__)
    __asm__("" : "+x"(x));
#else
    volatile double v = x; x = v;
#endif
    return x; }
#include "zurand_fdlibm.h"
int main(void) {
    char f; double x;
    while (scanf(" %c %la", &f, &x) == 2) {
        double y = f == 'e' ? zurand_exp(x) : zurand_log1p(x);
        printf("%a\n", y); fflush(stdout);
    }
    return 0;
}
