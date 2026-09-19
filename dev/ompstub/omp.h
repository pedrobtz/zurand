/* Minimal stub: enough for -fsyntax-only to parse zurand's OpenMP pragmas
 * and enforce default(none), on a machine with no libomp. Never linked. */
#ifndef OMP_STUB_H
#define OMP_STUB_H
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_thread_num(void)  { return 0; }
static inline int omp_get_num_threads(void) { return 1; }
#endif
