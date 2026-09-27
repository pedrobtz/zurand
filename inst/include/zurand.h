/* zurand C API, for packages that declare `LinkingTo: zurand` and
 * `Imports: zurand`.
 *
 * It exposes the same fills the R samplers run, writing into memory the
 * caller owns. Every function except zurand_key_get() calls no R API,
 * allocates nothing and changes no global state, so once the API table has
 * been fetched on the main thread its fills may be called from any thread
 * -- OpenMP, pthreads, TBB -- and return exactly what rng_uniform(),
 * rng_normal(), rng_integer() and rng_bits() return for the same key.
 *
 *   #include <zurand.h>
 *
 *   SEXP my_sim(SEXP keys) {
 *       const zurand_api *zr = zurand_get_api();   // main thread, once
 *       zurand_key k;
 *       zr->key_get(keys, 0, &k);                    // main thread
 *       double *buf = ...;                           // caller-owned
 *       // ... then, from any thread:
 *       if (zr->fill_normal(k, n, 0.0, 1.0, buf) != ZURAND_OK) ...
 *   }
 *
 * zurand.h includes <Rinternals.h>. With OpenMP, include <omp.h> first or
 * define R_NO_REMAP: R's `match` macro breaks the `declare variant` pragmas
 * in LLVM's omp.h.
 *
 * Keys come from R: create them with rng_key() or rng_fold() and pass the
 * vector down. zr->fold_int() derives a key from a whole number exactly
 * as rng_fold(key, i) does, without R.
 *
 * The table is versioned. Members are only ever appended; a package built
 * against a newer header than the installed zurand gets an R error from
 * zurand_get_api() instead of reading past the end of the table.
 */
#ifndef ZURAND_H
#define ZURAND_H

#include <stddef.h>
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ZURAND_API_VERSION 3

/* Engine codes, as recorded in zurand_key.engine. */
#define ZURAND_ENGINE_PHILOX4X64   0
#define ZURAND_ENGINE_THREEFRY4X64 1
#define ZURAND_ENGINE_XOSHIRO256PP 2

/* Normal methods, for fill_normal_method(): rng_normal(method = "ziggurat")
 * and rng_normal(method = "mcfarland"). Separate streams: the same key
 * gives unrelated values under each. */
#define ZURAND_NORMAL_ZIGGURAT  0
#define ZURAND_NORMAL_MCFARLAND 1

/* Return codes. */
#define ZURAND_OK      0
#define ZURAND_EINVAL  1   /* an argument is out of range, NaN or infinite */
#define ZURAND_EKEY    2   /* unknown engine, or a stream this zurand lacks */
#define ZURAND_STOPPED 3   /* a stream's consumer returned nonzero */

/* A stream's consumer. x[0 .. n-1] are positions start .. start + n - 1
 * (0-based) of the stream; n is the stream's chunk except for a shorter
 * last one. x is valid only during the call: zurand overwrites it for the
 * next chunk. Return 0 to continue, nonzero to stop the stream. */
typedef int (*zurand_consumer)(void *ctx, const double *x, size_t n,
                               uint64_t start);

/* One key: the 128 key bits, its engine and its stream version. Plain data;
 * copy it freely, including into worker threads. */
typedef struct zurand_key {
    uint64_t k0, k1;
    int engine;
    int stream;
} zurand_key;

typedef struct zurand_api {
    int version;   /* ZURAND_API_VERSION of the installed zurand */
    size_t size;   /* sizeof(zurand_api) in the installed zurand */

    /* Key i (0-based) of an rng_key vector. Main thread only: it reads the
     * R object, and raises an R error for anything that is not a valid key
     * vector. Returns ZURAND_EINVAL when i is out of range. */
    int (*key_get)(SEXP keys, R_xlen_t i, zurand_key *out);

    /* rng_fold(key, data) for one whole number `data`: the same key words
     * as rng_fold() with a length-one integer or double vector. */
    int (*fold_int)(zurand_key key, int64_t data, zurand_key *out);

    /* rng_uniform(key, n, min, max) into out[0 .. n-1]. */
    int (*fill_uniform)(zurand_key key, size_t n, double min, double max,
                        double *out);

    /* rng_normal(key, n, mean, sd) into out[0 .. n-1]. */
    int (*fill_normal)(zurand_key key, size_t n, double mean, double sd,
                       double *out);

    /* rng_integer(key, n, min, max) into out[0 .. n-1]; min and max are
     * inclusive and must not be NA_INTEGER. */
    int (*fill_integer)(zurand_key key, size_t n, int min, int max, int *out);

    /* rng_bits(key, n, bits = 64) as integers rather than hex strings. */
    int (*fill_bits64)(zurand_key key, size_t n, uint64_t *out);

    /* Version 2. */

    /* rng_normal(key, n, mean, sd, method) into out[0 .. n-1], with method
     * ZURAND_NORMAL_ZIGGURAT (what fill_normal() uses) or
     * ZURAND_NORMAL_MCFARLAND; any other method is ZURAND_EINVAL. */
    int (*fill_normal_method)(zurand_key key, size_t n, double mean,
                              double sd, int method, double *out);

    /* Version 3: positions other than 0, pulled or streamed. */

    /* rng_uniform(key, n, min, max, offset = start) into out[0 .. n-1]:
     * positions start .. start + n - 1, the values rng_uniform(key,
     * start + n)[(start + 1):(start + n)] has. Calling it for consecutive
     * ranges walks one stream in pieces of any size. */
    int (*fill_uniform_at)(zurand_key key, uint64_t start, size_t n,
                           double min, double max, double *out);

    /* rng_normal(key, n, mean, sd, method, offset = start), likewise. */
    int (*fill_normal_at)(zurand_key key, uint64_t start, size_t n,
                          double mean, double sd, int method, double *out);

    /* Positions 0 .. total - 1 of rng_uniform(key, total, min, max), handed
     * to fn in order, `chunk` at a time, through buf[0 .. chunk - 1], which
     * the caller owns and zurand reuses, so the values can be consumed from
     * cache without the whole stream ever being in memory. fn is called on
     * the calling thread only, sequentially. Returns ZURAND_OK once all
     * `total` values have been consumed, or ZURAND_STOPPED if fn stopped it. */
    int (*stream_uniform)(zurand_key key, uint64_t total, double min,
                          double max, double *buf, size_t chunk,
                          zurand_consumer fn, void *ctx);

    /* rng_normal(key, total, mean, sd, method), streamed likewise. */
    int (*stream_normal)(zurand_key key, uint64_t total, double mean,
                         double sd, int method, double *buf, size_t chunk,
                         zurand_consumer fn, void *ctx);
} zurand_api;

/* Fetch the API table. Call on the main thread; the pointer stays valid for
 * the session and may then be used from any thread. */
static inline const zurand_api *zurand_get_api(void) {
    typedef const zurand_api *(*getter)(void);
    getter get = (getter) R_GetCCallable("zurand", "zurand_api");
    const zurand_api *api = get();
    if (api->version < ZURAND_API_VERSION ||
        api->size < sizeof(zurand_api))
        Rf_error("the installed zurand provides C API version %d, but this "
                 "package was built against version %d; update zurand",
                 api->version, ZURAND_API_VERSION);
    return api;
}

#ifdef __cplusplus
}
#endif

#endif /* ZURAND_H */
