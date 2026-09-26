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

#define ZURAND_API_VERSION 1

/* Engine codes, as recorded in zurand_key.engine. */
#define ZURAND_ENGINE_PHILOX4X64   0
#define ZURAND_ENGINE_THREEFRY4X64 1
#define ZURAND_ENGINE_XOSHIRO256PP 2

/* Return codes. */
#define ZURAND_OK      0
#define ZURAND_EINVAL  1   /* an argument is out of range, NaN or infinite */
#define ZURAND_EKEY    2   /* unknown engine, or a stream this zurand lacks */

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
