#' Lazy random vectors
#'
#' `rng_lazy_uniform()` and `rng_lazy_normal()` return a numeric vector of
#' length `n` whose elements are computed only when read. Every value is a
#' pure function of the key and its position, so element `i` is exactly
#' element `i` of `rng_uniform(key, n, ...)` or `rng_normal(key, n, ...)`,
#' and `identical(x[], rng_normal(key, n))` holds. Creating one is instant
#' and takes no memory for the values, whatever `n` is: each lazy vector
#' holds its key, its arguments and a cache of one engine chunk, which is
#' `min(n, 5120)` doubles, so at most about 40 KB.
#'
#' @section Reading:
#' `x[i]`, `x[idx]`, `x[a:b]`, `sum(x)`, `mean(x)` and loops over
#' `x` compute only the values they touch; nothing is kept but the last
#' engine chunk read (at most 5120 values), so `sum()` over 1e9 lazy values
#' runs in the memory R itself uses. Full scans cost more than
#' `rng_normal()` followed by `sum()` would, because R reads lazy vectors
#' 512 values at a time on one thread. The gain is memory, and not paying
#' for values that are never read.
#'
#' @section Interrupts:
#' Scans that read a lazy vector in blocks -- `sum()`, `mean()`,
#' `max()`, `x[idx]` -- check for a user interrupt every \eqn{2^{22}} values,
#' so Ctrl-C stops `sum()` over 1e12 values. Reading single elements and
#' materialising do not check. For C and C++ code this means that
#' `REAL_ELT()` and `REAL()`/`REAL_RO()` on a lazy vector never signal an
#' interrupt, but `REAL_GET_REGION()` -- which cpp11's `doubles` iterator
#' uses -- can, by a longjmp that skips C++ destructors. Code that cannot be
#' unwound should take the data pointer first (which materialises the
#' vector) or read under `R_UnwindProtect()` (`cpp11::unwind_protect()`).
#'
#' @section Materialisation:
#' Anything that needs the whole vector as ordinary memory
#' -- most arithmetic such as `x * 2`, modifying `x`, or C code asking for
#' its data pointer -- fills it once, with the usual parallel fill, and
#' from then on `x` holds all `n` values like any vector.
#'
#' @section Saving:
#' [saveRDS()] and [serialize()] store only the key and the
#' arguments, so a lazy vector of 1e9 values saves as a few hundred bytes
#' and comes back lazy. A lazy vector that has been modified saves its
#' values. Reading is not modifying: after C code has had the data pointer
#' (`sd()`, `sort()` and `cumsum()` all take it), saving regenerates the
#' values and compares them with what the vector holds, which costs about
#' one more fill, and stores only the key when nothing changed. Restoring
#' checks the saved arguments and refuses ones it cannot have written.
#'
#' @section Engines:
#' The vector uses its key's engine, so its values are that
#' engine's. For scattered reads, `x[sample(n, k)]` or random indices in a
#' loop, key it with `rng_key(seed, engine = "philox4x64")`: Philox computes
#' any single value directly, several times faster per read than
#' `"xoshiro256pp"`, which reaches it through its 512-value block. For scans
#' and sequential loops the default is faster.
#'
#' @inheritParams rng_normal
#' @param key A single `rng_key`.
#' @return A double vector of length `n` (an ALTREP object).
#' @export
#' @examples
#' key <- rng_key(42L)
#' x <- rng_lazy_normal(key, 1e9)      # instant; no values yet
#' x[c(1, 5e8, 1e9)]
#' identical(rng_lazy_normal(key, 10)[], rng_normal(key, 10))
rng_lazy_normal <- function(
  key,
  n,
  mean = 0,
  sd = 1,
  method = c("ziggurat", "mcfarland")
) {
  method <- match.arg(method)
  # ZURAND_DIST_NORMAL_MCF = 2, ZURAND_DIST_NORMAL = 1 (src/zurand.c); 0 is
  # ZURAND_DIST_UNIFORM, below. Not the C API's method codes.
  .Call(C_rng_lazy, key, n, if (method == "mcfarland") 2L else 1L, mean, sd)
}

#' @rdname rng_lazy_normal
#' @param min,max Single finite numeric bounds.
#' @export
rng_lazy_uniform <- function(key, n, min = 0, max = 1) {
  .Call(C_rng_lazy, key, n, 0L, min, max)
}

# Whether a lazy vector has been materialised (NA for anything else); for
# tests and measurement.
rng_lazy_materialised <- function(x) .Call(C_rng_lazy_materialised, x)

# How a single read that misses the cache is served, for measurement: 0
# computes the value alone, 2 fills the engine chunk around it, 1 (the
# default) fills the chunk only when reads walk forward. Returns the
# previous mode.
rng_lazy_cache <- function(mode = NULL) {
  .Call(C_rng_lazy_cache, if (is.null(mode)) NULL else as.integer(mode))
}
