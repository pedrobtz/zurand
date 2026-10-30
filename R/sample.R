#' Draw uniform random values
#'
#' Stateless uniform sampling. `rng_uniform(key, n, min, max)` is a pure
#' function of its arguments and returns ordinary numeric vectors or matrices.
#'
#' With the default bounds every value lies strictly inside \eqn{(0, 1)}: it
#' is \eqn{(m + 1/2) 2^{-52}} for a 52-bit integer \eqn{m}, so neither
#' endpoint can occur. Other bounds are applied as `min + (max - min) * u`,
#' which, as with [stats::runif()], will not return either extreme value
#' unless `max == min` or `max - min` is small compared with `min`.
#'
#' @param key An `rng_key` vector.
#' @param n A single non-negative whole number. Defaults to one draw per key.
#' @param min,max Single finite numeric bounds.
#' @param offset A single non-negative whole number: the draws start at
#'   position `offset` of each key's stream, so the result holds positions
#'   `offset` to `offset + n - 1`. `f(key, n, offset = k)` equals
#'   `f(key, k + n)[(k + 1):(k + n)]` exactly, which lets a long simulation
#'   draw its next batch, or resume from a checkpoint, without holding the
#'   earlier values. `offset + n` must not exceed \eqn{2^{53}}.
#' @return A numeric vector of length `n` for one key, or an `n` by
#'   `length(key)` numeric matrix for multiple keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_uniform(key)
#' rng_uniform(key, 5L)
#' rng_uniform(rng_key(42L, n = 3L), 5L)
#' # the next five, on their own
#' identical(rng_uniform(key, 5L, offset = 5), rng_uniform(key, 10L)[6:10])
rng_uniform <- function(key, n = 1L, min = 0, max = 1, offset = 0) {
  .Call(C_rng_uniform, key, n, min, max, offset)
}

#' Draw normal random values
#'
#' Stateless normal sampling. `rng_normal(key, n, mean, sd)` is a pure function
#' of its arguments and returns ordinary numeric vectors or matrices.
#'
#' @inheritParams rng_uniform
#' @param mean A single finite numeric mean.
#' @param sd A single non-negative finite numeric standard deviation.
#' @param method The sampling algorithm:
#'   * `"ziggurat"` (the default): Marsaglia and Tsang's ziggurat on NumPy's
#'     256-layer tables.
#'   * `"mcfarland"`: McFarland's modified ziggurat (2016). Its common case
#'     needs no table comparison, so it is faster than the default on every
#'     CPU measured, by more on x86_64 than on arm64
#'     (`vignette("performance", package = "zurand")`).
#'
#'   Both are exact samplers of the standard normal (to double precision)
#'   and both are bit-reproducible across platforms. They are different
#'   streams: the same key gives unrelated values under the two methods, and
#'   changing the method changes every value. The default never changes.
#' @return A numeric vector of length `n` for one key, or an `n` by
#'   `length(key)` numeric matrix for multiple keys.
#' @references
#' Marsaglia, G. and Tsang, W. W. (2000). The ziggurat method for generating
#' random variables. *Journal of Statistical Software*, 5(8), 1-7.
#' \doi{10.18637/jss.v005.i08}
#'
#' McFarland, C. D. (2016). A modified ziggurat algorithm for generating
#' exponentially and normally distributed pseudorandom numbers. *Journal of
#' Statistical Computation and Simulation*, 86(7), 1281-1294.
#' \doi{10.1080/00949655.2015.1060234}
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_normal(key)
#' rng_normal(key, 5L)
#' rng_normal(key, 5L, method = "mcfarland")
#' identical(rng_normal(key, 5L, offset = 5), rng_normal(key, 10L)[6:10])
rng_normal <- function(key, n = 1L, mean = 0, sd = 1,
                       method = c("ziggurat", "mcfarland"), offset = 0) {
  method <- match.arg(method)
  .Call(C_rng_normal, key, n, mean, sd,
        # ZURAND_NORMAL_MCFARLAND = 1, ZURAND_NORMAL_ZIGGURAT = 0
        # (inst/include/zurand.h, shared with the C API).
        if (method == "mcfarland") 1L else 0L, offset)
}

#' Draw integer random values
#'
#' Stateless integer sampling over the inclusive range `[min, max]`.
#'
#' @inheritParams rng_uniform
#' @param min,max Single whole-number bounds in R's non-missing integer range.
#' @return An integer vector of length `n` for one key, or an `n` by
#'   `length(key)` integer matrix for multiple keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_integer(key, min = 1L, max = 6L)
#' rng_integer(key, 10L, min = 1L, max = 6L)
rng_integer <- function(key, n = 1L, min, max) {
  if (missing(min) || missing(max)) {
    stop("`min` and `max` are required", call. = FALSE)
  }
  .Call(C_rng_integer, key, n, min, max)
}

#' Draw low-level random bits
#'
#' `rng_bits()` exposes the raw generator stream. With `bits = 32`, it returns
#' unsigned 32-bit words stored exactly in a numeric vector. With `bits = 64`,
#' it returns fixed-width hexadecimal strings so no bits are lost to R's
#' numeric representation.
#'
#' @inheritParams rng_uniform
#' @param bits Either `32L` or `64L`.
#' @return A numeric vector or matrix for 32-bit words, or a character vector
#'   or matrix for 64-bit words.
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_bits(key)
#' rng_bits(key, 5L)
rng_bits <- function(key, n = 1L, bits = 32L) {
  .Call(C_rng_bits, key, n, bits)
}
