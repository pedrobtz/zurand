#' Create a stateless RNG key
#'
#' `rng_key()` creates immutable keys for the stateless generators in this
#' package. A key vector is an opaque S3 object with one key per row; pass it to
#' [rng_fold()] and the `rng_*()` samplers, but do not rely on its internal
#' integer representation.
#'
#' The samplers in this package never read or update `.Random.seed`. Calling a
#' sampler twice with the same key and arguments returns the same values.
#'
#' @param seed A single non-negative whole number.
#' @param n A single non-negative whole number: how many keys to create.
#' @param engine The engine. All three give the same guarantees: a value is
#'   a pure function of the key and the draw index, the same key gives the
#'   same values on every platform and thread count, and draw `i` does not
#'   depend on `n`. They differ in speed and in the stream they produce.
#'
#'   * `"xoshiro256pp"` (default) -- the fastest of the three, and on one
#'     thread about 1.6-2.5x dqrng on uniform and 1.9-2.7x RcppZiggurat on
#'     normal (x86_64 and Apple Silicon; see the performance article). Philox derives a fresh 256-bit xoshiro256++ state for each
#'     512-value chunk, then a cheap recurrence produces the chunk, so a
#'     chunk depends only on the key and its index. Reaching an arbitrary
#'     index costs at most 511 recurrence steps. Audited with PractRand to
#'     1 TB, including interleaved streams of sibling and folded keys.
#'   * `"philox4x64"` -- Philox4x64-10 from Random123, counter-based: each
#'     value is computed directly from its position.
#'   * `"threefry4x64"` -- Threefry4x64-13 from Random123, counter-based;
#'     about 1.1x Philox in practice.
#'
#'   A key records its engine, so a key made with one never produces the
#'   other's values. Changing engine changes every number you get.
#'
#' @section Stream version:
#' A key also records the version of the stream definition it was made
#' under, as `attr(key, "stream")`; every key made by this version of
#' zurand is stream `1L`. If a defect in a stream is ever found after
#' release, the fix becomes a new stream version for new keys, and existing
#' keys -- including saved ones -- keep producing exactly the values they
#' always did. A key saved before the attribute existed is treated as
#' stream 1. A key from a newer stream than this version of zurand
#' implements is an error rather than a silent reinterpretation.
#' @return An object of class `rng_key` containing `n` keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' keys <- rng_key(42L, n = 4L)
#' key
rng_key <- function(seed, n = 1L,
                    engine = c("xoshiro256pp", "philox4x64",
                               "threefry4x64")) {
  engine <- match.arg(engine)
  .Call(C_rng_key, seed, n, engine)
}

#' Create a key from R's RNG
#'
#' `rng_key_from_r()` is the one convenience function in this package that
#' intentionally consumes R's global RNG state. Use it when you want a fresh
#' stateless key from the current `.Random.seed`; all other functions in this
#' package leave `.Random.seed` untouched.
#'
#' @inheritParams rng_key
#' @return An object of class `rng_key` containing `n` keys.
#' @export
#' @examples
#' set.seed(1)
#' key <- rng_key_from_r()
#' rng_uniform(key, 3L)
#'
#' # The same seed gives the same key.
#' set.seed(1)
#' identical(rng_key_from_r(), key)
rng_key_from_r <- function(n = 1L,
                           engine = c("xoshiro256pp", "philox4x64",
                               "threefry4x64")) {
  engine <- match.arg(engine)
  .Call(C_rng_key_from_r, n, engine)
}

#' Fold data into a key
#'
#' `rng_fold()` deterministically derives a key from `key` and simple data.
#' This is useful for structured reproducibility, for example deriving a key
#' from names such as `"chain 3"` or counters such as `10L`.
#'
#' @param key An `rng_key` vector.
#' @param data A character, integer, double, logical or raw vector. Numeric
#'   values must be whole numbers; integer and double vectors holding the
#'   same values (`1L` and `1`) derive the same key. All other type
#'   distinctions matter: `rng_fold(key, "1")` and `rng_fold(key, 1)` differ.
#' @return An `rng_key` vector with the same length as `key`.
#' @export
#' @examples
#' key <- rng_key(42L)
#' layer_key <- rng_fold(key, "layer1")
rng_fold <- function(key, data) {
  .Call(C_rng_fold, key, data)
}

#' Draw uniform random values
#'
#' Stateless uniform sampling. `rng_uniform(key, n, min, max)` is a pure
#' function of its arguments and returns ordinary numeric vectors or matrices.
#'
#' With the default bounds every value lies strictly inside \eqn{(0, 1)}: it
#' is \eqn{(m + 1/2) 2^{-52}} for a 52-bit integer \eqn{m}, so neither
#' endpoint can occur. Other bounds are applied as `min + (max - min) * u`,
#' which, as with [stats::runif()], will not return either extreme value
#' unless `max = min` or `max - min` is small compared with `min`.
#'
#' @param key An `rng_key` vector.
#' @param n A single non-negative whole number. Defaults to one draw per key.
#' @param min,max Single finite numeric bounds.
#' @return A numeric vector of length `n` for one key, or an `n` by
#'   `length(key)` numeric matrix for multiple keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_uniform(key)
#' rng_uniform(key, 5L)
#' rng_uniform(rng_key(42L, n = 3L), 5L)
rng_uniform <- function(key, n = 1L, min = 0, max = 1) {
  .Call(C_rng_uniform, key, n, min, max)
}

#' Draw normal random values
#'
#' Stateless normal sampling. `rng_normal(key, n, mean, sd)` is a pure function
#' of its arguments and returns ordinary numeric vectors or matrices.
#'
#' @inheritParams rng_uniform
#' @param mean A single finite numeric mean.
#' @param sd A single non-negative finite numeric standard deviation.
#' @return A numeric vector of length `n` for one key, or an `n` by
#'   `length(key)` numeric matrix for multiple keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' rng_normal(key)
#' rng_normal(key, 5L)
rng_normal <- function(key, n = 1L, mean = 0, sd = 1) {
  .Call(C_rng_normal, key, n, mean, sd)
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

#' Control zurand's thread usage
#'
#' `rng_threads()` reports the maximum number of OpenMP threads the samplers
#' may use; `rng_threads(threads)` caps it for the rest of the session. The
#' cap applies only to zurand's own parallel fills, not to other OpenMP code
#' in the process. Results never depend on the thread count, so this is a
#' performance control, not a reproducibility one. In builds without OpenMP
#' the value is always 1 and setting a cap has no effect.
#'
#' Draws below an internal size threshold always run single-threaded.
#'
#' @param threads `NULL` to query the effective maximum, a single whole
#'   number of at least 1 to cap zurand's thread use, or `0` to remove the
#'   cap. Values above the machine's thread count are allowed and
#'   equivalent to no cap.
#' @return Querying returns the effective maximum thread count as an
#'   integer. Setting invisibly returns the previous *configured* cap,
#'   with `0L` meaning "no cap was set", so `old <- rng_threads(1L)`
#'   followed by `rng_threads(old)` restores the exact prior
#'   configuration -- including the uncapped state, which keeps tracking
#'   later changes to the process-wide OpenMP maximum.
#' @export
#' @examples
#' rng_threads()
#' old <- rng_threads(1L)
#' rng_threads(old)
rng_threads <- function(threads = NULL) {
  prev <- .Call(C_rng_threads, threads)
  if (is.null(threads)) prev else invisible(prev)
}

#' @rdname rng_key
#' @param x An `rng_key`.
#' @param ... Unused.
#' @export
format.rng_key <- function(x, ...) {
  .Call(C_rng_key_format, x)
}

#' @rdname rng_key
#' @export
print.rng_key <- function(x, ...) {
  n <- length(x)
  fp <- format(x)
  if (n == 1L) {
    cat("<rng_key> ", fp, "\n", sep = "")
  } else {
    cat("<rng_key[", n, "]>\n", sep = "")
    # "rng_key[" plus the first 12 of 32 hex digits
    fp <- paste0(substr(fp, 1L, 20L), "...]")
    show <- min(n, 5L)
    for (i in seq_len(show)) cat("[", i, "] ", fp[[i]], "\n", sep = "")
    if (n > show) cat("... ", n - show, " more\n", sep = "")
  }
  invisible(x)
}

#' @rdname rng_key
#' @export
length.rng_key <- function(x) {
  nrow(unclass(x))
}

#' @rdname rng_key
#' @param i Integer or logical index selecting keys. Missing values and
#'   fractional numbers are an error: a silently invented key would
#'   deterministically collide with every other key subset the same way.
#' @param drop Ignored; key subsetting always preserves the `rng_key` class.
#' @export
`[.rng_key` <- function(x, i, ..., drop = FALSE) {
  words <- unclass(x)
  if (missing(i)) i <- seq_len(nrow(words))
  if (anyNA(i)) {
    stop("`i` must not contain missing values", call. = FALSE)
  }
  if (is.numeric(i) && any(i != trunc(i))) {
    stop("`i` must contain whole numbers", call. = FALSE)
  }
  out <- words[i, , drop = FALSE]
  class(out) <- "rng_key"
  attr(out, "engine") <- attr(x, "engine", exact = TRUE)
  attr(out, "stream") <- key_stream(x)
  out
}

# A key without the attribute predates it, and is stream 1.
key_stream <- function(x) {
  s <- attr(x, "stream", exact = TRUE)
  if (is.null(s)) 1L else as.integer(s)
}

#' @rdname rng_key
#' @export
`[[.rng_key` <- function(x, i, ...) {
  out <- x[i]
  if (length(out) != 1L) {
    stop("`[[` must select exactly one key", call. = FALSE)
  }
  out
}

#' @rdname rng_key
#' @param recursive Ignored. Included for compatibility with [c()].
#' @export
c.rng_key <- function(..., recursive = FALSE) {
  dots <- Filter(Negate(is.null), list(...))
  if (!length(dots)) {
    return(rng_key(0, n = 0L))
  }

  ok <- vapply(dots, inherits, logical(1), "rng_key")
  if (!all(ok)) {
    stop("all inputs to `c.rng_key()` must be <rng_key> objects", call. = FALSE)
  }

  engines <- vapply(dots, attr, character(1), "engine", exact = TRUE)
  if (length(unique(engines)) != 1L) {
    stop("all keys must use the same engine", call. = FALSE)
  }
  streams <- vapply(dots, key_stream, integer(1))
  if (length(unique(streams)) != 1L) {
    stop("all keys must use the same stream version", call. = FALSE)
  }

  out <- do.call(rbind, lapply(dots, unclass))
  dimnames(out) <- NULL
  class(out) <- "rng_key"
  attr(out, "engine") <- engines[[1L]]
  attr(out, "stream") <- streams[[1L]]
  out
}

#' Report or disable the SIMD path
#'
#' `rng_simd()` reports which instruction path the `"xoshiro256pp"` engine
#' is using; `rng_simd(FALSE)` forces the portable one.
#'
#' This is a performance control and nothing else. Both paths emit the same
#' values: the vectorised one runs several independent 512-value sub-chunks
#' side by side -- four per AVX2 register on x86_64, two per NEON register
#' on arm64 -- each executing the same recurrence as the scalar path, so
#' output does not depend on the instruction set, the CPU, or this setting.
#' It is exposed because that claim is worth being able to check, and the
#' package's own tests check it by running both and comparing.
#'
#' Only `"xoshiro256pp"` has a vectorised path. The counter-based engines
#' are unaffected, as is every other function.
#'
#' @param enable `NULL` to query, `FALSE` to force the portable path, or
#'   `TRUE` to go back to using whatever the CPU supports.
#' @return A string: `"avx2"` (x86_64 with AVX2) or `"neon"` (arm64) if the
#'   vectorised path is active, otherwise `"none"`. When setting, the
#'   previous value, invisibly.
#' @export
#' @examples
#' rng_simd()
#' old <- rng_simd(FALSE)
#' identical(rng_uniform(rng_key(1L, engine = "xoshiro256pp"), 4L),
#'           { rng_simd(old != "none"); rng_uniform(rng_key(1L, engine = "xoshiro256pp"), 4L) })
rng_simd <- function(enable = NULL) {
  prev <- .Call(C_rng_simd, enable)
  if (is.null(enable)) prev else invisible(prev)
}
