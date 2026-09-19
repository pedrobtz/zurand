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
#' @param engine The counter-based engine. Both are Random123 generators, so
#'   both give the same guarantee: a value is a pure function of the key and
#'   the draw index, reachable directly without generating what precedes it.
#'   They differ only in speed and in the stream they produce.
#'
#'   * `"philox4x64"` (default) -- Philox4x64-10. The original engine; its
#'     stream is fixed and will not change.
#'   * `"threefry4x64"` -- Threefry4x64-13, roughly 1.8x faster per word.
#'
#'   A key records its engine, so a key made with one never produces the
#'   other's values. Changing engine changes every number you get.
#' @return An object of class `rng_key` containing `n` keys.
#' @export
#' @examples
#' key <- rng_key(42L)
#' keys <- rng_key(42L, n = 4L)
#' key
rng_key <- function(seed, n = 1L,
                    engine = c("philox4x64", "threefry4x64")) {
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
rng_key_from_r <- function(n = 1L,
                           engine = c("philox4x64", "threefry4x64")) {
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
  out
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

  out <- do.call(rbind, lapply(dots, unclass))
  dimnames(out) <- NULL
  class(out) <- "rng_key"
  attr(out, "engine") <- engines[[1L]]
  out
}

#' Allocate a buffer for the in-place fills
#'
#' `rng_buffer()` creates a numeric vector to be filled repeatedly by
#' [rng_fill_uniform()] and [rng_fill_normal()]. It is an ordinary numeric
#' vector carrying one marker attribute, which [rng_fill_uniform()] requires:
#' a vector that was not made for overwriting cannot be passed to them by
#' accident. Unlike `numeric(n)` it is not zero-filled, since the fill
#' overwrites it immediately.
#'
#' @param n A single non-negative whole number: the buffer length.
#' @return A numeric vector of length `n`.
#' @seealso [rng_fill_uniform()] for the hazard these buffers carry.
#' @export
#' @examples
#' buf <- rng_buffer(1000)
#' rng_fill_uniform(buf, rng_key(42L))
#' mean(buf)
rng_buffer <- function(n) {
  # Allocated in C: Rf_allocVector leaves the memory uninitialised, where
  # numeric(n) would zero-fill it first -- a second pass over memory that
  # the fill immediately overwrites.
  .Call(C_rng_buffer, n)
}

#' Fill a buffer in place
#'
#' These write directly into `buf` instead of returning a new vector, which
#' skips the allocation and garbage collection R otherwise pays for a fresh
#' result. That is worth about 15% at a million values and more than double
#' at fifty million, because the fill itself runs at a flat ~323 M values/s
#' at every size while allocation cost grows.
#'
#' @section The catch:
#'
#' **These break R's copy-on-write rule on purpose.** Everywhere else in R,
#' passing a vector to a function cannot change it. Here it does, and
#' anything else referring to the same vector changes with it:
#'
#' ```r
#' x <- rng_buffer(3); y <- x
#' rng_fill_uniform(x, rng_key(1L))
#' y   # also overwritten
#' ```
#'
#' The case that actually bites is storing the buffer while reusing it,
#' because every slot ends up holding the same vector:
#'
#' ```r
#' buf <- rng_buffer(1000); out <- list()
#' for (i in 1:3) {
#'   rng_fill_uniform(buf, keys[i])
#'   out[[i]] <- buf        # not a copy
#' }
#' # all three entries hold the last draw
#' ```
#'
#' **Nothing detects that loop.** An earlier version watched the buffer's
#' reference count for growth, which caught it; `REFCNT` is not part of R's
#' API, `R CMD check` reports it as a non-API call and CRAN does not permit
#' it. `MAYBE_SHARED` is permitted but is true for every realistic call, so
#' it cannot stand in. There is no API-legal way to tell a safe call from
#' an unsafe one, so requiring a buffer from [rng_buffer()] is the only
#' defence: it stops these functions being reached by accident, and does
#' nothing once you hand a buffer over deliberately.
#'
#' Use `out[[i]] <- buf[]` to copy results out, or [rng_uniform()]. Use `out[[i]] <- buf[]`
#' or plain [rng_uniform()] when you need to keep the results.
#'
#' Use these only for a buffer you allocated, fill, consume and refill.
#' When in doubt use [rng_uniform()] and [rng_normal()], which allocate and
#' are always safe.
#'
#' @param buf A numeric vector to overwrite, from [rng_buffer()].
#' @param key An `rng_key` vector of length one.
#' @param min,max Single finite numeric bounds.
#' @param mean A single finite numeric mean.
#' @param sd A single non-negative finite standard deviation.
#' @return `buf`, invisibly, filled. The return is for convenience; the
#'   argument has already been modified.
#' @export
#' @examples
#' buf <- rng_buffer(1000)
#' rng_fill_normal(buf, rng_key(42L), mean = 10, sd = 2)
#' round(mean(buf), 1)
rng_fill_uniform <- function(buf, key, min = 0, max = 1) {
  invisible(.Call(C_rng_fill_uniform, buf, key, min, max))
}

#' @rdname rng_fill_uniform
#' @export
rng_fill_normal <- function(buf, key, mean = 0, sd = 1) {
  invisible(.Call(C_rng_fill_normal, buf, key, mean, sd))
}
