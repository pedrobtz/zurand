#' Create a random number generator key
#'
#' An `rngat_key` identifies an independent random stream. Together with an
#' index (and optionally a domain) it fully determines every value returned
#' by [bits_at()], [runif_at()] and [rnorm_at()]: no state is kept or
#' advanced. Keys are 16-byte raw vectors and can be stored, compared and
#' serialized; the byte layout is platform independent.
#'
#' @param seed A single whole number.
#' @return A 128-bit key of class `rngat_key`.
#' @seealso [fold_in()] to derive new keys, [runif_at()] for drawing values.
#' @export
#' @examples
#' key <- rng_key(42)
#' key
rng_key <- function(seed) {
  .Call(C_rng_key, seed)
}

#' Derive a new key by folding data into an existing key
#'
#' Deterministically derives a new, statistically independent key from
#' `key` and `identity`. Use this to give substreams of a computation
#' (chains, workers, blocks) their own keys without coordinating between
#' them: `fold_in(key, i)` for distinct `i` yield unrelated streams.
#'
#' Key derivation is domain-separated from value draws: keys produced by
#' `fold_in(key, i)` share no bits with values drawn at index `i`.
#'
#' @param key An `rngat_key`, from [rng_key()] or [fold_in()].
#' @param identity A single whole number distinguishing the derived key.
#' @return A new `rngat_key`.
#' @export
#' @examples
#' key <- rng_key(42)
#' chain_keys <- lapply(1:4, function(i) fold_in(key, i))
fold_in <- function(key, identity) {
  .Call(C_fold_in, key, identity)
}

#' Draw random values at arbitrary indices
#'
#' Random access random numbers: each value is a pure function of
#' `(key, index, domain)`, computed with the Philox4x64-10 counter-based
#' generator from the Random123 library. Evaluating the same triple always
#' gives the same value, in any order, on any platform — there is no
#' generator state to set, advance or restore.
#'
#' * `bits_at()` returns raw random bits as a 32-bit unsigned integer
#'   value (stored in a double, range 0 to 2^32 - 1).
#' * `runif_at()` returns doubles strictly inside (0, 1), built from the
#'   top 53 bits of the generator output.
#' * `rnorm_at()` returns standard normal deviates via inversion, so
#'   `rnorm_at(k, i)` is identical to `qnorm(runif_at(k, i))`.
#'
#' @param key An `rngat_key`, from [rng_key()] or [fold_in()].
#' @param index Vector of whole numbers: the positions to evaluate.
#' @param domain Optional vector of whole numbers selecting an independent
#'   stream for the same indices (default `NULL` is domain 0). Either a
#'   single value or the same length as `index`.
#' @return A double vector as long as `index` (or `domain`, if `index`
#'   has length 1).
#' @export
#' @examples
#' key <- rng_key(42)
#' runif_at(key, 1:5)
#' runif_at(key, 3)              # element 3, without generating 1:2
#' runif_at(key, 1:5, domain = 1) # an unrelated stream at the same indices
#' rnorm_at(key, 1:5)
bits_at <- function(key, index, domain = NULL) {
  .Call(C_bits_at, key, index, if (is.null(domain)) 0 else domain)
}

#' @rdname bits_at
#' @export
runif_at <- function(key, index, domain = NULL) {
  .Call(C_runif_at, key, index, if (is.null(domain)) 0 else domain)
}

#' @rdname bits_at
#' @export
rnorm_at <- function(key, index, domain = NULL) {
  .Call(C_rnorm_at, key, index, if (is.null(domain)) 0 else domain)
}

#' Draw a contiguous run of random values
#'
#' Sequential fast path for the common case of a block of consecutive
#' indices: `runif_seq(key, start, len)` returns exactly
#' `runif_at(key, start + 0:(len - 1))` — bit for bit the same values —
#' but without building, reading or validating an index vector, and
#' emitting four values per Philox call throughout. Use it for bulk
#' generation; use the `*_at()` forms for scattered indices.
#'
#' There is deliberately no `by` argument: strides of 4 or more share no
#' Philox blocks, so a strided draw has no fast path — write it as
#' `runif_at(key, seq(start, by = ..., length.out = ...))` instead.
#'
#' @inheritParams bits_at
#' @param start A single whole number: the first index of the run.
#' @param len A single non-negative whole number: how many values.
#' @param domain Optional single whole number selecting an independent
#'   stream (default `NULL` is domain 0).
#' @return A double vector of length `len`.
#' @export
#' @examples
#' key <- rng_key(42)
#' identical(runif_seq(key, 1, 10), runif_at(key, 1:10))
#' runif_seq(key, 1e6, 5)   # 5 values starting at index 1e6
bits_seq <- function(key, start, len, domain = NULL) {
  .Call(C_bits_seq, key, start, len, if (is.null(domain)) 0 else domain)
}

#' @rdname bits_seq
#' @export
runif_seq <- function(key, start, len, domain = NULL) {
  .Call(C_runif_seq, key, start, len, if (is.null(domain)) 0 else domain)
}

#' @rdname bits_seq
#' @export
rnorm_seq <- function(key, start, len, domain = NULL) {
  .Call(C_rnorm_seq, key, start, len, if (is.null(domain)) 0 else domain)
}

#' @rdname rng_key
#' @param x An `rngat_key`.
#' @param ... Unused.
#' @export
format.rngat_key <- function(x, ...) {
  paste(format(unclass(x)), collapse = "")
}

#' @rdname rng_key
#' @export
print.rngat_key <- function(x, ...) {
  cat("<rngat_key> ", format(x), "\n", sep = "")
  invisible(x)
}
