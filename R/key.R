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
#'   * `"xoshiro256pp"` (default) -- the fastest of the three, and faster
#'     than dqrng and RcppZiggurat on one thread (x86_64 and Apple Silicon;
#'     `vignette("performance", package = "zurand")` has the numbers).
#'     Philox derives a fresh 256-bit xoshiro256++ state for each
#'     512-value chunk, then a cheap recurrence produces the chunk, so a
#'     chunk depends only on the key and its index. Reaching an arbitrary
#'     index costs at most 511 recurrence steps. Audited with PractRand to
#'     1 TB, including interleaved streams of sibling and folded keys.
#'   * `"philox4x64"` -- Philox4x64-10 from Random123, counter-based: each
#'     value is computed directly from its position.
#'   * `"threefry4x64"` -- Threefry4x64-13 from Random123, counter-based;
#'     close to Philox in speed.
#'
#'   A key records its engine, so a key made with one never produces the
#'   other's values. Changing engine changes every number you get.
#'
#' @section Key vectors:
#' A key vector behaves like a vector of opaque keys: `length()`, `[`,
#' `[[`, `c()`, `rev()`, `rep()`, `unique()`, `duplicated()` and
#' `as.list()` work key by key and keep the engine and stream version, so
#' `lapply(keys, f)`, `Map()` and the parallel map functions hand each task
#' one key. A key vector takes one index, as a vector does: `keys[1, 2]`
#' and `keys[, 1]` are an error, not a view of the words inside. `c()`
#' dispatches on its first argument, so combine key vectors with a key
#' first: `c(NULL, key)` returns the bare integer words, not keys, while
#' `c(key, NULL)` drops the `NULL` as `c()` does for any vector. Keys are
#' immutable: assigning into a key vector is an error, and so is sorting
#' one, since keys have no order; compare keys with [identical()].
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
#'   Only the storage type and the values are folded, so classed objects
#'   such as factors and dates are an error: a factor would fold as its
#'   level codes, whose meaning changes with the level set. Fold
#'   `as.character(f)` or `unclass(d)` instead. Names and dimensions are
#'   ignored: a matrix folds as its values in column order. Strings are
#'   folded as their UTF-8 bytes, whatever their declared encoding.
#' @details
#' The data are hashed to 64 bits before they derive the key, so two
#' different values can in principle fold to the same key. Among `m`
#' distinct values the chance of any such collision is about
#' `m^2 / 2^65`: under one in a million for a million values, but not
#' negligible near `2^32` of them. For that many independent streams,
#' request them up front with `rng_key(seed, n)`, which never collides.
#' @return An `rng_key` vector with the same length as `key`.
#' @export
#' @examples
#' key <- rng_key(42L)
#' layer_key <- rng_fold(key, "layer1")
rng_fold <- function(key, data) {
  if (is.object(data)) {
    stop("`data` must be a plain vector, not a <", class(data)[[1L]],
         ">; fold as.character() or unclass() of it instead", call. = FALSE)
  }
  .Call(C_rng_fold, key, data)
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
#' @param i Integer or logical index selecting keys. Missing values,
#'   fractional numbers, infinities and positions beyond the key vector are
#'   an error: R would turn such an index into a key of missing words, the
#'   same invented key whatever vector it came from, so unrelated subsets
#'   would silently share one stream.
#' @param drop Ignored; key subsetting always preserves the `rng_key` class.
#' @export
`[.rng_key` <- function(x, i, ..., drop = FALSE) {
  if (...length() > 0L) {
    stop("a key vector takes one index, `keys[i]`", call. = FALSE)
  }
  words <- unclass(x)
  if (missing(i)) i <- seq_len(nrow(words))
  check_key_index(i, nrow(words))
  out <- words[i, , drop = FALSE]
  class(out) <- "rng_key"
  attr(out, "engine") <- attr(x, "engine", exact = TRUE)
  attr(out, "stream") <- key_stream(x)
  out
}

# Matrix subsetting turns an index it cannot use -- an infinity, or a
# double beyond the integer range -- into NA, and a row of NA words is a
# valid opaque key: the same one from every vector. Refuse such indices
# here, before subsetting; the key words themselves cannot be checked,
# because any 32-bit pattern, NA_INTEGER's included, is a legitimate word.
check_key_index <- function(i, n) {
  if (anyNA(i)) {
    stop("`i` must not contain missing values", call. = FALSE)
  }
  if (is.numeric(i)) {
    if (!all(is.finite(i))) {
      stop("`i` must contain finite positions", call. = FALSE)
    }
    if (any(i != trunc(i))) {
      stop("`i` must contain whole numbers", call. = FALSE)
    }
    if (any(abs(i) > n)) {
      stop("`i` must select existing keys, positions 1 to ", n, call. = FALSE)
    }
  } else if (!is.logical(i)) {
    stop("`i` must be an integer, double or logical vector", call. = FALSE)
  }
}

#' @rdname rng_key
#' @param value Not used: keys are immutable, so assigning into a key vector
#'   is an error. Without these methods R would write the value into
#'   individual key words, silently producing different keys.
#' @export
`[<-.rng_key` <- function(x, i, ..., value) {
  stop_key_assignment()
}

#' @rdname rng_key
#' @export
`[[<-.rng_key` <- function(x, i, ..., value) {
  stop_key_assignment()
}

stop_key_assignment <- function() {
  stop("keys are immutable; build a new key vector with c() and `[` instead",
       call. = FALSE)
}

#' @rdname rng_key
#' @export
as.list.rng_key <- function(x, ...) {
  lapply(seq_len(length(x)), function(j) x[j])
}

#' @rdname rng_key
#' @export
rep.rng_key <- function(x, ...) {
  x[rep(seq_len(length(x)), ...)]
}

#' @rdname rng_key
#' @param incomparables Passed to [duplicated()].
#' @export
duplicated.rng_key <- function(x, incomparables = FALSE, ...) {
  # duplicated() on the word matrix compares rows, i.e. whole keys; it
  # returns a one-dimensional array, which as.vector() makes a plain logical.
  as.vector(duplicated(unclass(x), incomparables = incomparables, ...))
}

#' @rdname rng_key
#' @export
unique.rng_key <- function(x, incomparables = FALSE, ...) {
  x[!duplicated(x, incomparables = incomparables, ...)]
}

#' @rdname rng_key
#' @export
xtfrm.rng_key <- function(x) {
  stop("keys have no order; compare them with identical()", call. = FALSE)
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
  if (length(dots) == 0L) {
    return(rng_key(0, n = 0L))
  }

  ok <- vapply(dots, inherits, logical(1L), "rng_key")
  if (!all(ok)) {
    stop("all inputs to `c()` must be <rng_key> objects", call. = FALSE)
  }

  engines <- vapply(dots, attr, character(1L), "engine", exact = TRUE)
  if (length(unique(engines)) != 1L) {
    stop("all keys must use the same engine", call. = FALSE)
  }
  streams <- vapply(dots, key_stream, integer(1L))
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
