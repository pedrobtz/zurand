# Helpers for test-lazy.R, kept out of it so that the tests run shuffled.

engines <- c("xoshiro256pp", "philox4x64", "threefry4x64")

materialised <- function(x) zurand:::rng_lazy_materialised(x)

cache_mode <- function(mode = NULL) zurand:::rng_lazy_cache(mode)

# serialize(x) with its saved arguments replaced by `meta`, as a damaged
# file or a newer zurand would write them. `x` must be a lazy vector with
# the default scaling (mean 0, sd 1, or min 0, max 1), whose arguments are
# then c(length(x), dist, 0, 1): with their XDR length prefix, a byte
# pattern that occurs once in the stream.
tamper_lazy <- function(x, dist, meta) {
  xdr <- function(v) {
    c(writeBin(length(v), raw(), endian = "big"),
      writeBin(as.double(v), raw(), endian = "big"))
  }
  bytes <- serialize(x, NULL, xdr = TRUE)
  pat <- xdr(c(length(x), dist, 0, 1))
  starts <- seq_len(length(bytes) - length(pat) + 1L)
  at <- starts[vapply(
    starts,
    function(i) identical(bytes[i:(i + length(pat) - 1L)], pat),
    logical(1L)
  )]
  stopifnot(length(at) == 1L)
  c(bytes[seq_len(at - 1L)], xdr(meta), bytes[-seq_len(at + length(pat) - 1L)])
}
