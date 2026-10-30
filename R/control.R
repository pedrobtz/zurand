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
#' @section Page size on Linux:
#' R serves a large vector from fresh memory, so filling it is the first
#' touch of every page and pays one page fault per page. When the Linux
#' kernel's transparent huge pages are in `madvise` mode (the Ubuntu
#' default), the samplers ask for 2 MiB pages on output vectors of 4 MiB and
#' more before filling them, which about halves the time of large fills on
#' such systems; on kernels set to `always` the pages are large anyway, and
#' the request is not made. No value changes either way.
#' `options(zurand.hugepages = FALSE)` turns the request off.
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
#' key <- rng_key(1L, engine = "xoshiro256pp")
#' old <- rng_simd(FALSE)
#' portable <- rng_uniform(key, 4L)
#' rng_simd(old != "none")
#' identical(portable, rng_uniform(key, 4L))
rng_simd <- function(enable = NULL) {
  prev <- .Call(C_rng_simd, enable)
  if (is.null(enable)) prev else invisible(prev)
}
