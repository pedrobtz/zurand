# Package-global settings, changed for the rest of the calling test only,
# so that tests pass in any order (withr::defer restores on exit, error
# included).

local_threads <- function(n, .env = parent.frame()) {
  old <- rng_threads(n)
  withr::defer(rng_threads(old), envir = .env)
  invisible(old)
}

local_simd <- function(on, .env = parent.frame()) {
  old <- rng_simd(on)
  withr::defer(rng_simd(!identical(old, "none")), envir = .env)
  invisible(old)
}

local_cache_mode <- function(mode, .env = parent.frame()) {
  old <- cache_mode(mode)
  withr::defer(cache_mode(old), envir = .env)
  invisible(old)
}
