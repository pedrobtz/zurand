# Helpers for test-simd.R (the SIMD identity tests), kept out of it so that the
# tests run shuffled.

skip_if_no_simd <- function() {
  if (identical(rng_simd(), "none"))
    skip("no vectorised path on this CPU; both branches would be the scalar one")
}

with_simd <- function(on, code) {
  old <- rng_simd()
  on.exit(rng_simd(!identical(old, "none")), add = TRUE)
  rng_simd(on)
  force(code)
}
