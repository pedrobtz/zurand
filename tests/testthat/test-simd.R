# Dispatching on the instruction set is only acceptable if it cannot be
# observed in the output. That is not a property to assert in a comment --
# these run both paths and compare.

skip_if_no_simd <- function() {
  if (identical(rng_simd(), "none"))
    skip("no vectorised path on this CPU; both branches would be the scalar one")
}

with_simd <- function(on, code) {
  old <- rng_simd()
  on.exit(rng_simd(identical(old, "avx2")), add = TRUE)
  rng_simd(on)
  force(code)
}

test_that("the vectorised and portable paths produce identical values", {
  skip_if_no_simd()
  key <- rng_key(42L, engine = "xoshiro256pp")

  # Lengths chosen around the grouping: a group is 4 sub-chunks of 512, so
  # 2048 is one full group and the others leave a partial group, which is
  # the boundary the vectorised path declines and hands to the scalar one.
  for (n in c(1L, 7L, 511L, 512L, 513L, 2047L, 2048L, 2049L, 5000L)) {
    fast <- with_simd(TRUE,  rng_uniform(key, n))
    slow <- with_simd(FALSE, rng_uniform(key, n))
    expect_identical(fast, slow, info = paste("uniform, n =", n))
  }
  # Many full groups through the fused convert-and-store path, plus the
  # scaling pass that runs on its output.
  expect_identical(with_simd(TRUE,  rng_uniform(key, 100003L)),
                   with_simd(FALSE, rng_uniform(key, 100003L)))
  expect_identical(with_simd(TRUE,  rng_uniform(key, 5000L, -2, 7)),
                   with_simd(FALSE, rng_uniform(key, 5000L, -2, 7)))
  for (n in c(1L, 513L, 2048L, 5000L)) {
    expect_identical(with_simd(TRUE,  rng_normal(key, n)),
                     with_simd(FALSE, rng_normal(key, n)),
                     info = paste("normal, n =", n))
    expect_identical(with_simd(TRUE,  rng_integer(key, n, 1L, 6L)),
                     with_simd(FALSE, rng_integer(key, n, 1L, 6L)),
                     info = paste("integer, n =", n))
    expect_identical(with_simd(TRUE,  rng_bits(key, n)),
                     with_simd(FALSE, rng_bits(key, n)),
                     info = paste("bits, n =", n))
  }
  expect_identical(with_simd(TRUE,  rng_bits(key, 600L, bits = 64L)),
                   with_simd(FALSE, rng_bits(key, 600L, bits = 64L)))
  expect_identical(with_simd(TRUE,  format(rng_fold(key, "x"))),
                   with_simd(FALSE, format(rng_fold(key, "x"))))
})

test_that("multi-key and threaded fills are path-independent too", {
  skip_if_no_simd()
  keys <- rng_key(42L, n = 3L, engine = "xoshiro256pp")
  expect_identical(with_simd(TRUE,  rng_normal(keys, 3000L)),
                   with_simd(FALSE, rng_normal(keys, 3000L)))
  # 3 x 20000 values crosses the 32,768 threshold, so two threads really
  # split the fill; two is the most CRAN allows (see setup.R).
  raw <- rng_threads(1L)
  on.exit(rng_threads(raw), add = TRUE)
  one <- with_simd(TRUE, rng_uniform(keys, 20000L))
  rng_threads(2L)
  expect_identical(with_simd(TRUE, rng_uniform(keys, 20000L)), one)
})

test_that("rng_simd() validates and round-trips", {
  expect_true(rng_simd() %in% c("avx2", "none"))
  old <- rng_simd()
  expect_identical(rng_simd(FALSE), old)      # returns the previous value
  expect_identical(rng_simd(), "none")
  rng_simd(identical(old, "avx2"))
  expect_identical(rng_simd(), old)
  expect_error(rng_simd(NA), "TRUE or FALSE")
})
