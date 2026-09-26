test_that("keys record stream version 1", {
  expect_identical(attr(rng_key(1L), "stream"), 1L)
  expect_identical(attr(rng_key(1L, n = 3L, engine = "xoshiro256pp"), "stream"), 1L)
  expect_identical(attr(rng_key_from_r(2L), "stream"), 1L)
  expect_identical(attr(rng_fold(rng_key(1L), "a"), "stream"), 1L)
})

test_that("subsetting and c() keep the stream version", {
  keys <- rng_key(7L, n = 4L)
  expect_identical(attr(keys[2:3], "stream"), 1L)
  expect_identical(attr(keys[[1]], "stream"), 1L)
  expect_identical(attr(c(keys[1], keys[4]), "stream"), 1L)
})

test_that("a key saved before the stream attribute existed is stream 1", {
  key <- rng_key(42L)
  old <- key
  attr(old, "stream") <- NULL
  expect_identical(rng_uniform(old, 5L), rng_uniform(key, 5L))
  expect_identical(rng_normal(old, 5L), rng_normal(key, 5L))
  expect_identical(rng_integer(old, 5L, 1L, 6L), rng_integer(key, 5L, 1L, 6L))
  expect_identical(rng_bits(old, 5L), rng_bits(key, 5L))
  expect_identical(unclass(rng_fold(old, "x")), unclass(rng_fold(key, "x")))
  expect_identical(attr(old[1], "stream"), 1L)
  expect_identical(attr(c(old, key), "stream"), 1L)
})

test_that("a key from a stream this version does not implement is an error", {
  key <- rng_key(42L)
  attr(key, "stream") <- 2L
  expect_error(rng_uniform(key), "stream version 2")
  expect_error(rng_normal(key), "stream version 2")
  expect_error(rng_integer(key, 1L, 1L, 6L), "stream version 2")
  expect_error(rng_bits(key), "stream version 2")
  expect_error(rng_fold(key, "x"), "stream version 2")

  bad <- rng_key(42L)
  attr(bad, "stream") <- "one"
  expect_error(rng_uniform(bad), "invalid `stream`")
  attr(bad, "stream") <- 0L
  expect_error(rng_uniform(bad), "invalid `stream`")
})

test_that("c() refuses to mix stream versions", {
  a <- rng_key(1L)
  b <- rng_key(2L)
  attr(b, "stream") <- 2L
  expect_error(c(a, b), "same stream version")
})
