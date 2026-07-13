test_that("rng_key() is deterministic and opaque", {
  key <- rng_key(42L)
  keys <- rng_key(42L, n = 4L)

  expect_s3_class(key, "rng_key")
  expect_type(unclass(key), "integer")
  expect_equal(dim(unclass(key)), c(1L, 4L))
  expect_equal(dim(unclass(keys)), c(4L, 4L))
  expect_identical(length(key), 1L)
  expect_identical(length(keys), 4L)
  expect_identical(attr(key, "engine"), "philox4x64")
  expect_identical(key, rng_key(42L))
  expect_identical(keys[1L], key)
  expect_s3_class(keys[[2L]], "rng_key")
  expect_identical(length(keys[[2L]]), 1L)
  expect_false(identical(rng_key(1L), rng_key(2L)))

  expect_match(format(key), "^rng_key\\[[0-9a-f]{32}\\]$")
  expect_length(format(keys), 4L)
  expect_output(print(key), "<rng_key> rng_key\\[[0-9a-f]{32}\\]")
  expect_output(print(keys), "<rng_key\\[4\\]>")
})

test_that("rng_key vectors concatenate", {
  keys <- rng_key(42L, n = 4L)

  expect_identical(c(keys[1], keys[2:4]), keys)
  expect_identical(length(c(keys[1], keys[2])), 2L)
  expect_error(c(keys[1], "not a key"), "rng_key")
})

test_that("rng_fold() derives deterministic keys from typed data", {
  key <- rng_key(42L, n = 3L)

  expect_identical(rng_fold(key, "layer1"), rng_fold(key, "layer1"))
  expect_length(rng_fold(key, "layer1"), 3L)
  expect_false(identical(rng_fold(key, "layer1"), rng_fold(key, "layer2")))
  expect_false(identical(rng_fold(key, 1L), rng_fold(key, "1")))
  expect_false(identical(rng_fold(key, c("a", "b")), rng_fold(key, c("ab"))))
})

test_that("samplers are pure functions of key and arguments", {
  key <- rng_key(42L)
  keys <- rng_key(42L, n = 3L)

  expect_identical(rng_uniform(key, 10L), rng_uniform(key, 10L))
  expect_identical(rng_normal(key, 10L), rng_normal(key, 10L))
  expect_identical(rng_integer(key, 10L, 1L, 6L), rng_integer(key, 10L, 1L, 6L))
  expect_identical(rng_bits(key, 10L), rng_bits(key, 10L))

  expect_identical(rng_uniform(keys, 10L), rng_uniform(keys, 10L))
  expect_false(identical(rng_uniform(keys[1], 10L), rng_uniform(keys[2], 10L)))
})

test_that("samplers default to one draw per key", {
  key <- rng_key(42L)
  keys <- rng_key(42L, n = 3L)

  expect_identical(rng_uniform(key), rng_uniform(key, 1L))
  expect_identical(rng_normal(key), rng_normal(key, 1L))
  expect_identical(rng_integer(key, min = 1L, max = 6L), rng_integer(key, 1L, 1L, 6L))
  expect_identical(rng_bits(key), rng_bits(key, 1L))

  u <- rng_uniform(keys)
  z <- rng_normal(keys)
  i <- rng_integer(keys, min = 1L, max = 6L)
  b <- rng_bits(keys)
  b64 <- rng_bits(keys, bits = 64L)

  expect_equal(dim(u), c(1L, 3L))
  expect_equal(dim(z), c(1L, 3L))
  expect_equal(dim(i), c(1L, 3L))
  expect_equal(dim(b), c(1L, 3L))
  expect_equal(dim(b64), c(1L, 3L))
  expect_identical(u[, 2L], rng_uniform(keys[2]))
  expect_identical(z[, 2L], rng_normal(keys[2]))
  expect_identical(i[, 2L], rng_integer(keys[2], min = 1L, max = 6L))
  expect_identical(b[, 2L], rng_bits(keys[2]))
})

test_that("rng_uniform() returns values in the requested interval", {
  key <- rng_key(123L)
  keys <- rng_key(123L, n = 3L)
  u <- rng_uniform(key, 10000L, min = -2, max = 5)
  um <- rng_uniform(keys, 5L)

  expect_type(u, "double")
  expect_length(u, 10000L)
  expect_true(all(u > -2 & u < 5))
  expect_gt(mean(u), 1.4)
  expect_lt(mean(u), 1.6)
  expect_identical(rng_uniform(key, 3L, min = 2, max = 2), rep(2, 3L))
  expect_equal(dim(um), c(5L, 3L))
  expect_identical(um[, 1L], rng_uniform(keys[1], 5L))
})

test_that("rng_normal() returns shifted and scaled normal values", {
  key <- rng_key(123L)
  keys <- rng_key(123L, n = 2L)
  z <- rng_normal(key, 10000L, mean = 2, sd = 3)
  zm <- rng_normal(keys, 4L)

  expect_type(z, "double")
  expect_length(z, 10000L)
  expect_true(all(is.finite(z)))
  expect_gt(mean(z), 1.8)
  expect_lt(mean(z), 2.2)
  expect_identical(rng_normal(key, 3L, mean = 7, sd = 0), rep(7, 3L))
  expect_equal(dim(zm), c(4L, 2L))
  expect_identical(zm[, 2L], rng_normal(keys[2], 4L))
})

test_that("rng_integer() returns integers in the inclusive range", {
  key <- rng_key(7L)
  keys <- rng_key(7L, n = 2L)
  x <- rng_integer(key, 1000L, min = 1L, max = 6L)
  xm <- rng_integer(keys, 5L, min = 1L, max = 6L)

  expect_type(x, "integer")
  expect_length(x, 1000L)
  expect_true(all(x >= 1L & x <= 6L))
  expect_identical(rng_integer(key, 4L, min = -3L, max = -3L), rep(-3L, 4L))
  expect_equal(dim(xm), c(5L, 2L))
  expect_identical(xm[, 1L], rng_integer(keys[1], 5L, min = 1L, max = 6L))
})

test_that("rng_bits() returns exact 32-bit words or 64-bit hex words", {
  key <- rng_key(7L)
  keys <- rng_key(7L, n = 2L)
  b32 <- rng_bits(key, 1000L)
  b64 <- rng_bits(key, 5L, bits = 64L)
  b32m <- rng_bits(keys, 5L)
  b64m <- rng_bits(keys, 5L, bits = 64L)

  expect_type(b32, "double")
  expect_true(all(b32 == trunc(b32)))
  expect_true(all(b32 >= 0 & b32 < 2^32))
  expect_type(b64, "character")
  expect_match(b64, "^[0-9a-f]{16}$")
  expect_equal(dim(b32m), c(5L, 2L))
  expect_equal(dim(b64m), c(5L, 2L))
  expect_identical(b32m[, 2L], rng_bits(keys[2], 5L))
  expect_match(as.vector(b64m), "^[0-9a-f]{16}$")
})

test_that("stateless functions do not touch .Random.seed", {
  set.seed(123)
  before <- .Random.seed
  key <- rng_key(42L)

  rng_fold(key, "layer")
  rng_uniform(key, 10L)
  rng_normal(key, 10L)
  rng_integer(key, 10L, 1L, 6L)
  rng_bits(key, 10L)

  expect_identical(.Random.seed, before)
})

test_that("rng_key_from_r() consumes R's RNG state intentionally", {
  set.seed(123)
  key <- rng_key_from_r(n = 3L)
  after <- .Random.seed

  set.seed(123)
  expect_identical(rng_key_from_r(), key[1])
  set.seed(123)
  expect_identical(rng_key_from_r(n = 3L), key)
  expect_identical(.Random.seed, after)
})

test_that("invalid inputs error cleanly", {
  key <- rng_key(1L)

  expect_error(rng_uniform("not a key", 1L), "rng_key")
  expect_error(rng_key(-1L), "between 0")
  expect_error(rng_key(1.5), "whole number")
  expect_error(rng_key(2^53), "2\\^53")
  expect_error(rng_key(NA_integer_), "missing")
  expect_error(rng_key(1L, n = c(1L, 2L)), "single")
  expect_error(rng_key(1L, engine = "threefry"), "philox4x64")

  expect_error(rng_uniform(key, c(1L, 2L)), "single")
  expect_error(rng_uniform(key, "2"), "integer or double")
  expect_error(rng_uniform(key, -1L), "non-negative")
  expect_error(rng_uniform(key, 1.5), "non-negative whole")
  expect_error(rng_uniform(key, 1L, min = 2, max = 1), "less than or equal")

  expect_error(rng_normal(key, 1L, sd = -1), "non-negative")
  expect_error(rng_integer(key, 1L), "required")
  expect_error(rng_integer(key, 1L, min = 6L, max = 1L), "less than or equal")
  expect_error(rng_integer(key, 1L, min = NA_integer_, max = 6L), "missing")
  expect_error(rng_bits(key, 1L, bits = 16L), "32 or 64")

  expect_error(rng_fold(key, NA_character_), "missing")
  expect_error(rng_fold(key, 1.5), "whole numbers")
  expect_error(rng_fold(key, list("x")), "character, integer, double")
})
