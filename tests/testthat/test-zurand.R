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
  expect_identical(length(c(keys[0], keys[0])), 0L)
  expect_error(c(keys[1], "not a key"), "rng_key")
})

test_that("key subsetting rejects indices that would fabricate keys", {
  keys <- rng_key(42L, n = 4L)

  expect_error(keys[NA], "missing")
  expect_error(keys[NA_integer_], "missing")
  expect_error(keys[c(1L, NA, 3L)], "missing")
  expect_error(keys[1.5], "whole numbers")
  expect_error(keys[10L], "subscript out of bounds")

  expect_identical(keys[[2L]], keys[2L])
  expect_error(keys[[c(1L, 2L)]], "exactly one key")
  expect_error(keys[[TRUE]], "exactly one key")
  expect_error(keys[[integer()]], "exactly one key")
})

test_that("rng_fold() derives deterministic keys from typed data", {
  key <- rng_key(42L, n = 3L)

  expect_identical(rng_fold(key, "layer1"), rng_fold(key, "layer1"))
  expect_length(rng_fold(key, "layer1"), 3L)
  expect_false(identical(rng_fold(key, "layer1"), rng_fold(key, "layer2")))
  expect_false(identical(rng_fold(key, 1L), rng_fold(key, "1")))
  expect_false(identical(rng_fold(key, c("a", "b")), rng_fold(key, c("ab"))))
})

test_that("rng_fold() folds whole-number doubles and integers identically", {
  key <- rng_key(3L)

  expect_identical(rng_fold(key, 1), rng_fold(key, 1L))
  expect_identical(rng_fold(key, c(-2, 0, 7)), rng_fold(key, c(-2L, 0L, 7L)))
  expect_identical(rng_fold(key, 2^40), rng_fold(key, 2^40))
  expect_false(identical(rng_fold(key, 1), rng_fold(key, TRUE)))
  expect_false(identical(rng_fold(key, 1), rng_fold(key, as.raw(1))))
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
  zm <- rng_normal(keys, 5L)

  expect_type(z, "double")
  expect_length(z, 10000L)
  expect_true(all(is.finite(z)))
  expect_gt(mean(z), 1.8)
  expect_lt(mean(z), 2.2)
  expect_gt(sd(z), 2.8)
  expect_lt(sd(z), 3.2)
  expect_identical(rng_normal(key, 3L, mean = 7, sd = 0), rep(7, 3L))
  expect_equal(dim(zm), c(5L, 2L))
  expect_identical(zm[, 2L], rng_normal(keys[2], 5L))
})

test_that("scaling is the unfused mean + sd * z, for odd and even lengths", {
  # R evaluates `a + s * x` as two separately rounded vector operations, so
  # this pins the C scaling pass to the unfused definition. The pass works
  # two doubles at a time; odd lengths exercise its trailing element, and
  # multi-key output its full n x K extent.
  key <- rng_key(99L)
  keys <- rng_key(99L, n = 3L)
  for (n in c(1L, 2L, 7L, 1000L, 1001L)) {
    expect_identical(rng_normal(key, n, mean = 0.3, sd = 1.7),
                     0.3 + 1.7 * rng_normal(key, n))
    expect_identical(rng_uniform(key, n, min = -0.3, max = 2.9),
                     -0.3 + (2.9 - -0.3) * rng_uniform(key, n))
  }
  expect_identical(rng_normal(keys, 5L, mean = 0.3, sd = 1.7),
                   0.3 + 1.7 * rng_normal(keys, 5L))
})

test_that("rng_normal() draws from the normal distribution", {
  # deterministic given the key, so these are regression tests, not flaky
  # statistical ones; n exceeds the OpenMP threshold
  z <- rng_normal(rng_key(2024L), 50000L)

  ks <- suppressWarnings(stats::ks.test(z, "pnorm"))
  expect_gt(ks$p.value, 0.001)
  # the ziggurat fast path is capped at |z| < 3.6542; seeing larger values
  # proves the tail branch runs
  expect_gt(max(abs(z)), 3.6542)
  expect_lt(min(z), -3.6542)
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

test_that("large draws extend small draws bit-for-bit, including parallel paths", {
  # 70001 exceeds the OpenMP threshold (32768) and has an odd tail, so this
  # exercises the parallel block loops; the prefix comparisons pin the core
  # counter-mode invariant that draw i never depends on n.
  key <- rng_key(99L)
  keys <- rng_key(99L, n = 3L)
  n_big <- 70001L

  u <- rng_uniform(key, n_big)
  expect_identical(u[1:257], rng_uniform(key, 257L))
  z <- rng_normal(key, n_big)
  expect_identical(z[1:257], rng_normal(key, 257L))
  x <- rng_integer(key, n_big, -5L, 5L)
  expect_identical(x[1:257], rng_integer(key, 257L, -5L, 5L))
  b <- rng_bits(key, n_big)
  expect_identical(b[1:257], rng_bits(key, 257L))
  expect_identical(rng_bits(key, 9L, bits = 64L)[1:5], rng_bits(key, 5L, bits = 64L))

  # near-maximal range keeps the rejection path honest
  r <- rng_integer(key, n_big, 1L, 2000000000L)
  expect_identical(r[1:257], rng_integer(key, 257L, 1L, 2000000000L))

  # multi-key: columns of the parallel-over-columns fill match per-key draws
  um <- rng_uniform(keys, n_big)
  expect_identical(um[, 2L], rng_uniform(keys[2], n_big))
  expect_identical(um[1:9, 3L], rng_uniform(keys[3], 9L))
  zm <- rng_normal(keys, n_big)
  expect_identical(zm[, 1L], rng_normal(keys[1], n_big))
  xm <- rng_integer(keys, n_big, 0L, 9L)
  expect_identical(xm[, 3L], rng_integer(keys[3], n_big, 0L, 9L))
})

test_that("rng_threads() caps parallelism without changing results", {
  eff <- rng_threads()
  expect_gte(eff, 1L)

  key <- rng_key(31L)
  keys <- rng_key(31L, n = 3L)
  u <- rng_uniform(key, 70001L)
  z <- rng_normal(keys, 40001L)
  x <- rng_integer(key, 70001L, -9L, 9L)

  raw <- rng_threads(1L)
  on.exit(rng_threads(raw))
  expect_gte(raw, 0L)
  expect_identical(rng_threads(), 1L)
  expect_identical(rng_uniform(key, 70001L), u)
  expect_identical(rng_normal(keys, 40001L), z)
  expect_identical(rng_integer(key, 70001L, -9L, 9L), x)

  # setting returns the raw previous cap so restores round-trip exactly,
  # including back to the uncapped state
  prev <- rng_threads(0L)
  expect_identical(prev, 1L)
  expect_identical(rng_threads(), eff)
  expect_identical(rng_uniform(key, 70001L), u)

  # far above the machine count is allowed and equivalent to no cap
  rng_threads(1e9)
  expect_identical(rng_threads(), eff)
  rng_threads(0L)

  expect_error(rng_threads(-1L), "at least 1")
  expect_error(rng_threads(1.5), "at least 1")
  expect_error(rng_threads(c(1L, 2L)), "single value")
  expect_identical(rng_threads(), eff)
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
  expect_error(rng_key(1L, n = c(1L, 2L)), "single value")
  expect_error(rng_key(1L, engine = "mersenne"), "should be one of")

  expect_error(rng_uniform(key, c(1L, 2L)), "single value")
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
