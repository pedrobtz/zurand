test_that("rng_key() is deterministic and well-formed", {
  key <- rng_key(42)
  expect_s3_class(key, "rngat_key")
  expect_type(unclass(key), "raw")
  expect_length(unclass(key), 16)
  expect_identical(key, rng_key(42))
  expect_false(identical(rng_key(1), rng_key(2)))
})

test_that("draws are pure functions of (key, index, domain)", {
  key <- rng_key(42)
  expect_identical(runif_at(key, 1:10), runif_at(key, 1:10))
  # random access: single indices match the vectorized draw, in any order
  expect_identical(runif_at(key, 10:1), rev(runif_at(key, 1:10)))
  expect_identical(runif_at(key, 7), runif_at(key, 1:10)[7])
  # different keys, indices and domains give different values
  expect_false(runif_at(key, 1) == runif_at(key, 2))
  expect_false(runif_at(key, 1) == runif_at(rng_key(43), 1))
  expect_false(runif_at(key, 1) == runif_at(key, 1, domain = 1))
})

test_that("domain recycles against index", {
  key <- rng_key(1)
  expect_identical(
    runif_at(key, c(1, 1, 2), domain = c(0, 1, 1)),
    c(runif_at(key, 1), runif_at(key, 1, domain = 1), runif_at(key, 2, domain = 1))
  )
  # scalar index, vector domain
  expect_length(runif_at(key, 1, domain = 1:5), 5)
  expect_error(runif_at(key, 1:3, domain = 1:2), "length")
})

test_that("runif_at() lies strictly inside (0, 1)", {
  u <- runif_at(rng_key(123), 1:1e4)
  expect_true(all(u > 0 & u < 1))
  # crude uniformity sanity check
  expect_gt(mean(u), 0.45)
  expect_lt(mean(u), 0.55)
})

test_that("rnorm_at() is the inverse-CDF transform of runif_at()", {
  key <- rng_key(7)
  expect_identical(rnorm_at(key, 1:100), qnorm(runif_at(key, 1:100)))
})

test_that("bits_at() returns 32-bit unsigned integer values", {
  b <- bits_at(rng_key(7), 1:1e3)
  expect_true(all(b == trunc(b)))
  expect_true(all(b >= 0 & b < 2^32))
})

test_that("fold_in() derives distinct deterministic keys", {
  key <- rng_key(42)
  k1 <- fold_in(key, 1)
  expect_s3_class(k1, "rngat_key")
  expect_identical(k1, fold_in(key, 1))
  expect_false(identical(k1, key))
  expect_false(identical(k1, fold_in(key, 2)))
  # derived streams differ from the parent's
  expect_false(runif_at(k1, 1) == runif_at(key, 1))
})

test_that("indices beyond the integer range work", {
  key <- rng_key(1)
  expect_identical(runif_at(key, 2^40), runif_at(key, 2^40))
  expect_false(runif_at(key, 2^40) == runif_at(key, 2^40 + 1))
})

test_that("invalid inputs error cleanly", {
  key <- rng_key(1)
  expect_error(runif_at("not a key", 1), "rngat_key")
  expect_error(runif_at(structure(raw(16), class = "rngat_key"), 1), NA)
  expect_error(runif_at(key, 1.5), "whole number")
  expect_error(runif_at(key, NA_integer_), "missing")
  expect_error(rng_key(NULL), "integer or double")
  expect_error(rng_key(1:2), "single")
  expect_error(fold_in(key, "a"), "integer or double")
})

test_that("*_seq() is bit-identical to *_at() on the same run", {
  key <- rng_key(42)
  # sweep start offsets 0..5 and lengths 0..6 to cover every head/tail
  # alignment against the 4-word block boundary
  for (start in 0:5) {
    for (len in 0:6) {
      idx <- start + seq_len(len) - 1
      expect_identical(runif_seq(key, start, len), runif_at(key, idx))
      expect_identical(bits_seq(key, start, len), bits_at(key, idx))
    }
  }
  expect_identical(rnorm_seq(key, 3, 10), rnorm_at(key, 3:12))
  expect_identical(
    runif_seq(key, 5, 100, domain = 2),
    runif_at(key, 5:104, domain = 2)
  )
  # long run spanning many blocks, and a run far out in the index space
  expect_identical(runif_seq(key, 1, 1e4), runif_at(key, 1:1e4))
  expect_identical(runif_seq(key, 2^40, 8), runif_at(key, 2^40 + 0:7))
})

test_that("*_seq() validates inputs", {
  key <- rng_key(1)
  expect_identical(runif_seq(key, 1, 0), double(0))
  expect_error(runif_seq(key, 1, -1), "non-negative")
  expect_error(runif_seq(key, 1, 1.5), "non-negative whole")
  expect_error(runif_seq(key, 1.5, 10), "whole number")
  expect_error(runif_seq(key, 1, 10, domain = 1:2), "single")
  expect_error(runif_seq("not a key", 1, 10), "rngat_key")
})

test_that("keys print as hex", {
  out <- format(rng_key(42))
  expect_match(out, "^[0-9a-f]{32}$")
  expect_output(print(rng_key(42)), "<rngat_key> [0-9a-f]{32}")
})
