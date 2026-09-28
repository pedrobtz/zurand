# Lazy vectors (rng_lazy_uniform, rng_lazy_normal): every way of reading one
# must give the materialised sampler's values, on every engine, whatever
# state the one-chunk cache is in.

engines <- c("xoshiro256pp", "philox4x64", "threefry4x64")
materialised <- function(x) zurand:::rng_lazy_materialised(x)
cache_mode <- function(m) zurand:::rng_lazy_cache(m)

test_that("lazy vectors equal the samplers, read every way", {
  set.seed(1)
  for (engine in engines) {
    k <- rng_key(3L, engine = engine)
    # Lengths off the engine chunks (512; 4096 or 5120 for xoshiro256pp).
    for (n in c(1, 7, 513, 5121, 20001)) {
      cases <- list(
        list(function() rng_lazy_uniform(k, n), rng_uniform(k, n)),
        list(function() rng_lazy_uniform(k, n, -2, 5), rng_uniform(k, n, -2, 5)),
        list(function() rng_lazy_normal(k, n), rng_normal(k, n)),
        list(function() rng_lazy_normal(k, n, 1, 3, "mcfarland"),
             rng_normal(k, n, 1, 3, "mcfarland")),
        list(function() rng_lazy_normal(k, n, 4, 0), rng_normal(k, n, 4, 0)))
      for (cs in cases) {
        ref <- cs[[2]]
        idx <- c(n, 1, sample(n, min(n, 200), replace = TRUE),
                 seq_len(min(n, 40)) + max(0, n - 60))
        x <- cs[[1]]()
        expect_identical(length(x), as.integer(n))
        expect_identical(x[idx], ref[idx])
        expect_identical(vapply(idx, function(i) x[[i]], 0), ref[idx])
        expect_identical(x[seq_len(n)], ref)
        expect_identical(rev(x), rev(ref))
        expect_identical(sum(x), sum(ref))
        expect_identical(x[c(NA, n + 1, 1)], ref[c(NA, n + 1, 1)])
        expect_false(materialised(x))
        expect_identical(x + 0, ref + 0)     # materialises
        expect_true(materialised(x))
        expect_identical(x[idx], ref[idx])
      }
    }
  }
})

test_that("reads are right whatever the cache holds", {
  n <- 30000
  for (engine in engines) {
    k <- rng_key(9L, engine = engine)
    ref <- rng_normal(k, n, 2, 0.5)
    probe <- c(1, 511, 512, 513, 4096, 4097, 5120, 5121, 10240, 10241, 29999, n)
    region <- 5000:15500
    for (mode in 0:2) {
      old <- cache_mode(mode)
      prime <- list(
        empty = function(x) NULL,
        walked = function(x) for (i in 1:6000) x[[i]],
        scattered = function(x) { x[[17000]]; x[[3]]; x[[25000]] })
      for (p in names(prime)) {
        x <- rng_lazy_normal(k, n, 2, 0.5)
        prime[[p]](x)
        expect_identical(vapply(probe, function(i) x[[i]], 0), ref[probe],
                         info = paste(engine, mode, p))
        prime[[p]](x)
        expect_identical(x[probe], ref[probe], info = paste(engine, mode, p))
        prime[[p]](x)
        expect_identical(x[region], ref[region], info = paste(engine, mode, p))
        expect_false(materialised(x))
      }
      cache_mode(old)
    }
  }
})

test_that("long lazy vectors take double indices past 2^31", {
  skip_if(.Machine$sizeof.pointer < 8, "no long vectors")
  for (engine in c("xoshiro256pp", "philox4x64")) {
    k <- rng_key(5L, engine = engine)
    x <- rng_lazy_normal(k, 5e9)
    expect_identical(length(x), 5e9)
    at <- c(2^31 - 1, 2^31, 2^31 + 1, 2^32 + 7, 5e9)
    expect_identical(x[at], vapply(at, function(p) rng_normal(k, 1L, offset = p - 1), 0))
    # A run across 2^31, read as a region.
    across <- (2^31 - 40):(2^31 + 40)
    expect_identical(x[across], rng_normal(k, 81L, offset = 2^31 - 41))
    expect_false(materialised(x))
  }
})

test_that("lazy vectors serialise as their recipe, and modified ones as values", {
  k <- rng_key(8L)
  x <- rng_lazy_normal(k, 1e7, method = "mcfarland")
  bytes <- serialize(x, NULL)
  expect_lt(length(bytes), 1000)
  y <- unserialize(bytes)
  expect_false(materialised(y))
  expect_identical(y[c(1, 5e6, 1e7)], x[c(1, 5e6, 1e7)])

  z <- rng_lazy_uniform(k, 100)
  z[5] <- 0                                  # writes through the data pointer
  expect_identical(z[5], 0)
  w <- unserialize(serialize(z, NULL))
  expect_identical(w[5], 0)
  expect_identical(w[-5], rng_uniform(k, 100)[-5])
})

test_that("lazy uniforms with bounds more than DBL_MAX apart match the sampler", {
  big <- .Machine$double.xmax
  for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    k <- rng_key(3L, engine = engine)
    for (b in list(c(-1e308, 1e308), c(-big, big))) {
      ref <- rng_uniform(k, 70000L, b[1], b[2])
      expect_true(all(is.finite(ref)))
      x <- rng_lazy_uniform(k, 70000L, b[1], b[2])
      expect_identical(x[c(1, 2, 3, 600, 69999)], ref[c(1, 2, 3, 600, 69999)])
      expect_identical(vapply(1:40, function(i) x[[i]], 0), ref[1:40])
      expect_identical(x[5:40000], ref[5:40000])          # regions
      expect_identical(sum(x), sum(ref))
      y <- unserialize(serialize(x, NULL))
      expect_identical(y[c(7, 65000)], ref[c(7, 65000)])
      expect_identical(x[], ref)                          # materialised
    }
  }
})

test_that("materialisation happens once and copies stay independent", {
  k <- rng_key(4L)
  x <- rng_lazy_normal(k, 1000)
  y <- x
  y[1] <- 99                                 # copies y; x stays lazy
  expect_identical(y[1], 99)
  expect_identical(x[1], rng_normal(k, 1)[1])
  expect_false(materialised(x))
  x[2] <- -1                                 # materialises x itself
  expect_true(materialised(x))
  x[3] <- -2                                 # the same buffer, not a refill
  expect_identical(x[2:3], c(-1, -2))
  expect_identical(x[-(2:3)], rng_normal(k, 1000)[-(2:3)])
})

test_that("a written lazy vector no longer promises to have no NA", {
  x <- rng_lazy_normal(rng_key(1L), 1000L)
  expect_false(anyNA(x))
  x[5] <- NA
  expect_true(anyNA(x))
  expect_true(is.na(x[5]))
  y <- rng_lazy_uniform(rng_key(2L), 100L)
  y[3] <- NA
  expect_true(is.na(mean(y)))
  expect_identical(mean(y, na.rm = TRUE), mean(rng_uniform(rng_key(2L), 100L)[-3]))
})

test_that("real indices are range-checked before conversion", {
  x <- rng_lazy_normal(rng_key(1L), 10L)
  ref <- rng_normal(rng_key(1L), 10L)
  idx <- c(0.5, 1.5, 10.9, 11, 1e300, Inf, NA)
  expect_identical(x[idx], ref[idx])
  expect_true(all(is.na(x[c(11, 1e300, Inf)])))
})

test_that("large lazy sums match with threads on and off", {
  skip_on_cran()
  k <- rng_key(6L)
  ref <- sum(rng_uniform(k, 1e7, -1, 3))
  old <- rng_threads(1L)
  one <- sum(rng_lazy_uniform(k, 1e7, -1, 3))
  rng_threads(0L)
  all <- sum(rng_lazy_uniform(k, 1e7, -1, 3))
  rng_threads(old)
  expect_identical(one, ref)
  expect_identical(all, ref)
})

test_that("lazy vectors check their arguments", {
  expect_error(rng_lazy_normal(rng_key(1L, n = 2L), 10), "single key")
  expect_error(rng_lazy_normal(rng_key(1L), -1), "non-negative")
  expect_error(rng_lazy_normal(rng_key(1L), 10, sd = -1), "non-negative")
  expect_error(rng_lazy_uniform(rng_key(1L), 10, 2, 1), "less than or equal")
  expect_identical(length(rng_lazy_uniform(rng_key(1L), 0)), 0L)
})
