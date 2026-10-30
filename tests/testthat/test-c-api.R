# The C API (inst/include/zurand.h), exercised the way a downstream package
# uses it: a fixture package, tests/testthat/zurandclient, is built against
# the header, fetches the API table, and fills one key per column inside an
# OpenMP parallel loop. Everything it returns must equal the R samplers.
#
# Not on CRAN: it compiles a package, which takes seconds and needs a
# toolchain. Everywhere else a build failure is a test failure, not a skip.

test_that("the fixture is threaded wherever zurand is", {
  skip_on_cran()
  ns <- client()
  # zurand reports one thread when it was built without OpenMP; then the
  # fixture is serial too and there is nothing to compare (#60).
  skip_if(rng_threads() < 2L, "zurand was built without OpenMP")
  expect_gt(ns$client_threads(), 1L)
})

test_that("C API fills equal the R samplers, called from worker threads", {
  skip_on_cran()
  ns <- client()
  fill <- ns$client_fill

  for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    keys <- rng_key(20260927L, n = 6L, engine = engine)
    # 2600 values per key: several xoshiro chunks and a ragged tail.
    n <- 2600L
    expect_identical(fill(keys, n, "uniform"), rng_uniform(keys, n))
    expect_identical(fill(keys, n, "uniform", -3, 7.5),
                     rng_uniform(keys, n, min = -3, max = 7.5))
    expect_identical(fill(keys, n, "normal"), rng_normal(keys, n))
    expect_identical(fill(keys, n, "normal", 0.3, 1.7),
                     rng_normal(keys, n, mean = 0.3, sd = 1.7))
    expect_identical(fill(keys, n, "integer", -7, 1000),
                     rng_integer(keys, n, min = -7L, max = 1000L))
    expect_identical(fill(keys, n, "integer", -2e9, 2e9),
                     rng_integer(keys, n, min = -2000000000L, max = 2000000000L))
    expect_identical(fill(keys, n, "bits64"), rng_bits(keys, n, bits = 64L))
    # fill_normal_method(), API version 2.
    expect_identical(fill(keys, n, "normal_ziggurat"), rng_normal(keys, n))
    expect_identical(fill(keys, n, "normal_mcfarland"),
                     rng_normal(keys, n, method = "mcfarland"))
    expect_identical(fill(keys, n, "normal_mcfarland", 0.3, 1.7),
                     rng_normal(keys, n, mean = 0.3, sd = 1.7, method = "mcfarland"))
    # Bounds more than DBL_MAX apart, whose span overflows.
    big <- .Machine$double.xmax
    expect_identical(fill(keys, n, "uniform", -big, big),
                     rng_uniform(keys, n, -big, big))
    expect_identical(fill(keys, n, "uniform", -1e308, 1e308),
                     rng_uniform(keys, n, -1e308, 1e308))
    # Degenerate scales, handled before the fill.
    expect_identical(fill(keys, 5L, "uniform", 2, 2), rng_uniform(keys, 5L, 2, 2))
    expect_identical(fill(keys, 5L, "normal", 4, 0), rng_normal(keys, 5L, 4, 0))
    expect_identical(fill(keys, 5L, "normal_mcfarland", 4, 0),
                     rng_normal(keys, 5L, 4, 0, method = "mcfarland"))
  }
})

test_that("C API fold matches rng_fold() for whole numbers", {
  skip_on_cran()
  ns <- client()
  for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    key <- rng_key(5L, engine = engine)
    for (d in c(0, 1, 7, -3, 123456789)) {
      expect_identical(ns$client_fold(key, d), format(rng_fold(key, d)))
      expect_identical(ns$client_fold(key, d), format(rng_fold(key, as.integer(d))))
    }
  }
})

test_that("C API reports invalid arguments as return codes", {
  skip_on_cran()
  ns <- client()
  # version, then: key index out of range, min > max, sd < 0, infinite bound,
  # NA integer bound, unknown engine, unknown stream, n = 0 (fine); and for
  # fill_normal_method(): methods 2 and -1, sd < 0, unknown engine.
  # Version 3: positional fills (min > max, method, start + n overflowing)
  # and streams (chunk 0, no buffer, no consumer, engine, total 0 is fine).
  # Then a NULL output buffer for every fill and fold_int(), which used to
  # crash the session (#60), and a NULL buffer with n = 0, which is fine.
  # 0 ZURAND_OK, 1 ZURAND_EINVAL, 2 ZURAND_EKEY; "version" is the table's.
  want <- c(
    version = 3L, key_get_range = 1L, key_get_null_out = 1L,
    uniform_min_gt_max = 1L, uniform_inf = 1L, uniform_nan = 1L,
    uniform_too_long = 1L, normal_sd_neg = 1L, normal_mean_inf = 1L,
    normal_engine = 2L, normal_n0 = 0L, integer_na = 1L,
    integer_min_gt_max = 1L, bits64_stream = 2L, method_2 = 1L,
    method_neg = 1L, method_sd_neg = 1L, method_engine = 2L,
    uniform_at_min_gt_max = 1L, uniform_at_overflow = 1L,
    uniform_at_engine = 2L, normal_at_method = 1L, normal_at_overflow = 1L,
    normal_at_engine = 2L, stream_chunk_0 = 1L, stream_chunk_too_long = 1L,
    stream_no_buf = 1L, stream_no_fn = 1L, stream_uniform_engine = 2L,
    stream_normal_engine = 2L, stream_sd_neg = 1L, stream_method = 1L,
    stream_total_0 = 0L, fold_engine = 2L,
    null_uniform = 1L, null_uniform_constant = 1L, null_normal_method = 1L,
    null_integer = 1L, null_bits64 = 1L, null_normal_at = 1L,
    null_uniform_at = 1L, null_fold = 1L, null_uniform_n0 = 0L
  )
  expect_identical(ns$client_errors(rng_key(1L)), want)
  # key_get() is the one member that raises R errors, for a non-key.
  expect_error(ns$client_key_get(1:4), "must be an <rng_key>")
})

test_that("C API positional fills equal rng_*(offset =), from worker threads", {
  skip_on_cran()
  ns <- client()
  for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    keys <- rng_key(20260927L, n = 4L, engine = engine)
    # Offsets on and off the engine chunks (512; 4096 or 5120 for xoshiro).
    for (offset in c(0, 1, 511, 4097, 5121, 12345)) {
      expect_identical(ns$client_fill_at(keys, 3000L, "uniform", 0, 1, offset),
                       rng_uniform(keys, 3000L, offset = offset))
      expect_identical(ns$client_fill_at(keys, 3000L, "uniform", -2, 5, offset),
                       rng_uniform(keys, 3000L, -2, 5, offset = offset))
      expect_identical(
        ns$client_fill_at(keys, 3000L, "uniform", -1e308, 1e308, offset),
        rng_uniform(keys, 3000L, -1e308, 1e308, offset = offset))
      expect_identical(ns$client_fill_at(keys, 3000L, "normal", 0.5, 2, offset),
                       rng_normal(keys, 3000L, 0.5, 2, offset = offset))
      expect_identical(
        ns$client_fill_at(keys, 3000L, "normal", 0, 1, offset, "mcfarland"),
        rng_normal(keys, 3000L, method = "mcfarland", offset = offset))
    }
  }
})

test_that("C API streams deliver rng_*() in order, chunk by chunk", {
  skip_on_cran()
  ns <- client()
  total <- 23456L
  for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    key <- rng_key(7L, engine = engine)
    # Chunks below, at, between and above the engine chunks, and odd ones
    # (the scaling pass works in pairs).
    for (chunk in c(1L, 7L, 512L, 1000L, 4096L, 5120L, 5121L, 16384L, 30000L)) {
      s <- ns$client_stream(key, total, chunk, "normal", 0, 1)
      expect_identical(s[[1]], rng_normal(key, total))
      expect_identical(s[[4]], 0L)                     # ZURAND_OK
      expect_false(s[[5]])                             # in order
      n_full <- total %/% chunk
      expect_identical(s[[3]], as.double(c(rep(chunk, n_full),
                                           if (total %% chunk) total %% chunk)))
      expect_identical(s[[2]], as.double(chunk) * (seq_along(s[[3]]) - 1))
    }
    for (chunk in c(333L, 5121L)) {
      expect_identical(ns$client_stream(key, total, chunk, "normal", -1, 3,
                                        "mcfarland")[[1]],
                       rng_normal(key, total, -1, 3, method = "mcfarland"))
      expect_identical(ns$client_stream(key, total, chunk, "uniform", 0, 1)[[1]],
                       rng_uniform(key, total))
      expect_identical(ns$client_stream(key, total, chunk, "uniform", 10, 11)[[1]],
                       rng_uniform(key, total, 10, 11))
      # Degenerate scales stream the constant.
      expect_identical(ns$client_stream(key, 1001L, chunk, "normal", 4, 0)[[1]],
                       rng_normal(key, 1001L, 4, 0))
    }
  }
})

test_that("a C API stream stops when its consumer says so", {
  skip_on_cran()
  ns <- client()
  key <- rng_key(3L)
  s <- ns$client_stream(key, 10000L, 1000L, "uniform", 0, 1, stop_after = 3L)
  expect_identical(s[[4]], 3L)                         # ZURAND_STOPPED
  expect_identical(s[[3]], c(1000, 1000, 1000))
  expect_identical(s[[1]][1:3000], rng_uniform(key, 3000L))
  expect_true(all(is.na(s[[1]][3001:10000])))
})

test_that("C API streams run one per worker thread", {
  skip_on_cran()
  ns <- client()
  keys <- rng_key(11L, n = 8L)
  expect_identical(ns$client_stream_par(keys, 20000L, 3000L, "normal", 0, 1,
                                        "mcfarland"),
                   rng_normal(keys, 20000L, method = "mcfarland"))
  expect_identical(ns$client_stream_par(keys, 20000L, 4096L, "uniform", 0, 1),
                   rng_uniform(keys, 20000L))
  expect_identical(ns$client_stream_par(keys, 20000L, 333L, "uniform", -1e308, 1e308),
                   rng_uniform(keys, 20000L, -1e308, 1e308))
})
