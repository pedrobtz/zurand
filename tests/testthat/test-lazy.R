# Lazy vectors (rng_lazy_uniform, rng_lazy_normal): every way of reading one
# must give the materialised sampler's values, on every engine, whatever
# state the one-chunk cache is in.

test_that("lazy vectors equal the samplers, read every way", {
  withr::local_seed(1)
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
      for (ci in seq_along(cases)) {
        cs <- cases[[ci]]
        ref <- cs[[2]]
        idx <- c(n, 1, sample(n, min(n, 200), replace = TRUE),
                 seq_len(min(n, 40)) + max(0, n - 60))
        x <- cs[[1]]()
        info <- paste(engine, "n", n, "case", ci)
        expect_identical(length(x), as.integer(n), info = info)
        expect_identical(x[idx], ref[idx], info = info)
        expect_identical(vapply(idx, function(i) x[[i]], 0), ref[idx], info = info)
        expect_identical(x[seq_len(n)], ref, info = info)
        expect_identical(rev(x), rev(ref), info = info)
        expect_identical(sum(x), sum(ref), info = info)
        expect_identical(x[c(NA, n + 1, 1)], ref[c(NA, n + 1, 1)], info = info)
        expect_false(materialised(x), info = info)
        expect_identical(x + 0, ref + 0, info = info)  # materialises
        expect_true(materialised(x), info = info)
        expect_identical(x[idx], ref[idx], info = info)
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
      local_cache_mode(mode)
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
    }
  }
})

test_that("a short lazy vector carries a cache no longer than itself", {
  # #60: every lazy vector held a 5120-double (40 KB) cache, so 10,000 of
  # length 10 took about 400 MB. Now the cache is min(n, chunk) values.
  k <- rng_key(1L)
  used <- function() sum(gc(full = TRUE)[, 2L]) # Mb, Ncells + Vcells
  m0 <- used()
  xs <- lapply(1:2000, function(i) rng_lazy_normal(k, 10L))
  per_vector_kb <- (used() - m0) * 1024 / length(xs)
  expect_lt(per_vector_kb, 4)
  expect_identical(xs[[2000]][], rng_normal(k, 10L))
  expect_identical(rng_lazy_uniform(k, 0L)[], numeric())
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

test_that("long regions of lazy normals match the sampler at an offset", {
  # x[a:b] with b - a + 1 >= 32,768 goes through the parallel range fill
  # (ZURAND_LAZY_BIG_REGION); start it off every engine's chunk boundary.
  for (engine in engines) {
    k <- rng_key(12L, engine = engine)
    for (method in c("ziggurat", "mcfarland")) {
      x <- rng_lazy_normal(k, 1e5, mean = 2, sd = 0.5, method = method)
      expect_identical(
        x[777:(777 + 40000)],
        rng_normal(k, 40001, mean = 2, sd = 0.5, method = method, offset = 776),
        info = paste(engine, method)
      )
      expect_false(materialised(x))
    }
  }
})

test_that("a lazy vector that was only read still saves as its recipe", {
  k <- rng_key(8L)
  x <- rng_lazy_normal(k, 1e5)
  expect_identical(sd(x), sd(rng_normal(k, 1e5)))   # REAL(): writable pointer
  expect_true(materialised(x))
  expect_lt(length(serialize(x, NULL)), 1000)
  y <- unserialize(serialize(x, NULL))
  expect_false(materialised(y))
  expect_identical(y[], rng_normal(k, 1e5))

  x[7] <- x[7]                               # written, but unchanged
  expect_lt(length(serialize(x, NULL)), 1000)
  x[7] <- 0
  expect_gt(length(serialize(x, NULL)), 8e5)
  expect_identical(unserialize(serialize(x, NULL))[7], 0)

  z <- rng_lazy_uniform(k, 10, min = 3, max = 3)    # constant: s == 0
  invisible(sd(z))
  expect_lt(length(serialize(z, NULL)), 1000)
})

test_that("restoring a lazy vector refuses saved arguments it cannot have written", {
  x <- rng_lazy_normal(rng_key(1L), 10)
  u <- rng_lazy_uniform(rng_key(1L), 10)
  ok <- tamper_lazy(x, 1, c(10, 1, 0, 1))
  expect_identical(unserialize(ok)[], rng_normal(rng_key(1L), 10))
  bad <- list(
    distribution = c(10, 9, 0, 1),           # e.g. a later rng_lazy_exponential()
    distribution = c(10, 1.5, 0, 1),
    arguments = 10,
    arguments = c(10, 1, NaN, 1),
    arguments = c(10, 1, 0, Inf),
    length = c(-5, 1, 0, 1),
    length = c(1e300, 1, 0, 1),
    length = c(10.5, 1, 0, 1),
    scale = c(10, 1, 0, -1),
    scaling = c(10, 1, 0, 1, 1)              # the doubled scaling is uniform-only
  )
  for (i in seq_along(bad)) {
    expect_error(unserialize(tamper_lazy(x, 1, bad[[i]])),
                 paste("saved", names(bad)[i], "is invalid"),
                 info = deparse(bad[[i]]))
  }
  expect_error(unserialize(tamper_lazy(u, 0, c(10, 0, 0, 1, 2))),
               "saved scaling is invalid")
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
  # Base R itself errors on 1e300 as an index of a plain vector on 32-bit
  # builds, so the reference has nothing to compare against there.
  skip_if(.Machine$sizeof.pointer < 8)
  x <- rng_lazy_normal(rng_key(1L), 10L)
  ref <- rng_normal(rng_key(1L), 10L)
  idx <- c(0.5, 1.5, 10.9, 11, NA)
  expect_identical(x[idx], ref[idx])
  expect_true(is.na(x[11]))
  # Indices past R_XLEN_T_MAX give NA up to R 4.6 and "subscript too large"
  # from R-devel r90777 (PR#19261) on; a lazy vector does what a plain one
  # does on the running R.
  for (big in list(c(1, 1e300), c(1, Inf))) {
    want <- tryCatch(ref[big], error = function(e) e)
    if (inherits(want, "error")) {
      expect_error(x[big], conditionMessage(want), fixed = TRUE)
    } else {
      expect_identical(x[big], want)
    }
  }
})

test_that("large lazy sums match with threads on and off", {
  skip_on_cran()
  k <- rng_key(6L)
  ref <- sum(rng_uniform(k, 1e7, -1, 3))
  local_threads(1L)
  one <- sum(rng_lazy_uniform(k, 1e7, -1, 3))
  rng_threads(0L)
  all <- sum(rng_lazy_uniform(k, 1e7, -1, 3))
  expect_identical(one, ref)
  expect_identical(all, ref)
})

test_that("a huge lazy scan stops when interrupted", {
  # sum() over 1e12 lazy values runs for hours; before Get_region checked
  # for interrupts it ignored SIGINT and had to be killed.
  skip_on_cran()
  skip_on_os("windows")    # SIGINT
  skip_if_not_installed("processx")
  out <- interrupt_child(c(
    "x <- rng_lazy_normal(rng_key(1L), 1e12)",
    "writeLines('scanning', started)",
    "cat(tryCatch({ sum(x); 'finished' }, interrupt = function(e) 'interrupted'))"
  ), after = 1)
  expect_identical(out, "interrupted")
})

test_that("reading single elements never signals an interrupt", {
  # REAL_ELT is what C and C++ code calls in its own loops, which cannot
  # expect a longjmp (#77). is.na() reads a lazy vector through Elt alone;
  # interrupted while it runs (about 2 s here), the interrupt must wait for
  # it to return and land in the Sys.sleep() after it, with `y` assigned.
  skip_on_cran()
  skip_on_os("windows")    # SIGINT
  skip_if_not_installed("processx")
  out <- interrupt_child(c(
    "x <- rng_lazy_normal(rng_key(1L, engine = 'philox4x64'), 6e7)",
    "y <- NULL",
    "writeLines('reading', started)",
    "cat(tryCatch({ y <- is.na(x); Sys.sleep(60); 'finished' },",
    "  interrupt = function(e) if (is.null(y)) 'inside is.na()' else 'after is.na()'))"
  ), after = 0.3)
  expect_identical(out, "after is.na()")
})

test_that("lazy vectors check their arguments", {
  expect_error(rng_lazy_normal(rng_key(1L, n = 2L), 10), "single key")
  expect_error(rng_lazy_normal(rng_key(1L), -1), "non-negative")
  expect_error(rng_lazy_normal(rng_key(1L), 10, sd = -1), "non-negative")
  expect_error(rng_lazy_uniform(rng_key(1L), 10, 2, 1), "less than or equal")
  expect_identical(length(rng_lazy_uniform(rng_key(1L), 0)), 0L)
})
