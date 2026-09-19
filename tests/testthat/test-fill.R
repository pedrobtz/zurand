# The in-place fills deliberately break R's copy-on-write rule, so what
# needs pinning is not only that they produce the right numbers but that
# the guard around that choice behaves.

test_that("in-place fills match the allocating samplers exactly", {
  key <- rng_key(42L)
  tkey <- rng_key(42L, engine = "threefry4x64")

  b <- rng_buffer(1000L)
  rng_fill_uniform(b, key)
  expect_identical(as.vector(b), rng_uniform(key, 1000L))

  b <- rng_buffer(1000L)
  rng_fill_normal(b, key)
  expect_identical(as.vector(b), rng_normal(key, 1000L))

  # non-default arguments go through the scaling pass
  b <- rng_buffer(500L)
  rng_fill_uniform(b, key, min = -2, max = 5)
  expect_identical(as.vector(b), rng_uniform(key, 500L, min = -2, max = 5))

  b <- rng_buffer(500L)
  rng_fill_normal(b, key, mean = 2, sd = 3)
  expect_identical(as.vector(b), rng_normal(key, 500L, mean = 2, sd = 3))

  # sd = 0 short circuit, and a length that is not a multiple of 4 so the
  # tail block is exercised
  b <- rng_buffer(7L)
  rng_fill_normal(b, key, mean = 9, sd = 0)
  expect_identical(as.vector(b), rep(9, 7L))
  b <- rng_buffer(997L)
  rng_fill_uniform(b, key)
  expect_identical(as.vector(b), rng_uniform(key, 997L))

  # the second engine fills identically too
  b <- rng_buffer(256L)
  rng_fill_uniform(b, tkey)
  expect_identical(as.vector(b), rng_uniform(tkey, 256L))
})

test_that("only a vector from rng_buffer() can be overwritten", {
  key <- rng_key(42L)
  expect_error(rng_fill_uniform(numeric(10), key), "must come from rng_buffer")
  expect_error(rng_fill_normal(numeric(10), key), "must come from rng_buffer")
  expect_error(rng_fill_uniform(letters, key), "must be a numeric vector")
  # a buffer whose marker was stripped is no longer a buffer
  b <- rng_buffer(10L)
  attr(b, "zurand.buffer") <- NULL
  expect_error(rng_fill_uniform(b, key), "must come from rng_buffer")
})

test_that("refilling the same buffer many times is allowed", {
  # This is the whole point of the API and must not trip the guard.
  key <- rng_key(42L)
  b <- rng_buffer(64L)
  expect_silent(for (i in 1:200) rng_fill_uniform(b, key))
  expect_identical(as.vector(b), rng_uniform(key, 64L))
})

test_that("a buffer that is kept as well as refilled is refused", {
  # The dangerous pattern: every slot would end up holding the same vector,
  # silently, with output that looks random either way.
  key <- rng_key(42L)
  buf <- rng_buffer(64L)
  out <- list()
  expect_warning(
    for (i in 1:4) { rng_fill_uniform(buf, key); out[[i]] <- buf },
    "gained a reference")

  # copying out is the supported way to keep results, and must stay allowed
  buf2 <- rng_buffer(64L)
  kept <- list()
  expect_silent(for (i in 1:4) { rng_fill_uniform(buf2, key); kept[[i]] <- buf2[] })
  expect_identical(length(unique(lapply(kept, as.vector))), 1L)
})

test_that("in-place fills validate their arguments", {
  key <- rng_key(42L)
  # A fresh buffer per expectation: expect_error() itself takes a reference
  # each time, which is exactly what the sharing guard watches for.
  expect_error(rng_fill_uniform(rng_buffer(10L), rng_key(42L, n = 2L)),
               "exactly one key")
  expect_error(rng_fill_uniform(rng_buffer(10L), "not a key"), "rng_key")
  expect_error(rng_fill_uniform(rng_buffer(10L), key, min = 5, max = 1),
               "less than or equal")
  expect_error(rng_fill_normal(rng_buffer(10L), key, sd = -1), "non-negative")
  expect_silent(rng_fill_uniform(rng_buffer(0L), key))   # empty buffer
})
