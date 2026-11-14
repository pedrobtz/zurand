# The R layer's contract: findings of the 2026-11-14 bug hunt (#123-#131,
# #136). Every error is asserted by class (zurand_error), never by message.

test_that("every error zurand raises has class zurand_error", {
  key <- rng_key(1L)
  expect_error(rng_uniform(key, -1L), class = "zurand_error") # from C
  expect_error(rng_integer(key, 1L), class = "zurand_error") # from R
  expect_error(key[1, 2], class = "zurand_error")
  expect_error(rng_fold(key, factor("a")), class = "zurand_error")
  # The call still names the user's function, as Rf_error() did.
  e <- tryCatch(rng_uniform(key, -1L), error = identity)
  expect_identical(conditionCall(e)[[1L]], quote(rng_uniform))
})

test_that("a native non-ASCII string folds the same in every locale (#123)", {
  key <- rng_key(1L, engine = "philox4x64")
  utf8 <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xc3, 0xa9)))
  Encoding(utf8) <- "UTF-8"
  latin1 <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
  Encoding(latin1) <- "latin1"
  native <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xc3, 0xa9))) # unmarked
  ref <- rng_fold(key, utf8)
  escaped <- rng_fold(key, "caf<c3><a9>")
  for (loc in c("C", "en_US.UTF-8", "C.UTF-8")) {
    ok <- suppressWarnings(tryCatch(
      withr::with_locale(c(LC_CTYPE = loc), Sys.getlocale("LC_CTYPE") == loc),
      error = function(e) FALSE
    ))
    if (!isTRUE(ok)) {
      next
    }
    withr::with_locale(c(LC_CTYPE = loc), {
      expect_identical(rng_fold(key, latin1), ref, info = loc)
      expect_identical(rng_fold(key, native), ref, info = loc)
      expect_false(identical(rng_fold(key, native), escaped), info = loc)
    })
  }
  # In the C locale, native bytes that are not UTF-8 have no defined text.
  withr::with_locale(c(LC_CTYPE = "C"), {
    bad <- rawToChar(as.raw(c(0x63, 0xe9)))
    expect_error(rng_fold(key, bad), class = "zurand_error")
  })
  # ASCII is unaffected anywhere.
  withr::with_locale(c(LC_CTYPE = "C"), {
    expect_identical(
      rng_fold(key, "layer1"),
      withr::with_locale(c(LC_CTYPE = "en_US.UTF-8"), rng_fold(key, "layer1"))
    )
  })
})

test_that("zero keys give an n x 0 matrix or the dimension error (#124)", {
  ks <- rng_key(42, n = 5)[0]
  expect_identical(dim(rng_uniform(ks, 3L)), c(3L, 0L))
  big <- c(.Machine$integer.max + 1, 2^32 + 3)
  for (n in big) {
    expect_error(rng_uniform(ks, n), class = "zurand_error")
    expect_error(rng_normal(ks, n), class = "zurand_error")
    expect_error(rng_integer(ks, n, 1L, 2L), class = "zurand_error")
    expect_error(rng_bits(ks, n), class = "zurand_error")
    expect_error(rng_bits(ks, n, bits = 64L), class = "zurand_error")
  }
  expect_error(rng_uniform(ks, 2^40, offset = 5), class = "zurand_error")
})

test_that("arithmetic on keys is an error, comparison is by whole key (#125)", {
  keys <- rng_key(1, n = 5)
  for (f in list(
    function(k) k + 1L,
    function(k) -k,
    function(k) k %% 2L,
    function(k) k[1] - k[1],
    function(k) k * 2,
    function(k) !k,
    abs,
    round,
    sum,
    max,
    range,
    t
  )) {
    expect_error(f(keys), class = "zurand_error")
  }
  expect_identical(keys == keys, rep(TRUE, 5L))
  expect_identical(keys != keys[1], c(FALSE, rep(TRUE, 4L)))
  expect_identical(keys[2] == keys, c(FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(keys[0] == keys, logical())
  other <- rng_key(1, n = 5, engine = "philox4x64") # same words, other engine
  expect_identical(keys == other, rep(FALSE, 5L))
  expect_error(keys == 1L, class = "zurand_error")
  expect_error(keys[1:2] == keys[1:3], class = "zurand_error")
})

test_that("classed numbers are refused on every numeric argument (#126)", {
  key <- rng_key(1L)
  classed <- list(
    factor("5"),
    structure(2, class = "integer64"),
    as.Date("2020-01-01"),
    as.difftime(3, units = "secs")
  )
  for (x in classed) {
    expect_error(rng_key(x), class = "zurand_error")
    expect_error(rng_key(1L, n = x), class = "zurand_error")
    expect_error(rng_uniform(key, x), class = "zurand_error")
    expect_error(
      rng_uniform(key, 1L, min = x, max = 10),
      class = "zurand_error"
    )
    expect_error(rng_uniform(key, 1L, offset = x), class = "zurand_error")
    expect_error(rng_normal(key, 1L, mean = x), class = "zurand_error")
    expect_error(rng_normal(key, 1L, sd = x), class = "zurand_error")
    expect_error(rng_integer(key, 1L, x, 10L), class = "zurand_error")
    expect_error(rng_bits(key, 1L, bits = x), class = "zurand_error")
    expect_error(rng_lazy_uniform(key, x), class = "zurand_error")
    expect_error(rng_threads(x), class = "zurand_error")
  }
  expect_identical(rng_key(5), rng_key(5L)) # plain numbers still work
})

test_that("a key vector takes one index in [[ too (#128)", {
  keys <- rng_key(1, n = 5)
  expect_error(keys[[2, 3]], class = "zurand_error")
  expect_error(keys[[2, 99, "x"]], class = "zurand_error")
  expect_identical(keys[[2]], keys[2])
})

test_that("[ and c() keep the samplers' stream rule (#130)", {
  for (bad in list(1.5, "1", TRUE, c(1L, 1L), NA_integer_, 0L)) {
    x <- rng_key(1L)
    attr(x, "stream") <- bad
    expect_error(x[1], class = "zurand_error")
    expect_error(c(x, x), class = "zurand_error")
    expect_error(rng_uniform(x), class = "zurand_error")
  }
  x <- rng_key(1L)
  attr(x, "stream") <- 1 # a whole double is fine
  expect_identical(attr(x[1], "stream"), 1L)
})

test_that("rng_simd() takes NULL, TRUE, FALSE or its own return value (#131)", {
  old <- rng_simd()
  withr::defer(rng_simd(old))
  for (bad in list(c(TRUE, FALSE), 1, 0L, "yes", NA)) {
    expect_error(rng_simd(bad), class = "zurand_error")
  }
  expect_identical(rng_simd(), old) # a refused call changes nothing
  prev <- rng_simd(FALSE)
  expect_identical(prev, old)
  rng_simd(prev)
  expect_identical(rng_simd(), old)
})

test_that("key vectors work with head(), tail(), match() and as.character() (#129, #136)", {
  keys <- rng_key(1, n = 8)
  expect_identical(head(keys, 2), keys[1:2])
  expect_identical(tail(keys, 3), keys[6:8])
  expect_identical(head(keys[0]), keys[0])
  expect_identical(keys[NULL], keys[0])
  expect_identical(as.character(keys), format(keys))
  expect_identical(match(keys[c(3, 1)], keys), c(3L, 1L))
  expect_identical(keys[2:3] %in% keys[3:5], c(FALSE, TRUE))
  other <- rng_key(1, n = 8, engine = "philox4x64")
  expect_identical(match(other, keys), rep(NA_integer_, 8L))
  expect_error(names(keys) <- letters[1:8], class = "zurand_error")
})
