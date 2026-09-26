# The C API (inst/include/zurand.h), exercised the way a downstream package
# uses it: a fixture package, tests/testthat/zurandclient, is built against
# the header, fetches the API table, and fills one key per column inside an
# OpenMP parallel loop. Everything it returns must equal the R samplers.
#
# Not on CRAN: it compiles a package, which takes seconds and needs a
# toolchain. Everywhere else a build failure is a test failure, not a skip.

client_lib <- NULL

# Build once per file; returns the loaded fixture namespace.
client <- function() {
  if (is.null(client_lib)) client_lib <<- install_client()
  loadNamespace("zurandclient", lib.loc = client_lib)
}

install_client <- function() {
  inc <- system.file("include", package = "zurand")
  if (!nzchar(inc) || !file.exists(file.path(inc, "zurand.h")))
    stop("zurand.h not found under system.file(\"include\", package = \"zurand\")")
  src <- file.path(tempfile("zurandclient-src"), "zurandclient")
  dir.create(src, recursive = TRUE)
  file.copy(list.files(test_path("zurandclient"), full.names = TRUE), src,
            recursive = TRUE)
  writeLines(c(sprintf('PKG_CPPFLAGS = -I"%s"', normalizePath(inc, "/")),
               "PKG_CFLAGS = $(SHLIB_OPENMP_CFLAGS)",
               "PKG_LIBS = $(SHLIB_OPENMP_CFLAGS)"),
             file.path(src, "src", "Makevars"))
  lib <- tempfile("zurandclient-lib")
  dir.create(lib)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "--no-multiarch",
      paste0("--library=", shQuote(lib)), shQuote(src)),
    stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
    stop("building the C API fixture failed:\n", paste(out, collapse = "\n"))
  lib
}

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
    # Degenerate scales, handled before the fill.
    expect_identical(fill(keys, 5L, "uniform", 2, 2), rng_uniform(keys, 5L, 2, 2))
    expect_identical(fill(keys, 5L, "normal", 4, 0), rng_normal(keys, 5L, 4, 0))
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
  # NA integer bound, unknown engine, unknown stream, and n = 0 (fine).
  expect_identical(ns$client_errors(rng_key(1L)),
                   c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 0L))
})
