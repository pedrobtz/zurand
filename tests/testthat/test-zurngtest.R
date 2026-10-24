# zurand's streams tested with PractRand through zurngtest's C callable, the
# way a generator package written in C uses it: a fixture package,
# tests/testthat/zurandpract, is built against zurand.h and zurngtest.h.
#
# zurngtest is not a dependency of zurand: the test runs where it is
# installed and skips elsewhere, and never on CRAN, since it compiles a
# package and runs megabytes of tests.

install_pract_fixture <- function() {
  inc <- c(
    system.file("include", package = "zurand"),
    system.file("include", package = "zurngtest")
  )
  src <- file.path(tempfile("zurandpract-src"), "zurandpract")
  dir.create(src, recursive = TRUE)
  file.copy(list.files(test_path("zurandpract"), full.names = TRUE), src,
            recursive = TRUE)
  writeLines(
    sprintf("PKG_CPPFLAGS = %s",
            paste0('-I"', normalizePath(inc, "/"), '"', collapse = " ")),
    file.path(src, "src", "Makevars")
  )
  lib <- tempfile("zurandpract-lib")
  dir.create(lib)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "--no-multiarch",
      paste0("--library=", shQuote(lib)), shQuote(src)),
    stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
    stop("building the zurngtest fixture failed:\n", paste(out, collapse = "\n"))
  loadNamespace("zurandpract", lib.loc = lib)
}

test_that("zurand streams pass PractRand through zurngtest's C callable", {
  skip_on_cran()
  skip_if_not(requireNamespace("zurngtest", quietly = TRUE),
              "zurngtest is not installed")
  ns <- install_pract_fixture()
  for (engine in c("philox4x64", "threefry4x64", "xoshiro256pp")) {
    res <- ns$pract_zurand(
      rng_key(20261024L, engine = engine),
      list(min_bytes = 2^20, max_bytes = 2^24)
    )
    expect_s3_class(res, "zurngtest_result")
    expect_true(is.na(attr(res, "failed_at")), label = engine)
  }
})
