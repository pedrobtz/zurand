# zurand's streams tested with PractRand through zurngtest's C callable, the
# way a generator package written in C uses it: a fixture package,
# tests/testthat/zurandpract, is built against zurand.h and zurngtest.h.
#
# zurngtest is not a dependency of zurand: the test runs where it is
# installed and skips elsewhere, and never on CRAN, since it compiles a
# package and runs megabytes of tests.

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
