# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

# testthat is a suggested package, and CRAN runs a check flavour with
# suggested packages unavailable, where an unguarded library(testthat) is an
# ERROR. Running no tests there is the accepted outcome: that flavour checks
# that the package stands up without its suggested packages.
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(zurand)

  test_check("zurand")
}
