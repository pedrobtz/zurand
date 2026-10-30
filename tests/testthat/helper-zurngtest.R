# Helpers for test-zurngtest.R (the PractRand fixture), kept out of it so that the
# tests run shuffled.

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
