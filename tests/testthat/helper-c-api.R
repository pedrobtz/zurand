# Helpers for test-c-api.R (the C API fixture package, built once per session), kept out of it so that the
# tests run shuffled.

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
  # The OpenMP flags zurand's own configure chooses: R's SHLIB_OPENMP_CFLAGS,
  # except on macOS, where R leaves that empty and configure tries
  # -Xclang -fopenmp -lomp first, falling back when they do not build (#60).
  # Otherwise the fixture's loops would run on the main thread on macOS
  # even where zurand itself is threaded.
  flags <- list(c("$(SHLIB_OPENMP_CFLAGS)", "$(SHLIB_OPENMP_CFLAGS)"))
  if (Sys.info()[["sysname"]] == "Darwin")
    flags <- c(list(c("-Xclang -fopenmp", "-lomp")), flags)
  for (i in seq_along(flags)) {
    lib <- try_install_client(inc, flags[[i]][1], flags[[i]][2])
    if (is.character(lib)) return(lib)
  }
  stop("building the C API fixture failed:\n", paste(attr(lib, "log"), collapse = "\n"))
}

# The library path on success, otherwise FALSE with the build log attached.
try_install_client <- function(inc, cflags, libs) {
  src <- file.path(tempfile("zurandclient-src"), "zurandclient")
  dir.create(src, recursive = TRUE)
  file.copy(list.files(test_path("zurandclient"), full.names = TRUE), src,
            recursive = TRUE)
  writeLines(c(sprintf('PKG_CPPFLAGS = -I"%s"', normalizePath(inc, "/")),
               paste("PKG_CFLAGS =", cflags),
               paste("PKG_LIBS =", libs)),
             file.path(src, "src", "Makevars"))
  lib <- tempfile("zurandclient-lib")
  dir.create(lib)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "--no-multiarch",
      paste0("--library=", shQuote(lib)), shQuote(src)),
    stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
    return(structure(FALSE, log = out))
  lib
}
