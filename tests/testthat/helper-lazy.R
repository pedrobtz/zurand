# Helpers for test-lazy.R, kept out of it so that the tests run shuffled.

engines <- c("xoshiro256pp", "philox4x64", "threefry4x64")

materialised <- function(x) zurand:::rng_lazy_materialised(x)

cache_mode <- function(mode = NULL) zurand:::rng_lazy_cache(mode)

# serialize(x) with its saved arguments replaced by `meta`, as a damaged
# file or a newer zurand would write them. `x` must be a lazy vector with
# the default scaling (mean 0, sd 1, or min 0, max 1), whose arguments are
# then c(length(x), dist, 0, 1): with their XDR length prefix, a byte
# pattern that occurs once in the stream.
tamper_lazy <- function(x, dist, meta) {
  xdr <- function(v) {
    c(writeBin(length(v), raw(), endian = "big"),
      writeBin(as.double(v), raw(), endian = "big"))
  }
  bytes <- serialize(x, NULL, xdr = TRUE)
  pat <- xdr(c(length(x), dist, 0, 1))
  starts <- seq_len(length(bytes) - length(pat) + 1L)
  at <- starts[vapply(
    starts,
    function(i) identical(bytes[i:(i + length(pat) - 1L)], pat),
    logical(1L)
  )]
  stopifnot(length(at) == 1L)
  c(bytes[seq_len(at - 1L)], xdr(meta), bytes[-seq_len(at + length(pat) - 1L)])
}


# Runs `code` (lines of R) in a child Rscript with zurand loaded and sends
# it SIGINT `after` seconds after it writes "started" to the file named by
# the variable `started`; returns what the child printed. processx, not
# system2(wait = FALSE): a shell's background job starts with SIGINT ignored.
interrupt_child <- function(code, after) {
  path <- normalizePath(find.package("zurand"), "/")
  load <- if (dir.exists(file.path(path, "src"))) {
    # A source tree, as under devtools::test(): reuse its compiled library.
    sprintf("pkgload::load_all(%s, compile = FALSE, quiet = TRUE)", deparse(path))
  } else {
    sprintf("library(zurand, lib.loc = %s)", deparse(dirname(path)))
  }
  started <- tempfile()
  script <- tempfile(fileext = ".R")
  writeLines(c(load, sprintf("started <- %s", deparse(started)), code), script)
  child <- processx::process$new(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = "|", stderr = "|")
  on.exit({
    if (child$is_alive()) child$kill()
    unlink(c(started, script))
  })
  deadline <- Sys.time() + 60
  while (child$is_alive() && Sys.time() < deadline &&
         !(file.exists(started) && file.size(started) > 0)) {
    Sys.sleep(0.05)
  }
  testthat::expect_true(child$is_alive())
  Sys.sleep(after)
  child$interrupt()
  child$wait(15000)
  testthat::expect_false(child$is_alive())
  child$read_all_output()
}
