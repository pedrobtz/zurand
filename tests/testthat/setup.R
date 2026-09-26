# CRAN allows a package's checks at most two threads. zurand's fills go
# parallel above 32,768 values with every core the OpenMP runtime offers,
# so cap them at two for the whole test run there. Elsewhere (NOT_CRAN set,
# as in CI and devtools::test()) the machine's default stays, so the
# threaded paths are exercised at full width.
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  rng_threads(2L)
}
