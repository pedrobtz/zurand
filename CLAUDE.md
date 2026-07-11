# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project state

`rngat` is a newly scaffolded R package with a C backend. It was created with `usethis` and is essentially empty: `DESCRIPTION` still holds placeholder metadata, `NAMESPACE` exports nothing, `src/rngat.c` only includes R headers, and there are no R functions, tests, or docs yet. Expect to be building the package up from this scaffold rather than modifying existing behavior.

The package name (`rng` + `at`) and the sibling working directory point to the intended purpose: exposing the **Random123** counter-based RNG library (Philox, Threefry, ARS, AES) to R.

## Architecture

- **R ↔ C bridge.** [R/rngat-package.R](R/rngat-package.R) declares `@useDynLib rngat, .registration = TRUE`, so C routines are reached from R via `.Call()` with registered `R_CallMethodDef` entries. New C entry points go in [src/rngat.c](src/rngat.c) (using `R.h` / `Rinternals.h`, with `R_NO_REMAP` already set so all R API calls must use the `Rf_` prefix). When adding routines, register them with `R_registerRoutines` + `R_useDynamicSymbols(dll, FALSE)` and typically add an `R_init_rngat` function.
- **Reference library.** The Random123 headers live in the additional working directory `/Users/pbtz/Documents/repos/gh/public/random123` (upstream, header-only, `include/Random123/*.h`). This is the source material to bind, not part of this repo. To vendor it, headers would be copied under `src/` (e.g. `src/Random123/`) and referenced from `rngat.c`; the CBRNGs are stateless `result = CBRNGname(counter, key)` functions.
- **Generated files.** `NAMESPACE` and `man/*.Rd` are produced by roxygen2 (`Roxygen: list(markdown = TRUE)` in DESCRIPTION) — edit the roxygen comments above the code, never these files by hand. Compiled artifacts (`*.o`, `*.so`, `*.dll`) are gitignored under `src/`.

## Common commands

This is a devtools-based package (`PackageUseDevtools: Yes`). From an R session at the repo root:

```r
devtools::document()       # regenerate NAMESPACE + man/ from roxygen comments
devtools::load_all()       # compile src/ and load the package for interactive use
devtools::test()           # run testthat tests (once tests/ exists)
devtools::check()          # full R CMD check (build, install, tests, examples)
```

Single test file / single test (after `use_testthat()` sets up `tests/testthat/`):

```r
testthat::test_file("tests/testthat/test-<name>.R")
devtools::test(filter = "<name>")   # runs tests/testthat/test-<name>.R
```

Recompiling C code: `devtools::load_all()` (or `pkgbuild::compile_dll()`) rebuilds `src/`. Use `devtools::clean_dll()` if stale objects cause link errors.

From the shell instead of an R session:

```sh
R CMD build .        # produces rngat_<version>.tar.gz
R CMD check rngat_*.tar.gz --as-cran
```

## Before this is buildable

The scaffold placeholders in `DESCRIPTION` (`Title`, `Description`, `Authors@R`, `License`) are literal template text and will fail `R CMD check` until filled in — set a real license with e.g. `usethis::use_mit_license()`.
