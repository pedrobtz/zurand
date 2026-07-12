# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What rngat is

`rngat` provides random access ("at") random numbers for R, built on the Philox4x64-10 counter-based generator from the vendored **Random123** library. Every value is a pure function of `(key, index, domain)` — there is no generator state. Exported API: `rng_key(seed)`, `fold_in(key, identity)`, the vectorized draws `bits_at()` / `runif_at()` / `rnorm_at()`, and the contiguous-run fast path `bits_seq()` / `runif_seq()` / `rnorm_seq(key, start, len)` (bit-identical to the `_at` form over `start + 0:(len-1)`, enforced by tests; deliberately no `by` argument — strides share no Philox blocks, so a strided fast path can't exist).

## Architecture

- **Core construction** (all in [src/rngat.c](src/rngat.c)): counter = `{index >> 2, domain, purpose, 0}`, run through `philox4x64(ctr, key)`; the returned value is word `index & 3` of the 256-bit block, so four consecutive indices share one Philox call (the C loop caches the last block). The `purpose` word (0 = draws, 1 = `fold_in`) domain-separates key derivation from value draws — preserve this if adding new counter uses. Invariants relied on by tests (all relational, so the 4-per-block packing keeps them green): `rnorm_at == qnorm(runif_at)` (inversion from the same bits), `runif_at` strictly inside (0,1) (top 53 bits + 0.5, scaled by 2^-53), `bits_at` returns uint32-as-double. Note the value-vs-index mapping is *not* stable across this packing change — treat any change to the counter layout as breaking.
- **Keys** are 16-byte raw vectors of class `rngat_key`: the two 64-bit Philox key words serialized explicitly little-endian (see `load_le64`/`store_le64`), so keys and results are portable across architectures — don't replace this with `memcpy`.
- **Performance-critical inlining.** `rngat_block` must keep its `R123_FORCE_INLINE` forward declaration: plain `static inline` is only a hint, and clang declines it (the function is large after Philox's own forced inlining), leaving an out-of-line call per block whose 32-byte struct return goes through memory on arm64 — measured ~60% slower bulk draws. Likewise the `fill_seq_*` loops are macro-stamped per draw kind because a runtime `draw_kind` parameter did not get specialized by the inliner. If touching this, verify with `otool -tv rngat.so` that `draw_seq`/`draw_at` contain no `bl _rngat_block`. Benchmarks must use `R CMD INSTALL` (`-O2`), never `devtools::load_all()` (debug `-O0`, ~15× slower here). Reference numbers (M-series, 1e6 uniforms): `runif_seq` ≈ 2.6 ms, `runif_at(1:n)` ≈ 5.3 ms, randompack ≈ 2.9 ms.
- **Vendored code.** `src/Random123/` (philox.h, array.h, features/, LICENSE) is copied verbatim from upstream at `/Users/pbtz/Documents/repos/gh/public/random123` — don't edit it; re-copy from upstream to update. `src/Makevars` defines `-DR123_USE_MULHILO64_C99=1` as a portable fallback for platforms without `__uint128_t`; the fast path is still auto-selected elsewhere.
- **R ↔ C bridge.** [R/rngat-package.R](R/rngat-package.R) declares `@useDynLib rngat, .registration = TRUE`; entry points are registered in `R_init_rngat` ([src/rngat.c](src/rngat.c)) and called as `.Call(C_name, ...)` from [R/rngat.R](R/rngat.R). `R_NO_REMAP` is set, so all R API calls in C must use the `Rf_` prefix.
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

## Before CRAN / R CMD check

`Authors@R` and `License` in `DESCRIPTION` are still scaffold placeholders and will fail `R CMD check` until filled in (e.g. `usethis::use_mit_license()`). The vendored Random123 code is BSD 3-clause (D. E. Shaw Research) — its license text is kept at `src/Random123/LICENSE` and the copyright holder must be credited (cph role in `Authors@R` or `inst/COPYRIGHTS`) before release.
