# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What rngat is

`rngat` provides stateless random numbers for R, built on the Philox4x64-10 counter-based generator from the vendored **Random123** library. Every value is a pure function of an immutable key plus sampler arguments; ordinary samplers never read or mutate `.Random.seed`.

Exported API:

- `rng_key(seed, n = 1, engine = "philox4x64")` creates an opaque S3 vector of keys.
- `rng_fold(key, data)` derives one key per input key from typed data such as strings or whole-number counters.
- `rng_uniform()`, `rng_normal()`, `rng_integer()` and `rng_bits()` are pure samplers.
- `rng_key_from_r()` is the one convenience function that intentionally consumes R's global RNG state.

If two parts of a computation need independent random values, request a key vector up front with `rng_key(seed, n)` or fold the parent key by structured data.

## Architecture

- **Core construction** (all in [src/rngat.c](src/rngat.c)): counter = `{index >> 2, domain, purpose, 0}`, run through `philox4x64(ctr, key)`; the returned value is word `index & 3` of the 256-bit block, so four consecutive output positions share one Philox call. Public samplers use output positions `0:(n - 1)` as the index. The `purpose` word domain-separates bits, uniform, normal, integer and fold streams; preserve this if adding new counter uses.
- **Keys** are integer matrices with class `rng_key`, one key per row and four words per key: `{k0 low, k0 high, k1 low, k1 high}`. R integers are signed and reserve one bit pattern for `NA`, so the C code copies key words by exact 32-bit bit pattern with `memcpy()`. Do not replace this with ordinary integer coercion or NA validation; keys are opaque.
- **Sampling semantics.** Samplers default to `n = 1L`, i.e. one draw per key. `rng_uniform()` uses the top 53 bits plus 0.5 scaled by `2^-53`, so values are strictly inside `(min, max)` when `min < max`. `rng_normal()` uses Box-Muller pairs from its own domain-separated uniform stream. `rng_integer()` uses bounded 32-bit rejection sampling and returns ordinary R integer vectors. With one key, samplers return plain vectors; with multiple keys, samplers return an `n x length(key)` matrix, including a `1 x K` matrix when `n = 1`. `rng_bits(bits = 32)` returns uint32-as-double; `bits = 64` returns fixed-width hex strings to avoid losing bits in R doubles.
- **Optional OpenMP.** `src/Makevars` passes through R's `$(SHLIB_OPENMP_CFLAGS)`. When available, `rng_uniform()`, `rng_normal()` and `rng_integer()` parallelize over key columns above a size threshold, and within a single column over independent Philox blocks (bit-identical to serial output) when there is no column parallelism to exploit. Worker threads must not call R API; samplers hoist `INTEGER(key)` into a `const int *` before the parallel region and workers only read that pointer and write primitive numeric/integer output. `rng_bits(bits = 64)` remains single-threaded because it creates R strings. Apple clang has no bundled OpenMP, so on macOS `~/.R/Makevars` must define `SHLIB_OPENMP_CFLAGS` (`-Xclang -fopenmp` plus libomp include/lib flags); this machine is set up that way and links R 4.6's bundled `libomp.dylib`.
- **Performance-critical inlining.** `rngat_block` keeps the `R123_FORCE_INLINE` forward declaration: plain `static inline` is only a hint, and clang may decline it (the function is large after Philox's own forced inlining), leaving an out-of-line call per block whose 32-byte struct return goes through memory on arm64. If touching this, verify with `otool -tv rngat.so` that the hot sampler loops contain no call to `_rngat_block`. Benchmarks must use `R CMD INSTALL` (`-O2`), not `devtools::load_all()` debug builds.
- **Vendored code.** `src/Random123/` (philox.h, array.h, features/, LICENSE) is copied verbatim from upstream at `/Users/pbtz/Documents/repos/gh/public/random123` — don't edit it; re-copy from upstream to update. `src/Makevars` defines `-DR123_USE_MULHILO64_C99=1` as a portable fallback for platforms without `__uint128_t`; the fast path is still auto-selected elsewhere.
- **R ↔ C bridge.** [R/rngat-package.R](R/rngat-package.R) declares `@useDynLib rngat, .registration = TRUE`; entry points are registered in `R_init_rngat` ([src/rngat.c](src/rngat.c)) and called as `.Call(C_name, ...)` from [R/rngat.R](R/rngat.R). `R_NO_REMAP` is set, so all R API calls in C must use the `Rf_` prefix.
- **Generated files.** `NAMESPACE` and `man/*.Rd` are produced by roxygen2 (`Roxygen: list(markdown = TRUE)` in DESCRIPTION) — edit the roxygen comments above the code, then regenerate these files. Compiled artifacts (`*.o`, `*.so`, `*.dll`) are gitignored under `src/`.

## Common commands

This is a devtools-based package (`PackageUseDevtools: Yes`). From an R session at the repo root:

```r
devtools::document()       # regenerate NAMESPACE + man/ from roxygen comments
devtools::load_all()       # compile src/ and load the package for interactive use
devtools::test()           # run testthat tests
devtools::check()          # full R CMD check (build, install, tests, examples)
```

Single test file / single test:

```r
testthat::test_file("tests/testthat/test-rngat.R")
devtools::test(filter = "rngat")
```

Recompiling C code: `devtools::load_all()` (or `pkgbuild::compile_dll()`) rebuilds `src/`. Use `devtools::clean_dll()` if stale objects cause link errors.

From the shell instead of an R session:

```sh
R CMD build .        # produces rngat_<version>.tar.gz
R CMD check rngat_*.tar.gz --as-cran
```

## Before CRAN / R CMD check

`Authors@R` and `License` in `DESCRIPTION` are still scaffold placeholders and will fail or warn in `R CMD check` until filled in (e.g. `usethis::use_mit_license()`). The vendored Random123 code is BSD 3-clause (D. E. Shaw Research); its license text is kept at `src/Random123/LICENSE` and the copyright holder must be credited (cph role in `Authors@R` or `inst/COPYRIGHTS`) before release.
