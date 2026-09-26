# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What zurand is

`zurand` provides stateless random numbers for R, built on the vendored **Random123** library. Three engines: `philox4x64` (Philox4x64-10, the default), `threefry4x64` (Threefry4x64-13), and `xoshiro256pp`, a hybrid in which Philox derives a fresh xoshiro256++ state for each 512-value sub-chunk. Every value is a pure function of an immutable key plus sampler arguments; ordinary samplers never read or mutate `.Random.seed`.

Exported API:

- `rng_key(seed, n = 1, engine = c("philox4x64", "threefry4x64", "xoshiro256pp"))` creates an opaque S3 vector of keys; a key records its engine.
- `rng_fold(key, data)` derives one key per input key from typed data such as strings or whole-number counters.
- `rng_uniform()`, `rng_normal()`, `rng_integer()` and `rng_bits()` are pure samplers.
- `rng_key_from_r()` is the one convenience function that intentionally consumes R's global RNG state.
- `rng_threads(n)` caps OpenMP use; `rng_simd(enable)` reports or disables the AVX2 path for `xoshiro256pp`. Neither changes any output.

If two parts of a computation need independent random values, request a key vector up front with `rng_key(seed, n)` or fold the parent key by structured data.

## Architecture

- **Core construction** ([src/zurand.c](src/zurand.c) plus [src/zurand_engine.h](src/zurand_engine.h), which `zurand.c` includes once per engine with `ZE_*` macros defined, so each engine gets its own monomorphic copy of the sampler loops without a runtime branch): counter = `{index >> 2, domain, purpose, tag}`, run through `philox4x64(ctr, key)`; the returned value is word `index & 3` of the 256-bit block, so four consecutive output positions share one Philox call. Public samplers use output positions `0:(n - 1)` as the index. The `purpose` word domain-separates bits, uniform, normal, integer and fold streams; preserve this if adding new counter uses. Word 3 is the engine tag (`ZURAND_TAG_*`, `ZE_TAG` in the engine header): 0 for philox and threefry, 1 for every Philox call made on behalf of `xoshiro256pp` — its sub-chunk seeds `{s, 0, purpose, 1}` and its retry words. Without it the xoshiro stream was a function of the philox stream; any new engine that borrows Philox needs a tag of its own. `rng_fold()` deliberately uses tag 0 for xoshiro keys, so a seed folds to the same key words under both engines. Its AVX2 path runs four sub-chunks in four lanes, so it is bit-identical to the scalar path; `tests/testthat/test-simd.R` enforces this.
- **Keys** are integer matrices with class `rng_key`, one key per row and four words per key: `{k0 low, k0 high, k1 low, k1 high}`. R integers are signed and reserve one bit pattern for `NA`, so the C code copies key words by exact 32-bit bit pattern with `memcpy()`. Do not replace this with ordinary integer coercion or NA validation; keys are opaque.
- **Sampling semantics.** Samplers default to `n = 1L`, i.e. one draw per key. `rng_uniform()` stuffs the top 52 bits into the mantissa of a double in `[1, 2)` and subtracts `1 - 2^-53`, yielding exactly `(m + 0.5) * 2^-52` — strictly inside `(0, 1)`, with no int-to-float conversion in the hot loop. (Do not revert to the old `((bits >> 11) + 0.5) * 2^-53` form: the addition rounds for `bits >= 2^52` and can produce exactly 1.0.) `min`/`max` are applied in a separate pass only when non-default, like `mean`/`sd` for normals; `min + span * u` can round to exactly `min` or `max` when `|min|` is large relative to the span (about 12% of draws for `min = 1e15, max = 1e15 + 1`). `rng_normal()` uses a ziggurat sampler (NumPy's tables and layer layout, `src/numpyzig/`): the attempt-0 word from the shared block decides ~99% of draws with one table compare, and the rare wedge/tail paths draw more words at `(index, attempt >= 1)` in the domain slot so draw `i` stays a pure function of `(key, i)`; keep the cold `zig_normal_slow()` out of line so the fast path inlines into the fill loops. **No fused multiply-add, on any compiler:** every floating product that feeds an addition must go through `zurand_rounded()` (or `zurand_rounded2()` for the two-wide scaling pass, `affine_pass()`), because GCC ignores `#pragma STDC FP_CONTRACT OFF` and contracts by default on arm64, and one fused instruction changes output by an ulp. After touching arithmetic, compile `src/zurand.c` with `clang -O2 -mfma -ffp-contract=fast` and confirm the object has no `vfmadd`/`vfmsub`; the `ubuntu-24.04-arm` R-CMD-check leg is the end-to-end gate. `rng_integer()` uses bounded 32-bit rejection sampling and returns ordinary R integer vectors. With one key, samplers return plain vectors; with multiple keys, samplers return an `n x length(key)` matrix, including a `1 x K` matrix when `n = 1`. `rng_bits(bits = 32)` returns uint32-as-double; `bits = 64` returns fixed-width hex strings to avoid losing bits in R doubles.
- **Optional OpenMP.** `src/Makevars` passes through R's `$(SHLIB_OPENMP_CFLAGS)`. When available, `rng_uniform()`, `rng_normal()` and `rng_integer()` parallelize over key columns above a size threshold, and within a single column over independent Philox blocks (bit-identical to serial output) when there is no column parallelism to exploit. `rng_threads(n)` caps zurand's thread use per session via a package-local count passed as `num_threads()` on each pragma (not `omp_set_num_threads()`, which would be process-global); every parallel region must carry both the `if()` size gate and `num_threads(nt)`. Worker threads must not call R API; samplers hoist `INTEGER(key)` into a `const int *` before the parallel region and workers only read that pointer and write primitive numeric/integer output. `rng_bits(bits = 64)` remains single-threaded because it creates R strings. **On a machine without OpenMP the pragmas are never parsed**, so a `default(none)` region missing a data-sharing clause builds clean locally and fails every Linux leg in CI; run `sh dev/check-openmp.sh` before pushing, which syntax-checks them against a stub `omp.h`. Apple clang has no bundled OpenMP, so on macOS `~/.R/Makevars` must define `SHLIB_OPENMP_CFLAGS` (`-Xclang -fopenmp` plus libomp include/lib flags); this machine is set up that way and links R 4.6's bundled `libomp.dylib`.
- **Performance-critical inlining.** `zurand_block` keeps the `R123_FORCE_INLINE` forward declaration: plain `static inline` is only a hint, and clang may decline it (the function is large after Philox's own forced inlining), leaving an out-of-line call per block whose 32-byte struct return goes through memory on arm64. If touching this, verify with `otool -tv zurand.so` that the hot sampler loops contain no call to `_zurand_block`. Benchmarks must be run against an `-O2` build. `devtools::load_all()` defaults to `debug = TRUE`, which injects `-UNDEBUG -g -O0` and measures **11x slower** -- 17.7 vs 195 M normals/s. Pass `debug = FALSE` and it is equivalent to `R CMD INSTALL`: the compile line differs only by an extra `-Wall -pedantic`, and measured against a fixed in-session reference the two agree to 0.2% (1.058 vs 1.056 over three alternating rounds). Either is fine; `load_all(debug = FALSE)` iterates faster. If `devtools::load_all(debug = FALSE)` warns "Arguments in `...` must be used", nothing was recompiled and the existing objects were reused with whatever flags built them -- run `devtools::clean_dll()` first when switching from a debug build. Note `R CMD INSTALL .` builds in-place, so `rm src/*.o` first if you want to be certain nothing is reused.
- **Vendored code.** `src/Random123/` (philox.h, array.h, features/, LICENSE) is copied verbatim from upstream at `/Users/pbtz/Documents/repos/gh/public/random123` — don't edit it; re-copy from upstream to update. `src/Makevars` defines `-DR123_USE_MULHILO64_C99=1` as a portable fallback for platforms without `__uint128_t`; the fast path is still auto-selected elsewhere. `src/numpyzig/` (ziggurat_constants.h, LICENSE) is NumPy's ziggurat table header, also verbatim (BSD 3-clause); only the normal-double tables are used and the include site suppresses `-Wunused-const-variable` for the rest. `src/zurand_fdlibm.h` is **derived, not verbatim**: fdlibm's `exp` and `log1p` (Sun, permissive notice kept in the file), which the ziggurat's slow path uses instead of the platform libm, because Apple and glibc/musl disagree often enough to change `rng_normal()` streams (the whole-stream digests in `tests/testthat/test-digest.R` catch it). It is edited — no-fusion barriers, memcpy word access — so it must be maintained by hand, never re-copied. `src/zigbounds.h` is generated — not vendored — by `tools/generate-zig-bounds.R` from the NumPy tables: fixed-point brackets of the wedge test that let `zig_normal_slow()` resolve most wedge decisions without `exp()`; rerun the script instead of editing it.
- **R ↔ C bridge.** [R/zurand-package.R](R/zurand-package.R) declares `@useDynLib zurand, .registration = TRUE`; entry points are registered in `R_init_zurand` ([src/zurand.c](src/zurand.c)) and called as `.Call(C_name, ...)` from [R/zurand.R](R/zurand.R). `R_NO_REMAP` is set, so all R API calls in C must use the `Rf_` prefix.
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
testthat::test_file("tests/testthat/test-zurand.R")
devtools::test(filter = "zurand")
```

Recompiling C code: `devtools::load_all()` (or `pkgbuild::compile_dll()`) rebuilds `src/`. Use `devtools::clean_dll()` if stale objects cause link errors.

From the shell instead of an R session:

```sh
R CMD build .        # produces zurand_<version>.tar.gz
R CMD check zurand_*.tar.gz --as-cran
```

## Before CRAN / R CMD check

`Authors@R`, `License` (MIT) and `inst/COPYRIGHTS` are filled in, and they credit both vendored components: Random123 (D. E. Shaw Research) and NumPy's ziggurat tables (NumPy Developers). What is still open is tracked in `dev/roadmap.md` (the version is still `0.0.0.9000`, and `cran-comments.md` does not exist yet).
