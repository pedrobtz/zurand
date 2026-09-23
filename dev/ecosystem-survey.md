# What other RNG ecosystems do, and what zurand could take from them

Survey date: 2026-09-23. Companion to `dev/review-2026-09-23.md`, which
reviews the package and the roadmap against the goal; this file is the
standing reference for the field and should be updated as the field moves.

Covers NumPy, JAX, Julia, MATLAB, Fortran/HPC (gfortran, Intel MKL, SPRNG)
and randompack, the closest R competitor. Ideas are sorted by whether they
must land before the stream freeze (they change the API surface or the
output) or can be added afterwards.

Facts marked *verified* were checked against the project's own README or
paper on the survey date; the rest is working knowledge of those libraries.
Sources are listed at the end.

## 1. The landscape in one table

| ecosystem | default engine | Gaussian method | parallel-stream model | notable |
|---|---|---|---|---|
| NumPy `Generator` | PCG64 (DXSM recommended) | ziggurat | `SeedSequence.spawn()` tree | `out=` in-place fill, `dtype=float32`, BitGenerator C API |
| JAX | threefry2x32 | **inversion** (`erfinv`) | `split()` / `fold_in()` | same key model as zurand; sharding-invariant output |
| Julia `Random` | xoshiro256++ | ziggurat | per-task splitting (DotMix), reproducible regardless of scheduling (*verified*) | `rand!`/`randn!` in place; VectorizedRNG.jl: SIMD **Box-Muller**, 2.1x Base `randn!` on AVX2, 8 xoshiro lanes (*verified*) |
| MATLAB `RandStream` | MT19937 | ziggurat, **user-selectable** `NormalTransform` = Ziggurat / Polar / Inversion | `Substream`, `StreamIndex`; philox and threefry engines | `Antithetic`, `FullPrecision` properties |
| Fortran / HPC | gfortran: xoshiro256**; MKL VSL | MKL: vectorised Box-Muller and **ICDF** | skip-ahead, leapfrog, SPRNG parameterised streams | MKL ARS5 (AES-NI) and Philox engines; Sobol |
| randompack (R, CRAN) | 13 engines: xoshiro256++, PCG64, Philox, ChaCha20, **squares64**, SIMD xoshiro/sfc64 with AVX2/AVX-512/**NEON** (*verified*) | not stated | jump-ahead, spawn keys, seed_seq_fe128 | bit-identical U/N/Exp/int across platforms **and languages** (R, Python, Julia, Fortran); 14 distributions; permutations; float and double |

Two things stand out. Every SIMD-first Gaussian in the field (JAX, MKL,
VectorizedRNG) abandons the ziggurat for a branch-free transform, so the
roadmap's "research problem" has three shipped answers. And randompack is
the closest competitor on positioning, not dqrng: it already ships
cross-language identical streams, NEON, and a squares64 counter-based
engine. zurand's benchmark ran randompack's `philox` engine, not its SIMD
xoshiro one; re-run the field with randompack's fastest engine before
quoting ratios against it.

## 2. Before the freeze (API surface or output)

| # | idea | from | why now |
|--:|---|---|---|
| A1 | **`offset` argument on every sampler**: `rng_uniform(key, n, offset = 0)` draws positions `offset:(offset+n-1)`. | NumPy `advance()`, MKL skip-ahead, MATLAB `Substream` | The API has no way to draw the *next* n values after a batch; users must fold by batch index. For the counter engines this is O(1) and is the one property philox keeps that nothing exposes; for xoshiro it costs at most 511 recurrence steps to align. Signature change, so before 0.1.0. |
| A2 | **Vector `mean`/`sd`/`min`/`max`, recycled to `n`**. | base R `rnorm(n, mu_vec, sd_vec)`, MATLAB `normrnd(mu, sigma)` arrays, NumPy broadcasting | The scaling pass already exists; it just reads scalars. Without this zurand is not a drop-in for the hierarchical-model idiom. Also interacts with 1.3 (fused scaling): design the fused transform to take per-element scale from the start. |
| A3 | **A stream-format version on keys**, next to `engine`. | NumPy's `RandomState` vs `Generator` stability split; StableRNGs.jl | Costs nothing now; if a bug ever forces a stream change after 0.1.0, old keys keep producing old values. Cannot be added later without breaking every saved key. |
| A4 | **Whole-stream digest test for `rng_normal()`** across CI platforms (hash of `rng_normal(key, 1e7)` bit patterns). | Julia gets cross-platform identical `randn` because it ships its own libm | The ziggurat tail calls the platform `log1p()` and the ambiguous wedge band calls `exp()`. Eight pinned tail values do not show that glibc, musl, macOS libm and MinGW agree on the ~2,600 tail draws in a 1e7 fill. If the digest disagrees anywhere, vendor deterministic `exp`/`log1p` (FDLIBM/openlibm, as musl and Julia do). This is the one item here that could invalidate the core promise, and it is a 30 ms test. |
| A5 | `rng_split(key, n)` as a named alias for folding by `1:n`; `rng_key()` with no seed drawing from OS entropy. | JAX `split`, NumPy `default_rng()` | Sugar, but users arriving from JAX and NumPy look for exactly these names. Optional. |

## 3. After the freeze (additive)

| # | idea | from | notes |
|--:|---|---|---|
| B1 | **`rng_normal(method = "inversion")`**: vectorised Wichura AS241 (the algorithm behind R's `qnorm`) on the 52-bit uniform. | MATLAB `NormalTransform='Inversion'`, JAX, MKL ICDF | The central branch (85% of draws) is a rational polynomial with no branches or gathers, so it auto-vectorises under a `target("avx2")` attribute with no intrinsics. Gives a monotone map from uniform to normal: common random numbers, antithetic variates and quasi-Monte Carlo all need this, and none work with a ziggurat. Also a second, independent check of the ziggurat in the audit. Tail branches use `log`/`sqrt`, so vendor them (item A4 above). Own purpose value, so additive. |
| B2 | **SIMD Box-Muller** as the pure speed play. | VectorizedRNG.jl 2.1x on AVX2; MKL; cuRAND | Needs vectorised `log` and `sincos` polynomials (SLEEF-style), which are deterministic by construction because they never touch libm. Higher effort than B1 and a different stream, so a third method name, not a replacement. Do B1 first; do B2 only if the AVX2 ziggurat prototype in the review's §3.4 loses. |
| B3 | **squares64 as a counter-based engine.** | Widynski 2020-22; randompack ships it; PractRand to 32 TB reported by its author | Keeps `value = f(key, index)` at a fraction of Philox's cost (rounds of 64-bit squaring, no 128-bit multiply). If a microbenchmark shows it within 10% of xoshiro in the fill, the default could stay counter-based instead of switching to the hybrid, which is a cleaner story. Measure before the freeze, since it bears on the default-engine decision; ship as opt-in regardless until it has independent scrutiny. |
| B4 | **ARS (AES-based) engine with runtime dispatch**: AES-NI on x86, ARMv8 crypto extensions on Apple Silicon, software AES fallback for identical output elsewhere. | MKL ARS5; zurand's own 3.4x-vs-philox measurement | The fastest counter-based option measured. Only worth it if B3 disappoints. |
| B5 | **`RNGkind("user-supplied")` bridge**: a buffered adapter that installs zurand behind base `runif()`, `rnorm()` and every `r*` function via R's `user_unif_rand`/`user_norm_rand` hooks. | precedent: randtoolbox | Stateful by construction, so opt-in and clearly labelled. It is the adoption on-ramp for users who will not change code: base `rnorm` at 24 M/s would go through the ziggurat instead of inversion. |
| B6 | **`rng_permutation(key, n)` and `rng_sample()`.** | dqrng `dqsample` (its most used function), NumPy, Julia, MATLAB, randompack | Stateless and parallel: generate 64-bit keys, sort. Exactly uniform, O(n log n). A Kensler-style keyed bijection gives O(1) per element but is not exactly uniform; document or avoid. |
| B7 | **`rng_exponential()`**: NumPy's exponential ziggurat tables are already vendored. | everyone | ~40 lines mirroring the normal sampler. |
| B8 | **Published stream spec plus test vectors** (`dev/kat/*.py` is most of a Python reference already). | randompack's cross-language claim | Cheap, but no longer a differentiator; randompack has it. Do it for the documentation value, not as a headline. |

## 4. Deliberately not

- **Fourteen distributions.** That is randompack's turf and the scope
  creep would dilute the one claim zurand can own: fastest uniform and
  Gaussian, stateless, thread-invariant.
- **float32 output.** R has no native type; only relevant to the C API,
  where a float fill halves memory traffic for users of the `float`
  package. Add there if asked, never at the R level.
- **A stateful default API.** B5 is a bridge, not the interface.

## Sources

Read on the survey date:

- VectorizedRNG.jl README: <https://github.com/JuliaSIMD/VectorizedRNG.jl>
  (SIMD Box-Muller, eight xoshiro256++ lanes, 2.1x Base `randn!` on AVX2)
- randompack README: <https://github.com/jonasson2/randompack>; paper
  "Randompack: Cross-Platform Reproducible Random Number Generation and
  Distribution Sampling", arXiv:2605.05099
- Widynski, "Squares: A Fast Counter-Based RNG", arXiv:2004.06278
- Julia `Random` documentation: <https://docs.julialang.org/en/v1/stdlib/Random/>

Not read (blocked from the review environment), so the corresponding rows
are working knowledge and worth re-checking before being cited:

- MathWorks, `RandStream` reference (NormalTransform, Antithetic,
  FullPrecision, Substream, engine list)
- NumPy, `numpy.random.Generator` reference (`out=`, `dtype`,
  `SeedSequence.spawn`)
- Intel MKL VSL documentation (Gaussian methods, skip-ahead, leapfrog,
  ARS5 and Philox basic generators)
