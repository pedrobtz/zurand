# zurand roadmap

**Goal:** the fastest uniform and Gaussian random vectors in R, where every
value is a pure function of `(key, position)`, so output is bit-identical
across thread counts, SIMD paths, platforms and call order.

- `dev/design.md` says **what** zurand is at v0.1.0: the stream spec to be
  frozen, the R and C surface, the reproducibility contract, and the rule
  deciding which changes may land when.
- This file says **in what order**, with an exit criterion per milestone.
- `dev/log.md` keeps the history: the original phases, both amendments,
  the memory-traffic and SIMD measurements, and the adopted revision.
- `dev/benchmark-baseline.md` holds the numbers; `tools/benchmark.R`
  regenerates them. `dev/review-2026-09-23.md` and
  `dev/ecosystem-survey.md` are the evidence behind this revision.

Rewritten 2026-09-25 from the review. The old principle "the default
stream never changes" was adopted at 0.0.0.9000 with no users, and it cost
the package its headline, because the default was the slowest engine. It
is replaced by one output-changing window that closes at v0.1.0
(design §2).

---

## Measurement protocol

Unchanged, and binding on every performance claim:

```
R CMD INSTALL . && Rscript tools/benchmark.R
```

- **`-O2` build only.** `devtools::load_all()` defaults to `debug = TRUE`
  (`-O0`, 11x slower). Use `R CMD INSTALL` or `load_all(debug = FALSE)`
  after `devtools::clean_dll()`. `debug` is read only when something
  recompiles, so a warm debug build is otherwise reused silently.
- **Ratios, not absolutes.** Measure zurand and the reference in the same
  `bench::mark()` call, over several rounds. This machine's absolutes move
  13% under background load.
- **A platform header on every result**, or the result is discarded.
- **Metrics, in priority order:** Gaussian and uniform at one thread and
  n = 1e7; the same at all threads; per-call overhead at n in {1, 100,
  1000}, which must not regress.
- **Field:** base R, dqrng, RcppZiggurat (MT, LZLLV, GSL, QL), randompack
  on its fastest engine, sitmo, rTRNG.
- **Acceptance:** the ratio interval separates from the baseline's, no
  other metric regresses, and output is bit-identical unless the change is
  listed in milestone A.

## Status (2026-09-26)

n = 1e7, M/s, from the `benchmark` workflow on GitHub runners, each package
in the same `bench::mark()` call:

| | EPYC 7763, 1 thread | EPYC, 4 threads | Apple M1, 1 thread |
|---|---:|---:|---:|
| uniform: zurand `xoshiro256pp` | **864** | **1801** | **1364** |
| uniform: randompack `x256++simd` | 688 | 684 | 1087 |
| uniform: dqrng | 368 | 367 | 230 |
| Gaussian: zurand `xoshiro256pp` | **366** | **919** | **624** |
| Gaussian: randompack `x256++simd` | 355 | 355 | 503 |
| Gaussian: RcppZiggurat MT | 125 | 124 | 193 |

After the fused AVX2 uniform (#34), the NEON path (#35, #38), `DC ZVA` on
Apple (#39), the signed-table ziggurat (#40), huge pages (#42, no effect
on these runners, whose kernels use them anyway) and the slow path's
codegen (#43). zurand leads randompack on every row; x86 Gaussian by 3%
on this Zen 3 and 8% on a Xeon 6973P (the A/B behind #43). The M1 runner
is a noisy VM: the same code measured 1837/661 M/s in the morning and
1364/624 here; same-machine A/Bs on three M1 runners gave 1.05-1.15x
(uniform) and 1.05-1.20x (Gaussian). Both packages' normals are bounded
by the scalar ziggurat's table lookups; the ziggurat implementation is
frozen (below), its compiler attributes are not.

Done: KAT against Random123 for both counter engines; golden values in hex
floats; threads = serial and SIMD = scalar identity tests; two-pass uniform
fill; the threefry and xoshiro engines; AVX2 dispatch; the audit tooling;
the whole-field benchmark; the full sanitizer/valgrind/rchk/LTO/gctorture
CI. As of PR #17 the golden tests also run, and pass, on 32-bit i386 (x87)
and musl.

---

## Milestone A -- stream freeze

Everything that changes existing output, done once, deliberately. The
order matters: evidence first, then the changes, then one regeneration of
the golden files.

| # | task | why / evidence | verify |
|--:|---|---|---|
| A1 | ~~Add an `ubuntu-24.04-arm` leg to R-CMD-check~~ **DONE** #19 | GCC contracts `a + b * c` into an FMA by default on aarch64 and ignores the `FP_CONTRACT` pragma, so the golden tests are expected to **fail** there today. No leg covers GCC on aarch64 | the leg runs and its failure (or pass) is recorded |
| A2 | ~~Forbid contraction on GCC~~ **DONE** #19, with a code-level barrier (design §3.4 item 4) | a bug fix: it restores the contract on builds that break it today | A1 green |
| A3 | ~~Whole-stream digest test~~ **DONE** #20: a hash of `rng_normal(key, 1e7)` (and uniform, integer) bit patterns, checked on every CI leg | the tail calls libm `log1p()` and the wedge fallback calls `exp()`; eight pinned tail values do not show that glibc, musl, macOS and UCRT agree on the ~thousands of libm calls in a 1e7 fill | digest identical on every leg, i386 included |
| A4 | ~~Deterministic `exp`/`log1p`~~ **DONE** #20: A3 found Apple/Windows against glibc/musl; fdlibm port in `src/zurand_fdlibm.h` | deterministic C under the no-contraction rule, so identical everywhere | A3 green |
| A5 | ~~`rng_integer()` signed overflow~~ **DONE** #21: `(int)((int64_t)min + offset)`, plus a golden case at `min = -2e9, max = 2e9` | UB for ranges above 2^31 (UBSan-verified in the review); output-neutral | sanitizer leg covers the new case |
| A6 | **Engine tag in counter word 3** for every Philox call made on behalf of `xoshiro256pp` (design §3.2) | today xoshiro's stream for a seed is a function of philox's (verified: first word = `rotl(s0+s3,23)+s0` of philox block 0) | new xoshiro KAT; philox and threefry KATs unchanged |
| A7 | **`stream = 1L` attribute on keys**; samplers reject unknown versions; a missing attribute means 1 | makes a post-release fix possible without breaking saved keys | tests for missing, 1, and 2 (error) |
| A8 | ~~Audit `xoshiro256pp`~~ **DONE** #24: clean to 1 TB (bits, uniform, both cross-key modes) and 512 GB (normal). Originally: add it to `statistical-audit.yml`'s engine choice, add the two cross-key modes (interleaved `rng_key(s, 64)` and interleaved `rng_fold(key, 1:64)`), and run at least 1 TB of bits and 256 GB of `pnorm(rng_normal())` through PractRand | reseeding xoshiro from Philox every 512 words is a novel construction, and 512 MB is a smoke test; the package is used as many short keyed streams, so cross-key correlation is the likeliest weakness | results recorded in `dev/statistical-audit.md`; nothing worse than "unusual" |
| A9 | ~~Make `xoshiro256pp` the default~~ **DONE** in `rng_key()` and `rng_key_from_r()`. Rewrite `?rng_key` (it still says "Both are Random123 generators" above three engines) and DESCRIPTION's Description (it names only Philox) | conditions met: A6 done, A8 passed, and no counter-based engine within 10% (squares64 measured at 0.62x, design §3.4) | golden files regenerated **once**, in a commit that says so |
| A10 | **Tag `stream-1`** | the freeze itself | tag on the commit after A9 |

**Exit criterion:** every CI leg green, including i386, musl and GCC on
aarch64, with the regenerated golden values and the digest. The audit is
recorded, and `stream-1` is tagged. After this, design §2 applies without
exception.

---

## Milestone B -- v0.1.0 on CRAN

No output changes. The work is adoption and packaging.

| # | task | why |
|--:|---|---|
| B1 | ~~C API~~ **DONE** (design §5): `inst/include/zurand.h`, `R_RegisterCCallable()` for reentrant fill functions, R samplers reimplemented as thin wrappers over them | the in-place fill PR #4 attempted, delivered where it is safe; lets simulation, MCMC and bootstrap packages draw from their own threads. Test it by building a tiny fixture package in CI that `LinkingTo`s zurand and compares against the R functions |
| B2 | **Docs**: state `rng_uniform()`'s interval the way `?runif` does; README benchmark table with a platform header, one thread and all threads, randompack on its fastest engine; vignettes "why stateless" and "performance"; NEWS.md; pkgdown | the case for the package currently lives in `dev/` |
| B3 | ~~Threading on CRAN~~ **DONE**: confirm the tests stay within CRAN's two-core limit (cap with `rng_threads(2L)` in `tests/testthat/setup.R` if `OMP_THREAD_LIMIT` is not honoured), and say plainly in the README that CRAN's macOS binaries have no OpenMP | a common reason for a CRAN bounce; and on macOS the one-thread number is the number |
| B4 | ~~Release hygiene~~ **DONE**: version 0.1.0; `cran-comments.md` explaining the `-Wunused-const-variable` pragma and any FP-contract mechanism from A2; the upstream Random123 URL instead of the local path in CLAUDE.md; a pass with the `cran-extrachecks` skill | the old Phase 4 items |
| B5 | Submit | |

**Exit criterion:** zurand 0.1.0 is on CRAN.

---

## Milestone C -- after release (neutral or additive only)

Ordered by value per effort. Each item either leaves every bit unchanged
or claims a new argument, purpose value or engine name.

| # | task | kind | note |
|--:|---|---|---|
| C1 | `offset = 0` on every sampler (and in the C API from B1) | additive | counter engines: O(1); xoshiro: at most 511 steps |
| C2 | Vector `mean`/`sd`/`min`/`max`, recycled | additive | same per-element formula as the scalar case |
| C3 | Fuse the `mean`/`sd` and `min`/`max` scaling into the per-chunk transform; ~~packed `{ki, wi}` ziggurat table~~ | neutral | old item 1.3; the packed table (old 1.2) was remeasured at the R level on x86 and M1 and gave nothing -- both tables sit in L1 -- so it is dropped with the ziggurat freeze below |
| C4 | ~~NEON path for `xoshiro256pp`~~ **DONE** #35: 8 sub-chunks in four 2-lane registers; M1 uniform 632 -> 1208 M/s, level with randompack | neutral | **highest priority after release**: on Apple M1 randompack's SIMD engine is 2x zurand on uniform (1309 vs 638 M/s) and CRAN's macOS binary has no threads to make up for it; the M-series baseline now exists (2026-09-26 CI benchmark) |
| C5 | `rng_exponential()` | additive, purpose 5 | NumPy's exponential tables are already vendored |
| C6 | ~~AVX2 ziggurat fast path~~ **TRIED, NOT MERGED** #36: gathers 1.8x slower on Downfall-patched Intel; plain loads +0-9% xoshiro, philox -6..+7% across four CPUs; fails the no-regression rule | neutral | closed by the ziggurat freeze below: in isolation every SIMD fast path was slower than the scalar loop |
| C7 | `rng_normal(method = "inversion")` | additive, purpose 6 | monotone in `u`, for CRN, antithetics and QMC |
| C8 | `rng_permutation()`, `rng_sample()` | additive | dqrng's most-used function |
| C9 | Split `src/zurand.c` along engine and sampler lines | neutral | 1,050 lines and three engines |
| C10 | Make `src/zigbounds.h` platform-independent (exact rationals or `Rmpfr`) | build only | old 0.4 |
| C11 | `RNGkind("user-supplied")` bridge | additive, opt-in | only if users ask; stateful by nature |
| C13 | ~~`rng_normal(method = "mcfarland")`~~ **DONE**: McFarland's modified ziggurat, purpose 7; tables from `tools/generate-mcfarland-tables.py` (mpmath) | additive | transform -21..-39% per value on eight CPUs (`dev/simd/mcfarland.c`); R level, one thread, xoshiro, 1e7: 1.17x the default on Zen 3, 1.25x Zen 4, 1.10x Neoverse N2, 1.05-1.07x M1; 1.2-1.7x randompack; 1e9-draw chi-square clean on two engines |
| C14 | ~~`rng_normal(method = "boxmuller")`, vectorised Box-Muller (VectorizedRNG.jl)~~ **TRIED, NOT ADDED** | additive | deterministic version built (`dev/simd/boxmuller.c`: fdlibm log, fitted sine polynomial, no FMA; AVX2/AVX-512/NEON bit-identical to scalar), ns per normal with words in cache vs McFarland: Xeon 8573C AVX-512 2.41 vs 0.61, Zen 3/4 AVX2 3.45 vs 0.80-0.85, M1 NEON 5.69 vs 0.53, Neoverse N2 NEON 7.96 vs 0.63. With FMA, as VectorizedRNG uses: 1.81 / 2.76-2.92 / 4.27 / 8.22 -- 0.97-1.33x faster, still 3-13x slower, and 7% of values change by up to 9e-16, i.e. no cross-platform identity. Divide, sqrt and two polynomials cost more than a table lookup; VectorizedRNG's 2x is against Julia's buffered scalar randn |
| C12 | AVX-512 path for `xoshiro256pp` (8 lanes, `vprolq`) | neutral | randompack has one; +12-15% uniform on Xeon 6973P, Zen 4 inconsistent in one sample (`dev/simd/avx512.c`); needs three Zen 4 and three Intel samples before deciding |

---

### Also done after the freeze

- Fused uniform conversion in the AVX2 path (#34): EPYC one thread
  569 -> 887 M/s, four threads 1370 -> 1928. Non-temporal stores on top
  were measured and rejected (slower at the R level on Intel and AMD).
- `dev/simd/ab.R`: same-machine A/B of two builds. Runner-to-runner
  comparisons are not enough, since GitHub hands out Intel and AMD CPUs at
  random.

### Ziggurat implementation frozen (2026-09-26)

The normal sampler's *implementation*, not only its stream, is frozen at
#40. Every angle below was measured bit-identical to the shipped output;
only the two in bold were merged. Reopen only for a different method
(say, inversion or Box-Muller with deterministic vector `log`/`sin`,
which is additive: C7), never to tune this one again.

| angle | result | where |
|---|---|---|
| **signed `wi` table + `rabs - 1 < ki - 1` accept test** | -15% per draw on x86 (four instructions fewer) | #40, N5 |
| **fast-path loop unrolled by 4** | -18-22% per draw on M1, level on x86 | #40, N6 |
| SIMD fast path: AVX2 gathers, AVX2 plain loads, NEON lane loads | slower than scalar on every CPU; gathers 1.8x slower under Intel's Downfall microcode | #36, `dev/simd/zig_compare.c` N3 |
| branch-free with a fix-up pass, plain and unrolled | 1.4-1.6x slower: the ~1% reject branch is nearly free | N1, N2 |
| randompack's structure (words in the output, converted in place, backwards) | slower on x86 and M1 | N4 |
| packed `{ki, wi}` table | no gain: both tables sit in L1 | microbench |
| BMI2 field extraction (`bextr`) | inconsistent across CPUs | N7 |
| batched slow path: rejects recorded, their retry Philox/Threefry blocks computed four at a time | i5, 10 alternating rounds, one thread: xoshiro +2.3-2.7%, philox +0.7%, **threefry -1.6-2.6% (0 of 20 rounds faster)**; M1 +2-4% at 1e7, mixed at 1e6; x86 runners within noise. Capped by the ~1.15% reject rate; fails the no-regression rule | not merged |
| `DC ZVA` before large normal fills (Apple) | level or slower on three M1 runners (1e7: 1.05-1.13 -> 0.87-1.02 of randompack). Normals are one sequential store stream per chunk, which Apple cores already write without reading; the uniform path's ten interleaved 4 KiB blocks are what `DC ZVA` fixes (#39) | not merged |

What remains is the table lookups, the same in every ziggurat, and memory.

### x86 analysis, second pass (2026-09-26)

Asked: with the ziggurat frozen, is anything left to beat randompack
comfortably on every architecture, x86 in particular? Everything below was
measured, bit-identical, on the machines named; harnesses in `dev/simd/`.

| idea | result | outcome |
|---|---|---|
| **Transparent huge pages** for large outputs (`madvise(MADV_HUGEPAGE)`, Linux) | a fresh 80 MB vector on 4 KiB pages costs 2.6-3.9 ns/value in faults, more than the generation; on the ubuntu-24.04-arm runner (THP `madvise`, the Ubuntu default) the hint took uniform 268 -> 532 M/s and normal 214 -> 344, against randompack's 271 / 200; in `always` mode the explicit request costs 6-9% (compaction), so it is made only in `madvise` mode | **#42** |
| **Slow path codegen**: fdlibm `exp`/`log1p` as calls, no `cold` | the slow path compiled badly, not computed slowly: stubs of every kind "gained" 15-19% on the i5, and so did simply not inlining `zurand_exp()`; `exp()` runs on ~6% of wedge tests (`zig_slow.c`); GCC compiles `cold` for size (2-5% on Zen 3) | **#43**: +7-9% normals on Xeon 6973P, +15-19% with Apple clang, level on Neoverse |
| fused generator + transform (no word buffer) | slower everywhere: i5 2.04 vs 1.81 ns, Zen 3 2.82 vs 1.75, Zen 4 2.99 vs 1.53, Xeon 2.01 vs 1.42 (`fuse_normal.c`) | rejected |
| strip-mined 2 KiB word buffer | level to 10% slower | rejected |
| non-temporal stores on fresh memory (128-bit, R's alignment) | slower everywhere, 1.3-3.3x on fresh pages (`fuse_normal.c` UN) | rejected, closes the earlier question |
| `MADV_POPULATE_WRITE` pre-faulting | no gain (`pagefault.c`) | rejected |
| BMI2 `bextr` field extraction | +3-7% on the transform on Zen 3/4 only, level on Intel | not adopted (AMD-only, needs dispatch) |
| batched retry blocks | +2.5% xoshiro, -2% threefry, measured before #43 | not adopted |
| **AVX-512 generator** (8 lanes, `vprolq`, 8x8 transpose; `avx512.c`) | Intel Xeon 6973P: uniform +12-15%, normal +2%; Zen 4 (one sample): uniform slower in cache (0.65 vs 0.51 ns), normal faster (0.96 vs 1.27) -- inconsistent, needs more samples. randompack ships an AVX-512 path (`-mavx512f -mavx512dq` via configure) | **open**, C12 |
| other Gaussian methods (vectorised Box-Muller or inversion) | ~5-7 cycles/value on AVX2 against the scalar ziggurat's 3.7; only AVX-512 with full-width units comes close; would need deterministic vector `log`/`sincos` and a new purpose value | not pursued |
| where the time goes now, one thread, EPYC 7763, n = 1e7 (`fuse_normal.c`, `pagefault.c`) | words 0.54 ns, uniform conversion+store 0.13, ziggurat transform 1.03, fresh pages 0.35-0.5 (2 MiB) or 2.6-3.9 (4 KiB); the R level matches the harness for uniform exactly | -- |

Note on the i5-8500B: small (<10%) effects there are unreliable until
compared on a second machine -- the three slow-path stubs all measured
+17% for reasons that turned out to be codegen, not the work removed.

## Retired

| item | why |
|---|---|
| 1.4, 1.5 (M1 chunk sweep and ILP retest) | folded into C4; tuning for one machine is not a milestone |
| 2.3 `philox4x32-10`, 2.4 AVX-512 Philox | the default no longer runs on Philox after A9 (Amendment 2) |
| 3.1 R-level in-place fill | PR #4 closed for the copy-on-write hazard; replaced by B1 |
| 3.2 threading headline | now part of B2 and B3 |
| 3.3 `configure` for OpenMP | needed only if A2 picks the `configure` route; otherwise open as a robustness item |
| squares64 and ARS as engines | squares64: 0.62x scalar xoshiro in the fill and no AVX2 path; ARS does not compile without AES-NI and has no NEON path |

## Decisions log

| date | decision | rationale |
|---|---|---|
| 2026-09-19 | API keeps the `rng_*` prefix through the rename | the prefix is the API's, not the package's |
| 2026-09-19 | `src/zigbounds.h` kept as generated on the M1 | regeneration drifts by libm ulps; the guard absorbs it |
| 2026-09-19 | `-DR123_USE_MULHILO64_C99=1` stays in `Makevars` | harmless on 64-bit; the 32-bit fallback, exercised by `arch` |
| 2026-09-19 | `philox4x64-7` rejected | 0.95x, slower than 10 rounds |
| 2026-09-19 | no R-level in-place fill (PR #4 closed) | a copy-on-write hazard in a permanent public API |
| 2026-09-19 | SIMD vectorises across sub-chunks, never within one | same bits as scalar, so no stream depends on the ISA |
| 2026-09-25 | **"the default stream never changes" replaced by one output-changing window, closing at v0.1.0** | the principle predated any users and cost the default its speed; design §2 |
| 2026-09-25 | squares64 not adopted | 0.62x scalar xoshiro in the two-pass fill; no AVX2 path |
| 2026-09-26 | no fused multiply-add, enforced in code, not by compiler flag | flags can be overridden by user `CFLAGS` and draw CRAN notes; #19 |
| 2026-09-26 | fdlibm `exp`/`log1p` instead of the platform libm | platforms disagreed; deterministic beats correctly rounded here; #20 |
| 2026-09-26 | 32-bit x87 outside the reproducibility contract | double rounding on every add; no CRAN platform since R 4.2.0 |
| 2026-09-26 | **default engine becomes `xoshiro256pp`** (was proposed 2026-09-25) | 2.4x dqrng on uniform versus 0.96x today; conditional on A6 and A8 |
| 2026-09-26 | huge pages requested only in THP `madvise` mode, never in `always` | 2x on large fills where it applies; 6-9% compaction cost where it does not; #42 |
| 2026-09-26 | the ziggurat's slow path stays plain `noinline`, never `cold`; fdlibm `exp`/`log1p` stay calls | attributes are not the algorithm; #43 |
| 2026-09-26 | **ziggurat implementation frozen at #40** | every remaining angle measured level or worse on some engine or CPU; see "Ziggurat implementation frozen" |

## Open questions

- C API shape: static-inline shims in the header over `R_GetCCallable()`
  (dqrng's pattern) or a struct of function pointers behind one
  `zurand_api()` call? The latter versions more cleanly.
- Should C11 exist at all? Decide from user requests after release.

## Non-goals

- `-march=native` or any machine-specific output.
- Changing any `(engine, stream)` output after v0.1.0.
- float32 at the R level; distributions beyond uniform, normal,
  exponential and integer; a stateful default API.
