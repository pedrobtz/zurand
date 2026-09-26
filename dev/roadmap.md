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

## Status (2026-09-25)

x86_64 i5-8500B, one thread, n = 1e7, from `dev/benchmark-baseline.md`:

| engine | uniform | Gaussian | note |
|---|---:|---:|---|
| `xoshiro256pp` + AVX2 | **776 M/s**, 2.4x dqrng | **376 M/s**, 1.9x RcppZiggurat MT | fastest in the field on both |
| `philox4x64` (current default) | 314, 0.96x dqrng | 215, 1.09x | |
| `threefry4x64` | ~1.1x philox | ~1.15x philox | |

No M-series baseline exists yet, and no all-threads baseline has been
published.

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
| B4 | **Release hygiene**: version 0.1.0; `cran-comments.md` explaining the `-Wunused-const-variable` pragma and any FP-contract mechanism from A2; the upstream Random123 URL instead of the local path in CLAUDE.md; a pass with the `cran-extrachecks` skill | the old Phase 4 items |
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
| C3 | Fuse the `mean`/`sd` and `min`/`max` scaling into the per-chunk transform; packed `{ki, wi}` ziggurat table | neutral | old items 1.3 and 1.2; 1.2 measured 1.10x on the transform |
| C4 | NEON path for `xoshiro256pp` (2 sub-chunks per register) and an M-series baseline from a `macos-14` runner | neutral | old 5.2 and the missing primary-machine baseline; matters more to Mac users than threads do |
| C5 | `rng_exponential()` | additive, purpose 5 | NumPy's exponential tables are already vendored |
| C6 | AVX2 ziggurat fast path: two gathers, a compare, a scalar loop for rejected lanes | neutral | land only if the ratio interval separates; realistic target 500 M/s |
| C7 | `rng_normal(method = "inversion")` | additive, purpose 6 | monotone in `u`, for CRN, antithetics and QMC |
| C8 | `rng_permutation()`, `rng_sample()` | additive | dqrng's most-used function |
| C9 | Split `src/zurand.c` along engine and sampler lines | neutral | 1,050 lines and three engines |
| C10 | Make `src/zigbounds.h` platform-independent (exact rationals or `Rmpfr`) | build only | old 0.4 |
| C11 | `RNGkind("user-supplied")` bridge | additive, opt-in | only if users ask; stateful by nature |

---

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
