# zurand roadmap

**Goal:** the fastest Gaussian and uniform random *vector* generators in the
R ecosystem, without giving up the property that makes zurand different --
every value is a pure function of `(key, index)`, so output is bit-identical
across thread counts, platforms and call order.

This document is the plan. `dev/benchmark-baseline.md` holds the numbers it
is measured against; `tools/benchmark.R` regenerates them.

---

## 0. What "fastest" means here

A claim nobody can check is not a claim. Every number in this project is
produced by one protocol:

```
R CMD INSTALL . && Rscript tools/benchmark.R
```

- **Installed build only.** `load_all()` is a debug build and understates
  zurand by a wide margin.
- **Primary machine: the Apple Silicon M1** the package is tuned on. The
  x86_64 Intel machine is a secondary datapoint. The script prints a
  platform header; a result without one is discarded.
- **Four metrics**, in priority order:
  1. Gaussian, single thread, n = 1e7 -- algorithm vs algorithm.
  2. Uniform, single thread, n = 1e7.
  3. Gaussian and uniform, all threads, n = 1e7 -- what a laptop user gets.
  4. Per-call overhead at n in {1, 100, 1000} -- already at parity with
     base R and ahead of dqrng; must not regress.
- **Competitors:** base R, dqrng, RcppZiggurat (MT, LZLLV, GSL, QL),
  randompack, sitmo, rTRNG. Add any package that enters the field.

**Acceptance for a change:** at least +2% on the metric it targets, no
regression >1% on any other, and `rng_*()` output bit-identical to before
unless the change is a new opt-in engine (see §3).

### Targets (proposed, revise after the first M1 run)

| metric | today (x86_64) | target |
|---|---|---|
| Gaussian, 1 thread | 213 M/s, 1.07x next best | **>= 1.25x next best** |
| Uniform, 1 thread | 189 M/s, 0.57x best (dqrng) | **within 10% of best** |
| Gaussian, all threads | unmeasured locally | **>= 3x any competitor** |
| n = 1000 overhead | 9.4 us (0.75x dqrng) | no regression |

---

## 1. Where we are

Measured 2026-09-19 on x86_64 (Intel i5-8500B, no OpenMP). See the
baseline document for the full tables.

**Gaussian already leads, narrowly.** 213.5 M/s vs RcppZiggurat MT 199.1,
randompack 184.3, dqrng 148.8, base R 22.8.

**Uniform is third.** 189.1 M/s vs dqrng 333.1, randompack 255.2.

**The cost split is the whole story:**

| component | ns/value | share |
|---|--:|--:|
| Philox4x64-10 (one 64-bit word per double) | 2.74 | 58% |
| ziggurat transform + chunk buffer | 1.94 | 42% |

Philox's **365 Mword/s is a hard ceiling** for any path built on it. That
ceiling is barely above dqrng's *current* uniform throughput, and it means
Gaussian cannot exceed ~1.7x today's number even if the transform were free.

**Conclusion:** the sampler side (ziggurat, uniform conversion, Lemire
integers) is close to optimal for scalar code. The generator is the lever.
Everything in §2 buys single-digit percentages; §3 is where the decisive
gains are.

---

## 2. Guiding principles

1. **The default stream never changes.** A key created today must produce
   the same values forever. Faster generators are *new engines*, opted into
   via `rng_key(engine = ...)`. The abstraction exists already -- `engine`
   is validated by `match.arg()` in R and `check_engine()` in C and stored
   on every key -- it simply has one implementation.
2. **Measure on the M1 before deciding.** Microarchitectural conclusions do
   not transfer between the two machines (128- vs 64-byte cache lines,
   `mul`+`umulh` vs `mulq`, 31 vs 16 registers, OpenMP present vs absent).
3. **Every speed change ships with its benchmark delta** in the commit
   message, from the protocol above.
4. **Every speed change is bit-identity tested.** A golden-value test
   (§Phase 0) is the gate. If it fails, either the change is wrong or it is
   a new engine and must be labelled as one.
5. **Vendored code stays verbatim.** Derived tables are generated from it
   at build time, never edited in place.

---

## Phase 0 -- Trust (before any optimisation)

Speed claims invite scrutiny. The suite currently proves the generator is
*deterministic*, not that it is *correct*. Close that first so every later
change has a gate.

| # | task | why | verify | effort |
|--:|---|---|---|---|
| 0.1 | **Known-answer test against Random123's reference vectors.** Upstream ships `tests/kat_vectors` with philox4x64 entries. Compare `rng_bits(key, n, bits = 64L)` for the reference keys/counters. | Nothing today proves the Philox is Philox. | new test file; ~15 lines | 1 h |
| 0.2 | **Golden-value test for the full pipeline.** Hardcode ~256 values each of `rng_uniform`, `rng_normal`, `rng_integer` for a fixed key, including values known to hit the ziggurat wedge and tail paths. | CI runs on arm64, x86_64, i386 and musl; nothing asserts they agree. This is the reproducibility guarantee stateless generation exists to give, and the gate every later change passes through. | new test file | 2 h |
| 0.3 | **One-time statistical audit, recorded in `dev/`.** PractRand (or TestU01 SmallCrush) on `rng_bits()`; the same on `pnorm(rng_normal())` to exercise the custom wedge shortcut and tail. Not CI -- a document. | The uniform conversion and wedge brackets are custom code. One KS test at n=50k is a smoke test, not evidence. | `dev/statistical-audit.md` | 1 day |
| 0.4 | **Make `src/zigbounds.h` platform-independent.** Generated on the M1, it drifts by +/-1 low bit when regenerated on x86_64 (libm ulp differences). `ZURAND_ZIG_GUARD` absorbs it, so it is not a correctness bug, but "rerun the script" yields spurious diffs off the M1. Either compute the brackets with `Rmpfr`/exact rationals, or document "regenerate on arm64 only" in the header itself. | Reproducibility of the build, not of the output. | regenerate on both machines, `diff` empty | 2 h |

**Exit criterion:** 0.1 and 0.2 green on all CI legs. 0.3 and 0.4 can trail.

---

## Phase 1 -- Free wins (non-breaking, default stream unchanged)

Each is measured; each is small. Together perhaps +5% Gaussian, +20% uniform.

| # | task | expected | evidence | effort |
|--:|---|---|---|---|
| 1.1 | **Two-pass fill for uniform.** `fill_normal_column` buffers a chunk of Philox output then transforms it; `fill_uniform_column` still fuses generate-and-transform. | **~+20% uniform** | uniform costs 5.29 ns/value vs normal's 4.68 despite less arithmetic -- the fused loop cannot pipeline across Philox's 10-round chain | 1 h |
| 1.2 | **Packed ziggurat table.** `ki_double[idx]` and `wi_double[idx]` are separate arrays: two cache lines per draw. Derive a `{uint64_t ki; double wi;}` array from the vendored header at build time (a generated `src/zigtable.h`, like `zigbounds.h`). | **~+4% Gaussian** on x86; re-measure on M1's 128-byte lines | standalone microbenchmark: 1.910 -> 1.737 ns/draw, 1.10x | 2 h |
| 1.3 | **Fuse `mean`/`sd` into the transform.** `C_rng_normal` fills standard normals then sweeps the whole output again with `mean + sd * x`. The transform already ends in a multiply: `x = s * (wi[idx] * sd) + mean` is one FMA with `wi*sd` hoisted per call. | removes a full memory-bound pass whenever `mean`/`sd` are non-default | code reading; measure with `mean = 1, sd = 2` | 2 h |
| 1.4 | **`ZURAND_CHUNK_BLOCKS` sweep.** Currently 128 blocks (4 KiB). Try 32, 64, 256, 512 on the M1. | unknown, likely neutral | none yet | 30 min |
| 1.5 | **Re-test manual ILP on the M1.** Interleaving 2/4/8 independent Philox blocks was 0.79-0.87x on x86 -- register pressure with 16 GPRs. arm64 has 31. | unknown; do not carry the x86 conclusion across | `/tmp/philox_bench.c` from the baseline work; rerun | 30 min |

**Exit criterion:** baseline re-recorded on the M1 with 1.1-1.3 landed.

---

## Phase 2 -- The generator (decisive gains, via new engines)

The only route to the targets in §0. Each candidate is an *additional*
engine selected by `rng_key(engine = ...)`; `"philox4x64"` stays the
default and stays bit-identical.

| # | candidate | Philox cost | why it might win | why it might not | effort |
|--:|---|---|---|---|---|
| 2.1 | **`philox4x64-7`** -- same generator, 7 rounds | ~0.7x | Random123's authors report 7 rounds passes BigCrush; 10 is their safety margin. Simplest change, largest certain gain: **~+21% Gaussian, ~+40% uniform**. | Thinner margin above the failure threshold (BigCrush fails at 6). Needs its own Phase 0.3 audit. | 1 day incl. tests |
| 2.2 | **Threefry4x64** | ? | Add/rotate/xor only. Vectorises 4-wide on AVX2 and NEON, where Philox's 64x64->128 multiply cannot. | Slower than Philox *scalar*; the win exists only with SIMD. Measure before committing. | 2 days to prototype |
| 2.3 | **`philox4x32-10`** | ? | 32x32->64 multiplies vectorise (8 lanes AVX2, 4 lanes NEON). | A double needs two 32-bit words; per-double cost is unclear and the ziggurat wants a 64-bit word. | 2 days to prototype |
| 2.4 | **SIMD Philox4x64 on AVX-512** | ~2-4x on capable CPUs | `_mm512_mullo_epi64` + high-half emulation. | Neither dev machine has AVX-512; CRAN cannot ship `-march=native`. Runtime dispatch only. Park unless a target machine appears. | -- |

**Method:** prototype each as a standalone C microbenchmark first (as was
done for ILP and the table layout), on both machines. Only the winner
becomes an engine. Then: `engine_<name>.c`, the `engine` attribute drives
dispatch, KAT + golden values for the new engine, its own statistical
audit, and benchmark tables gain a row.

**Decision needed before starting:** does a 7-round engine's smaller
statistical margin fit a package that positions itself on reproducibility?
Recommendation: yes, *as opt-in* -- the user who chooses
`engine = "philox4x64-7"` is choosing throughput, and the documentation
says so.

---

## Phase 3 -- Beyond the vector

| # | task | why | effort |
|--:|---|---|---|
| 3.1 | **In-place fill API.** `rng_fill_normal(buf, key)` writing into a caller-supplied vector. | First-touch page faults on a fresh 80 MB output are ~17% of the n = 1e7 time. Every R generator returning a new vector pays it; only an in-place API skips it. Additive, no existing behaviour changes. | 1 day |
| 3.2 | **Make the threading story the headline.** Benchmark all-threads on the M1 and publish it. | Stateless generation gives bit-identical multi-threaded output; no competitor does. This is where the lead can be 3x, not 7%. | 2 h |
| 3.3 | **`configure` script for OpenMP.** Actually compile-and-link a test program; emit `src/Makevars` from `Makevars.in`. | Forwarding `SHLIB_OPENMP_CFLAGS` blindly broke the UBSan build (`-fopenmp` with no libomp). `__has_include` guards the compile; nothing guards the link. Standard practice for optional OpenMP on CRAN. | half day |
| 3.4 | **Split `src/zurand.c`** into `engine_*.c` / `sampler_*.c` once a second engine lands. | 950 lines is fine for one engine; the purpose-word design makes the seam obvious. | with 2.x |

---

## Phase 4 -- Release

| # | task |
|--:|---|
| 4.1 | Version off `0.0.0.9000` (the one remaining `R CMD check` NOTE under your control) |
| 4.2 | `cran-comments.md` justifying the `-Wunused-const-variable` pragma around NumPy's unused tables |
| 4.3 | README benchmark section generated from `tools/benchmark.R` output on the M1, with the platform header |
| 4.4 | Replace the local path in `CLAUDE.md` (`/Users/pbtz/...random123`) with the upstream URL now the repo is public |

---

## Decisions log

| date | decision | rationale |
|---|---|---|
| 2026-09-19 | Exported API keeps the `rng_*` prefix through the rename to zurand | The prefix is the API's, not the package's; renaming it is a separate breaking change |
| 2026-09-19 | `src/zigbounds.h` numeric tables kept as committed, not regenerated on x86_64 | Regeneration drifts by libm ulps; the guard absorbs it; the committed values are the M1's |
| 2026-09-19 | Faster generators are new engines, never a changed default | Reproducibility is the product; throughput is opt-in |
| 2026-09-19 | `-DR123_USE_MULHILO64_C99=1` stays in `Makevars` | Verified harmless: `philox.h` selects `__uint128_t` first. It is the fallback for 32-bit, which CI's `arch` job exercises weekly |

## Open questions

- Is the 7-round engine acceptable as opt-in? (Phase 2 gate.)
- Does the M1 benchmark change the competitive ranking? (Phase 1 exit.)
- Should the in-place API be exported, or kept `@keywords internal` for
  packages that Import zurand? (Phase 3.1.)
- Should `rng_normal()` gain a `method = c("ziggurat", ...)` argument if a
  SIMD-friendly transform ever beats the ziggurat? Not before Phase 2 says
  the generator is no longer the bottleneck.

## Explicit non-goals

- Beating dqrng on uniform *with the default engine*. Philox4x64-10's
  ceiling makes that a coin flip at best; the honest path is a new engine.
- `-march=native` or any build that produces machine-specific output.
- Changing what any existing key produces. Ever.
