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

- **`-O2` build only.** `devtools::load_all()` defaults to `debug = TRUE`,
  which compiles `-O0` and measures 11x slower. Use `R CMD INSTALL` or
  `load_all(debug = FALSE)` -- verified equivalent, agreeing to 0.2% against
  a fixed in-session reference. `R CMD INSTALL .` builds in-place, so clear
  `src/*.o` first if a stale object is possible.
- **Primary machine: the Apple Silicon M1** the package is tuned on. The
  x86_64 Intel machine is a secondary datapoint. The script prints a
  platform header; a result without one is discarded.
- **Ratios, not absolutes.** This machine's absolute throughput moved 13%
  across six identical rounds under background load, and the ratio against a
  competitor moved 12.6% -- both larger than the ~5% effect being chased.
  Measure zurand and the reference in the *same* `bench::mark` call and
  report the ratio, repeated over several rounds, never a single median from
  a run made an hour earlier.
- **Four metrics**, in priority order:
  1. Gaussian, single thread, n = 1e7 -- algorithm vs algorithm.
  2. Uniform, single thread, n = 1e7.
  3. Gaussian and uniform, all threads, n = 1e7 -- what a laptop user gets.
  4. Per-call overhead at n in {1, 100, 1000} -- already at parity with
     base R and ahead of dqrng; must not regress.
- **Competitors:** base R, dqrng, RcppZiggurat (MT, LZLLV, GSL, QL),
  randompack, sitmo, rTRNG. Add any package that enters the field.

**Acceptance for a change:** the ratio interval over repeated rounds must
separate from the baseline's -- a bare +2% median is not evidence at the
noise levels measured here -- with no regression on any other metric, and `rng_*()` output bit-identical to before
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
| 0.1 | ~~KAT against Random123's reference vectors~~ **DONE, both engines** Upstream ships `tests/kat_vectors` with philox4x64 entries. Compare `rng_bits(key, n, bits = 64L)` for the reference keys/counters. | Nothing today proves the Philox is Philox. | new test file; ~15 lines | 1 h |
| 0.2 | **Golden-value test for the full pipeline.** Hardcode ~256 values each of `rng_uniform`, `rng_normal`, `rng_integer` for a fixed key, including values known to hit the ziggurat wedge and tail paths. | CI runs on arm64, x86_64, i386 and musl; nothing asserts they agree. This is the reproducibility guarantee stateless generation exists to give, and the gate every later change passes through. | new test file | 2 h |
| 0.3 | ~~One-time statistical audit~~ **TOOLING DONE** -- `tools/statistical-audit.R` plus an on-demand workflow; still needs a long run and its results recorded. Originally: PractRand (or TestU01 SmallCrush) on `rng_bits()`; the same on `pnorm(rng_normal())` to exercise the custom wedge shortcut and tail. Not CI -- a document. | The uniform conversion and wedge brackets are custom code. One KS test at n=50k is a smoke test, not evidence. | `dev/statistical-audit.md` | 1 day |
| 0.4 | **Make `src/zigbounds.h` platform-independent.** Generated on the M1, it drifts by +/-1 low bit when regenerated on x86_64 (libm ulp differences). `ZURAND_ZIG_GUARD` absorbs it, so it is not a correctness bug, but "rerun the script" yields spurious diffs off the M1. Either compute the brackets with `Rmpfr`/exact rationals, or document "regenerate on arm64 only" in the header itself. | Reproducibility of the build, not of the output. | regenerate on both machines, `diff` empty | 2 h |

**Exit criterion:** 0.1 and 0.2 green on all CI legs. 0.3 and 0.4 can trail.

---

## Phase 1 -- Free wins (non-breaking, default stream unchanged)

Each is measured; each is small. Together perhaps +5% Gaussian, +20% uniform.

| # | task | expected | evidence | effort |
|--:|---|---|---|---|
| 1.1 | ~~Two-pass fill for uniform~~ **DONE** | **+56% uniform** (189 -> 295 M/s); ratio to dqrng 0.67 -> 0.92 | delivered; far exceeded the +20% estimate | done |
| 1.2 | **Packed ziggurat table.** `ki_double[idx]` and `wi_double[idx]` are separate arrays: two cache lines per draw. Derive a `{uint64_t ki; double wi;}` array from the vendored header at build time (a generated `src/zigtable.h`, like `zigbounds.h`). | **~+4% Gaussian** on x86; re-measure on M1's 128-byte lines | standalone microbenchmark: 1.910 -> 1.737 ns/draw, 1.10x | 2 h |
| 1.3 | **Fuse `mean`/`sd` into the transform.** `C_rng_normal` fills standard normals then sweeps the whole output again with `mean + sd * x`. The transform already ends in a multiply: `x = s * (wi[idx] * sd) + mean` is one FMA with `wi*sd` hoisted per call. | removes a full memory-bound pass whenever `mean`/`sd` are non-default | code reading; measure with `mean = 1, sd = 2` | 2 h |
| 1.4 | ~~`ZURAND_CHUNK_BLOCKS` sweep~~ **inconclusive on x86** | none measurable | swept 16/32/64/128/256: uniform ratios 0.93-0.99, normal 0.99-1.31, i.e. pure noise on a loaded machine. Kept 128. Redo on a quiet M1. | done (retry) |
| 1.5 | **Re-test manual ILP on the M1.** Interleaving 2/4/8 independent Philox blocks was 0.79-0.87x on x86 -- register pressure with 16 GPRs. arm64 has 31. | unknown; do not carry the x86 conclusion across | `/tmp/philox_bench.c` from the baseline work; rerun | 30 min |

**Exit criterion:** baseline re-recorded on the M1 with 1.1-1.3 landed.

---

## Phase 2 -- The generator (decisive gains, via new engines)

The only route to the targets in §0. Each candidate is an *additional*
engine selected by `rng_key(engine = ...)`; `"philox4x64"` stays the
default and stays bit-identical.

| # | candidate | Philox cost | why it might win | why it might not | effort |
|--:|---|---|---|---|---|
| 2.1 | ~~`philox4x64-7`~~ **REJECTED** | **0.95x -- slower than 10 rounds**, reproducibly, over four runs | would have broken every stream for a regression | closed |
| 2.2 | ~~Threefry4x64-13~~ **SHIPPED** as `engine = "threefry4x64"` | 1.78x raw generator, but only **1.09x uniform / 1.15x normal** end to end | past L1 the fill is not generator-bound; see the note below | done |
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

## Amendment, after Phase 2 (2026-09-19)

Faster generators are not the lever this roadmap assumed.

threefry4x64-13 is 1.78x philox4x64-10 per word in a register-only
microbenchmark, and delivers 1.09x on uniform and 1.15x on normal at
n = 1e7. Sweeping the working set shows why: the speedup is 1.53x at
n = 1000, where the output fits L1, and is gone by n = 10000. A real fill
stores to memory, and past L1 the generator is no longer the constraint.

This is consistent with the one big win so far: the two-pass uniform fill,
a pure memory-layout change, bought +56%, while a 1.78x faster generator
bought +9%.

So Phase 3's memory-traffic items outrank any further engine work, and
should be measured before they are built:

- Is the chunk buffer's store-and-reload actually costing what was assumed?
  The in-place variant measured ~12% slower for normal, but that predates
  the two-pass uniform change and deserves a re-test.
- Is the 17% first-touch page-fault figure real? It was one median on a
  loaded machine, and single medians have been wrong three times today.

ARS remains the fastest thing measured (3.40x) and cannot be a default:
without AES-NI `ars4x32_ctr_t` is not declared, and neither ars.h nor aes.h
has a NEON path, so it does not compile on the M1 at all.

## Memory-traffic measurements (2026-09-19)

The amendment above said to measure Phase 3's assumptions before building
them. Done, in C so the allocator and the fill can be separated, then
confirmed at the R level.

### The chunk buffer is not costing us -- item closed

Three layouts, identical output, one allocation reused across reps so
page-fault cost is excluded:

| layout | n=1e5 | n=1e6 | n=1e7 |
|---|---:|---:|---:|
| **(a) two-pass, stack buffer (current)** | **351 M/s** | **334** | **305** |
| (b) fused, no buffer | 200 (0.57x) | 198 (0.59x) | 189 (0.62x) |
| (c) raw words into out[], transform in place | 340 (0.97x) | 293 (0.88x) | 278 (0.91x) |

The current design wins at every size. The stack buffer stays in L1 while
the output is written once and streamed; (c) turns that into
write-read-write over the full output, and (b) makes every store wait on a
10-round Philox chain. Nothing to win here -- the store-and-reload
hypothesis is wrong.

### Allocation is the real cost, and it is R-wide

Same fill, reused buffer vs a fresh malloc per call:

| n | reused | fresh alloc | penalty |
|---|---:|---:|---:|
| 1e6 (7.6 MB) | 300 M/s | 296 | +1.3% |
| 1e7 (76 MB) | 290 M/s | 253 | **+14.6%** |
| 5e7 (381 MB) | 289 M/s | 134 | **+116%** |

R pays it, and `bench` reports a GC in every iteration at the larger sizes:

| n | zurand | dqrng | ratio |
|---|---:|---:|---:|
| 1e6 | 251 M/s | 275 | 0.91 |
| 1e7 | 242 M/s | 267 | 0.91 |
| 5e7 | 142 M/s | 139 | **1.03** |

Two things follow. The tax is **not zurand-specific** -- dqrng degrades
identically, and at 5e7 the allocation cost so dominates that zurand edges
ahead. And it is **large**: at 5e7 R gets 142 M/s where the same fill into
a reused buffer gets 289.

So the earlier "17% at n=1e7" figure was about right, and understated what
happens above it.

### What this means for item 3.1

An in-place fill is the only lever that touches this, and no competitor
offers one, so it is a real differentiator rather than a micro-optimisation.
The risk noted earlier stands and has a standard answer: do not mutate a
caller's vector unconditionally. Take the vector, check `MAYBE_REFERENCED`,
fill in place only when unshared and duplicate otherwise, and return it --
the same contract R's own subassignment uses. A caller writing
`x <- rng_fill_uniform(x, key)` in a loop then gets the benefit safely, and
a caller who aliased the vector gets correct results instead of corruption.

Still unbuilt, and still a permanent public API commitment, so it wants an
explicit decision rather than being taken as implied by the numbers.

## Amendment 2: the counter-based constraint was self-imposed (2026-09-19)

The roadmap assumed every engine had to be a Random123 counter-based
generator, and concluded from measurement that uniform could not beat
dqrng. The measurement was right and the conclusion was too broad: it held
only under that assumption, and the assumption was never required.

`value = f(key, index)` for an arbitrary index is not something the API
exposes -- there is no `rng_at(key, i)`. What the API promises is that a
key reproduces its stream, that draw i does not depend on n, and that a
threaded fill matches a serial one. All three survive seeding a cheap
recurrence per chunk.

`xoshiro256pp` does that: one Philox call per 512 values derives an
xoshiro256++ state, which then costs about five operations per word
instead of the 2.5-3 ns a counter-based engine needs to recompute from its
counter. Measured through the R API at n=1e7, single-threaded:

  uniform    249 -> 581 M/s    1.79x dqrng   (philox was 0.88x)
  normal     183 -> 316 M/s    1.61x RcppZiggurat MT (philox was 0.97x)

philox4x64 remains the default and keeps the indexed-access property open
for anyone who needs it.

Remaining headroom, roughly: SIMD xoshiro (4-8 lanes) could approach the
~1200 M/s memory write floor on uniform, but CRAN cannot ship -march and
it needs runtime dispatch. The Gaussian path is separately bounded by the
scalar ziggurat at about 1.9 ns/draw, so beating ~400 M/s there needs a
SIMD-friendly transform, which is a research problem rather than an
optimisation.

## Phase 5 -- SIMD (not started)

Recorded while the context is fresh. `xoshiro256pp` closed the gap to the
field; SIMD is what is left, and it is a portability project rather than
an optimisation.

### The ceiling it is chasing

Measured on x86_64, single thread, reused buffer:

| | M values/s |
|---|---:|
| philox4x64-10 fill | ~320 |
| xoshiro256pp fill (current best) | ~580 uniform, ~320 normal |
| `numeric(n)` -- allocate and zero-fill | ~1180 |

That last row is roughly the single-core write bandwidth, about 9.4 GB/s,
and nothing returning a fresh vector can pass it. So uniform has maybe
2x of headroom left and Gaussian rather more, but Gaussian is bounded
first by the scalar ziggurat at ~1.9 ns/draw.

### 5.1 SIMD xoshiro for uniform

xoshiro256++ is add/shift/rotate/xor and vectorises cleanly: four
independent lanes in AVX2 (4 x 64-bit), eight in AVX-512, two per NEON
register on arm64 -- though the M1 has four NEON units, so it issues
wider than that suggests. `u01_open` vectorises too and already does at
SSE2 width.

**The constraint that matters, and it is easy to miss.** Output must not
depend on whether SIMD was available at build or run time, or the golden
tests break and the package's central promise with them. That means the
lane assignment is part of the *engine definition*, not an implementation
detail: pick it once -- lane L takes words L, L+4, L+8, ... of a chunk,
or contiguous quarters -- and make the scalar path reproduce exactly that
interleaving. Writing the scalar version first and vectorising it later
will produce two different streams.

This also argues for deciding it now, before `xoshiro256pp` has users: if
the lane layout is baked into the engine from the start, a later SIMD
implementation is a pure speedup with no stream change. As it stands the
current engine is defined as a single sequential recurrence per chunk, so
adding SIMD later *would* change its output and need a fourth engine name.
Worth fixing while nothing depends on it.

### 5.2 SIMD Gaussian

Harder, and possibly not worth it. The ziggurat's per-draw table lookup is
the obstacle: AVX2 has `vgatherqpd` but it is slow on older cores (~12
cycles on Coffee Lake), which can cost more than it saves. Options, none
free:

- Vectorise the fast path with gather, keep a scalar fallback for the ~1%
  that reject. Needs measuring before it is believed.
- A transform without table lookups -- Box-Muller vectorises well but
  needs `sin`/`cos`/`log`, so it wants polynomial approximations and is
  slower scalar. It would also be a different stream, i.e. another engine.
- Accept ~400 M/s as the Gaussian ceiling and put the effort elsewhere.

### 5.3 Runtime dispatch

CRAN does not accept `-march=native`, so the SIMD path cannot simply be
compiled in. The standard shape is a separate translation unit built with
`-mavx2`, selected at run time behind a CPUID check, with a baseline
build always present. `configure` (item 3.3) would grow the job of
deciding whether the compiler accepts the flags at all.

Combined with 5.1's constraint, the dispatch must be a *performance*
choice only: both paths must emit identical bits, and a golden test
should run under each.

### Order

1. Settle the lane layout for `xoshiro256pp` before it has users (5.1a).
2. Prototype SIMD xoshiro as a standalone microbenchmark on both x86_64
   and the M1 -- the pattern that has worked all along, and the one that
   stopped philox4x64-7 and ARS from being built.
3. Only then decide whether the dispatch machinery earns its complexity.
