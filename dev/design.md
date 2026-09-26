# zurand design, v1 stream

Status: proposed 2026-09-25. Supersedes the implicit design spread across
`dev/roadmap.md`, its amendments and `dev/philox.md`; builds on
`dev/review-2026-09-23.md` and `dev/ecosystem-survey.md`.

This document fixes what zurand *is* at v0.1.0: the output contract, the
stream specification that will be frozen, the public surface in R and C,
and the rule that decides whether a future change is allowed. `dev/roadmap.md`
says in what order it gets built.

---

## 1. Position

**The fastest uniform and Gaussian vectors in R, where every value is a pure
function of `(key, position)`.** That one property is what the speed is
spent on and what no competitor offers together with it: output is
bit-identical across thread counts, SIMD paths, platforms, call order and
chunking, and draw `i` does not depend on `n`.

What zurand competes on, in order:

1. **Reproducible parallelism.** Threaded fills equal serial fills, and the
   C API is reentrant, so other packages can draw from their own worker
   threads. dqrng and base R cannot offer the second; randompack's streams
   are stateful.
2. **Throughput.** Today (x86_64, one thread, n = 1e7) `xoshiro256pp` with
   AVX2 is 776 M/s uniform (2.4x dqrng) and 376 M/s Gaussian (1.9x
   RcppZiggurat). The default engine gets 314 and 215.
3. **A small, stable surface.** Uniform, normal, integer, bits, and an
   exponential. Fourteen distributions is randompack's niche, not this one.

## 2. The rule for changes

Every change is one of three kinds, and the kind decides when it may land.

| kind | test | when |
|---|---|---|
| **output-neutral** | every existing call returns identical bits (SIMD paths, threading, fusing a pass, faster code) | any time |
| **additive** | new argument with a default that reproduces today, new sampler with a fresh purpose value, new engine name | any time |
| **output-changing** | any existing call returns different bits | only before the freeze; afterwards only as a new engine name or a new stream version |

Sorting the review's pre-freeze list by this rule shrinks it. Vector
`mean`/`sd`/`min`/`max`, an `offset` argument and fusing the
`mean + sd * z` scaling into the per-chunk transform are all **additive or
neutral**. The fused form keeps the same formula and only moves it from a
second whole-array pass into the chunk stage, while the chunk is still in
L1. None of them needs to block v0.1.0. What does is in §3.4.

## 3. The stream specification (frozen at v0.1.0 as stream 1)

### 3.1 Keys

A key is 128 bits `{k0, k1}` plus two attributes:

- `engine`: one of `"xoshiro256pp"` (**default from v0.1.0**),
  `"philox4x64"`, `"threefry4x64"`.
- `stream`: integer, `1L`. Samplers reject a stream version they do not
  implement. A key without the attribute is stream 1, so keys saved before
  it existed keep working. This field makes a post-release fix possible
  without breaking saved keys: a bug found after the freeze becomes
  stream 2, and stream-1 keys keep producing stream-1 values.

`rng_key(seed, n)` takes key `i` from consecutive `splitmix64` outputs
seeded by `seed`. The key words do not depend on the engine: the same seed
gives the same `{k0, k1}` under every engine, and engines are separated in
the counter instead (§3.2). `rng_fold(key, data)` derives
`{k0, k1} = philox4x64({h, h ^ 0x9e3779b97f4a7c15, PURPOSE_FOLD, 0}, key)[0:1]`,
with `h` the 64-bit hash of the typed data, for
philox- and xoshiro-engine keys, and the threefry analogue for threefry
keys.

### 3.2 Counters

Every Philox or Threefry evaluation uses a 256-bit counter with one meaning
per word:

| word | meaning |
|---|---|
| 0 | block index (`position >> 2` for the counter engines; sub-chunk index for xoshiro seeding) |
| 1 | domain: 0 on the fast path; attempt number `g >= 1` for ziggurat and integer retry words; `h ^ 0x9e37...` for fold (whose word 0 is `h`) |
| 2 | purpose: 0 bits, 1 uniform, 2 normal, 3 integer, 4 fold; 5 exponential and 6 normal-by-inversion reserved; 7+ free |
| 3 | **engine tag**: 0 for the philox engine's own output, **1 for every Philox evaluation made on behalf of `xoshiro256pp`** (sub-chunk seeds and ziggurat retry words) |

Word 3 is the one change from today. At present `xoshiro256pp` seeds
sub-chunk `s` from `{s, 0, purpose, 0}`, the same counter as the philox
engine's block `s`. That makes the two engines' streams for the same seed
functions of each other: the first xoshiro word of `rng_key(42)` is
`rotl(s0 + s3, 23) + s0` of philox's first block, verified 2026-09-25.

### 3.3 Samplers

Each sampler is defined on the standard draw at a position, and each
per-element formula is evaluated in IEEE binary64 **without contraction**
(no FMA) and without excess precision.

| sampler | per position | parameters |
|---|---|---|
| bits | the 64-bit word; `bits = 32` takes the low half | -- |
| uniform | `u = ((w >> 12) \| 0x3ff0...) as double - (1 - 2^-53)`, exactly `(m + 0.5) 2^-52`, so `u` is in (0, 1) | `min + (max - min) * u` when non-default |
| normal | ziggurat on NumPy's 256-layer tables: word -> layer, sign, 52-bit `rabs`; accept `rabs < ki` gives `z = ±rabs * wi`; wedge and tail draw words at `(position, g >= 1)` | `mean + sd * z` when non-default |
| integer | Lemire bounded 32-bit on the low word, retries at `(position, g >= 1)` | inclusive `[min, max]` |
| exponential (additive) | ziggurat on NumPy's exponential tables, purpose 5 | `rate` |

`uniform`'s interval is stated the way base R states `runif`'s: strictly
inside `(min, max)` for the default bounds, and able to return a bound only
when `max - min` is small relative to `|min|`. That is a documentation
change, not a clamp. A clamp would change output and buys nothing base R
users expect.

### 3.4 What must be settled before the freeze

These change existing output, so they happen before v0.1.0 or never:

1. **Engine tag in counter word 3** (§3.2). Changes `xoshiro256pp` only.
2. **Default engine becomes `xoshiro256pp`.** Conditions, in order:
   - Separation lands first (item 1).
   - A long audit passes (§6).
   - No counter-based engine is within 10% of it in the fill.
     **Settled 2026-09-25:** squares64 (Widynski), in zurand's own
     two-pass fill on the i5-8500B at n = 1e7, runs 428 M/s against scalar
     xoshiro's 692, i.e. **0.62x**, stable over three rounds. It cannot
     use AVX2 (no 64-bit multiply-low below AVX-512), while xoshiro gains
     another 1.29-1.44x there. Philox measured 0.49x in the same harness.

   Philox's one structural advantage, O(1) access to any position, is not
   lost: xoshiro reaches position `p` by seeding sub-chunk `p / 512` and
   stepping at most 511 times, about a microsecond. `offset` (§4) works
   for every engine.
3. **libm independence of `rng_normal()`.** *Done (#20).* Whole-stream
   digests split CI in two, Apple and Windows against glibc and musl,
   because a one-ulp difference in the tail's `log1p()` or the wedge's
   `exp()` flips accept/reject decisions. Neither side is correctly
   rounded: Apple's `log1p` is off by one ulp on 6% of the calls tested.
   The ziggurat now uses its own port of fdlibm's `exp` and `log1p`
   (`src/zurand_fdlibm.h`), a fixed sequence of IEEE operations and the
   same choice Java's StrictMath made. The resulting stream equals what
   glibc and musl already produced, so only macOS and Windows output
   moved.
4. **Contraction is forbidden on every compiler.** *Done (#19).* GCC
   ignores `#pragma STDC FP_CONTRACT OFF` and fuses by default on arm64;
   the new `ubuntu-24.04-arm` leg failed the golden tests by one ulp
   before the fix. A compiler flag was rejected: R appends the user's
   `CFLAGS` after the package's, so `-ffp-contract=fast` in a
   `~/.R/Makevars` would undo it, and `R CMD check` flags any `-f` flag
   in `src/Makevars`. Instead every product that feeds an addition goes
   through `zurand_rounded()`, an empty `asm` no compiler sees through,
   two lanes at a time in the scaling passes so they stay SIMD. No
   measurable cost. `-ffast-math` stays outside the contract.
5. **The `stream` attribute on keys** (§3.1). Strictly additive, but cheap
   now and awkward to explain later.

Everything else in the review's pre-freeze list is additive or neutral by
§2 and moves after the release.

## 4. The R surface at v0.1.0

```r
rng_key(seed, n = 1, engine = c("xoshiro256pp", "philox4x64", "threefry4x64"))
rng_key_from_r(n = 1, engine = ...)
rng_fold(key, data)

rng_uniform(key, n = 1L, min = 0, max = 1)
rng_normal (key, n = 1L, mean = 0, sd = 1)
rng_integer(key, n = 1L, min, max)
rng_bits   (key, n = 1L, bits = 32L)

rng_threads(n = NULL); rng_simd(enable = NULL)
```

After v0.1.0, additively:

- `offset = 0` on every sampler: positions `offset:(offset + n - 1)`. It
  draws the next batch without folding, and it is how a caller checkpoints
  a long simulation.
- Vector `mean`, `sd`, `min` and `max`, recycled to `n` (the
  hierarchical-model idiom). Each element uses the same formula as the
  scalar case, so a constant vector gives the scalar result.
- `rng_exponential(key, n, rate = 1)`.
- `rng_normal(method = "inversion")`: Wichura AS241 on the uniform,
  purpose 6. It is monotone in `u`, which common random numbers,
  antithetic variates and QMC all need.
- `rng_permutation(key, n)` and `rng_sample()`: sort by 64-bit keys,
  exactly uniform.

Deliberately absent: an R-level in-place fill (PR #4 was closed for the
copy-on-write hazard; the C API replaces it) and a stateful default. A
`RNGkind("user-supplied")` bridge may come later as a clearly labelled
opt-in.

## 5. The C API (for packages that `LinkingTo: zurand`)

This is the adoption multiplier and the safe home of the in-place fill.

```c
/* inst/include/zurand.h -- header-only shims over R_GetCCallable() */
typedef struct { uint64_t k0, k1; int engine; int stream; } zurand_key;

int  zurand_api_version(void);                 /* bumped on ABI change */
int  zurand_key_get(SEXP keys, R_xlen_t i, zurand_key *out); /* main thread */
zurand_key zurand_fold_u64(zurand_key k, uint64_t data);

void zurand_fill_uniform(zurand_key k, uint64_t offset, size_t n,
                         double min, double max, double *out);
void zurand_fill_normal (zurand_key k, uint64_t offset, size_t n,
                         double mean, double sd, double *out);
void zurand_fill_integer(zurand_key k, uint64_t offset, size_t n,
                         int min, int max, int *out);
void zurand_fill_bits64 (zurand_key k, uint64_t offset, size_t n,
                         uint64_t *out);
```

Contract: everything except `zurand_key_get` is **reentrant**. The
functions call no R API, allocate nothing and touch no mutable global
state. The SIMD dispatch flag is read-only and output-neutral. So a caller
may invoke them from its own OpenMP or pthreads workers and get the same
values as the R functions. The caller owns `out`, so there is no
copy-on-write question. The fill functions are the same code the R
samplers run, since the R entry points become thin wrappers over them.

## 6. Evidence required before the freeze

| gate | exists | still needed |
|---|---|---|
| KAT against Random123 vectors, both counter engines | yes | -- |
| golden values, hex floats, all samplers | yes | regenerate once at the freeze |
| SIMD = scalar, threads = serial | yes | -- |
| whole-stream digests, 1e6 draws per case, on every CI leg | **yes** (#20) | -- |
| golden tests on i386 and musl | **pass**, first run 2026-09-25 (PR #17; before it, the `arch` legs died at `library(testthat)` and the tests had never run there) | keep the weekly leg green |
| golden tests on GCC + aarch64 | **yes** (#19) | -- |
| long audit, `xoshiro256pp` | 512 MB smoke run only | at least 1 TB of `rng_bits` and 256 GB of `pnorm(rng_normal())` through PractRand, plus the two cross-key tests below |

The audit must also test **across keys**, because that is how the package
is used: many short streams, not one long one. Two extra PractRand inputs:
round-robin interleaving of the streams of `rng_key(s, 64)`, and of
`rng_fold(key, 1:64)`. A weakness in the splitmix64 key derivation, the
fold, or the engine tag shows up there and nowhere else.

## 7. Reproducibility contract (what the README will promise)

Bit-identical output for the same key, sampler, arguments and positions:

- across thread counts, SIMD on or off, `n`, and call order: **guaranteed
  and tested**.
- across x86_64 and aarch64 on Linux (glibc and musl), macOS and Windows:
  **guaranteed**, and tested on every CI leg by golden values and
  whole-stream digests.
- on 32-bit x86 with x87 arithmetic: **not guaranteed**. x87 rounds
  each addition to 64 bits and then to 53 on store, and the digests
  confirmed the difference (scaled uniform, which uses no libm). No CRAN
  platform has been 32-bit since R 4.2.0, and Debian's i386 R is the last
  common source. The `arch` leg still runs there for crashes and
  undefined behaviour; the digest tests skip.
- across zurand versions: guaranteed per `(engine, stream)` from v0.1.0.
  Output before v0.1.0 is not preserved.

## 8. Non-goals

- `-march=native` or any machine-specific output.
- float32 at the R level (R has no type for it). Possible in the C API if
  asked.
- Distributions beyond uniform, normal, exponential and integer. Anything
  else composes from these.
- A stateful default API.
