# Gaussian throughput baseline

Reference point for the "fastest Gaussian generator for R" work. Regenerate
with `R CMD INSTALL . && Rscript tools/benchmark.R`, which now prints the
machine it ran on as a header -- these numbers mean nothing without it.

> **This run is on x86_64, which is NOT the machine zurand was tuned on.**
> The package was developed and tuned on an Apple Silicon (M1) Mac, and
> several decisions in CLAUDE.md are explicitly arm64: the
> `R123_FORCE_INLINE` requirement is justified by a 32-byte struct return
> going through memory *on arm64*, and the verification step is `otool -tv`.
> Read the numbers below as a secondary-platform datapoint. An M1 run is
> the one that should drive tuning decisions.
>
> Three things differ enough to change conclusions:
> - **Cache line**: 128 bytes on M1, 64 on Intel -- directly changes how
>   much the ki/wi table split costs.
> - **Philox's 64x64->128 multiply**: one instruction on x86 (`mulq`),
>   two on arm64 (`mul` + `umulh`).
> - **OpenMP**: present on the M1 dev machine, absent here.

## Environment (this run)

| | |
|---|---|
| CPU | Intel Core i5-8500B @ 3.0 GHz (Coffee Lake, 6 physical / 6 logical) |
| SIMD | AVX2, FMA, BMI1/2, SSE4.2 — **no AVX-512** |
| R | 4.5.2, x86_64-apple-darwin20 |
| Compiler | Apple clang 17.0.0, `-O2 -falign-functions=64` |
| OpenMP | **absent on this machine** — no libomp, no `~/.R/Makevars`, so `rng_threads()` is 1 and every number below is single-threaded. CLAUDE.md's claim that the machine is configured for OpenMP is correct *for the M1 dev machine*, not for this one. |

zurand must be measured after `R CMD INSTALL` (`-O2`). A `load_all()` debug
build understates it.

## Normal, n = 1e7, single thread

| rank | method | M values/s | ns/value | vs zurand |
|---:|---|---:|---:|---:|
| 1 | **zurand** | **213.5** | 4.68 | — |
| 2 | RcppZiggurat MT | 199.1 | 5.02 | 1.07x |
| 3 | randompack | 184.3 | 5.43 | 1.16x |
| 4 | RcppZiggurat LZLLV | 158.8 | 6.30 | 1.34x |
| 5 | dqrng | 148.8 | 6.72 | 1.43x |
| 6 | RcppZiggurat QL | 72.8 | 13.74 | 2.93x |
| 7 | RcppZiggurat GSL | 67.2 | 14.89 | 3.18x |
| 8 | base R (inversion) | 22.8 | 43.90 | 9.37x |
| 9 | rTRNG | 16.2 | 61.64 | 13.16x |

zurand already leads, but only by 7% over the next generator.

## Uniform, n = 1e7, single thread

| rank | method | M values/s | ns/value |
|---:|---|---:|---:|
| 1 | dqrng | 333.1 | 3.00 |
| 2 | randompack | 255.2 | 3.92 |
| 3 | **zurand** | **189.1** | 5.29 |
| 4 | sitmo | 183.0 | 5.46 |
| 5 | rTRNG | 171.9 | 5.82 |
| 6 | base R | 78.8 | 12.69 |

## Cost decomposition

Measured, not assumed. Raw Philox from a standalone C microbenchmark; the
rest by difference against the R-level timings.

| component | ns/value | share of normal |
|---|---:|---:|
| Philox4x64-10 (1 word per draw) | 2.74 | 58% |
| ziggurat transform + buffer round trip | 1.94 | 42% |
| **total** | **4.68** | 100% |

Philox itself is the floor: **365 Mword/s / 2.74 ns per 64-bit word**. No
Gaussian path built on philox4x64-10 can exceed that, so 213 M/s today is
58% of the attainable ceiling.

A standalone microbenchmark of the ziggurat fast path alone measures
1.91 ns/draw, independently confirming the 1.94 ns above.

## Findings

**Ruled out on x86_64 — manual ILP.** Interleaving 2, 4 or 8 independent
Philox blocks per iteration to hide the 10-round dependency chain is
*slower* than the plain loop (0.79x, 0.83x, 0.87x). clang already
software-pipelines the independent iterations; explicit unrolling only adds
register pressure. **Worth re-testing on M1** before ruling it out there:
arm64 has 31 general-purpose registers to x86-64's 16, so the register
pressure that sinks it here may not apply.

**Ruled out — a suspected mulhilo bug.** `src/Makevars` defines
`-DR123_USE_MULHILO64_C99=1`, which looks like it forces the portable
4-multiply fallback. It does not: `philox.h` tests `R123_USE_GNU_UINT128`
first, so x86_64 still takes the single 128-bit multiply. Verified by
`#warning` probe and by timing both ways — 10.96 vs 10.97 ns/block.

**Confirmed win on x86_64 — interleaved ziggurat tables.** `ki_double[idx]`
and `wi_double[idx]` are separate arrays, so each draw touches two cache
lines. Packing them into one `{uint64_t ki; double wi;}` array is **1.10x**
on the transform (1.910 -> 1.737 ns/draw), worth roughly +3.7% end to end.
The vendored NumPy header must stay verbatim, so the combined table has to
be derived rather than edited in place. The gain should differ on M1, whose
128-byte lines hold twice as many table entries per line.

**Confirmed gap — uniform does not use the two-pass fill.**
`fill_normal_column` buffers a chunk of Philox output and transforms it in a
second pass, which is why its transform pipelines well.
`fill_uniform_column` still fuses generate-and-transform. That is why
uniform costs *more* per value (5.29 ns) than normal (4.68 ns) despite doing
strictly less arithmetic. Porting the two-pass fill to uniform should close
most of that gap.

## Ranked next steps

| # | change | est. gain on normal | breaking? |
|---:|---|---|---|
| 1 | Philox rounds 10 -> 7 | **~+21%** (258 M/s) | **yes** — changes every value |
| 2 | Interleaved ki/wi table | +3.7% | no |
| 3 | `ZURAND_CHUNK_BLOCKS` sweep (currently 128) | unknown, cheap to test | no |
| 4 | Two-pass fill for uniform | none on normal; ~+20% on uniform | no |
| 5 | SIMD Philox | poor fit: AVX2 has no 64x64->128 multiply, and this CPU has no AVX-512 | no |

Item 1 is the only one that moves the needle decisively, and it is a design
decision rather than an optimization: Random123's authors report philox4x64-7
passing BigCrush, but dropping to 7 rounds changes the output of every
existing key.

## Note: `src/zigbounds.h` is platform-stamped

Regenerating the header on this x86_64 machine rewrites the gap table by +/-1
in the low bits relative to the committed values, which were generated on the
M1. `tools/generate-zig-bounds.R` goes through R's `exp()`/`log()`, and arm64
and x86_64 libm differ by an ulp on some inputs.

This is not rot, and it is not a correctness problem: `ZURAND_ZIG_GUARD`
widens the ambiguous band precisely so a bracket that is off by a few
fixed-point units defers to the `exp()` fallback rather than deciding wrongly.

It does mean the committed header records *the machine that generated it*, so
"rerun the script instead of editing it" produces a spurious diff on any other
machine. Regenerate on the M1, or the diff is noise.
