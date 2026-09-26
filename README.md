# zurand

<!-- badges: start -->
[![R-CMD-check](https://github.com/pedrobtz/zurand/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/pedrobtz/zurand/actions/workflows/R-CMD-check.yml)
[![native-checks](https://github.com/pedrobtz/zurand/actions/workflows/native-checks.yml/badge.svg)](https://github.com/pedrobtz/zurand/actions/workflows/native-checks.yml)
![Coverage](https://github.com/pedrobtz/zurand/raw/main/.github/badges/coverage.svg)
<!-- badges: end -->

Fast, stateless random numbers for R. The default engine is xoshiro256++,
reseeded every 512 values from the Philox4x64-10 counter-based generator of
[Random123](https://github.com/DEShawResearch/random123); Philox and Threefry
are available directly.

Every value is a pure function of an immutable **key** plus the sampler's
arguments. Samplers never read or mutate `.Random.seed`, so results do not
depend on call order, on how many values you drew earlier, or on how many
threads the machine has.

## Installation

From CRAN, once released:

```r
install.packages("zurand")
```

The development version from GitHub:

```r
# install.packages("pak")
pak::pak("pedrobtz/zurand")
```

## Quick start

Create a key, then draw from it. The same key and arguments always give the
same values:

```r
library(zurand)

key <- rng_key(42L)
key
#> <rng_key> rng_key[2feb6e95bdd73226b266f10328efe333]

rng_uniform(key, 5L)
#> [1] 0.6436165 0.4170741 0.1279353 0.8522212 0.3491436

rng_uniform(key, 5L)   # same key, same draws
#> [1] 0.6436165 0.4170741 0.1279353 0.8522212 0.3491436
```

Drawing more values extends the sequence rather than changing it — draw `i` is
a function of `(key, i)` alone:

```r
identical(rng_uniform(key, 100L)[1:5], rng_uniform(key, 5L))
#> [1] TRUE
```

## Independent substreams

There are two ways to get values that are independent of each other. Ask for a
key *vector* up front:

```r
keys <- rng_key(42L, n = 3L)
keys
#> <rng_key[3]>
#> [1] rng_key[2feb6e95bdd7...]
#> [2] rng_key[130f9f524752...]
#> [3] rng_key[244823f209bc...]
```

With several keys, a sampler returns an `n` x `length(key)` matrix — one
column per key:

```r
rng_uniform(keys, 4L)
#>           [,1]       [,2]      [,3]
#> [1,] 0.6436165 0.85879755 0.3794245
#> [2,] 0.4170741 0.60044715 0.6979140
#> [3,] 0.1279353 0.09813149 0.6872125
#> [4,] 0.8522212 0.77553624 0.4570548
```

Or derive a key from structured data with `rng_fold()`, which is useful when
the pieces of a computation are named rather than numbered:

```r
rng_fold(key, "layer1")
#> <rng_key> rng_key[8dbacb9f407ebc3ae10ed1561177db53]

rng_uniform(rng_fold(key, "chain 3"), 3L)
#> [1] 0.53081894 0.07105536 0.40177051
```

Folding is typed: `rng_fold(key, 1)` and `rng_fold(key, "1")` give different
keys. Integer and double vectors holding the same whole numbers give the same
key, so `1L` and `1` agree.

## Samplers

| Function | Returns |
|---|---|
| `rng_uniform(key, n, min, max)` | doubles in `(min, max)` |
| `rng_normal(key, n, mean, sd)` | normal deviates, via a ziggurat sampler |
| `rng_integer(key, n, min, max)` | integers in `[min, max]`, inclusive |
| `rng_bits(key, n, bits)` | raw generator words |

```r
rng_normal(key, 5L)
#> [1]  1.6306140 -1.1186343 -0.6639268  0.2393690 -0.1899009

rng_integer(key, 10L, min = 1L, max = 6L)
#>  [1] 3 1 4 2 5 3 3 1 1 5
```

`rng_bits()` exposes the generator stream directly. With `bits = 32` it
returns unsigned 32-bit words stored exactly in a double; with `bits = 64` it
returns fixed-width hex strings, because a double cannot hold 64 bits without
loss:

```r
rng_bits(key, 3L)
#> [1] 2206215613  289400228  609344871

rng_bits(key, 3L, bits = 64L)
#> [1] "43c5f4cd83802dbd" "78da9ead113fe5a4" "299a3bb52451dd67"
```

## Engines

A key records its engine, chosen when it is created:

| `engine` | what it is |
|---|---|
| `"xoshiro256pp"` (default) | xoshiro256++, reseeded from Philox every 512 values; zurand's fastest engine, with an AVX2 path that gives the same bits as the scalar one |
| `"philox4x64"` | Philox4x64-10, counter-based: every value computed directly from its position |
| `"threefry4x64"` | Threefry4x64-13, counter-based |

Every engine gives the same guarantees: the same key gives the same values
on every platform and thread count, and draw `i` does not depend on how
many values you ask for. The default passed a PractRand audit to 1 TB,
including interleaved streams of sibling and folded keys
([dev/statistical-audit.md](https://github.com/pedrobtz/zurand/blob/main/dev/statistical-audit.md)). Changing engine
changes every value, so name it explicitly if you need the stream to
survive a change of default:

```r
key <- rng_key(42L, engine = "xoshiro256pp")
```

## Performance

Millions of values per second for `n = 1e7`, median of five rounds, from the
`benchmark` workflow on GitHub-hosted runners (2026-09-26). Every package is
measured in the same `bench::mark()` call on the same machine.

| | Linux x86_64, 1 thread | Linux x86_64, 4 threads | macOS arm64, 1 thread |
|---|---:|---:|---:|
| **uniform** | | | |
| zurand (`xoshiro256pp`) | 569 | **1370** | 638 |
| zurand (`philox4x64`) | 258 | 587 | 516 |
| randompack (`x256++simd`) | **617** | 656 | **1309** |
| dqrng | 350 | 352 | 256 |
| base R `runif()` | 99 | 103 | 185 |
| **normal** | | | |
| zurand (`xoshiro256pp`) | 329 | **822** | 437 |
| zurand (`philox4x64`) | 197 | 468 | 363 |
| randompack (`x256++simd`) | **333** | 324 | **594** |
| RcppZiggurat (MT) | 123 | 121 | 224 |
| dqrng | 140 | 140 | 229 |
| base R `rnorm()` | 28 | 28 | 43 |

Linux: AMD EPYC 7763, AVX2, OpenMP. macOS: Apple M1, no OpenMP, as in the
CRAN binary. On one thread zurand is ahead of dqrng and RcppZiggurat on
both platforms, and level with randompack's SIMD engine on x86_64;
randompack is faster on Apple Silicon, where zurand has no vector path yet.
With threads, zurand's output is still bit-identical to one thread, and
none of the others parallelizes at the R level. The
[performance article](https://github.com/pedrobtz/zurand/blob/main/vignettes/performance.Rmd)
has the method and how to reproduce these numbers.

The same comparison on your own machine, one call per distribution:

```r
library(zurand)
library(bench)

n <- 1e7
key <- rng_key(42L)                    # zurand's default engine, xoshiro256pp
rp <- randompack::randompack_rng()     # randompack's fastest engine
dqrng::dqset.seed(42L)

bench::mark(
  zurand     = rng_uniform(key, n),
  dqrng      = dqrng::dqrunif(n),
  randompack = rp$unif(len = n),
  base       = runif(n),
  check = FALSE
)
#> # A tibble: 4 × 13
#>   expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>
#> 1 zurand       12.7ms   12.8ms     77.5     76.3MB    96.9     12    15      155ms <NULL>
#> 2 dqrng        30.1ms   31.3ms     32.2     76.3MB    32.2      7     7      218ms <NULL>
#> 3 randompack   21.9ms     22ms     45.5     76.3MB    40.9     10     9      220ms <NULL>
#> 4 base        125.7ms  126.1ms      7.93    76.3MB     7.93     2     2      252ms <NULL>
#> # ℹ 3 more variables: memory <list>, time <list>, gc <list>

bench::mark(
  zurand       = rng_normal(key, n),
  dqrng        = dqrng::dqrnorm(n),
  RcppZiggurat = RcppZiggurat::zrnormMT(n),
  randompack   = rp$normal(len = n),
  base         = rnorm(n),
  check = FALSE
)
#> # A tibble: 5 × 13
#>   expression        min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result
#>   <bch:expr>   <bch:tm> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>
#> 1 zurand           26ms  26.2ms     38.0     76.3MB    38.0      8     8      211ms <NULL>
#> 2 dqrng          66.9ms    67ms     14.5     76.3MB    19.4      3     4      206ms <NULL>
#> 3 RcppZiggurat   50.6ms  51.1ms     19.4     76.5MB    19.4      5     5      258ms <NULL>
#> 4 randompack     37.3ms  37.6ms     26.5     76.3MB    26.5      6     6      226ms <NULL>
#> 5 base          418.7ms 418.7ms      2.39    76.3MB     2.39     1     1      419ms <NULL>
#> # ℹ 3 more variables: memory <list>, time <list>, gc <list>
```

That output is from an Intel Core i5-8500B (x86_64, macOS, one thread, R
4.5.2), where zurand leads both calls. Which package is fastest on one
thread depends on the processor, as the table above shows.

## Using zurand from C

Packages can draw from zurand in C, into buffers they own and from their
own worker threads, and get exactly the values the R functions return.
Declare `LinkingTo: zurand` and `Imports: zurand`, then:

```c
#include <zurand.h>

const zurand_api *zr = zurand_get_api();     /* main thread, once */
zurand_key k;
zr->key_get(keys, 0, &k);                     /* keys: an rng_key from R */

/* then from any thread: */
if (zr->fill_normal(k, n, 0.0, 1.0, buf) != ZURAND_OK) { /* bad arguments */ }
```

The header, [inst/include/zurand.h](https://github.com/pedrobtz/zurand/blob/main/inst/include/zurand.h), documents the
rest: uniform, integer and 64-bit fills, and `fold_int()`, which derives a
key exactly as `rng_fold(key, i)` does.

## Interop with R's RNG

`rng_key_from_r()` is the one function here that deliberately consumes R's
global RNG state, for when you want a fresh key seeded from `set.seed()`:

```r
set.seed(123)
k <- rng_key_from_r()
```

Everything else leaves `.Random.seed` untouched.

## Parallelism

When the package is compiled with OpenMP, the uniform, normal and integer
samplers parallelize large fills. Output is bit-identical regardless of the
thread count — this is a performance control only, never a reproducibility
one.

```r
rng_threads()        # effective maximum
old <- rng_threads(1L)
rng_threads(old)     # restore
```

On macOS, Apple clang ships without OpenMP, and R's macOS configuration
leaves `SHLIB_OPENMP_CFLAGS` empty, so the CRAN binary for macOS is
single-threaded. On a Mac the single-thread speed is the speed. A source
build gets threads if `~/.R/Makevars` defines `SHLIB_OPENMP_CFLAGS` against
a `libomp` installation.

## License

MIT, see [LICENSE](https://github.com/pedrobtz/zurand/blob/main/LICENSE).

This package bundles three third-party components, credited in
[inst/COPYRIGHTS](https://github.com/pedrobtz/zurand/blob/main/inst/COPYRIGHTS): the Random123 library (D. E. Shaw
Research, BSD 3-clause), NumPy's ziggurat constant tables (NumPy Developers,
BSD 3-clause), and fdlibm's `exp` and `log1p` (Sun Microsystems, permissive
notice), which keep normal draws identical across platforms.
