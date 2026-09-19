# zurand

<!-- badges: start -->
[![R-CMD-check](https://github.com/pedrobtz/zurand/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/pedrobtz/zurand/actions/workflows/R-CMD-check.yml)
[![native-checks](https://github.com/pedrobtz/zurand/actions/workflows/native-checks.yml/badge.svg)](https://github.com/pedrobtz/zurand/actions/workflows/native-checks.yml)
![Coverage](https://github.com/pedrobtz/zurand/raw/main/.github/badges/coverage.svg)
<!-- badges: end -->

Stateless random numbers for R, built on the Philox4x64-10 counter-based
generator from [Random123](https://github.com/DEShawResearch/random123).

Every value is a pure function of an immutable **key** plus the sampler's
arguments. Samplers never read or mutate `.Random.seed`, so results do not
depend on call order, on how many values you drew earlier, or on how many
threads the machine has.

## Installation

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
#> [1] 0.7901089 0.7292057 0.5056636 0.1168211 0.1449297

rng_uniform(key, 5L)   # same key, same draws
#> [1] 0.7901089 0.7292057 0.5056636 0.1168211 0.1449297
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
#>           [,1]      [,2]      [,3]
#> [1,] 0.7901089 0.1729620 0.7636065
#> [2,] 0.7292057 0.9580783 0.8891381
#> [3,] 0.5056636 0.9048795 0.1261067
#> [4,] 0.1168211 0.5774249 0.4911186
```

Or derive a key from structured data with `rng_fold()`, which is useful when
the pieces of a computation are named rather than numbered:

```r
rng_fold(key, "layer1")
#> <rng_key> rng_key[8dbacb9f407ebc3ae10ed1561177db53]

rng_uniform(rng_fold(key, "chain 3"), 3L)
#> [1] 0.2692752 0.3617336 0.7737181
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
#> [1] 0.2397225 0.6828912 0.2552504 1.5763296 0.7653972

rng_integer(key, 10L, min = 1L, max = 6L)
#>  [1] 2 2 1 5 5 2 3 3 2 2
```

`rng_bits()` exposes the generator stream directly. With `bits = 32` it
returns unsigned 32-bit words stored exactly in a double; with `bits = 64` it
returns fixed-width hex strings, because a double cannot hold 64 bits without
loss:

```r
rng_bits(key, 3L)
#> [1]  601333868 2746897238 3549456219

rng_bits(key, 3L, bits = 64L)
#> [1] "c0e6592123d7a06c" "0c89edc8a3ba5356" "6eb27f65d390675b"
```

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

On macOS, Apple clang ships without OpenMP, so a build there is
single-threaded unless `~/.R/Makevars` defines `SHLIB_OPENMP_CFLAGS` against
a `libomp` installation.

## License

MIT, see [LICENSE](LICENSE).

This package bundles two BSD 3-clause components, credited in
[inst/COPYRIGHTS](inst/COPYRIGHTS): the Random123 library (D. E. Shaw
Research) and NumPy's ziggurat constant tables (NumPy Developers).
