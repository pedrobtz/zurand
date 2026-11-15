# zurand

zurand is a random number generator for R in which every value is a pure
function of a key and its position, so drawing never reads or changes
`.Random.seed`.

- A key gives the same values on Linux, macOS and Windows, on x86_64 and
  ARM, and at any thread count.
- Drawing more values, or drawing in another order, leaves every other
  result unchanged, so parallel tasks need no seeds.
- On one thread it is 8 to 17 times faster than
  [`runif()`](https://rdrr.io/r/stats/Uniform.html) and
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html), and large draws use
  several threads without changing a value.

It builds on the counter-based generators of
[Random123](https://github.com/DEShawResearch/random123) and on
xoshiro256++.

## Installation

Install the released version from CRAN:

``` r

install.packages("zurand")
```

Or the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("pedrobtz/zurand")
```

## Usage

A key gives the same numbers every time you use it, and asking for more
extends the sequence rather than changing it:

``` r

library(zurand)

key <- rng_key(2026)
rng_normal(key, 3)
#> [1]  0.006686381 -1.233412249  0.591067060
rng_normal(key, 5)
#> [1]  0.006686381 -1.233412249  0.591067060  0.835300390  0.249764157
```

For several independent streams, such as the replicates of a simulation,
make a vector of keys. Each column is the stream of one key:

``` r

keys <- rng_key(2026, n = 3)
rng_normal(keys, 4)
#>              [,1]       [,2]        [,3]
#> [1,]  0.006686381 -2.2711494  0.42616283
#> [2,] -1.233412249 -0.2468333  0.05688975
#> [3,]  0.591067060  0.9256739 -1.02992948
#> [4,]  0.835300390  1.1623784  0.68884687
```

zurand also draws uniform values and integers, derives keys from names,
and makes lazy vectors of any length; Get started shows how.

## Learn more

- [Get started](https://pedrobtz.github.io/zurand/articles/zurand.html):
  keys, samplers and folding, through one small bootstrap.
- [Run a reproducible simulation
  study](https://pedrobtz.github.io/zurand/articles/simulations.html):
  scenarios, parallel workers, and stopping and resuming.
- [Reproducibility
  guarantees](https://pedrobtz.github.io/zurand/articles/reproducibility.html):
  what stays the same across platforms, engines and versions.
- [Performance](https://pedrobtz.github.io/zurand/articles/performance.html):
  measured against base R, dqrng and randompack.
- The [reference index](https://pedrobtz.github.io/zurand/reference/)
  lists every function.

The licences of the code zurand bundles, from Random123, NumPy and
fdlibm, are in
[inst/COPYRIGHTS](https://github.com/pedrobtz/zurand/blob/main/inst/COPYRIGHTS).
