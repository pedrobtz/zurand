# zurand

## Overview

Random numbers in R come from a hidden, global state.
[`set.seed()`](https://rdrr.io/r/base/Random.html) fixes where it
starts, and every draw moves it on. So the numbers you get depend on
everything drawn before: add one
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) early in a script and
every later result changes, and split a loop across workers and the
results depend on how the work was divided.

zurand has no hidden state. You make a **key**, and each value is worked
out from the key and the value’s position, so the same key always gives
the same numbers.

- **The same numbers everywhere.** A key gives the same values on Linux,
  macOS and Windows, on Intel and ARM, and whatever the number of
  threads.
- **Nothing to keep track of.** Drawing never changes anything, so
  drawing more, or drawing in a different order, leaves every other
  result alone.
- **Parallel work without worry.** Give each task its own key, and the
  results no longer depend on how the tasks are scheduled.
- **Fast.** Several times faster than
  [`runif()`](https://rdrr.io/r/stats/Uniform.html) and
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) on one thread, and it
  uses more threads for large draws without changing a single value.
- **Huge vectors for free.** A lazy vector of a billion values takes
  almost no memory: its values are computed only when you read them.

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

``` r

library(zurand)
```

A key gives the same numbers every time you use it:

``` r

key <- rng_key(2026)
rng_normal(key, 3)
#> [1]  0.006686381 -1.233412249  0.591067060
rng_normal(key, 3)
#> [1]  0.006686381 -1.233412249  0.591067060
```

Asking for more values extends the sequence rather than changing it:

``` r

rng_normal(key, 5)
#> [1]  0.006686381 -1.233412249  0.591067060  0.835300390  0.249764157
```

For a simulation with several replicates, make one key per replicate.
Each column holds the stream of one replicate:

``` r

keys <- rng_key(2026, n = 3)
rng_uniform(keys, 4)
#>            [,1]       [,2]      [,3]
#> [1,] 0.15037421 0.02841554 0.4769039
#> [2,] 0.91830293 0.65836513 0.1586149
#> [3,] 0.06342998 0.48395427 0.8484967
#> [4,] 0.74915105 0.79447021 0.9430020
```

Replicate 3 gets the same numbers whether or not replicates 1 and 2 ever
run, on this machine or another:

``` r

rng_uniform(keys[3], 4)
#> [1] 0.4769039 0.1586149 0.8484967 0.9430020
```

When the parts of a computation have names rather than numbers, derive a
key from the name:

``` r

rng_integer(rng_fold(key, "treatment"), 6, min = 1, max = 6)
#> [1] 6 3 5 6 4 3
rng_integer(rng_fold(key, "control"), 6, min = 1, max = 6)
#> [1] 6 6 6 3 5 5
```

A lazy vector is created instantly, whatever its length, and computes
only the values you read. They are the same values
[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md)
gave above:

``` r

x <- rng_lazy_normal(key, 1e9)
length(x)
#> [1] 1000000000
x[c(1, 2, 3)]
#> [1]  0.006686381 -1.233412249  0.591067060
```

Every function is named `rng_*()`:
[`rng_key()`](https://pedrobtz.github.io/zurand/reference/rng_key.md)
and
[`rng_fold()`](https://pedrobtz.github.io/zurand/reference/rng_fold.md)
make keys, and
[`rng_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_uniform.md),
[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md),
[`rng_integer()`](https://pedrobtz.github.io/zurand/reference/rng_integer.md)
and
[`rng_bits()`](https://pedrobtz.github.io/zurand/reference/rng_bits.md)
draw from them.

## Learn more

- [Stateless random
  numbers](https://pedrobtz.github.io/zurand/articles/zurand.html): key
  vectors, folding and engines.
- [Reproducible simulation
  studies](https://pedrobtz.github.io/zurand/articles/simulations.html):
  studies that stay reproducible as they grow, run in parallel, and stop
  and resume.
- [Lazy random
  vectors](https://pedrobtz.github.io/zurand/articles/lazy.html): random
  vectors too big for memory.
- [Performance](https://pedrobtz.github.io/zurand/articles/performance.html):
  how fast zurand is against other generators, and how threads are used.
- [Using zurand from
  C](https://pedrobtz.github.io/zurand/articles/c-api.html): the same
  values from C, in other packages.
- The [reference index](https://pedrobtz.github.io/zurand/reference/)
  lists every function, and the
  [website](https://pedrobtz.github.io/zurand/) has all of it.

## Acknowledgements

zurand builds on the Philox and Threefry generators of
[Random123](https://github.com/DEShawResearch/random123) (Salmon,
Moraes, Dror and Shaw, 2011), the xoshiro256++ generator of Blackman and
Vigna (2021), the ziggurat method of Marsaglia and Tsang (2000) with
NumPy’s tables, McFarland’s (2016) modified ziggurat, and fdlibm’s
[`exp()`](https://rdrr.io/r/base/Log.html) and
[`log1p()`](https://rdrr.io/r/base/Log.html). The licences of the code
it bundles are in
[inst/COPYRIGHTS](https://github.com/pedrobtz/zurand/blob/main/inst/COPYRIGHTS).
zurand itself is MIT licensed.
