# Draw uniform random values

Stateless uniform sampling. `rng_uniform(key, n, min, max)` is a pure
function of its arguments and returns ordinary numeric vectors or
matrices.

## Usage

``` r
rng_uniform(key, n = 1L, min = 0, max = 1, offset = 0)
```

## Arguments

- key:

  An `rng_key` vector.

- n:

  A single non-negative whole number. Defaults to one draw per key.

- min, max:

  Single finite numeric bounds.

- offset:

  A single non-negative whole number: the draws start at position
  `offset` of each key's stream, so the result holds positions `offset`
  to `offset + n - 1`. `f(key, n, offset = k)` equals
  `f(key, k + n)[(k + 1):(k + n)]` exactly, which lets a long simulation
  draw its next batch, or resume from a checkpoint, without holding the
  earlier values. `offset + n` must not exceed \\2^{53}\\.

## Value

A numeric vector of length `n` for one key, or an `n` by `length(key)`
numeric matrix for multiple keys.

## Details

With the default bounds every value lies strictly inside \\(0, 1)\\: it
is \\(m + 1/2) 2^{-52}\\ for a 52-bit integer \\m\\, so neither endpoint
can occur. Other bounds are applied as `min + (max - min) * u`, which,
as with [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html), will
not return either extreme value unless `max == min` or `max - min` is
small compared with `min`.

## Examples

``` r
key <- rng_key(42L)
rng_uniform(key)
#> [1] 0.6436165
rng_uniform(key, 5L)
#> [1] 0.6436165 0.4170741 0.1279353 0.8522212 0.3491436
rng_uniform(rng_key(42L, n = 3L), 5L)
#>           [,1]       [,2]      [,3]
#> [1,] 0.6436165 0.85879755 0.3794245
#> [2,] 0.4170741 0.60044715 0.6979140
#> [3,] 0.1279353 0.09813149 0.6872125
#> [4,] 0.8522212 0.77553624 0.4570548
#> [5,] 0.3491436 0.07162715 0.9324355
# the next five, on their own
identical(rng_uniform(key, 5L, offset = 5), rng_uniform(key, 10L)[6:10])
#> [1] TRUE
```
