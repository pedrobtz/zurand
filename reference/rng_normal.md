# Draw normal random values

Stateless normal sampling. `rng_normal(key, n, mean, sd)` is a pure
function of its arguments and returns ordinary numeric vectors or
matrices.

## Usage

``` r
rng_normal(
  key,
  n = 1L,
  mean = 0,
  sd = 1,
  method = c("ziggurat", "mcfarland"),
  offset = 0
)
```

## Arguments

- key:

  An `rng_key` vector.

- n:

  A single non-negative whole number. Defaults to one draw per key.

- mean:

  A single finite numeric mean.

- sd:

  A single non-negative finite numeric standard deviation.

- method:

  The sampling algorithm:

  - `"ziggurat"` (the default): Marsaglia and Tsang's ziggurat on
    NumPy's 256-layer tables.

  - `"mcfarland"`: McFarland's modified ziggurat (2016). Its common case
    needs no table comparison, so it is faster than the default on every
    CPU measured, by more on x86_64 than on arm64
    ([`vignette("performance", package = "zurand")`](https://pedrobtz.github.io/zurand/articles/performance.md)).

  Both are exact samplers of the standard normal (to double precision)
  and both are bit-reproducible across platforms. They are different
  streams: the same key gives unrelated values under the two methods,
  and changing the method changes every value. The default never
  changes.

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

## References

Marsaglia, G. and Tsang, W. W. (2000). The ziggurat method for
generating random variables. *Journal of Statistical Software*, 5(8),
1-7. [doi:10.18637/jss.v005.i08](https://doi.org/10.18637/jss.v005.i08)

McFarland, C. D. (2016). A modified ziggurat algorithm for generating
exponentially and normally distributed pseudorandom numbers. *Journal of
Statistical Computation and Simulation*, 86(7), 1281-1294.
[doi:10.1080/00949655.2015.1060234](https://doi.org/10.1080/00949655.2015.1060234)

## Examples

``` r
key <- rng_key(42L)
rng_normal(key)
#> [1] 1.630614
rng_normal(key, 5L)
#> [1]  1.6306140 -1.1186343 -0.6639268  0.2393690 -0.1899009
rng_normal(key, 5L, method = "mcfarland")
#> [1] -0.111325527  0.434353516  1.234857312 -0.557763925  0.008950617
identical(rng_normal(key, 5L, offset = 5), rng_normal(key, 10L)[6:10])
#> [1] TRUE
```
