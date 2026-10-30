# Draw integer random values

Stateless integer sampling over the inclusive range `[min, max]`.

## Usage

``` r
rng_integer(key, n = 1L, min, max)
```

## Arguments

- key:

  An `rng_key` vector.

- n:

  A single non-negative whole number. Defaults to one draw per key.

- min, max:

  Single whole-number bounds in R's non-missing integer range.

## Value

An integer vector of length `n` for one key, or an `n` by `length(key)`
integer matrix for multiple keys.

## Examples

``` r
key <- rng_key(42L)
rng_integer(key, min = 1L, max = 6L)
#> [1] 3
rng_integer(key, 10L, min = 1L, max = 6L)
#>  [1] 3 1 4 2 5 3 3 1 1 5
```
