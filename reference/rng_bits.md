# Draw low-level random bits

`rng_bits()` exposes the raw generator stream. With `bits = 32`, it
returns unsigned 32-bit words stored exactly in a numeric vector. With
`bits = 64`, it returns fixed-width hexadecimal strings so no bits are
lost to R's numeric representation.

## Usage

``` r
rng_bits(key, n = 1L, bits = 32L)
```

## Arguments

- key:

  An `rng_key` vector.

- n:

  A single non-negative whole number. Defaults to one draw per key.

- bits:

  Either `32L` or `64L`.

## Value

A numeric vector or matrix for 32-bit words, or a character vector or
matrix for 64-bit words.

## Examples

``` r
key <- rng_key(42L)
rng_bits(key)
#> [1] 2206215613
rng_bits(key, 5L)
#> [1] 2206215613  289400228  609344871  751839911  715752478
```
