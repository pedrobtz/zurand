# Report or disable the SIMD path

`rng_simd()` reports which instruction path the `"xoshiro256pp"` engine
is using; `rng_simd(FALSE)` forces the portable one.

## Usage

``` r
rng_simd(enable = NULL)
```

## Arguments

- enable:

  `NULL` to query, `FALSE` to force the portable path, or `TRUE` to go
  back to using whatever the CPU supports.

## Value

A string: `"avx2"` (x86_64 with AVX2) or `"neon"` (arm64) if the
vectorised path is active, otherwise `"none"`. When setting, the
previous value, invisibly.

## Details

This is a performance control and nothing else. Both paths emit the same
values: the vectorised one runs several independent 512-value sub-chunks
side by side – four per AVX2 register on x86_64, two per NEON register
on arm64 – each executing the same recurrence as the scalar path, so
output does not depend on the instruction set, the CPU, or this setting.
It is exposed because that claim is worth being able to check, and the
package's own tests check it by running both and comparing.

Only `"xoshiro256pp"` has a vectorised path. The counter-based engines
are unaffected, as is every other function.

## Examples

``` r
rng_simd()
#> [1] "avx2"
key <- rng_key(1L, engine = "xoshiro256pp")
old <- rng_simd(FALSE)
portable <- rng_uniform(key, 4L)
rng_simd(old != "none")
identical(portable, rng_uniform(key, 4L))
#> [1] TRUE
```
