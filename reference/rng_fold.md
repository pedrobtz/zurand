# Fold data into a key

`rng_fold()` deterministically derives a key from `key` and simple data.
This is useful for structured reproducibility, for example deriving a
key from names such as `"chain 3"` or counters such as `10L`.

## Usage

``` r
rng_fold(key, data)
```

## Arguments

- key:

  An `rng_key` vector.

- data:

  A character, integer, double, logical or raw vector. Numeric values
  must be whole numbers; integer and double vectors holding the same
  values (`1L` and `1`) derive the same key. All other type distinctions
  matter: `rng_fold(key, "1")` and `rng_fold(key, 1)` differ. Only the
  storage type and the values are folded, so classed objects such as
  factors and dates are an error: a factor would fold as its level
  codes, whose meaning changes with the level set. Fold
  `as.character(f)` or `unclass(d)` instead. Names and dimensions are
  ignored: a matrix folds as its values in column order. Strings are
  folded as their UTF-8 bytes, whatever their declared encoding.

## Value

An `rng_key` vector with the same length as `key`.

## Details

The data are hashed to 64 bits before they derive the key, so two
different values can in principle fold to the same key. Among `m`
distinct values the chance of any such collision is about `m^2 / 2^65`:
under one in a million for a million values, but not negligible near
`2^32` of them. For that many independent streams, request them up front
with `rng_key(seed, n)`, which never collides.

## Examples

``` r
key <- rng_key(42L)
layer_key <- rng_fold(key, "layer1")
```
