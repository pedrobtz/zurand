# Create a key from R's RNG

`rng_key_from_r()` is the one convenience function in this package that
intentionally consumes R's global RNG state. Use it when you want a
fresh stateless key from the current `.Random.seed`; all other functions
in this package leave `.Random.seed` untouched.

## Usage

``` r
rng_key_from_r(
  n = 1L,
  engine = c("xoshiro256pp", "philox4x64", "threefry4x64")
)
```

## Arguments

- n:

  A single non-negative whole number: how many keys to create.

- engine:

  The engine. All three give the same guarantees: a value is a pure
  function of the key and the draw index, the same key gives the same
  values on every platform and thread count, and draw `i` does not
  depend on `n`. They differ in speed and in the stream they produce.

  - `"xoshiro256pp"` (default) – the fastest of the three, and faster
    than dqrng and RcppZiggurat on one thread (x86_64 and Apple Silicon;
    [`vignette("performance", package = "zurand")`](https://pedrobtz.github.io/zurand/articles/performance.md)
    has the numbers). Philox derives a fresh 256-bit xoshiro256++ state
    for each 512-value chunk, then a cheap recurrence produces the
    chunk, so a chunk depends only on the key and its index. Reaching an
    arbitrary index costs at most 511 recurrence steps. Audited with
    PractRand to 1 TB, including interleaved streams of sibling and
    folded keys.

  - `"philox4x64"` – Philox4x64-10 from Random123, counter-based: each
    value is computed directly from its position.

  - `"threefry4x64"` – Threefry4x64-13 from Random123, counter-based;
    close to Philox in speed.

  A key records its engine, so a key made with one never produces the
  other's values. Changing engine changes every number you get.

## Value

An object of class `rng_key` containing `n` keys.

## Examples

``` r
set.seed(1)
key <- rng_key_from_r()
rng_uniform(key, 3L)
#> [1] 0.2330005 0.8687607 0.5942335

# The same seed gives the same key.
set.seed(1)
identical(rng_key_from_r(), key)
#> [1] TRUE
```
