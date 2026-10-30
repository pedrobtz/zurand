# Lazy random vectors

`rng_lazy_uniform()` and `rng_lazy_normal()` return a numeric vector of
length `n` whose elements are computed only when read. Every value is a
pure function of the key and its position, so element `i` is exactly
element `i` of `rng_uniform(key, n, ...)` or `rng_normal(key, n, ...)`,
and `identical(x[], rng_normal(key, n))` holds. Creating one is instant
and takes no memory for the values, whatever `n` is: each lazy vector
holds its key, its arguments and a cache of one engine chunk, which is
`min(n, 5120)` doubles, so at most about 40 KB.

## Usage

``` r
rng_lazy_normal(key, n, mean = 0, sd = 1, method = c("ziggurat", "mcfarland"))

rng_lazy_uniform(key, n, min = 0, max = 1)
```

## Arguments

- key:

  A single `rng_key`.

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

- min, max:

  Single finite numeric bounds.

## Value

A double vector of length `n` (an ALTREP object).

## Reading

`x[i]`, `x[idx]`, `x[a:b]`, `sum(x)`, `mean(x)` and loops over `x`
compute only the values they touch; nothing is kept but the last engine
chunk read (at most 5120 values), so
[`sum()`](https://rdrr.io/r/base/sum.html) over 1e9 lazy values runs in
the memory R itself uses. Full scans cost more than
[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md)
followed by [`sum()`](https://rdrr.io/r/base/sum.html) would, because R
reads lazy vectors 512 values at a time on one thread. The gain is
memory, and not paying for values that are never read.

## Interrupts

Scans that read a lazy vector in blocks –
[`sum()`](https://rdrr.io/r/base/sum.html),
[`mean()`](https://rdrr.io/r/base/mean.html),
[`max()`](https://rdrr.io/r/base/Extremes.html), `x[idx]` – check for a
user interrupt every \\2^{22}\\ values, so Ctrl-C stops
[`sum()`](https://rdrr.io/r/base/sum.html) over 1e12 values. Reading
single elements and materialising do not check. For C and C++ code this
means that `REAL_ELT()` and `REAL()`/`REAL_RO()` on a lazy vector never
signal an interrupt, but `REAL_GET_REGION()` – which cpp11's `doubles`
iterator uses – can, by a longjmp that skips C++ destructors. Code that
cannot be unwound should take the data pointer first (which materialises
the vector) or read under `R_UnwindProtect()`
(`cpp11::unwind_protect()`).

## Materialisation

Anything that needs the whole vector as ordinary memory – most
arithmetic such as `x * 2`, modifying `x`, or C code asking for its data
pointer – fills it once, with the usual parallel fill, and from then on
`x` holds all `n` values like any vector.

## Saving

[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`serialize()`](https://rdrr.io/r/base/serialize.html) store only the
key and the arguments, so a lazy vector of 1e9 values saves as a few
hundred bytes and comes back lazy. A lazy vector that has been modified
saves its values. Reading is not modifying: after C code has had the
data pointer ([`sd()`](https://rdrr.io/r/stats/sd.html),
[`sort()`](https://rdrr.io/r/base/sort.html) and
[`cumsum()`](https://rdrr.io/r/base/cumsum.html) all take it), saving
regenerates the values and compares them with what the vector holds,
which costs about one more fill, and stores only the key when nothing
changed. Restoring checks the saved arguments and refuses ones it cannot
have written.

## Engines

The vector uses its key's engine, so its values are that engine's. For
scattered reads, `x[sample(n, k)]` or random indices in a loop, key it
with `rng_key(seed, engine = "philox4x64")`: Philox computes any single
value directly, several times faster per read than `"xoshiro256pp"`,
which reaches it through its 512-value block. For scans and sequential
loops the default is faster.

## Examples

``` r
key <- rng_key(42L)
x <- rng_lazy_normal(key, 1e9)      # instant; no values yet
x[c(1, 5e8, 1e9)]
#> [1]  1.6306140 -0.1649847  1.0731953
identical(rng_lazy_normal(key, 10)[], rng_normal(key, 10))
#> [1] TRUE
```
