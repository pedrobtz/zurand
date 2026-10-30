# Create a stateless RNG key

`rng_key()` creates immutable keys for the stateless generators in this
package. A key vector is an opaque S3 object with one key per row; pass
it to
[`rng_fold()`](https://pedrobtz.github.io/zurand/reference/rng_fold.md)
and the `rng_*()` samplers, but do not rely on its internal integer
representation.

## Usage

``` r
rng_key(seed, n = 1L, engine = c("xoshiro256pp", "philox4x64", "threefry4x64"))

# S3 method for class 'rng_key'
format(x, ...)

# S3 method for class 'rng_key'
print(x, ...)

# S3 method for class 'rng_key'
length(x)

# S3 method for class 'rng_key'
x[i, ..., drop = FALSE]

# S3 method for class 'rng_key'
x[i, ...] <- value

# S3 method for class 'rng_key'
x[[i, ...]] <- value

# S3 method for class 'rng_key'
as.list(x, ...)

# S3 method for class 'rng_key'
rep(x, ...)

# S3 method for class 'rng_key'
duplicated(x, incomparables = FALSE, ...)

# S3 method for class 'rng_key'
unique(x, incomparables = FALSE, ...)

# S3 method for class 'rng_key'
xtfrm(x)

# S3 method for class 'rng_key'
x[[i, ...]]

# S3 method for class 'rng_key'
c(..., recursive = FALSE)
```

## Arguments

- seed:

  A single non-negative whole number.

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

- x:

  An `rng_key`.

- ...:

  Unused.

- i:

  Integer or logical index selecting keys. Missing values, fractional
  numbers, infinities and positions beyond the key vector are an error:
  R would turn such an index into a key of missing words, the same
  invented key whatever vector it came from, so unrelated subsets would
  silently share one stream.

- drop:

  Ignored; key subsetting always preserves the `rng_key` class.

- value:

  Not used: keys are immutable, so assigning into a key vector is an
  error. Without these methods R would write the value into individual
  key words, silently producing different keys.

- incomparables:

  Passed to [`duplicated()`](https://rdrr.io/r/base/duplicated.html).

- recursive:

  Ignored. Included for compatibility with
  [`c()`](https://rdrr.io/r/base/c.html).

## Value

An object of class `rng_key` containing `n` keys.

## Details

The samplers in this package never read or update `.Random.seed`.
Calling a sampler twice with the same key and arguments returns the same
values.

## Key vectors

A key vector behaves like a vector of opaque keys:
[`length()`](https://rdrr.io/r/base/length.html), `[`, `[[`,
[`c()`](https://rdrr.io/r/base/c.html),
[`rev()`](https://rdrr.io/r/base/rev.html),
[`rep()`](https://rdrr.io/r/base/rep.html),
[`unique()`](https://rdrr.io/r/base/unique.html),
[`duplicated()`](https://rdrr.io/r/base/duplicated.html) and
[`as.list()`](https://rdrr.io/r/base/list.html) work key by key and keep
the engine and stream version, so `lapply(keys, f)`,
[`Map()`](https://rdrr.io/r/base/funprog.html) and the parallel map
functions hand each task one key. A key vector takes one index, as a
vector does: `keys[1, 2]` and `keys[, 1]` are an error, not a view of
the words inside. [`c()`](https://rdrr.io/r/base/c.html) dispatches on
its first argument, so combine key vectors with a key first:
`c(NULL, key)` returns the bare integer words, not keys, while
`c(key, NULL)` drops the `NULL` as
[`c()`](https://rdrr.io/r/base/c.html) does for any vector. Keys are
immutable: assigning into a key vector is an error, and so is sorting
one, since keys have no order; compare keys with
[`identical()`](https://rdrr.io/r/base/identical.html).

## Stream version

A key also records the version of the stream definition it was made
under, as `attr(key, "stream")`; every key made by this version of
zurand is stream `1L`. If a defect in a stream is ever found after
release, the fix becomes a new stream version for new keys, and existing
keys – including saved ones – keep producing exactly the values they
always did. A key saved before the attribute existed is treated as
stream 1. A key from a newer stream than this version of zurand
implements is an error rather than a silent reinterpretation.

## Examples

``` r
key <- rng_key(42L)
keys <- rng_key(42L, n = 4L)
key
#> <rng_key> rng_key[2feb6e95bdd73226b266f10328efe333]
```
