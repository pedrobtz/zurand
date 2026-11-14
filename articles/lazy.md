# Lazy random vectors

``` r

library(zurand)
```

In zurand every random value is a function of its key and its position.
So a vector of random values does not have to exist before it is read:
its values can be worked out when they are needed. That is a lazy
vector.

## Created instantly, whatever its length

A lazy vector of a billion normal values is created at once. It holds
its key, its arguments and a small cache, at most about 40 KB, where an
ordinary vector of that length would need 8 GB:

``` r

key <- rng_key(2026)
x <- rng_lazy_normal(key, 1e9)
length(x)
#> [1] 1000000000
```

([`object.size()`](https://rdrr.io/r/utils/object.size.html) reports the
8 GB an ordinary vector would take: it counts the length, and cannot see
that the values are not there.)

It behaves like any numeric vector. Reading part of it computes only
that part:

``` r

x[1:5]
#> [1]  0.006686381 -1.233412249  0.591067060  0.835300390  0.249764157
x[c(1e3, 1e6, 1e9)]
#> [1] -0.06916703  0.23749869  1.72024893
```

## The same values as an ordinary draw

Element `i` of a lazy vector is exactly element `i` of the matching
draw, so the two can be used interchangeably:

``` r

y <- rng_lazy_normal(key, 1e5)
identical(y[], rng_normal(key, 1e5))
#> [1] TRUE
identical(sum(y), sum(rng_normal(key, 1e5)))
#> [1] TRUE
```

[`rng_lazy_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_lazy_normal.md)
does the same for uniform values:

``` r

u <- rng_lazy_uniform(key, 1e9, min = 10, max = 20)
u[1:3]
#> [1] 11.50374 19.18303 10.63430
identical(u[1:3], rng_uniform(key, 3, min = 10, max = 20))
#> [1] TRUE
```

## Small files

Saving a lazy vector stores only its key and arguments, a few hundred
bytes, and it comes back lazy:

``` r

f <- tempfile(fileext = ".rds")
saveRDS(x, f)
file.size(f)
#> [1] 198
identical(readRDS(f)[1:5], x[1:5])
#> [1] TRUE
```

The file needs zurand to be read back.

## Sums and scans in constant memory

[`sum()`](https://rdrr.io/r/base/sum.html),
[`mean()`](https://rdrr.io/r/base/mean.html) and
[`max()`](https://rdrr.io/r/base/Extremes.html) read a lazy vector in
pieces, so they work on vectors far larger than memory. They take as
long as generating the values, so this example uses ten million:

``` r

z <- rng_lazy_normal(key, 1e7)
mean(z)
#> [1] -0.0002918956
```

## When the vector is filled in

Anything that needs the whole vector as ordinary memory fills it in
once: arithmetic such as `z * 2`, changing an element, or C code asking
for its data. From then on it is an ordinary vector of that length, so
do this only with vectors that fit in memory:

``` r

w <- rng_lazy_normal(key, 10)
w2 <- w * 2             # fills in w, then doubles it
w[3] <- 0               # an ordinary change to an ordinary vector
w[1:4]
#> [1]  0.006686381 -1.233412249  0.000000000  0.835300390
```

A changed lazy vector saves its values rather than its recipe, so saving
never loses a change.

## Scattered reads

For reading values at random positions, such as `x[sample(n, k)]`, use a
key with the `"philox4x64"` engine. It computes any single value
directly, which makes scattered reads several times faster than with the
default engine:

``` r

p <- rng_lazy_normal(rng_key(2026, engine = "philox4x64"), 1e9)
p[c(7, 7e8, 42)]
#> [1] -0.1846474 -1.3280065 -2.6424239
```

For scans and reading in order, keep the default engine.

## When to use them

Lazy vectors help when a vector is too big to hold, when only part of it
will be read, or when it is saved and shared. For vectors that fit in
memory and are read in full, an ordinary draw with
[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md)
or
[`rng_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_uniform.md)
is simpler and faster. See
[`?rng_lazy_normal`](https://pedrobtz.github.io/zurand/reference/rng_lazy_normal.md)
for the details.
