# Get started with zurand

``` r

library(zurand)
```

R’s own generator keeps a hidden state.
[`set.seed()`](https://rdrr.io/r/base/Random.html) fixes where it
starts, and every call to
[`runif()`](https://rdrr.io/r/stats/Uniform.html) or
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) moves it on, so each
number depends on everything drawn before it. zurand has no such state.
In this vignette we bootstrap the mean of a small sample, and along the
way meet keys, the samplers, key vectors and folding. Why the values are
the same on every machine, and how fast they are, is left to the
articles.

## Make a key

A **key** is an immutable value made from a seed. Every random value is
a pure function of a key and the value’s **position**, so the same key
always gives the same numbers:

``` r

key <- rng_key(2026)
key
#> <rng_key> rng_key[91948d23db9c5598ed35455d78bc927d]
rng_uniform(key, 3)
#> [1] 0.15037421 0.91830293 0.06342998
rng_uniform(key, 3)
#> [1] 0.15037421 0.91830293 0.06342998
```

Asking for more values extends the sequence rather than changing it. The
first three of six are the three above:

``` r

rng_uniform(key, 6)
#> [1] 0.15037421 0.91830293 0.06342998 0.74915105 0.60528900 0.30633444
```

Drawing uses nothing up: the key is the same afterwards, so there is no
state to pass between functions, save or restore.

## Draw from a key

The samplers are
[`rng_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_uniform.md),
[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md),
[`rng_integer()`](https://pedrobtz.github.io/zurand/reference/rng_integer.md)
and
[`rng_bits()`](https://pedrobtz.github.io/zurand/reference/rng_bits.md).
Each takes a key and a count, then the parameters of its distribution.
Here is the sample we will bootstrap, and one resample of it, drawn as
positions between 1 and 8:

``` r

x <- c(2.1, 3.4, 1.9, 5.6, 4.2, 3.3, 2.8, 4.9)
i <- rng_integer(key, length(x), min = 1, max = length(x))
i
#> [1] 3 7 2 7 2 8 5 2
mean(x[i])
#> [1] 3.35
```

[`rng_normal()`](https://pedrobtz.github.io/zurand/reference/rng_normal.md)
takes `mean` and `sd`, and
[`rng_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_uniform.md)
takes `min` and `max`, in the same way.

## One key per replicate

A bootstrap needs many resamples, each from its own stream.
[`rng_key()`](https://pedrobtz.github.io/zurand/reference/rng_key.md)
makes a **key vector** when asked for `n` keys, and a sampler given a
key vector returns a matrix with one column per key:

``` r

keys <- rng_key(2026, n = 200)
length(keys)
#> [1] 200
idx <- rng_integer(keys, length(x), min = 1, max = length(x))
dim(idx)
#> [1]   8 200
```

Each column is one resample, so the bootstrap distribution of the mean
is one line away:

``` r

boot <- apply(idx, 2, function(i) mean(x[i]))
quantile(boot, c(0.025, 0.975))
#>     2.5%    97.5% 
#> 2.799062 4.451250
```

A key vector can be indexed like any vector. Replicate 7 alone gives the
same resample as column 7 of the matrix, so a replicate can be rerun
without the others:

``` r

identical(rng_integer(keys[7], length(x), min = 1, max = length(x)), idx[, 7])
#> [1] TRUE
```

## Name the parts with folding

Counting keys works while a computation has one part. When it has
several, derive each part’s key from what the part *is*.
[`rng_fold()`](https://pedrobtz.github.io/zurand/reference/rng_fold.md)
takes a key and some data, a name or a number, and returns a new key;
the same data always gives the same key:

``` r

boot_key <- rng_fold(key, "bootstrap")
replicate_key <- function(b) rng_fold(boot_key, b)
rng_uniform(replicate_key(3), 2)
#> [1] 0.3819716 0.5885736
rng_uniform(replicate_key(3), 2)
#> [1] 0.3819716 0.5885736
```

Because replicate `b` gets its numbers from its own key, the bootstrap
gives the same answer whatever order the replicates run in. Here they
run once in order and once shuffled, as parallel workers might run them
([`sample()`](https://rdrr.io/r/base/sample.html) only picks the order):

``` r

boot_mean <- function(b) {
  i <- rng_integer(replicate_key(b), length(x), min = 1, max = length(x))
  mean(x[i])
}
in_order <- vapply(1:200, boot_mean, numeric(1))

shuffled <- numeric(200)
run_order <- sample(200)
shuffled[run_order] <- vapply(run_order, boot_mean, numeric(1))
identical(in_order, shuffled)
#> [1] TRUE
```

A new part of the analysis, folded under another name, cannot change
these values. The same holds for zurand’s own threads:
[`rng_threads()`](https://pedrobtz.github.io/zurand/reference/rng_threads.md)
changes how fast a large draw is made, never what it contains.

## Vectors of any length

[`rng_lazy_normal()`](https://pedrobtz.github.io/zurand/reference/rng_lazy_normal.md)
and
[`rng_lazy_uniform()`](https://pedrobtz.github.io/zurand/reference/rng_lazy_normal.md)
return a **lazy vector**: it holds only its key and arguments, and
computes a value when it is read. A billion values are created at once,
and they are the values the ordinary sampler gives:

``` r

z <- rng_lazy_normal(key, 1e9)
length(z)
#> [1] 1000000000
z[1:3]
#> [1]  0.006686381 -1.233412249  0.591067060
rng_normal(key, 3)
#> [1]  0.006686381 -1.233412249  0.591067060
```

## A common mistake

The samplers take a key, not a seed. Passing the seed itself is an error
that says what to do instead:

``` r

rng_normal(2026, 3)
#> Error in `rng_normal()`:
#> ! `key` must be an <rng_key>; create one with rng_key()
```

## Where next

- [Run a reproducible simulation
  study](https://pedrobtz.github.io/zurand/articles/simulations.html):
  scenarios, parallel workers, common random numbers, and stopping and
  resuming with `offset`.
- [Work with vectors too big for
  memory](https://pedrobtz.github.io/zurand/articles/lazy.html): lazy
  vectors, saving them, and when they are filled in.
- [Reproducibility
  guarantees](https://pedrobtz.github.io/zurand/articles/reproducibility.html):
  the engines, what stays bit-identical, and stream versions.
- [Performance](https://pedrobtz.github.io/zurand/articles/performance.html)
  and [Use zurand from
  C](https://pedrobtz.github.io/zurand/articles/c-api.html).
- The [reference index](https://pedrobtz.github.io/zurand/reference/)
  has every function and argument.
