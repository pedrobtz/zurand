# Reproducible simulation studies

``` r

library(zurand)
```

A simulation study runs the same experiment many times, under several
scenarios, and compares the results. This article shows how zurand keeps
such a study reproducible while it grows, while it runs in parallel and
while it is stopped and restarted.

## The problem with one global seed

With R’s own generator, every number depends on how many were drawn
before it. Here an analysis draws ten values, and then someone adds a
single extra draw at the start:

``` r

set.seed(1)
before <- rnorm(10)

set.seed(1)
extra <- runif(1)          # one new draw, early in the script
after <- rnorm(10)

identical(before, after)
#> [1] FALSE
```

Every later result has changed. With zurand, each piece of work draws
from its own key, so adding work elsewhere cannot disturb it:

``` r

study <- rng_key(2026)
before <- rng_normal(rng_fold(study, "analysis"), 10)

extra <- rng_uniform(rng_fold(study, "new step"), 1)
after <- rng_normal(rng_fold(study, "analysis"), 10)

identical(before, after)
#> [1] TRUE
```

## A key for every part of the study

A study has a natural shape: scenarios, and replicates within each. Fold
the study key by the scenario’s name, then by the number of the
replicate, and every replicate of every scenario has its own stream:

``` r

replicate_key <- function(scenario, r) {
  rng_fold(rng_fold(study, scenario), r)
}

one_run <- function(scenario, r, n = 30) {
  key <- replicate_key(scenario, r)
  shift <- c(null = 0, small = 0.3, large = 0.8)[[scenario]]
  x <- rng_normal(key, n, mean = shift)
  t.test(x)$p.value < 0.05
}

power <- function(scenario, reps = 200) {
  mean(vapply(seq_len(reps), function(r) one_run(scenario, r), logical(1)))
}

results <- c(null = power("null"), small = power("small"))
results
#>  null small 
#> 0.025 0.360
```

Adding a scenario later leaves the earlier ones exactly as they were,
because no scenario’s numbers depend on another’s:

``` r

results2 <- c(null = power("null"), small = power("small"),
              large = power("large"))
identical(results, results2[c("null", "small")])
#> [1] TRUE
results2
#>  null small large 
#> 0.025 0.360 0.990
```

The same holds for more replicates: replicates 1 to 200 are the same
whether the study runs 200 or 1,000 of them.

## Running in parallel

Because the numbers of each replicate come from its own key, it does not
matter which worker runs it, or in which order. No worker needs a seed:

``` r

cores <- if (.Platform$OS.type == "windows") 1L else 2L
serial <- vapply(1:100, function(r) one_run("small", r), logical(1))
parallel <- unlist(parallel::mclapply(
  1:100, function(r) one_run("small", r), mc.cores = cores
))
identical(serial, parallel)
#> [1] TRUE
```

`future.apply::future_lapply()` or any other parallel map gives the same
result. zurand’s own threads work the same way:
[`rng_threads()`](https://pedrobtz.github.io/zurand/reference/rng_threads.md)
changes how fast a large draw is made, never its values.

## Common random numbers

To compare two methods fairly, give them the same random inputs. The
difference between them then reflects the methods, not the luck of the
draw. With keys this is just a matter of using the same key for both.

Here two estimators of a centre, the mean and the median, are compared
on data with an occasional outlier:

``` r

make_data <- function(key, n = 25) {
  x <- rng_normal(key, n)
  outlier <- rng_uniform(rng_fold(key, "outlier"), n) < 0.1
  x + outlier * 8
}

compare <- function(r, shared) {
  k1 <- replicate_key("estimators", r)
  k2 <- if (shared) k1 else rng_fold(k1, "second")
  abs(mean(make_data(k1))) - abs(median(make_data(k2)))
}

diff_shared <- vapply(1:300, compare, numeric(1), shared = TRUE)
diff_indep <- vapply(1:300, compare, numeric(1), shared = FALSE)

rbind(
  shared = c(mean = mean(diff_shared), sd = sd(diff_shared)),
  independent = c(mean = mean(diff_indep), sd = sd(diff_indep))
)
#>                  mean        sd
#> shared      0.5335896 0.4433732
#> independent 0.5308661 0.5043226
```

Both designs estimate the same average difference, the extra error of
the mean from the outliers, but with shared keys the difference in each
replicate is less noisy (a smaller `sd`), so the same precision needs
fewer replicates.

## Stopping and resuming

A long simulation can draw its random numbers in pieces with `offset`.
Position `k` of a key’s stream is the same however it is reached, so
drawing in batches gives exactly what one long draw would:

``` r

key <- rng_fold(study, "long run")
all_at_once <- rng_normal(key, 1000)

batches <- lapply(0:3, function(b) rng_normal(key, 250, offset = 250 * b))
identical(all_at_once, unlist(batches))
#> [1] TRUE
```

To resume after a stop, all a checkpoint needs to record is how far the
run got; the random numbers for the rest follow from the key:

``` r

done <- 600                                     # saved at the checkpoint
rest <- rng_normal(key, 1000 - done, offset = done)
identical(rest, all_at_once[(done + 1):1000])
#> [1] TRUE
```

## Saving keys

A key is a small R object. Save it with the results, and anyone can
regenerate the exact random numbers later, on any machine:

``` r

f <- tempfile(fileext = ".rds")
saveRDS(study, f)
identical(rng_normal(readRDS(f), 3), rng_normal(study, 3))
#> [1] TRUE
```
