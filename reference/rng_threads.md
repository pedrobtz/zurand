# Control zurand's thread usage

`rng_threads()` reports the maximum number of OpenMP threads the
samplers may use; `rng_threads(threads)` caps it for the rest of the
session. The cap applies only to zurand's own parallel fills, not to
other OpenMP code in the process. Results never depend on the thread
count, so this is a performance control, not a reproducibility one. In
builds without OpenMP the value is always 1 and setting a cap has no
effect.

## Usage

``` r
rng_threads(threads = NULL)
```

## Arguments

- threads:

  `NULL` to query the effective maximum, a single whole number of at
  least 1 to cap zurand's thread use, or `0` to remove the cap. Values
  above the machine's thread count are allowed and equivalent to no cap.

## Value

Querying returns the effective maximum thread count as an integer.
Setting invisibly returns the previous *configured* cap, with `0L`
meaning "no cap was set", so `old <- rng_threads(1L)` followed by
`rng_threads(old)` restores the exact prior configuration – including
the uncapped state, which keeps tracking later changes to the
process-wide OpenMP maximum.

## Details

Draws below an internal size threshold always run single-threaded.

## Page size on Linux

R serves a large vector from fresh memory, so filling it is the first
touch of every page and pays one page fault per page. When the Linux
kernel's transparent huge pages are in `madvise` mode (the Ubuntu
default), the samplers ask for 2 MiB pages on output vectors of 4 MiB
and more before filling them, which about halves the time of large fills
on such systems; on kernels set to `always` the pages are large anyway,
and the request is not made. No value changes either way.
`options(zurand.hugepages = FALSE)` turns the request off.

## Examples

``` r
rng_threads()
#> [1] 4
old <- rng_threads(1L)
rng_threads(old)
```
