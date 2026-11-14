source("/w/hunt-2/common.R"); library(parallel)
key <- rng_key(5L); v <- commandArgs(TRUE)[1]
p <- rng_normal(key, 32768)   # parent enters a parallel region (smallest threaded size)
lz <- rng_lazy_normal(key, 1e6)
f <- switch(v,
  eager_small = function(i) length(rng_normal(key, 32767)),
  eager_min = function(i) length(rng_normal(key, 32768)),
  lazy_sum = function(i) sum(lz),
  lazy_region = function(i) length(lz[1:32768]),
  lazy_mat = function(i) length(lz + 0),
  lazy_save_after_ptr = function(i) { z <- rng_lazy_normal(key, 40000); length(z) },
  cap1 = function(i) { rng_threads(1L); length(rng_normal(key, 1e6)) },
  integer = function(i) length(rng_integer(key, 1e5, 1L, 10L)))
r <- mclapply(1:2, f, mc.cores = 2); cat(v, "OK", unlist(r), "\n")
