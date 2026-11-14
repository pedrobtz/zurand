source("/w/hunt-2/common.R"); library(parallel)
key <- rng_key(5L); v <- commandArgs(TRUE)[1]
x <- rng_lazy_normal(key, 1e5)
s <- sd(x)                       # materialises in the parent with a threaded fill; pointer went out writable
cat("parent materialised:", mat(x), "\n")
f <- switch(v,
  ret = function(i) x,                         # worker only returns the vector: mclapply serializes it
  read = function(i) sum(x),                   # worker only reads materialised data
  ser = function(i) length(serialize(x, NULL)))
r <- mclapply(1:2, f, mc.cores = 2); cat(v, "OK\n"); str(r)
