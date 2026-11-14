.libPaths("/w/linuxlib"); library(zurand); key <- rng_key(3L)
if (commandArgs(TRUE)[1] == "warm") p <- rng_normal(key, 32768)
r <- parallel::mclapply(1:2, function(i) rng_normal(key, 32768), mc.cores = 2)
cat("ok", identical(r[[1]], rng_normal(key, 32768)), "\n")
