.libPaths("HUNT/rerun/lib"); library(zurand)
x <- rng_lazy_normal(rng_key(1L), 64); y <- structure(x, unit = "cm")
sum(y); saveRDS(y, tempfile()); print(try(sum(y)))
