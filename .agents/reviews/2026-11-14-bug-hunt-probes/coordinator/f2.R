.libPaths("HUNT/rerun/lib"); library(zurand); ks <- rng_key(42, n = 5)[0]
print(dim(rng_uniform(ks, 2^32 + 3))); print(dim(rng_uniform(ks, 2^40, offset = 5))); print(try(dim(rng_uniform(ks, 2^31)), silent=TRUE)); print(try(rng_uniform(rng_key(42,n=5), 2^32+3), silent=TRUE))
keys <- rng_key(1, n = 5); print(identical(keys[[2, 3]], keys[2]))
print(class(keys + 1L)); print(identical(rng_key(factor("5")), rng_key(1)))
print(identical(rng_key(753666613934570, n = 527), rng_key(0, n = 6000)[5474:6000]))
x <- rng_key(1); attr(x, "stream") <- 1.5; print(inherits(try(rng_uniform(x), silent=TRUE), "try-error")); print(length(rng_uniform(x[1])))
print(try(rng_simd(c(TRUE, FALSE)))); print(try(head(keys, 2), silent=TRUE))
