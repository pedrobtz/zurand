source("helper.R")
a <- rng_key(0, n = 6000); b <- rng_key(753666613934570, n = 527); b2 <- rng_key(7138832632025594, n=100)
cat("rng_key(753666613934570, n=527) == rng_key(0, n=6000)[5474:6000]:", identical(b, a[5474:6000]), "\n")
cat("rng_key(7138832632025594, n=100) == rng_key(0, n=6000)[2890:2989]:", identical(b2, a[2890:2989]), "\n")
cat("dups in c(a,b):", sum(duplicated(c(a, b))), "\n")
cat("same draws:", identical(rng_normal(b[1], 5), rng_normal(a[5474], 5)), "\n")
# small seeds: adjacent seeds, any shared keys among first 1e5?
x <- c(rng_key(1, n=1e5), rng_key(2, n=1e5)); cat("dups seeds 1,2:", sum(duplicated(x)), "\n")
