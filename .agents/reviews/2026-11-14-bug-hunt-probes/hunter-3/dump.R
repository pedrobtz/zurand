args <- commandArgs(TRUE); lib <- args[1]; out <- args[2]
.libPaths(lib); library(zurand)
N <- 300000L
ranges <- list(c(1,6), c(-1, 2147483647), c(-2147483647, 2147483647), c(5,5), c(-1610612736, 1610612736), c(0, 2^31-2))
for (e in c("philox4x64","threefry4x64","xoshiro256pp")) {
  k <- rng_key(42, 3, engine = e)[2]
  f <- function(nm) file.path(out, paste0(e, "_", nm))
  writeLines(rng_bits(k, 50000, 64), f("bits"))
  writeBin(rng_uniform(k, N), f("unif"))
  writeBin(rng_normal(k, N), f("norm"))
  writeBin(rng_normal(k, N, method = "mcfarland"), f("mcf"))
  for (i in seq_along(ranges)) writeBin(rng_integer(k, N, ranges[[i]][1], ranges[[i]][2]), f(paste0("int", i)))
  offs <- c(511, 4095, 2^32 - 7, 2^53 - 20000)
  for (i in seq_along(offs)) {
    writeBin(rng_uniform(k, 20000, offset = offs[i]), f(paste0("unifoff", i)))
    writeBin(rng_normal(k, 20000, offset = offs[i]), f(paste0("normoff", i)))
    writeBin(rng_normal(k, 20000, offset = offs[i], method = "mcfarland"), f(paste0("mcfoff", i)))
  }
}
cat(format(rng_key(42, 3)), "\n")
