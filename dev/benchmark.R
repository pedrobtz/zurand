# Benchmarks must run against the installed package (R CMD INSTALL, -O2);
# pkgload::load_all() builds understate rngat — see CLAUDE.md.
library(rngat)

x <- randompack::randompack_rng("philox")
y1 <- rng_key(123)
y10 <- rng_key(123, n = 10)

# NOTE: rng_*(y10, n = 1e6) draws 1e6 values per key = 1e7 values total, so
# multi-key comparisons need len = 1e7 on the randompack side.

res_unif <- bench::mark(
  check = F,
  rngat = rng_uniform(y1, n = 1e6),
  randompack = x$unif(len = 1e6),
)

res_unif_multikey <- bench::mark(
  check = F,
  rngat = rng_uniform(y10, n = 1e6),
  randompack = x$unif(len = 1e6 * 10),
)

res_norm <- bench::mark(
  check = F,
  rngat = rng_normal(y1, n = 1e6),
  randompack = x$normal(len = 1e6),
)

res_norm_multikey <- bench::mark(
  check = F,
  rngat = rng_normal(y10, n = 1e6),
  randompack = x$normal(len = 1e6 * 10),
)
