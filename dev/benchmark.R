pkgload::load_all(compile = T, debug = F)

x <- randompack::randompack_rng("philox")
y <- rngat::rng_key(seed = 123, n = 10)


res <- bench::mark(
  check = F,
  rngat = rngat::rng_uniform(y, n = 1e6),
  randompack = x$unif(len = 1e6 * 10),
)

res1 <- bench::mark(
  check = F,
  rngat = rngat::rng_uniform(y, n = 1e6),
  randompack = x$unif(len = 1e6),
)

res2 <- bench::mark(
  check = F,
  rngat = rngat::rng_normal(y, n = 1e6),
  randompack = x$normal(len = 1e6),
)
