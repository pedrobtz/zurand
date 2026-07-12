pkgload::load_all(compile = T, debug = F)

x <- randompack::randompack_rng("philox")
y <- rngat::rng_key(123)


res <- bench::mark(
  check = F,
  rngat = rngat::runif_seq(y, start = 1, len = 1e6),
  randompack = x$unif(len = 1e6),
)
