# Extracted from test-c-api.R:155

# test -------------------------------------------------------------------------
skip_on_cran()
ns <- client()
n <- 3000L
for (engine in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    keys <- rng_key(11L, n = 2L, engine = engine)
    for (base in c(2^60, 2^63 + 2^40, 2^64 - 2^16)) {
      for (step in c(256, 4096)) {
        info <- paste(engine, format(base, digits = 20), step)
        tail <- step + seq_len(n)
        for (method in c("ziggurat", "mcfarland")) {
          long <- ns$client_fill_at(keys, n + step, "normal", 1, 2, base, method)
          expect_identical(
            ns$client_fill_at(keys, n, "normal", 1, 2, base + step, method),
            long[tail, , drop = FALSE], info = paste(info, method)
          )
        }
        long <- ns$client_fill_at(keys, n + step, "uniform", -1, 3, base)
        expect_identical(ns$client_fill_at(keys, n, "uniform", -1, 3, base + step),
                         long[tail, , drop = FALSE], info = info)
      }
    }
  }
