args <- commandArgs(TRUE)
.libPaths(args[1]); library(zurand)
bad <- 0L
chk <- function(ok, ...) { if (!isTRUE(ok)) { bad <<- bad + 1L; cat("FAIL:", ..., "\n") } }
err <- function(expr) tryCatch({ force(expr); "OK" }, error = function(e) paste("ERR:", conditionMessage(e)))
M <- .Machine$double.xmax; tiny <- 2^-1074; eps <- 2^-52
nx <- function(x, up = TRUE) {
  if (x == 0) return(if (up) tiny else -tiny)
  e <- floor(log2(abs(x))); step <- max(2^(e - 52), tiny)
  away <- (x > 0) == up
  if (!away && abs(x) == 2^e && step / 2 >= tiny) step <- step / 2
  x + (if (up) step else -step)
}
stopifnot(nx(1) == 1 + eps, nx(1, FALSE) == 1 - eps / 2, nx(-1) == -1 + eps / 2)
bounds <- list(c(-M, M), c(-M, 2^970), c(-2^970, M), c(-M, 0), c(0, M), c(-M, 2^969), c(-M / 2, M), c(-M, M / 2), c(-2^1023, 2^1023), c(-2^1023, M),
  c(-M, nx(-M)), c(nx(M, FALSE), M), c(1, 1 + eps), c(0, tiny), c(-tiny, tiny), c(tiny, 2 * tiny), c(-0, 0), c(0, -0), c(-tiny, 0), c(2^-1022, nx(2^-1022)),
  c(-1, nx(-1)), c(1e308, nx(1e308)), c(-1, 1.5 * 2^-53), c(-1, 3), c(1 - eps / 2, 1), c(3, 3), c(-M, -M), c(M, M), c(-2^1023, nx(2^1023, FALSE)), c(-M, nx(2^970)), c(-M, nx(2^970, FALSE)),
  c(-5e-324, 1e-320), c(2^52, 2^52 + 1), c(2^53, 2^53 + 2), c(-1e-310, 1e-310), c(-0.1, 0.7), c(1e300, 1.5e300))
for (e in c("xoshiro256pp", "philox4x64")) { k <- rng_key(1, engine = e); n <- 200000; u <- rng_uniform(k, n)
  chk(all(u > 0 & u < 1) && all((u * 2^53) %% 2 == 1), "default u form")
  for (b in bounds) {
    x <- rng_uniform(k, n, b[1], b[2])
    span <- b[2] - b[1]
    want <- if (span == 0) rep(b[1], n) else if (is.finite(span)) b[1] + span * u else 2 * (b[1] / 2 + (b[2] / 2 - b[1] / 2) * u)
    chk(identical(x, want), "uniform formula", e, sprintf("%a %a", b[1], b[2]))
    chk(all(is.finite(x)) && all(x >= b[1] & x <= b[2]), "uniform range", e, sprintf("%a %a", b[1], b[2]), range(x))
    if (span != 0) cat(sprintf("%s [%a, %a]: hits min %d, hits max %d\n", e, b[1], b[2], sum(x == b[1]), sum(x == b[2])))[e == "philox4x64"]
  }
  z <- rng_normal(k, n); zm <- rng_normal(k, n, method = "mcfarland")
  for (p in list(c(0, 1), c(1, 2), c(-0, 1), c(M, 1), c(M, M), c(-M, M), c(0, M), c(0, tiny), c(tiny, tiny), c(1, eps), c(3, 0), c(-0, 0), c(0, -0), c(1e308, 1e292), c(0, 2^-1022), c(-0.0, 2))) {
    x <- rng_normal(k, n, p[1], p[2]); xm <- rng_normal(k, n, p[1], p[2], "mcfarland")
    want <- if (p[2] == 0) rep(p[1], n) else if (p[1] == 0 && p[2] == 1) z else p[1] + p[2] * z
    wantm <- if (p[2] == 0) rep(p[1], n) else if (p[1] == 0 && p[2] == 1) zm else p[1] + p[2] * zm
    chk(identical(x, want) && identical(xm, wantm), "normal formula", e, sprintf("%a %a", p[1], p[2]))
    if (!all(is.finite(x))) cat(sprintf("%s normal mean=%a sd=%a: %d non-finite (Inf %d, NaN %d)\n", e, p[1], p[2], sum(!is.finite(x)), sum(is.infinite(x)), sum(is.nan(x))))
  }
  cat("signed zeros: mean=-0, sd=2 any -0?", any(1 / rng_normal(k, n, -0.0, 2) == -Inf), " sd=0 mean=-0 gives", 1 / rng_normal(k, 1, -0.0, 0), " unif(-0,0):", 1 / rng_uniform(k, 1, -0.0, 0), " unif(0,-0):", 1 / rng_uniform(k, 1, 0, -0.0), "\n")
}
k <- rng_key(1)
calls <- list(
  quote(rng_uniform(k, NaN)), quote(rng_uniform(k, Inf)), quote(rng_uniform(k, -0.0)), quote(rng_uniform(k, -1)), quote(rng_uniform(k, 0.5)), quote(rng_uniform(k, 2^53)), quote(rng_uniform(k, 2^52)),
  quote(rng_uniform(k, NA_integer_)), quote(rng_uniform(k, 1:2)), quote(rng_uniform(k, "1")), quote(rng_uniform(k, TRUE)), quote(rng_uniform(k, 1, NaN, 1)), quote(rng_uniform(k, 1, 0, Inf)), quote(rng_uniform(k, 1, -Inf, 0)), quote(rng_uniform(k, 1, 1, 0)), quote(rng_uniform(k, 1, NA, 1)), quote(rng_uniform(k, 1, NA_real_, 1)),
  quote(rng_uniform(k, 1, offset = 2^53)), quote(rng_uniform(k, 1, offset = 2^53 - 1)), quote(rng_uniform(k, 0, offset = 2^53)), quote(rng_uniform(k, 2, offset = 2^53 - 1)), quote(rng_uniform(k, 1, offset = -1)), quote(rng_uniform(k, 1, offset = -0.0)), quote(rng_uniform(k, 1, offset = NaN)), quote(rng_uniform(k, 1, offset = 0.5)), quote(rng_uniform(k, 1, offset = Inf)), quote(rng_uniform(k, 1, offset = 2^60)),
  quote(rng_normal(k, 1, 0, -1)), quote(rng_normal(k, 1, 0, -0.0)), quote(rng_normal(k, 1, 0, NaN)), quote(rng_normal(k, 1, Inf, 1)), quote(rng_normal(k, 1, 0, Inf)), quote(rng_normal(k, 1, 0, -tiny)), quote(rng_normal(k, 1, method = "zig")), quote(rng_normal(k, 1, method = c("mcfarland", "ziggurat"))), quote(rng_normal(k, 1, offset = 2^53 - 1)), quote(rng_normal(k, 1, offset = 2^53)),
  quote(rng_integer(k, 1, -2^31, 0)), quote(rng_integer(k, 1, -2^31 + 1, 2^31 - 1)), quote(rng_integer(k, 1, 0, 2^31)), quote(rng_integer(k, 1, 0.5, 2)), quote(rng_integer(k, 1, 2, 1)), quote(rng_integer(k, 1, NA_integer_, 1L)), quote(rng_integer(k, 1, -0.0, 0)), quote(rng_integer(k, 1, NaN, 1)), quote(rng_integer(k, 1, 1, Inf)), quote(rng_integer(k, 1, 1)), quote(rng_integer(k, 3, 2^31 - 1, 2^31 - 1)),
  quote(rng_bits(k, 1, 32.0)), quote(rng_bits(k, 1, 64)), quote(rng_bits(k, 1, 33)), quote(rng_bits(k, 1, NaN)), quote(rng_bits(k, 1, 16)),
  quote(rng_key(2^53)), quote(rng_key(2^53 - 1)), quote(rng_key(-0.0)), quote(rng_key(-1)), quote(rng_key(0.5)), quote(rng_key(NaN)), quote(rng_key(1, 0)), quote(rng_key(1, -1)), quote(rng_key(1, 2^31)), quote(rng_key(1, engine = "x")),
  quote(rng_threads(0)), quote(rng_threads(-1)), quote(rng_threads(2^40)), quote(rng_threads(NaN)), quote(rng_threads(1.5)), quote(rng_simd(NA)), quote(rng_simd("a")), quote(rng_simd(c(TRUE, FALSE))), quote(rng_simd(0)))
for (cl in calls) cat(sprintf("%-55s %s\n", deparse(cl), err(eval(cl))))
rng_threads(0); rng_simd(TRUE)
chk(identical(rng_uniform(k, 3, offset = -0.0), rng_uniform(k, 3)), "offset -0")
chk(identical(rng_key(-0.0), rng_key(0)), "seed -0"); chk(identical(rng_key(2^53 - 1), rng_key(2^53 - 1)), "seed max")
chk(identical(rng_uniform(k, 5, offset = 2^53 - 5), rng_uniform(k, 5, offset = 2^53 - 5)), "offset max")
chk(identical(rng_integer(k, 1000, 5, 5), rep(5L, 1000)), "range 1"); chk(identical(rng_integer(k, 10, -0.0, 0), rep(0L, 10)), "range -0")
x <- rng_integer(k, 1e6, -2147483647, 2147483647); chk(!anyNA(x), "full range no NA"); cat("full range: range", range(x), "\n")
x <- rng_integer(k, 1e6, 2147483646, 2147483647); chk(all(x %in% c(2147483646L, 2147483647L)) && abs(mean(x == 2147483647L) - .5) < .003, "top range")
cat("numbers BAD =", bad, "\n")
