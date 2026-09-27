# The normal methods side by side, one build, one thread: zurand's default
# ziggurat, rng_normal(method = "mcfarland") and randompack's fastest engine
# in the same bench::mark() call, alternating processes. Reports M/s and
# McFarland's speed-up over the default and over randompack.
args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench)})
  rng_threads(1L)
  rp <- if (requireNamespace("randompack", quietly = TRUE)) randompack::randompack_rng()
  for (eng in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    key <- rng_key(7L, engine = eng)
    for (n in c(1e5, 1e6, 1e7)) {
      it <- if (n < 1e7) 60 else 15
      exprs <- list(z = quote(rng_normal(key, n)),
                    m = quote(rng_normal(key, n, method = "mcfarland")))
      if (!is.null(rp)) exprs$r <- quote(rp$normal(len = n))
      b <- bench::mark(exprs = exprs, check = FALSE, min_iterations = it, filter_gc = FALSE)
      t <- as.numeric(b$median)
      cat(sprintf("%s\t%g\t%.0f\t%.0f\t%.0f\n", eng, n, n / t[1] / 1e6, n / t[2] / 1e6,
                  if (length(t) > 2) n / t[3] / 1e6 else NA))
    }
  }
  quit(save = "no")
}
res <- character()
for (r in 1:5) res <- c(res, system2("Rscript", c("dev/simd/ab_methods.R", "child"), stdout = TRUE))
d <- read.delim(text = res, header = FALSE, col.names = c("eng", "n", "ziggurat", "mcfarland", "randompack"))
d$vs_zig <- d$mcfarland / d$ziggurat; d$vs_rp <- d$mcfarland / d$randompack
f <- function(x) sprintf("%.2f [%.2f-%.2f]", median(x), min(x), max(x))
g <- function(x) sprintf("%.0f", median(x))
out <- aggregate(cbind(ziggurat, mcfarland, randompack) ~ eng + n, d, g)
out2 <- aggregate(cbind(vs_zig, vs_rp) ~ eng + n, d, f)
print(merge(out, out2), row.names = FALSE)
