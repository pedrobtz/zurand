# Same-machine A/B against randompack: `main` in lib-a, a candidate in lib-b,
# alternating processes; each times zurand and randompack's fastest engine
# in the same bench::mark() call and reports zurand's speed-up over it.
args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench)})
  rp <- randompack::randompack_rng()
  key <- rng_key(7L)
  for (n in c(1e5, 1e6, 1e7)) {
    it <- if (n < 1e7) 60 else 12
    u <- bench::mark(z = rng_uniform(key, n), r = rp$unif(len = n), check = FALSE,
                     min_iterations = it, filter_gc = FALSE)
    g <- bench::mark(z = rng_normal(key, n), r = rp$normal(len = n), check = FALSE,
                     min_iterations = it, filter_gc = FALSE)
    cat(sprintf("%s\t%g\t%.4f\t%.4f\t%.0f\t%.0f\n", args[2], n,
                as.numeric(u$median[2]) / as.numeric(u$median[1]),
                as.numeric(g$median[2]) / as.numeric(g$median[1]),
                n / as.numeric(u$median[1]) / 1e6, n / as.numeric(g$median[1]) / 1e6))
  }
  quit(save = "no")
}
res <- character()
for (r in 1:5) for (v in c("a", "b")) {
  Sys.setenv(R_LIBS = normalizePath(paste0("lib-", v)))
  res <- c(res, system2("Rscript", c("dev/simd/ab_randompack.R", "child", v), stdout = TRUE))
}
d <- read.delim(text = res, header = FALSE,
                col.names = c("v", "n", "unif_vs_rp", "norm_vs_rp", "unif_Mps", "norm_Mps"))
f <- function(x) sprintf("%.2f [%.2f-%.2f]", median(x), min(x), max(x))
out <- aggregate(cbind(unif_vs_rp, norm_vs_rp) ~ v + n, d, f)
out$unif_Mps <- aggregate(unif_Mps ~ v + n, d, median)$unif_Mps
out$norm_Mps <- aggregate(norm_Mps ~ v + n, d, median)$norm_Mps
print(out, row.names = FALSE)
