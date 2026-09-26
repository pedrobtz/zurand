# Page size at the R level: the same build run three ways through the
# ZURAND_THP switch (off = no madvise, no = 4 KiB pages, huge = 2 MiB), one
# thread, randompack timed in the same call as the reference that does not
# ask for anything. On a runner with transparent_hugepage "always", `off`
# and `huge` coincide and `no` shows what a "madvise" system pays today.
args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench)})
  rng_threads(1L); key <- rng_key(7L); rp <- randompack::randompack_rng()
  for (n in c(1e6, 4e6, 1e7)) {
    it <- if (n < 1e7) 30 else 12
    u <- bench::mark(z = rng_uniform(key, n), r = rp$unif(len = n), check = FALSE, min_iterations = it, filter_gc = FALSE)
    g <- bench::mark(z = rng_normal(key, n), r = rp$normal(len = n), check = FALSE, min_iterations = it, filter_gc = FALSE)
    cat(sprintf("%s\t%g\t%.0f\t%.0f\t%.0f\t%.0f\n", args[2], n,
                n / as.numeric(u$median[1]) / 1e6, n / as.numeric(g$median[1]) / 1e6,
                n / as.numeric(u$median[2]) / 1e6, n / as.numeric(g$median[2]) / 1e6))
  }
  quit(save = "no")
}
res <- character()
for (r in 1:4) for (m in c("off", "no", "huge")) {
  Sys.setenv(ZURAND_THP = m)
  res <- c(res, system2("Rscript", c("dev/simd/ab_thp.R", "child", m), stdout = TRUE))
}
d <- read.delim(text = res, header = FALSE, col.names = c("mode", "n", "zurand_unif", "zurand_norm", "rp_unif", "rp_norm"))
f <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), min(x), max(x))
print(aggregate(cbind(zurand_unif, zurand_norm, rp_unif, rp_norm) ~ mode + n, d, f), row.names = FALSE)
