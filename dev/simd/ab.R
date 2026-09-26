# Same-machine A/B of two zurand builds. Install the baseline into lib-a and
# the candidate into lib-b in one directory, then run this from there. It
# alternates processes, each timing zurand against RcppZiggurat (normals)
# and dqrng (uniforms) in the same bench::mark() call. Runner-to-runner
# comparisons are not enough: GitHub hands out Intel and AMD CPUs at random.args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench); library(RcppZiggurat); library(dqrng)})
  for (eng in c("xoshiro256pp", "philox4x64")) {
    key <- rng_key(7L, engine = eng)
    for (n in c(1e6, 1e7)) {
      g <- bench::mark(z = rng_normal(key, n), r = zrnormMT(n), check = FALSE,
                       min_iterations = if (n < 1e7) 40 else 8, filter_gc = FALSE)
      u <- bench::mark(z = rng_uniform(key, n), r = dqrunif(n), check = FALSE,
                       min_iterations = if (n < 1e7) 40 else 8, filter_gc = FALSE)
      cat(sprintf("%s\t%s\t%g\t%.4f\t%.4f\n", args[2], eng, n,
                  as.numeric(g$median[2]) / as.numeric(g$median[1]),
                  as.numeric(u$median[2]) / as.numeric(u$median[1])))
    }
  }
  quit(save = "no")
}
res <- character()
for (r in 1:5) for (v in c("a", "b")) {
  Sys.setenv(R_LIBS = normalizePath(paste0("lib-", v)))
  res <- c(res, system2("Rscript", c("dev/simd/ab.R", "child", v), stdout = TRUE))
}
d <- read.delim(text = res, header = FALSE, col.names = c("v", "eng", "n", "normal", "unif"))
f <- function(x) sprintf("%.2f [%.2f-%.2f]", median(x), min(x), max(x))
print(aggregate(cbind(normal, unif) ~ v + eng + n, d, f), row.names = FALSE)
