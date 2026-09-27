# Same-machine A/B of the existing fills, for a change that must not move
# them (feat/stream factored each fill's loop body into a forced-inline
# per-chunk function). Install the baseline into lib-a and the candidate
# into lib-b in one directory and run this from there. Processes alternate;
# each times zurand against dqrng (uniform) or RcppZiggurat (normal) in the
# same bench::mark() call, so the reported ratio b/a is of zurand relative
# to that in-process reference, which cancels drift between processes.
args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench); library(dqrng)
                    library(RcppZiggurat)})
  for (th in c(1L, 0L)) {
    rng_threads(th)
    for (eng in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
      key <- rng_key(7L, engine = eng)
      for (n in c(1e5, 1e6, 1e7)) {
        it <- if (n < 1e7) 60 else 12
        u <- bench::mark(z = rng_uniform(key, n), r = dqrunif(n),
                         check = FALSE, min_iterations = it, filter_gc = FALSE)
        g <- bench::mark(z = rng_normal(key, n), r = zrnormMT(n),
                         check = FALSE, min_iterations = it, filter_gc = FALSE)
        f <- bench::mark(z = rng_normal(key, n, method = "mcfarland"),
                         r = zrnormMT(n), check = FALSE, min_iterations = it,
                         filter_gc = FALSE)
        rel <- function(b) as.numeric(b$median[2]) / as.numeric(b$median[1])
        cat(sprintf("%s\t%d\t%s\t%g\t%.4f\t%.4f\t%.4f\n", args[2],
                    if (th == 0L) 8L else 1L, eng, n, rel(u), rel(g), rel(f)))
      }
    }
  }
  quit(save = "no")
}
res <- character()
for (r in 1:5) for (v in c("a", "b")) {
  Sys.setenv(R_LIBS = normalizePath(paste0("lib-", v)))
  res <- c(res, system2("Rscript", c(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE))), "child", v),
    stdout = TRUE))
}
d <- read.delim(text = res, header = FALSE,
                col.names = c("v", "threads", "eng", "n", "uniform", "normal",
                              "mcfarland"))
key <- c("threads", "eng", "n")
a <- aggregate(cbind(uniform, normal, mcfarland) ~ threads + eng + n,
               d[d$v == "a", ], median)
b <- aggregate(cbind(uniform, normal, mcfarland) ~ threads + eng + n,
               d[d$v == "b", ], median)
out <- merge(a, b, by = key, suffixes = c(".a", ".b"))
for (s in c("uniform", "normal", "mcfarland"))
  out[[paste0(s, " b/a")]] <- round(out[[paste0(s, ".b")]] / out[[paste0(s, ".a")]], 3)
print(out[order(out$threads, out$eng, out$n),
          c(key, "uniform b/a", "normal b/a", "mcfarland b/a")], row.names = FALSE)
