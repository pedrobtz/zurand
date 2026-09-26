# Normal-only same-machine A/B, one thread: lib-a vs lib-b, ten alternating
# rounds, three engines, two sizes; reports b/a with its range.
args <- commandArgs(TRUE)
if (length(args) && args[1] == "child") {
  suppressMessages({library(zurand); library(bench)})
  for (eng in c("xoshiro256pp", "philox4x64", "threefry4x64")) {
    key <- rng_key(7L, engine = eng); rng_threads(1L)
    for (n in c(1e6, 1e7)) {
      b <- bench::mark(rng_normal(key, n), min_iterations = if (n < 1e7) 200 else 30, filter_gc = FALSE)
      cat(sprintf("%s\t%s\t%g\t%.1f\n", args[2], eng, n, n / as.numeric(b$median) / 1e6))
    }
  }
  quit(save = "no")
}
res <- character()
for (r in 1:10) for (v in c("a", "b")) {
  Sys.setenv(R_LIBS = normalizePath(paste0("lib-", v)))
  res <- c(res, system2("Rscript", c("dev/simd/norm_ab.R", "child", v), stdout = TRUE))
}
d <- read.delim(text = res, header = FALSE, col.names = c("v", "eng", "n", "Mps"))
w <- reshape(transform(d, r = ave(Mps, v, eng, n, FUN = seq_along)), idvar = c("eng", "n", "r"), timevar = "v", direction = "wide")
w$ratio <- w$Mps.b / w$Mps.a
print(do.call(rbind, lapply(split(w, list(w$eng, w$n)), function(x)
  data.frame(eng = x$eng[1], n = x$n[1], a = median(x$Mps.a), b = median(x$Mps.b),
             ratio = median(x$ratio), lo = min(x$ratio), hi = max(x$ratio),
             b_wins = sum(x$ratio > 1)))), row.names = FALSE)
