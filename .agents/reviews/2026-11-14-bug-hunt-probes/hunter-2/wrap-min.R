source(if (file.exists("/w/hunt-2/common.R")) "/w/hunt-2/common.R" else "common.R")
cat(R.version.string, "\n")
for (n in c(63, 64)) { x <- rng_lazy_normal(rng_key(1L), n)
  y <- structure(x, unit = "cm")
  r <- tryCatch({ s1 <- sum(y); saveRDS(y, tf <- tempfile()); s2 <- sum(y); paste("ok", identical(s1, s2)) }, error = function(e) conditionMessage(e))
  cat("n =", n, ": sum; saveRDS; sum ->", r, "\n")
  r <- tryCatch({ m <- mean(y); "ok" }, error = function(e) conditionMessage(e)); cat("   later mean(y) ->", r, "\n")
  cat("   file restores right:", identical(readRDS(tf), structure(rng_normal(rng_key(1L), n), unit = "cm")), "\n") }
