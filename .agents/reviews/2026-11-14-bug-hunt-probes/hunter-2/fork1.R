source("/w/hunt-2/common.R"); library(parallel)
cat("threads:", rng_threads(), "\n")
key <- rng_key(5L)
arg <- commandArgs(TRUE)
e1 <- rng_threads(1L); ref <- rng_normal(key, 1e6); rng_threads(e1)   # serial reference (no parallel region: nt=1)
if ("warm" %in% arg) { p <- rng_normal(key, 1e6); chk("parent threaded == serial", identical(p, ref)) }
t0 <- Sys.time()
r <- mclapply(1:4, function(i) { x <- rng_normal(key, 1e6); l <- rng_lazy_normal(key, 1e6); c(identical(x, ref), identical(sum(l), sum(ref)), identical(l + 0, ref)) }, mc.cores = 2)
print(Sys.time() - t0); str(r)
chk("children ok", all(unlist(r)))
# mcparallel
j <- mcparallel({ x <- rng_lazy_uniform(key, 1e6); invisible(cumsum(x)); list(identical(x[], rng_uniform(key, 1e6)), serialize(x, NULL)) })
res <- mccollect(j, wait = TRUE)[[1]]
chk("mcparallel", res[[1]] && length(res[[2]]) < 2000 && identical(unserialize(res[[2]])[], rng_uniform(key, 1e6)))
cat("failures:", nfail, "\n")
