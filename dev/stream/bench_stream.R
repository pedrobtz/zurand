# Does streaming keep the values out of RAM, and does it pay?
#
#   Rscript dev/stream/bench_stream.R [lib]
#
# Run from the repository root against an -O2 install of zurand (R CMD
# INSTALL, or pass the library it went into). Two parts:
#
# 1. C level (dev/stream/bench_stream.c): K keys x `total` uniforms or
#    normals, summed as squares, one key per thread, at 1, 2, 4 and 8 threads. A whole-stream
#    buffer (reused, or fresh per key as R allocates) against stream_normal()
#    through a cache-sized chunk. Identical sums are checked.
# 2. R level: a loop over chunks of one big vector against the same loop
#    pulling each chunk with rng_normal(offset =). Same FUN, same values,
#    identical result; time and memory allocated.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressMessages({library(zurand); library(bench)})

# ---- 1. C level ----
src <- normalizePath("dev/stream/bench_stream.c")
inc <- normalizePath("inst/include")
dir <- tempfile("bench_stream"); dir.create(dir)
invisible(file.copy(src, dir))
omp <- if (Sys.info()[["sysname"]] == "Darwin") {
  c("-Xclang -fopenmp", "-lomp")
} else {
  c("-fopenmp", "-fopenmp")
}
# No contraction: clang otherwise fuses the consumer's x * x + acc where it
# inlines it and not where it is called through a pointer, and the sums
# then differ in the last bit although the values are the same.
Sys.setenv(PKG_CPPFLAGS = sprintf('-I"%s"', inc),
           PKG_CFLAGS = paste(omp[1], "-ffp-contract=off"), PKG_LIBS = omp[2])
old <- setwd(dir)
stopifnot(system2(file.path(R.home("bin"), "R"),
                  c("CMD", "SHLIB", "-o", "bench_stream.so", "bench_stream.c"),
                  stdout = FALSE) == 0)
setwd(old)
dyn.load(file.path(dir, "bench_stream.so"))
run <- function(keys, total, nt, mode, dist, chunk = 16384L)
  .Call("bs_run", keys, as.double(total), as.integer(nt), as.integer(mode),
        as.integer(chunk), as.integer(dist))

K <- 16L; total <- 1e7
keys <- rng_key(20260927L, n = K)
for (dist in 0:1) {
  ref <- run(keys, total, 1L, 0L, dist)
  stopifnot(identical(run(keys, total, 1L, 2L, dist), ref),
            identical(run(keys, total, 8L, 2L, dist), ref),
            identical(run(keys, total, 8L, 1L, dist), ref))
}
cat(sprintf("C level: %d keys x %g values, sum of squares; identical sums: yes\n",
            K, total))
rows <- list()
for (r in 1:3) for (dist in 0:1) for (nt in c(1L, 2L, 4L, 8L)) for (mode in 0:2) {
  t <- as.numeric(bench::mark(run(keys, total, nt, mode, dist),
                              min_iterations = 3, check = FALSE,
                              filter_gc = FALSE)$median)
  rows[[length(rows) + 1]] <- data.frame(round = r,
    dist = c("uniform", "normal")[dist + 1], threads = nt,
    mode = c("big buffer", "fresh buffer", "stream 16k")[mode + 1],
    Mps = K * total / t / 1e6)
}
d <- do.call(rbind, rows)
a <- aggregate(Mps ~ dist + threads + mode, d, median)
w <- reshape(a, idvar = c("dist", "threads"), timevar = "mode", direction = "wide")
names(w) <- sub("Mps.", "", names(w))
w$`stream vs big` <- round(w$`stream 16k` / w$`big buffer`, 2)
w$`stream vs fresh` <- round(w$`stream 16k` / w$`fresh buffer`, 2)
for (m in c("big buffer", "fresh buffer", "stream 16k")) w[[m]] <- round(w[[m]])
print(w[order(-xtfrm(w$dist), w$threads), ], row.names = FALSE)

# ---- 2. R level ----
cat("\nR level: the same FUN over 16384-value chunks, one key\n")
key <- rng_key(42L)
chunk <- 16384L
big_loop <- function(N, FUN) {
  x <- rng_normal(key, N)
  acc <- 0
  for (s in seq(0, N - 1, by = chunk))
    acc <- acc + FUN(x[(s + 1):min(s + chunk, N)])
  acc
}
pull_loop <- function(N, FUN) {
  acc <- 0
  for (s in seq(0, N - 1, by = chunk))
    acc <- acc + FUN(rng_normal(key, min(chunk, N - s), offset = s))
  acc
}
funs <- list(`sum(x)` = function(x) sum(x), `sum(x^2)` = function(x) sum(x^2))
for (th in c(1L, 0L)) {
  rng_threads(th)
  for (N in c(1e7, 1e8)) for (fn in names(funs)) {
    stopifnot(identical(big_loop(1e5, funs[[fn]]), pull_loop(1e5, funs[[fn]])))
    b <- bench::mark(big = big_loop(N, funs[[fn]]), pull = pull_loop(N, funs[[fn]]),
                     min_iterations = if (N < 1e8) 5 else 2, check = TRUE,
                     filter_gc = FALSE, memory = TRUE)
    cat(sprintf("threads %s  N %g  %-9s  big %6.0f ms %8s   pull %6.0f ms %8s   pull/big time %.2f\n",
                if (th == 0L) "all" else "1", N, fn,
                1e3 * as.numeric(b$median[1]), format(b$mem_alloc[1]),
                1e3 * as.numeric(b$median[2]), format(b$mem_alloc[2]),
                as.numeric(b$median[2]) / as.numeric(b$median[1])))
  }
}
rng_threads(0L)
