# Lazy vectors (rng_lazy_normal) against the materialised vector, per engine.
#
#   Rscript dev/lazy/bench_lazy.R [lib]
#
# From the repository root, against an -O2 install. n = 1e8 normals, all
# threads, xoshiro256pp and philox4x64 keys. "eager" is what a user would
# otherwise write: rng_normal() and then the operation; "plain op" is the
# operation alone on an existing vector.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressMessages({library(zurand); library(bench)})
cache <- zurand:::rng_lazy_cache
ms <- function(b) sprintf("%.3g", 1e3 * as.numeric(b$median))
n <- 1e8
set.seed(1)
idx <- sample(n, 1e5)
res <- list()
add <- function(eng, what, lazy, eager = NA, plain = NA, mem = NA)
  res[[length(res) + 1]] <<- data.frame(engine = eng, operation = what,
    lazy_ms = lazy, eager_ms = eager, plain_ms = plain, lazy_mem = mem)

for (eng in c("xoshiro256pp", "philox4x64")) {
  key <- rng_key(42L, engine = eng)
  m <- rng_normal(key, n)
  lz <- rng_lazy_normal(key, n)
  stopifnot(identical(sum(lz), sum(m)))

  b <- bench::mark(rng_lazy_normal(key, n), min_iterations = 50)
  add(eng, "create", ms(b), ms(bench::mark(rng_normal(key, n), min_iterations = 3)))

  for (op in c("sum", "mean", "head1e6", "times2")) {
    f <- switch(op, sum = sum, mean = mean, head1e6 = function(x) x[1:1e6],
                times2 = function(x) x * 2)
    bl <- bench::mark(f(rng_lazy_normal(key, n)), min_iterations = 3,
                      filter_gc = FALSE)
    be <- bench::mark(f(rng_normal(key, n)), min_iterations = 3, filter_gc = FALSE)
    bp <- bench::mark(f(m), min_iterations = 3, filter_gc = FALSE)
    add(eng, op, ms(bl), ms(be), ms(bp), format(bl$mem_alloc))
  }

  for (mode in 0:2) {
    cache(mode)
    tag <- c("value alone", "adaptive", "chunk always")[mode + 1]
    x <- rng_lazy_normal(key, n)
    b <- bench::mark(x[idx], min_iterations = 5)
    add(eng, paste("x[sample(n, 1e5)],", tag), ms(b), NA, ms(bench::mark(m[idx])))
    loop <- function(x, ii) { s <- 0; for (i in ii) s <- s + x[i]; s }
    b <- bench::mark(loop(x, idx), min_iterations = 3)
    add(eng, paste("for over 1e5 random i,", tag), ms(b), NA,
        ms(bench::mark(loop(m, idx), min_iterations = 3)))
    b <- bench::mark(loop(x, 1:1e6), min_iterations = 3)
    add(eng, paste("for (i in 1:1e6),", tag), ms(b), NA,
        ms(bench::mark(loop(m, 1:1e6), min_iterations = 3)))
    stopifnot(!zurand:::rng_lazy_materialised(x))
  }
  cache(1L)
  add(eng, "serialize bytes", length(serialize(lz, NULL)), length(serialize(m, NULL)))
  rm(m, lz); invisible(gc())
}

# sum at 1e9 in fresh processes, peak resident memory from /usr/bin/time:
# loading zurand alone, then summing a lazy 1e9 vector.
peak <- function(expr) {
  lib <- if (length(args)) sprintf(".libPaths(c('%s', .libPaths())); ", args[1]) else ""
  code <- sprintf("%ssuppressMessages(library(zurand)); t <- system.time(v <- %s); cat(t[['elapsed']], v, '\\n')",
                  lib, expr)
  flag <- if (Sys.info()[["sysname"]] == "Darwin") "-l" else "-v"
  out <- system2("/usr/bin/time", c(flag, file.path(R.home("bin"), "Rscript"),
                                    "-e", shQuote(code)), stdout = TRUE, stderr = TRUE)
  kb <- if (flag == "-l") as.numeric(sub("^ *([0-9]+) +maximum resident.*", "\\1",
                                         grep("maximum resident", out, value = TRUE))) / 1024
        else as.numeric(sub(".*: ", "", grep("Maximum resident", out, value = TRUE)))
  list(mb = kb / 1024, out = out[1])
}
p0 <- peak("0")
p1 <- peak("sum(rng_lazy_normal(rng_key(42L), 1e9))")
cat(sprintf("sum(rng_lazy_normal(key, 1e9)): %s s; peak RSS %.0f MB vs %.0f MB for R + zurand alone (a materialised 1e9 vector is 7629 MB)\n",
            strsplit(p1$out, " ")[[1]][1], p1$mb, p0$mb))

d <- do.call(rbind, res)
print(d, row.names = FALSE)
