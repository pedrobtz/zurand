args <- commandArgs(TRUE)
.libPaths(c(args[2], args[1])); library(zurand); library(zc3)
outdir <- args[3]
bad <- 0L
chk <- function(ok, ...) { if (!isTRUE(ok)) { bad <<- bad + 1L; cat("FAIL:", ..., "\n") } }
engs <- c("philox4x64", "threefry4x64", "xoshiro256pp")
K <- function(e, seed = 7) { k <- rng_key(seed, engine = e); z <- zc_key_get(k, 0); attr(z, "rc") <- NULL; list(r = k, c = z) }
# 1 key_get
for (e in engs) {
  k <- rng_key(3, 5, engine = e)
  for (i in 0:4) { z <- zc_key_get(k, i); f <- format(k[i + 1])
    # format prints words k0lo k0hi k1lo k1hi
    h <- sub("rng_key\\[(.*)\\]", "\\1", f); w <- substring(h, c(1, 9, 17, 25), c(8, 16, 24, 32))
    chk(identical(z$k0, paste0("0x", w[2], w[1])) && identical(z$k1, paste0("0x", w[4], w[3])) && z$stream == 1L && z$engine == match(e, engs) - 1L, "key_get", e, i) }
  chk(attr(zc_key_get(k, 5), "rc") == 1L, "key_get oob hi"); chk(attr(zc_key_get(k, -1), "rc") == 1L, "key_get oob lo")
  chk(attr(zc_key_get(k, 2^40), "rc") == 1L, "key_get oob huge")
}
ns <- c(0, 1, 3, 4, 5, 511, 512, 513, 2047, 2048, 2049, 4095, 4096, 4097, 5119, 5120, 5121, 8191, 8192, 8193, 10239, 10240, 10241, 20481)
for (e in engs) { k <- K(e)
  for (n in ns) {
    r <- zc_fill(k$c, n, 0); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_uniform(k$r, n)), "fill_uniform", e, n)
    r <- zc_fill(k$c, n, 0, -3, 5.5); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_uniform(k$r, n, -3, 5.5)), "fill_uniform scaled", e, n)
    r <- zc_fill(k$c, n, 0, -.Machine$double.xmax, .Machine$double.xmax); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_uniform(k$r, n, -.Machine$double.xmax, .Machine$double.xmax)) && all(is.finite(r$x)), "fill_uniform twice", e, n)
    r <- zc_fill(k$c, n, 1); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n)), "fill_normal", e, n)
    r <- zc_fill(k$c, n, 1, 2, 3); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n, 2, 3)), "fill_normal scaled", e, n)
    r <- zc_fill(k$c, n, 4, 0, 1, 1L); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n, method = "mcfarland")), "fill_normal mcf", e, n)
    r <- zc_fill(k$c, n, 4, 0, 1, 0L); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n)), "fill_normal_method zig", e, n)
    r <- zc_fill(k$c, n, 2, imin = -1L, imax = 2147483647L); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_integer(k$r, n, -1L, 2147483647L)), "fill_integer", e, n)
    r <- zc_fill(k$c, n, 2, imin = -2147483647L, imax = 2147483647L); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_integer(k$r, n, -2147483647L, 2147483647L)), "fill_integer full", e, n)
    r <- zc_fill(k$c, n, 3); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_bits(k$r, n, 64)), "fill_bits64", e, n)
  }
  # positional
  for (st in c(0, 1, 2, 3, 5, 511, 512, 513, 4095, 4096, 4097, 5119, 5120, 5121, 2^32 - 3, 2^53 - 6000))
    for (n in c(0, 1, 2, 7, 511, 512, 513, 4096, 5120, 5121, 6000)) {
      sh <- sprintf("%.0f", st)
      r <- zc_fill_at(k$c, sh, n, 0); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_uniform(k$r, n, offset = st)), "fill_uniform_at", e, st, n)
      r <- zc_fill_at(k$c, sh, n, 0, 1, 1 + 2^-52); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_uniform(k$r, n, 1, 1 + 2^-52, offset = st)), "fill_uniform_at adj", e, st, n)
      r <- zc_fill_at(k$c, sh, n, 1); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n, offset = st)), "fill_normal_at", e, st, n)
      r <- zc_fill_at(k$c, sh, n, 1, -1, 0.5, 1L); chk(r$rc == 0 && r$guard_ok && identical(r$x, rng_normal(k$r, n, -1, 0.5, "mcfarland", offset = st)), "fill_normal_at mcf", e, st, n)
    }
  # pieces vs one
  full <- rng_normal(k$r, 30000)
  set.seed(1); cuts <- sort(c(0, sample(29999, 40), 30000)); acc <- numeric()
  for (i in seq_len(length(cuts) - 1)) acc <- c(acc, zc_fill_at(k$c, cuts[i], cuts[i + 1] - cuts[i], 1)$x)
  chk(identical(acc, full), "pieces normal", e)
  r <- zc_par(k$c, "0", 30000, 37, 1, 0L); chk(r$rc == 0 && identical(r$x, full), "par pieces", e)
  r <- zc_par(k$c, "100", 30000, 64, 0, 0L); chk(r$rc == 0 && identical(r$x, rng_uniform(k$r, 30000, offset = 100)), "par pieces unif", e)
  # near 2^64: dump for the python oracle
  for (st in c("0xffffffffffffe000", "0xfffffffffffffbfd", "0x8000000000000000", "0x0020000000000000", "0x001ffffffffffff0")) {
    n <- if (st == "0xfffffffffffffbfd") 1026 else 3000
    u <- zc_fill_at(k$c, st, n, 0); z <- zc_fill_at(k$c, st, n, 1); m <- zc_fill_at(k$c, st, n, 1, 0, 1, 1L)
    chk(u$rc == 0 && z$rc == 0 && m$rc == 0 && u$guard_ok && z$guard_ok && m$guard_ok, "high start rc", e, st)
    writeBin(c(u$x, z$x, m$x), file.path(outdir, paste0("hi_", e, "_", st)))
  }
  # wrap
  r <- zc_fill_at(k$c, "0xffffffffffffffff", 1, 0); cat(e, "start=2^64-1 n=1 rc", r$rc, "untouched", is.nan(r$x), "\n")
  r <- zc_fill_at(k$c, "0xffffffffffffffff", 0, 0); cat(e, "start=2^64-1 n=0 rc", r$rc, "\n")
  r <- zc_fill_at(k$c, "0xfffffffffffffffe", 1, 0); cat(e, "start=2^64-2 n=1 rc", r$rc, r$x, "\n")
  r <- zc_fill_at(k$c, "0xfffffffffffffffe", 2, 0); cat(e, "start=2^64-2 n=2 rc", r$rc, "untouched", all(is.nan(r$x)), "\n")
  r <- zc_fill_at(k$c, "0xfffffffffffffff0", 17, 1); cat(e, "start=2^64-16 n=17 rc", r$rc, "untouched", all(is.nan(r$x)), "\n")
  r <- zc_fill_at(k$c, "0xfffffffffffffff0", 17, 1, 5, 0); cat(e, "start=2^64-16 n=17 sd=0 rc", r$rc, "untouched", all(is.nan(r$x)), "\n")
  # invalid arguments: out must be untouched
  inv <- list(
    list(0, NaN, 1), list(0, 0, NaN), list(0, -Inf, 1), list(0, 0, Inf), list(0, 2, 1), list(0, 1, 1 - 2^-53),
    list(1, NaN, 1), list(1, 0, NaN), list(1, 0, -1), list(1, 0, -2^-1074), list(1, Inf, 1), list(1, 0, Inf), list(1, 0, -0.0))
  for (a in inv) {
    r1 <- zc_fill(k$c, 5, a[[1]], a[[2]], a[[3]]); r2 <- zc_fill_at(k$c, "9", 5, a[[1]], a[[2]], a[[3]])
    s <- zc_stream(k$c, "20", 5, a[[1]], a[[2]], a[[3]])
    exp_rc <- if (identical(a, list(1, 0, -0.0))) 0L else 1L
    chk(r1$rc == exp_rc && r2$rc == exp_rc && s$rc == exp_rc, "invalid rc", e, unlist(a), r1$rc, r2$rc, s$rc)
    if (exp_rc == 1L) chk(all(is.nan(r1$x)) && all(is.nan(r2$x)) && s$calls == 0, "invalid touched out", e, unlist(a))
    else chk(identical(r1$x, rep(0, 5)) && identical(r2$x, rep(0, 5)), "sd=-0 values", e)
  }
  for (m in c(-1L, 2L, 255L, 256L, .Machine$integer.max)) {
    r1 <- zc_fill(k$c, 5, 4, 0, 1, m); r2 <- zc_fill_at(k$c, "9", 5, 1, 0, 1, m); s <- zc_stream(k$c, "20", 5, 1, 0, 1, m)
    chk(r1$rc == 1 && r2$rc == 1 && s$rc == 1 && all(is.nan(r1$x)) && all(is.nan(r2$x)) && s$calls == 0, "invalid method", e, m)
  }
  for (bk in list(c(3L, 1L), c(-1L, 1L), c(match(e, engs) - 1L, 0L), c(match(e, engs) - 1L, 2L), c(match(e, engs) - 1L, -1L), c(256L, 1L))) {
    kk <- k$c; kk$engine <- bk[1]; kk$stream <- bk[2]
    rcs <- c(zc_fill(kk, 5, 0)$rc, zc_fill(kk, 5, 1)$rc, zc_fill(kk, 5, 2, imin = 1L, imax = 6L)$rc, zc_fill(kk, 5, 3)$rc, zc_fill(kk, 5, 4)$rc,
             zc_fill_at(kk, "1", 5, 0)$rc, zc_fill_at(kk, "1", 5, 1)$rc, zc_stream(kk, "9", 4, 0)$rc, zc_stream(kk, "9", 4, 1)$rc, attr(zc_fold(kk, "1"), "rc"),
             zc_fill(kk, 5, 0, NaN, 1)$rc, zc_fill(kk, 5, 4, 0, 1, 9L)$rc)
    chk(all(rcs == 2L), "bad key rc", e, bk, ":", rcs)
    f <- zc_fold(kk, "1"); chk(f$engine == -7L, "fold wrote out on bad key")
  }
  r <- zc_fill(k$c, 5, 2, imin = NA_integer_, imax = 5L); chk(r$rc == 1 && all(r$x == -777777777L), "int NA min")
  r <- zc_fill(k$c, 5, 2, imin = 1L, imax = NA_integer_); chk(r$rc == 1 && all(r$x == -777777777L), "int NA max")
  r <- zc_fill(k$c, 5, 2, imin = 6L, imax = 1L); chk(r$rc == 1 && all(r$x == -777777777L), "int min>max")
  r <- zc_fill(k$c, 5, 2, imin = 6L, imax = 6L); chk(r$rc == 0 && all(r$x == 6L), "int min==max")
  r <- zc_fill(k$c, 9, 2, imin = -2147483647L, imax = -2147483647L); chk(r$rc == 0 && all(r$x == -2147483647L), "int min==max low")
  nl <- zc_null(k$c); cat(e, "null/edge rcs:\n"); print(nl)
  # fold
  for (d in c(0, 1, -1, 2, 2^31 - 1, -2^31, 2^31, 2^53 - 1, -(2^53 - 1))) {
    f <- zc_fold(k$c, sprintf("%.0f", d)); g <- zc_key_get(rng_fold(k$r, d), 0)
    chk(attr(f, "rc") == 0 && identical(f$k0, g$k0) && identical(f$k1, g$k1) && f$engine == g$engine && f$stream == g$stream, "fold_int", e, d)
    if (abs(d) < 2^31) { g <- zc_key_get(rng_fold(k$r, as.integer(d)), 0); chk(identical(f$k0, g$k0) && identical(f$k1, g$k1), "fold_int vs int", e, d) }
  }
  # streams
  for (tot in c(0, 1, 5, 511, 512, 513, 4096, 5120, 5121, 12000))
    for (ch in c(1, 2, 7, 511, 512, 513, 1000, 4096, 5120, 5121, 20000)) {
      if (tot / ch > 4000) next
      for (cfg in list(list(0, 0, 1, 0L), list(0, -2, 9, 0L), list(1, 0, 1, 0L), list(1, 1, 2, 1L), list(0, 4, 4, 0L), list(1, 4, 0, 0L), list(0, -.Machine$double.xmax, .Machine$double.xmax, 0L))) {
        s <- zc_stream(k$c, sprintf("%.0f", tot), ch, cfg[[1]], cfg[[2]], cfg[[3]], cfg[[4]])
        ref <- if (cfg[[1]] == 0) rng_uniform(k$r, tot, cfg[[2]], cfg[[3]]) else rng_normal(k$r, tot, cfg[[2]], cfg[[3]], if (cfg[[4]] == 1L) "mcfarland" else "ziggurat")
        ncall <- ceiling(tot / ch)
        en <- if (tot == 0) numeric() else c(rep(ch, ncall - 1), tot - ch * (ncall - 1))
        ok <- s$rc == 0 && s$guard_ok && identical(s$x, ref) && s$calls == ncall && identical(s$n, as.numeric(en)) &&
              identical(s$start, if (tot == 0) numeric() else (seq_len(ncall) - 1) * ch) && s$bufbad == 0 && s$seqbad == 0
        chk(ok, "stream", e, tot, ch, unlist(cfg), "rc", s$rc, "calls", s$calls, ncall)
      }
    }
  for (sa in c(1L, 2L, 3L)) {
    s <- zc_stream(k$c, "12000", 1000, 1, 0, 1, 0L, stop_after = sa)
    chk(s$rc == 3 && s$calls == sa && identical(s$x, rng_normal(k$r, 1000 * sa)), "stream stop", e, sa)
    s <- zc_stream(k$c, "12000", 1000, 0, 5, 5, 0L, stop_after = sa); chk(s$rc == 3 && s$calls == sa, "stream stop const", e, sa)
  }
  s <- zc_stream(k$c, "1000", 1000, 1, 0, 1, 0L, stop_after = 1L); cat(e, "stop on the last chunk: rc", s$rc, "\n")
  # huge total, stop early
  s <- zc_stream(k$c, "0xffffffffffffffff", 5120, 0, 0, 1, 0L, stop_after = 3L); chk(s$rc == 3 && identical(s$x, rng_uniform(k$r, 15360)), "stream huge total", e)
  s <- zc_stream(k$c, "0xffffffffffffffff", 3, 0, 4, 4, 0L, stop_after = 3L); chk(s$rc == 3 && identical(s$x, rep(4, 9)), "stream huge total const", e)
}
cat("capi BAD =", bad, "\n")
