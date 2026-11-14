source(if (file.exists("/w/hunt-2/common.R")) "/w/hunt-2/common.R" else "common.R")
if (LINUX) cat("THP:", readLines("/sys/kernel/mm/transparent_hugepage/enabled"), "\n")
cat("threads", rng_threads(), "simd", rng_simd(), "\n")
# rng_threads contract
t0 <- rng_threads(); old <- rng_threads(1L); chk("old 0", identical(old, 0L)); chk("now 1", rng_threads() == 1L)
o2 <- rng_threads(old); chk("o2 1", identical(o2, 1L)); chk("restored", rng_threads() == t0)
o <- rng_threads(1e10); chk("big cap", rng_threads() == t0); o <- rng_threads(o); chk("ret big", is.integer(o) && !is.na(o)); 
o <- rng_threads(2); chk("cap 2", rng_threads() == min(2L, t0)); rng_threads(0)
for (bad in list(NA, NA_integer_, -1, 1.5, "2", c(1,2), NULL, NaN, Inf, TRUE, integer(0))) { r <- tryCatch({rng_threads(bad); "ok"}, error = function(e) "err"); cat("rng_threads(", deparse(bad), ") ->", r, " now", rng_threads(), "\n") }
rng_threads(0)
for (bad in list(NA, 1, "a", c(TRUE,FALSE), 0L)) { r <- tryCatch({rng_simd(bad); "ok"}, error = function(e) "err"); cat("rng_simd(", deparse(bad), ") ->", r, " now", rng_simd(), "\n") }
rng_simd(TRUE)
# failed setter must not change state
rng_threads(2L); try(rng_threads(-3), silent = TRUE); chk("failed set keeps cap", identical(rng_threads(0L), 2L))
# hugepage straddle x threads
for (eng in engines) { key <- rng_key(9L, engine = eng)
 for (n in c(524287, 524288, 524289, 1048576+3)) {
  rng_threads(1L); options(zurand.hugepages = FALSE); ref_u <- rng_uniform(key, n); ref_n <- rng_normal(key, n); ref_m <- rng_normal(key, n, 1, 2, method = "mcfarland")
  for (hp in list(TRUE, FALSE, NULL, NA, "no", 0L)) for (th in c(0L, 1L, 2L, 3L)) { options(zurand.hugepages = hp); rng_threads(th); lab <- paste(eng, n, deparse(hp), th)
    chk(paste(lab,"u"), identical(rng_uniform(key, n), ref_u)); chk(paste(lab,"n"), identical(rng_normal(key, n), ref_n))
    x <- rng_lazy_normal(key, n, 1, 2, method = "mcfarland"); chk(paste(lab,"lazy mat"), identical(x + 0, ref_m)); chk(paste(lab, "lazy ser"), length(serialize(x, NULL)) < 2000)
    x <- rng_lazy_uniform(key, n); chk(paste(lab,"lazy region"), identical(x[2:n], ref_u[-1])) }
 } }
options(zurand.hugepages = NULL); rng_threads(0L)
# rng_bits 64 interrupt
pid <- Sys.getpid(); system(sprintf("(sleep 1; kill -INT %d)", pid), wait = FALSE); t0 <- proc.time()[[3]]
r <- tryCatch({ b <- rng_bits(rng_key(1L), 3e8, bits = 64); "completed" }, interrupt = function(c) "interrupted"); cat("bits64:", r, proc.time()[[3]] - t0, "\n")
chk("bits64 interrupted", r == "interrupted")
b1 <- rng_bits(rng_key(1L), 6000, bits = 64); b2 <- rng_bits(rng_key(1L), 6000, bits = 64); chk("bits after intr", identical(b1, b2) && !anyNA(b1) && all(nchar(b1) == 16))
cat("failures:", nfail, "\n")
