source("helper.R")
k <- rng_key(1L, engine="philox4x64")
E <- "philox4x64"
f <- function(...) format(rng_key(..., engine=E))
# seed matrix
vals <- list(int=5L, dbl=5, big=123456789.0, g31=2^31+1, g53m=2^53-1, g53=2^53, g53p=2^53+2, neg=-1, frac=1.5, negzero=-0, NAi=NA_integer_, NAr=NA_real_, NA_=NA, NaN_=NaN, Inf_=Inf, lgl=TRUE, chr="5", fac=factor("5"), fac2=factor(c("a","b"))[2], date=as.Date("1970-01-06"), posix=as.POSIXct(5, origin="1970-01-01", tz="UTC"), difft=as.difftime(5, units="secs"), lst=list(5), null=NULL, len0=integer(0), len2=c(5L,6L), named=c(a=5), cplx=5+0i, raw=as.raw(5), i64=structure(2, class="integer64"), mat=matrix(5), arr=array(5, c(1,1,1)), nv=numeric_version("5"), ts=ts(5), tiny=5+1e-12, e=5e0, hex=0x5)
for (nm in names(vals)) { cat(sprintf("seed %-8s: ", nm)); r <- tryCatch(f(vals[[nm]]), error=function(e) paste("ERR", conditionMessage(e))); cat(r, "\n") }
cat("ref 5: ", f(5L), " ref 2: ", f(2), "\n")
for (nm in names(vals)) { cat(sprintf("n    %-8s: ", nm)); r <- tryCatch(length(rng_key(1, n=vals[[nm]], engine=E)), error=function(e) paste("ERR", conditionMessage(e))); cat(r, "\n") }
for (nm in names(vals)) { cat(sprintf("unif n %-8s: ", nm)); r <- tryCatch(length(rng_uniform(k, vals[[nm]])), error=function(e) paste("ERR", conditionMessage(e))); cat(r, "\n") }
for (nm in names(vals)) { cat(sprintf("offset %-8s: ", nm)); r <- tryCatch(rng_uniform(k, 1, offset=vals[[nm]]), error=function(e) paste("ERR", conditionMessage(e))); cat(format(r, digits=17), "\n") }
cat("ref offset 5:", format(rng_uniform(k, 6)[6], digits=17), "\n")
for (nm in names(vals)) { cat(sprintf("min    %-8s: ", nm)); r <- tryCatch(rng_uniform(k, 1, min=vals[[nm]], max=2^60), error=function(e) paste("ERR", conditionMessage(e))); cat(format(r, digits=17), "\n") }
for (nm in names(vals)) { cat(sprintf("imax   %-8s: ", nm)); r <- tryCatch(rng_integer(k, 1, min=0, max=vals[[nm]]), error=function(e) paste("ERR", conditionMessage(e))); cat(format(r, digits=17), "\n") }
for (nm in names(vals)) { cat(sprintf("sd     %-8s: ", nm)); r <- tryCatch(rng_normal(k, 1, sd=vals[[nm]]), error=function(e) paste("ERR", conditionMessage(e))); cat(format(r, digits=17), "\n") }
for (nm in names(vals)) { cat(sprintf("bits   %-8s: ", nm)); r <- tryCatch(rng_bits(k, 1, bits=vals[[nm]]), error=function(e) paste("ERR", conditionMessage(e))); cat(format(r, digits=17), "\n") }
for (b in list(32L, 32, 64L, 64, 32.5, "32", TRUE, c(32,64), NULL, NA, factor(32))) { r <- tryCatch(rng_bits(k, 1, bits=b), error=function(e) paste("ERR", conditionMessage(e))); cat("bits", deparse(b), ":", format(r), "\n") }
for (nm in names(vals)) { cat(sprintf("thr    %-8s: ", nm)); r <- tryCatch({o <- rng_threads(vals[[nm]]); c(o, rng_threads())}, error=function(e) paste("ERR", conditionMessage(e))); cat(format(r), "\n"); rng_threads(0) }
for (nm in names(vals)) { cat(sprintf("fold   %-8s: ", nm)); r <- tryCatch(format(rng_fold(k, vals[[nm]])), error=function(e) paste("ERR", conditionMessage(e))); cat(r, "\n") }
cat("fold ref 5L", format(rng_fold(k, 5L)), "\n")
# engine
for (e in list("philox", "x", "philox4x64 ", c("philox4x64","threefry4x64"), NULL, NA, 1, character(0), factor("philox4x64"), "PHILOX4X64", "")) { r <- tryCatch(attr(rng_key(1, engine=e), "engine"), error=function(e) paste("ERR", conditionMessage(e))); cat("engine", deparse(e), ":", r, "\n") }
for (m in list("mc", "z", "", NULL, c("mcfarland","ziggurat"), NA, 1L, "Mcfarland")) { r <- tryCatch(rng_normal(k, 1, method=m), error=function(e) paste("ERR", conditionMessage(e))); cat("method", deparse(m), ":", format(r, digits=17), "\n") }
cat("zig", format(rng_normal(k,1), digits=17), "mcf", format(rng_normal(k,1,method="mcfarland"), digits=17), "\n")
# simd
for (s in list(NULL, TRUE, FALSE, NA, 1L, 0, "avx2", "none", c(TRUE,FALSE), logical(0))) { r <- tryCatch({o <- rng_simd(s); paste(o, "->", rng_simd())}, error=function(e) paste("ERR", conditionMessage(e))); cat("simd", deparse(s), ":", r, "\n"); rng_simd(TRUE) }
