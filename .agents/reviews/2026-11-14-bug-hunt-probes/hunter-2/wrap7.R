source(if (file.exists("/w/hunt-2/common.R")) "/w/hunt-2/common.R" else "common.R")
key <- rng_key(3L); n <- 100; e <- rng_normal(key, n, 1, 2)
steps <- list(A = quote(chk("v", identical(w[c(1, n, 7:40)], e[c(1, n, 7:40)]))), B = quote(chk("s", identical(sum(w), sum(e)))), C = quote(chk("m", !mat(x))), D = quote(s <- serialize(w, NULL)), E = quote(chk("u", identical(unserialize(s)[], e))))
run <- function(which, mkx) { x <<- mkx(); w <<- .Internal(wrap_meta(x, 0L, TRUE)); for (st in steps[which]) eval(st, globalenv()); rm(x, envir = globalenv())
  tryCatch({ eval(quote(w[3] <- 7), globalenv()); if (identical(w[], replace(e, 3, 7))) "ok" else "WRONG" }, error = function(e) conditionMessage(e)) }
lz <- function() rng_lazy_normal(key, n, 1, 2); pl <- function() rng_normal(key, n, 1, 2)
subs <- unlist(lapply(2:4, function(k) combn(c("A","B","C","D","E"), k, simplify = FALSE)), recursive = FALSE); for (wh in subs) { r <- run(wh, lz); if (r != "ok") cat(paste(wh, collapse = ""), r, "\n") }
