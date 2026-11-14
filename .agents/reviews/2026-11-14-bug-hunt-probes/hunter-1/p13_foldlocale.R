source("helper.R")
k <- rng_key(1L, engine="philox4x64")
nm <- readLines("names.txt")
cat(Sys.getlocale("LC_CTYPE"), "| Encoding:", Encoding(nm), "| fold(readLines):", format(rng_fold(k, nm)), "| fold('caf<c3><a9>'):", format(rng_fold(k, "caf<c3><a9>")), "| draw:", format(rng_uniform(rng_fold(k, nm)), digits=17), "\n")
