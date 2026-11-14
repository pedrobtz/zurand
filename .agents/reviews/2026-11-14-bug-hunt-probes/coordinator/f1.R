.libPaths("HUNT/rerun/lib"); library(zurand); k <- rng_key(1L, engine="philox4x64")
s <- readLines("HUNT/rerun/names.txt"); cat(Sys.getlocale("LC_CTYPE"), Encoding(s), format(rng_fold(k, s)), format(rng_fold(k, "caf<c3><a9>")), "\n")
