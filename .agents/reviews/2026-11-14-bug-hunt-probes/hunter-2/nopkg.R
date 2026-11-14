source("/w/hunt-2/common.R")
x <- rng_lazy_normal(rng_key(5L), 1e5); saveRDS(x, "/tmp/x.rds"); saveRDS(list(a = 1, x = x), "/tmp/l.rds")
out <- system2("Rscript", c("--vanilla", "-e", shQuote("z <- readRDS('/tmp/x.rds'); cat('class', class(z), 'length', length(z), '\n'); l <- readRDS('/tmp/l.rds'); str(l)")), stdout = TRUE, stderr = TRUE); cat(out, sep = "\n")
