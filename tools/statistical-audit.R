#!/usr/bin/env Rscript
# Write a raw binary stream from one zurand sampler to stdout, for a
# statistical test battery (PractRand, dieharder, TestU01) to consume.
#
#   Rscript tools/statistical-audit.R --sampler=normal --engine=philox4x64 \
#           --gb=4 | RNG_test stdin32
#
# Why each sampler is worth testing separately:
#
#   bits     The generator as zurand drives it. Philox and Threefry are
#            well studied and test-kat.R already proves this is really
#            them, so what this exercises is the counter layout -- that
#            {index, domain, purpose, 0} does not accidentally correlate
#            streams.
#   uniform  Adds u01_open(), the mantissa-stuffing conversion.
#   normal   Adds the ziggurat, including the fixed-point wedge shortcut
#            and the tail, which are this package's own code and the least
#            scrutinised thing in it. Normals are mapped back to uniforms
#            with pnorm() before testing, so a battery built for uniform
#            input can see a distributional error.
#
# Output is little-endian uint32. Chunked, so memory stays bounded however
# large --gb is.

suppressMessages(library(zurand))

args <- commandArgs(TRUE)
getarg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1]])
}

sampler <- match.arg(getarg("sampler", "normal"), c("normal", "uniform", "bits"))
engine  <- match.arg(getarg("engine", "philox4x64"), c("philox4x64", "threefry4x64", "xoshiro256pp"))
gb      <- as.numeric(getarg("gb", "1"))
seed    <- as.integer(getarg("seed", "20260919"))
chunk   <- 1e7                                   # values per block, ~80 MB

total_words <- ceiling(gb * 2^30 / 4)            # one uint32 per value
key <- rng_key(seed, engine = engine)
# Destination. --out=- streams to stdout, anything else is a file path.
#
# Two traps here, both of which fail by writing zero bytes and exiting 0:
# file("stdout", "wb") does not work under Rscript (pipe("cat", "wb") is
# the working equivalent), and on.exit() does nothing at script top level,
# so the connection has to be closed explicitly or the buffer is discarded.
out <- getarg("out", "-")
con <- if (out == "-") pipe("cat", "wb") else file(out, "wb")

# uint32 -> the same 32 bits as a signed R integer, which writeBin emits
# verbatim. Values at or above 2^31 wrap to negative; the bit pattern is
# unchanged, which is all the test battery sees.
as_u32 <- function(x) as.integer(x - 4294967296 * (x >= 2147483648))

done <- 0
block <- 0
while (done < total_words) {
  n <- as.integer(min(chunk, total_words - done))
  # A fresh key per block via rng_fold keeps memory flat while still
  # producing one continuous stream: folding by block index is exactly the
  # substream mechanism the package documents.
  k <- rng_fold(key, block)
  w <- switch(sampler,
    bits    = rng_bits(k, n, bits = 32L),
    uniform = floor(rng_uniform(k, n) * 4294967296),
    # pnorm() maps N(0,1) back to U(0,1); a wrong ziggurat shows up as a
    # non-uniform result here.
    normal  = floor(pnorm(rng_normal(k, n)) * 4294967296))
  w[w > 4294967295] <- 4294967295                # guard the open endpoint
  writeBin(as_u32(w), con, size = 4L, endian = "little")
  done <- done + n
  block <- block + 1
}
close(con)
if (out != "-")
  message(sprintf("wrote %s: %.2f GB of %s from %s",
                  out, file.size(out) / 2^30, sampler, engine))
