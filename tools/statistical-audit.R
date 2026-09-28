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
#            input can see a distributional error. --method=mcfarland
#            audits McFarland's modified ziggurat instead of the default.
#            (tests/testthat/test-shortcuts.R checks both samplers' edge
#            shortcuts decision by decision against exp().)
#
# Output is little-endian uint32. Chunked, so memory stays bounded however
# large --gb is.
#
# --mode picks what the stream is made of:
#
#   stream   (default) a new key every 1e7 values, rng_fold(key, block),
#            each read from its first position: a concatenation of the
#            first 1e7 values of many folded keys, not one key's stream.
#   continuous  one key advanced without re-keying, via offset =: block
#            b holds positions b * 1e7 onwards, so the output is exactly
#            the key's stream from position 0, however long the run. Tests
#            a sampler far along a single key. uniform and normal only
#            (the positional samplers).
#   keys     cross-key: each block takes the next 64 consecutive sibling
#            keys of one rng_key(seed, n) vector and interleaves their first
#            --words values (position 0 of all 64, then position 1, ...).
#            A correlation between sibling keys -- the splitmix64 key
#            derivation, or an engine's seeding -- is invisible in any
#            single stream and shows up here.
#   folds    cross-key: the same, with the 64 keys rng_fold(key, i) for
#            consecutive integers i. Tests the fold.
#
# The cross-key modes use rng_bits(), the raw words, and interleave short
# streams because that is how the package is used: many keys, each drawing
# a modest number of values. At the default --words = 1e5, each key spans
# about 200 of xoshiro256pp's 512-word sub-chunk seedings.

suppressMessages(library(zurand))

args <- commandArgs(TRUE)
getarg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1]])
}

sampler <- match.arg(getarg("sampler", "normal"), c("normal", "uniform", "bits"))
mode    <- match.arg(getarg("mode", "stream"), c("stream", "continuous", "keys", "folds"))
method  <- match.arg(getarg("method", "ziggurat"), c("ziggurat", "mcfarland"))
if (mode == "continuous" && sampler == "bits")
  stop("--mode=continuous needs a positional sampler: uniform or normal")
words   <- as.integer(getarg("words", "100000"))  # per key, cross-key modes
nkeys   <- 64L
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

# Cross-key modes: one block = nkeys keys x `words` positions, interleaved
# by position. Returns the block's uint32 values, or NULL when done.
cross_block <- function(block) {
  ids <- block * nkeys + seq_len(nkeys) - 1L          # 0-based, consecutive
  keys <- if (mode == "keys") {
    all_keys[ids + 1L]
  } else {
    do.call(c, lapply(ids, function(i) rng_fold(key, i)))
  }
  m <- rng_bits(keys, words, bits = 32L)             # words x nkeys
  as.vector(t(m))                                    # interleave by position
}
if (mode == "keys") {
  # Every sibling key the run will use, from one seed, generated once.
  nblocks <- ceiling(total_words / (as.numeric(nkeys) * words))
  all_keys <- rng_key(seed, n = nblocks * nkeys, engine = engine)
}

done <- 0
block <- 0
while (mode %in% c("keys", "folds") && done < total_words) {
  w <- cross_block(block)
  w <- w[seq_len(min(length(w), total_words - done))]
  writeBin(as_u32(w), con, size = 4L, endian = "little")
  done <- done + length(w)
  block <- block + 1
}
while (mode %in% c("stream", "continuous") && done < total_words) {
  n <- as.integer(min(chunk, total_words - done))
  # stream: a fresh key per block via rng_fold(), each from position 0 --
  # the documented substream mechanism, but not one key's long stream.
  # continuous: the same key throughout, positions done .. done + n - 1.
  k <- if (mode == "stream") rng_fold(key, block) else key
  at <- if (mode == "stream") 0 else done
  w <- switch(sampler,
    bits    = rng_bits(k, n, bits = 32L),
    uniform = floor(rng_uniform(k, n, offset = at) * 4294967296),
    # pnorm() maps N(0,1) back to U(0,1); a wrong ziggurat shows up as a
    # non-uniform result here.
    normal  = floor(pnorm(rng_normal(k, n, method = method, offset = at)) * 4294967296))
  w[w > 4294967295] <- 4294967295                # guard the open endpoint
  writeBin(as_u32(w), con, size = 4L, endian = "little")
  done <- done + n
  block <- block + 1
}
close(con)
if (out != "-")
  message(sprintf("wrote %s: %.2f GB of %s%s from %s, mode %s",
                  out, file.size(out) / 2^30, sampler,
                  if (sampler == "normal") paste0(" (", method, ")") else "",
                  engine, mode))
