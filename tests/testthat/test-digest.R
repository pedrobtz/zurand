# Whole-stream digests: every bit of a long draw, on every platform CI
# reaches. The golden tests pin a few hundred values; these pin millions.
#
# What they are for is the platform's libm. The ziggurat's tail calls
# log1p() and the ambiguous band of its wedge calls exp(); per 1e6 normals
# that is about 500 log1p() and 930 exp() calls for each engine. glibc,
# musl, macOS and Windows do not promise correctly rounded results, and a
# one-ulp disagreement at the wrong point flips an accept/reject decision,
# which changes that draw completely. If any CI leg disagrees here, the fix
# is to vendor a deterministic exp/log1p (roadmap A4), not to update the
# digest.
#
# Uniform and integer draws make no libm calls; they are here because a
# digest is cheap and the golden tests sample only their first few values.
#
# Every key names its engine explicitly, so changing the default engine
# does not touch these. Changing a stream does, and then the digests are
# regenerated once, in the commit that says so.
stream_md5 <- function(x) {
  # tools::md5sum(bytes = ) needs R >= 4.5; CI's i386 leg runs 4.2.
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(x, f, endian = "little")
  unname(tools::md5sum(f))
}

# Computed 2026-09-26 on macOS x86_64 (Apple libm), R 4.5.2.
digest_normal <- c(
  xoshiro256pp = "cb7b0bd1f39e5027bccada544db7245a",
  philox4x64   = "a8fbaa4906aab3013560d36bebb5494a",
  threefry4x64 = "257bbe55425f91e66d2302cd6b3660fe"
)
digest_normal_scaled  <- "48326dd8120757d4b30fec0b22ed7391"
digest_uniform_scaled <- "5a74357222df42c83f569102ead0c79e"
digest_integer        <- "0c68ddc22108c826117675ace3faa71c"

test_that("rng_normal() streams agree bit for bit, libm paths included", {
  n <- 1e6
  for (engine in names(digest_normal)) {
    key <- rng_key(20260926L, engine = engine)
    expect_identical(stream_md5(rng_normal(key, n)), digest_normal[[engine]],
                     label = paste("rng_normal digest,", engine))
  }
  key <- rng_key(20260926L, engine = "xoshiro256pp")
  expect_identical(stream_md5(rng_normal(key, n, mean = 0.3, sd = 1.7)),
                   digest_normal_scaled)
})

test_that("rng_uniform() and rng_integer() streams agree bit for bit", {
  key <- rng_key(20260926L, engine = "xoshiro256pp")
  expect_identical(stream_md5(rng_uniform(key, 1e6, min = -0.3, max = 2.9)),
                   digest_uniform_scaled)
  expect_identical(stream_md5(rng_integer(key, 1e6, min = -7L, max = 1000L)),
                   digest_integer)
})
