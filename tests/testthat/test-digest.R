# Whole-stream digests: every bit of a long draw, on every platform CI
# reaches. The golden tests pin a few hundred values; these pin millions.
#
# What they were written for is libm. The ziggurat's tail needs log1p()
# and the ambiguous band of its wedge needs exp(); per 1e6 normals that is
# about 500 log1p() and 930 exp() calls for each engine. Against the
# platform's libm these digests split CI in two -- Apple and Windows one
# way, glibc and musl the other -- because a one-ulp difference can flip an
# accept/reject decision. zurand now computes both with its own fdlibm port
# (src/zurand_fdlibm.h), so a leg that disagrees here has broken IEEE
# double evaluation somewhere -- a fused multiply-add, or x87 excess
# precision -- and the fix is in the code, not the digest.
#
# Uniform and integer draws make no libm calls; they are here because a
# digest is cheap and the golden tests sample only their first few values.
#
# Every key names its engine explicitly, so changing the default engine
# does not touch these. Changing a stream does, and then the digests are
# regenerated once, in the commit that says so.
# 32-bit x86 with x87 arithmetic rounds each addition twice (to 64 bits, then
# to 53 on store), and is outside the reproducibility contract
# (dev/design.md, section 7). Its CI leg still runs everything else.
skip_if_x87 <- function() {
  if (grepl("^i[3-6]86$", R.version$arch))
    skip("32-bit x86 (x87): outside the bit-reproducibility contract")
}

stream_md5 <- function(x) {
  # tools::md5sum(bytes = ) needs R >= 4.5; CI's i386 leg runs 4.2.
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(x, f, endian = "little")
  unname(tools::md5sum(f))
}

# Computed 2026-09-26 on macOS x86_64 with the fdlibm port. The philox and
# threefry digests equal what every glibc and musl leg produced with the
# system libm before the port; the xoshiro-keyed ones were regenerated when
# the xoshiro engine got its own counter tag (dev/design.md, section 3.2).
digest_normal <- c(
  xoshiro256pp = "47dbe72c26c885915260b168f01bae0a",
  philox4x64   = "de08c301e6229eb344d9c570dcf23535",
  threefry4x64 = "257bbe55425f91e66d2302cd6b3660fe"
)
digest_normal_scaled  <- "3db5888105cfe7d7d34d1c01491dd4fc"
digest_uniform_scaled <- "685b24f007ecea32b544a1cb9495699f"
digest_integer        <- "e76d991a8fbdce2645b8f69bf06f6852"

test_that("rng_normal() streams agree bit for bit, libm paths included", {
  skip_if_x87()
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
  skip_if_x87()
  key <- rng_key(20260926L, engine = "xoshiro256pp")
  expect_identical(stream_md5(rng_uniform(key, 1e6, min = -0.3, max = 2.9)),
                   digest_uniform_scaled)
  expect_identical(stream_md5(rng_integer(key, 1e6, min = -7L, max = 1000L)),
                   digest_integer)
})
