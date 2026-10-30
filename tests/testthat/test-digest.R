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

test_that("rng_normal(method = \"mcfarland\") streams agree bit for bit", {
  skip_if_x87()
  n <- 1e6
  for (engine in names(digest_mcfarland)) {
    key <- rng_key(20260926L, engine = engine)
    expect_identical(stream_md5(rng_normal(key, n, method = "mcfarland")),
                     digest_mcfarland[[engine]],
                     label = paste("rng_normal(method = \"mcfarland\") digest,", engine))
  }
  key <- rng_key(20260926L, engine = "xoshiro256pp")
  expect_identical(stream_md5(rng_normal(key, n, mean = 0.3, sd = 1.7, method = "mcfarland")),
                   digest_mcfarland_scaled)
})

test_that("rng_uniform() and rng_integer() streams agree bit for bit", {
  skip_if_x87()
  key <- rng_key(20260926L, engine = "xoshiro256pp")
  expect_identical(stream_md5(rng_uniform(key, 1e6, min = -0.3, max = 2.9)),
                   digest_uniform_scaled)
  expect_identical(stream_md5(rng_integer(key, 1e6, min = -7L, max = 1000L)),
                   digest_integer)
})
