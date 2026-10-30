# Helpers for test-digest.R (the whole-stream digests), kept out of it so that the
# tests run shuffled.

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
# rng_normal(method = "mcfarland"), computed 2026-09-27 on macOS x86_64,
# after the edge began reusing the fast-path word (#46).
# Per 1e6 draws its edge takes ~11,700 of them through every overhang kind
# and ~280 through the tail (log1p), and calls exp() about 3,000 times.
digest_mcfarland <- c(
  xoshiro256pp = "282ba89bad7e711f730832366eef95c7",
  philox4x64   = "7234891b18858dc38e3c31dc5aabe2a9",
  threefry4x64 = "4ca018a6a95d2d07ecb49fbf351f75df"
)
digest_mcfarland_scaled <- "d12db748d5fd5a07d0650f7815151ebe"
digest_uniform_scaled <- "685b24f007ecea32b544a1cb9495699f"
digest_integer        <- "e76d991a8fbdce2645b8f69bf06f6852"
