# The normal samplers settle most edge draws with integer shortcuts and
# leave exp() a band around each chord. Golden and digest tests prove the
# streams are reproducible; these prove the shortcuts are right: every
# decision a shortcut makes must agree with the plain acceptance test,
# y < exp(-x^2 / 2), computed here in R from the defining tables alone.
# Draws are built with zurand itself, so the tests are deterministic and
# leave .Random.seed alone.

tables <- function() .Call(C_zurand_normal_tables)

shortcut <- function(method, layer, a, b) {
  .Call(C_zurand_shortcut, method, as.integer(layer),
        a$hi, a$lo, b$hi, b$lo)
}

# A word from its high and low 32-bit halves, and back.
halves <- function(hi, lo) list(hi = hi, lo = lo)
split_word <- function(w) {        # w: a whole number below 2^53
  hi <- floor(w / 2^32)
  halves(hi, w - hi * 2^32)
}

# The plain test, with a tolerance of a few ulps of the curve: a decision
# closer to the curve than rounding can resolve is not a contradiction.
contradictions <- function(decision, x, y) {
  f <- exp(-x * x / 2)
  inside <- y < f
  decided <- decision != 0L
  clear <- abs(y - f) > 8 * .Machine$double.eps * f
  sum(decided & clear & (decision > 0L) != inside)
}

test_that("ziggurat wedge shortcuts agree with the exp() test", {
  t <- tables()
  ki <- t$ki
  wi <- t$wi
  fi <- t$fi
  layers <- rep(1:255, each = 400L)
  n <- length(layers)
  key <- rng_key(20260928L, engine = "philox4x64")
  u <- rng_uniform(key, 4L * n)

  # A magnitude anywhere in the wedge, rabs in [ki, 2^52).
  lo_edge <- ki[layers + 1L]
  rabs <- pmin(lo_edge + floor(u[seq_len(n)] * (2^52 - lo_edge)), 2^52 - 1)
  # The chord's height at rabs, as a wedge fraction: R / L in the C code.
  chord <- (2^52 - rabs) / (2^52 - lo_edge)
  # Half the draws uniform over the wedge's height, half within 1e-3 to
  # 1e-12 of the chord on either side, where the bands and guard act.
  near <- seq_len(n) %% 2L == 0L
  scale <- 10^-c(3, 6, 9, 12)[(seq_len(n) %% 4L) + 1L]
  sign <- ifelse(u[n + seq_len(n)] < 0.5, -1, 1)
  frac <- ifelse(near, chord + sign * scale * u[2L * n + seq_len(n)],
                 u[3L * n + seq_len(n)])
  frac <- pmin(pmax(frac, 0), 1 - 2^-53)
  # Y from the fraction: top 32 bits, then the next 32.
  y_hi <- floor(frac * 2^32)
  y_lo <- floor((frac * 2^32 - y_hi) * 2^32)

  decision <- shortcut(0L, layers, split_word(rabs), halves(y_hi, y_lo))

  # Independently: u = (Y >> 11) * 2^-53, the height fi[idx] + u * rise,
  # the abscissa rabs * wi[idx].
  u_y <- (y_hi * 2^21 + floor(y_lo / 2^11)) * 2^-53
  x <- rabs * wi[layers + 1L]
  y <- fi[layers + 1L] + (fi[layers] - fi[layers + 1L]) * u_y
  expect_identical(contradictions(decision, x, y), 0L)

  # The shortcuts do the work: both outcomes, and exp() only rarely.
  expect_gt(sum(decision == 1L), n / 10)
  expect_gt(sum(decision == -1L), n / 10)
  expect_lt(mean(decision == 0L), 0.25)
  expect_gt(sum(decision == 0L & near), 0)
  # All three regimes: below, at and above the inflection layer.
  for (side in list(layers < t$zig_inflection, layers == t$zig_inflection,
                    layers > t$zig_inflection)) {
    expect_true(any(decision[side] == 1L) && any(decision[side] == -1L))
  }
})

test_that("McFarland edge shortcuts agree with the exp() test", {
  t <- tables()
  mx <- t$mcf_x
  my <- t$mcf_y
  overhangs <- rep(seq_len(length(mx) - 1L), each = 400L)  # j = 1 .. 253
  n <- length(overhangs)
  key <- rng_key(20260928L, engine = "threefry4x64")
  bits <- rng_bits(key, 4L * n)                           # uint32 as double
  u1 <- halves(bits[seq_len(n)] %% 2^31, bits[n + seq_len(n)])
  u2 <- halves(bits[2L * n + seq_len(n)] %% 2^31, bits[3L * n + seq_len(n)])
  value <- function(w) w$hi * 2^32 + w$lo                 # near-exact, 63 bits

  # The concave overhangs see the attempt ordered u1 <= u2, as the sampler
  # arranges it by reflection.
  concave <- overhangs < t$mcf_inflection
  swap <- concave & value(u2) < value(u1)
  a <- halves(ifelse(swap, u2$hi, u1$hi), ifelse(swap, u2$lo, u1$lo))
  b <- halves(ifelse(swap, u1$hi, u2$hi), ifelse(swap, u1$lo, u2$lo))

  decision <- shortcut(1L, overhangs, a, b)

  # Independently: x and y interpolate the layer edges (the tables are
  # scaled by 2^-63), and the draw is inside when y < exp(-x^2 / 2).
  j <- overhangs + 1L                                     # 1-based X_j
  x <- mx[j] * 2^63 + value(a) * (mx[j - 1L] - mx[j])
  y <- my[j] * 2^63 + value(b) * (my[j - 1L] - my[j])
  expect_identical(contradictions(decision, x, y), 0L)

  convex <- overhangs > t$mcf_inflection
  expect_gt(sum(decision[convex] == 1L), 0)
  expect_gt(sum(decision[convex] == -1L), 0)
  expect_gt(sum(decision[concave] == 1L), 0)
  expect_true(all(decision[overhangs == t$mcf_inflection] == 0L))
  expect_true(all(decision[concave] >= 0L))
})

test_that("the shortcut entry point refuses draws outside the edges", {
  one <- halves(0, 0)
  expect_error(shortcut(0L, 0L, one, one), "ziggurat wedge")
  expect_error(shortcut(0L, 128L, one, one), "ziggurat wedge")  # rabs < ki
  expect_error(shortcut(1L, 254L, one, one), "McFarland overhang")
  expect_error(shortcut(1L, 5L, halves(2^31, 0), one), "63 bits")
  expect_error(shortcut(0L, 5L, halves(2^20, 0), one), "52 bits")
})
