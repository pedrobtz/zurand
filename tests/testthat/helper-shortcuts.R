# Helpers for test-shortcuts.R (the ziggurat shortcut tests), kept out of it so that the
# tests run shuffled.

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
