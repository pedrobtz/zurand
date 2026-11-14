# Every error zurand raises has the root class `zurand_error`, so callers
# and tests can catch zurand's errors by class rather than by message
# (#62, R9). R code calls zurand_abort(); C code calls zurand_error() in
# src/zurand.c, which evaluates zurand_abort_c() in this namespace.

zurand_abort <- function(..., call = NULL) {
  stop(zurand_condition(paste0(...), call))
}

# Called from C. sys.call(-1) is the R function whose .Call() raised the
# error, so the message reads "Error in rng_uniform(...)" as Rf_error()'s
# did.
zurand_abort_c <- function(message) {
  stop(zurand_condition(message, sys.call(-1L)))
}

zurand_condition <- function(message, call) {
  structure(
    class = c("zurand_error", "error", "condition"),
    list(message = message, call = call)
  )
}
