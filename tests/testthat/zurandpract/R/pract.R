# Tests key 1 of `keys` with zurngtest's pract_test(); `options` is a named
# list of its other arguments.
pract_zurand <- function(keys, options = NULL) {
  .Call(pract_zurand_c, keys, options)
}
