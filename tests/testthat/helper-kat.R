# Helpers for test-kat.R (the known-answer tests), kept out of it so that the
# tests run shuffled.

# A key is an opaque S3 object;;these tests reach past that on purpose, to
# set the exact key words a reference vector requires. Package code must
# not do this.
raw_key <- function(words) {
  structure(matrix(as.integer(words), nrow = 1L),
            class = "rng_key", engine = "philox4x64")
}
