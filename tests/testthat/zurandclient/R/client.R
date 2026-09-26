client_fill <- function(keys, n, what, a = 0, b = 1) {
  .Call(C_client_fill, keys, as.integer(n),
        match(what, c("uniform", "normal", "integer", "bits64")) - 1L,
        as.double(a), as.double(b))
}
client_fold <- function(key, data) .Call(C_client_fold, key, as.double(data))
client_errors <- function(key) .Call(C_client_errors, key)
