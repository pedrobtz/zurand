client_fill <- function(keys, n, what, a = 0, b = 1) {
  .Call(C_client_fill, keys, as.integer(n),
        match(what, c("uniform", "normal", "integer", "bits64",
                    "normal_ziggurat", "normal_mcfarland")) - 1L,
        as.double(a), as.double(b))
}
client_fold <- function(key, data) .Call(C_client_fold, key, as.double(data))
client_errors <- function(key) .Call(C_client_errors, key)

dist_code <- function(what) match(what, c("uniform", "normal")) - 1L
method_code <- function(method) if (method == "mcfarland") 1L else 0L
client_fill_at <- function(keys, n, what, a, b, offset, method = "ziggurat") {
  .Call(C_client_fill_at, keys, as.integer(n), dist_code(what), as.double(a),
        as.double(b), as.double(offset), method_code(method))
}
client_stream <- function(key, total, chunk, what, a, b, method = "ziggurat",
                          stop_after = 0L) {
  .Call(C_client_stream, key, as.integer(total), as.integer(chunk),
        dist_code(what), as.double(a), as.double(b), method_code(method),
        as.integer(stop_after))
}
client_stream_par <- function(keys, total, chunk, what, a, b,
                              method = "ziggurat") {
  .Call(C_client_stream_par, keys, as.integer(total), as.integer(chunk),
        dist_code(what), as.double(a), as.double(b), method_code(method))
}
client_threads <- function() .Call(C_client_threads)
