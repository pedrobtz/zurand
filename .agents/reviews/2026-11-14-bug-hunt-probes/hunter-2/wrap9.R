source(if (file.exists("/w/hunt-2/common.R")) "/w/hunt-2/common.R" else "common.R")
key <- rng_key(3L)
for (mk in list(function() rnorm(1000), function() rng_lazy_normal(key, 1000))) {
 x <- mk(); y <- x; attr(y, "u") <- 1; .Internal(inspect(y, 1))
 x <- mk(); y <- x; names(y) <- NULL; class(y) <- "foo"; .Internal(inspect(y, 1))
 x <- mk(); y <- x; dim(y) <- c(10, 100); .Internal(inspect(y, 1))
 x <- mk(); y <- structure(x, class = "bar"); .Internal(inspect(y, 1))
 x <- mk(); y <- x; levels(y) <- "a"; .Internal(inspect(y, 1))
 x <- mk(); y <- x; attributes(y) <- list(a = 1); .Internal(inspect(y, 1))
 x <- mk(); l <- list(x = x); l$x[2] <- 0; .Internal(inspect(l$x, 1))
}
