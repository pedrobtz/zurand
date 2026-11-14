# GNU OpenMP's thread pool does not survive fork(): a forked child that
# enters a parallel region after the parent has used one waits for threads
# that no longer exist (#122). parallel::mclapply() is exactly that, so
# zurand runs its fills on one thread in a forked child. The child's work
# runs under mcparallel() and is collected with a timeout, so a regression
# fails this test instead of hanging the suite.

fork_collect <- function(expr, timeout = 60) {
  job <- parallel::mcparallel(expr)
  out <- parallel::mccollect(job, wait = FALSE, timeout = timeout)
  if (is.null(out)) {
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = TRUE)
    return(NULL)
  }
  out[[1]]
}

test_that("fills in a forked child finish after the parent has run threaded", {
  skip_on_os("windows")
  skip_if(rng_threads() < 2L, "needs OpenMP and at least two threads")
  key <- rng_key(3L)
  n <- 70001L                          # above the OpenMP threshold
  lz <- rng_lazy_normal(key, 1e5)
  ref <- list(
    normal = rng_normal(key, n),       # the parent enters a parallel region
    integer = rng_integer(key, n, 1L, 10L),
    region = lz[1:40000],
    lazy = lz + 0                      # materialised, threaded
  )

  got <- fork_collect(list(
    threads = rng_threads(),
    normal = rng_normal(key, n),
    integer = rng_integer(key, n, 1L, 10L),
    region = lz[1:40000],
    lazy = lz + 0
  ))

  expect_false(is.null(got), label = "the forked child finished")
  expect_identical(got$threads, 1L)
  expect_identical(got[names(ref)], ref)
})

test_that("a forked child can return a lazy vector the parent materialised", {
  skip_on_os("windows")
  skip_if(rng_threads() < 2L, "needs OpenMP and at least two threads")
  x <- rng_lazy_normal(rng_key(5L), 1e5)
  s <- sd(x)                           # materialises x with threads
  # Returning x serialises it, which compares it with a threaded refill.
  got <- fork_collect(x)
  expect_false(is.null(got), label = "the forked child finished")
  expect_identical(got[], rng_normal(rng_key(5L), 1e5))
})
