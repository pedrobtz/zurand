# Run expr in a forked child (parallel::mcparallel()) and collect its value,
# or NULL if it has not finished within timeout seconds -- in which case the
# child is killed, so a hang fails the calling test instead of the suite.
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
