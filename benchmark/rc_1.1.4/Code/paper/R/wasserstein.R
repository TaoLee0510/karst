# Wasserstein distance between two transport::wpp objects that cannot stall.
#
# transport's default solver (networkflow) can run without end on some
# inputs: on ABM_varyMisseg, 111 parent karyotypes against a KARST forecast on
# 13,692 did not finish in an hour, while revsimplex found the same optimum
# (0.30297) in 3 s (docs/fig2d_misseg_6ae2a24.md). Problems with at least 1e5
# point pairs therefore run networkflow in a forked child and switch to
# revsimplex after WASSERSTEIN_TIMEOUT seconds (default 60). Every other call,
# and every call that finishes in time, returns what transport::wasserstein()
# returns; an error in the child is raised again.
wasserstein_safe <- function(a, b, timeout = as.numeric(Sys.getenv("WASSERSTEIN_TIMEOUT", "60"))) {
  if (as.numeric(a$N) * as.numeric(b$N) < 1e5 || .Platform$OS.type != "unix") return(transport::wasserstein(a, b))
  job <- parallel::mcparallel(transport::wasserstein(a, b), silent = TRUE)
  res <- parallel::mccollect(job, wait = FALSE, timeout = timeout)
  if (!is.null(res)) {
    v <- res[[1]]
    if (inherits(v, "try-error")) stop(attr(v, "condition"))
    if (is.numeric(v)) return(v)
  } else {
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(job))
  }
  message("wasserstein_safe: networkflow did not finish in ", timeout, " s (", a$N, " x ", b$N, "); revsimplex")
  transport::wasserstein(a, b, method = "revsimplex")
}
