test_that("every fit records the objective at the posterior mode it reports", {
  sim <- simulate_sparse_joint_counts()
  one <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = 300))
  value <- one$diagnostics$posterior_objective
  expect_true(is.finite(value))
  # the same value from a search at the fitted scales started at the stored
  # mode (the rebuilt objective; the rebuild itself takes the stored mode and
  # reports this value)
  rebuilt <- sparse_joint_rebuild(one)
  expect_equal(rebuilt$solver$value, value)
  expect_equal(as.numeric(rebuilt$objective$fn(numeric(0))), value, tolerance = 1e-6)
})

test_that("candidate histories are compared on the objective at the posterior mode", {
  sim <- simulate_sparse_joint_counts()
  history <- c(sim$dt, 300)
  fit_at <- function(h) suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta,
                                                          major_frequency = 0.05, max_cn = 4, history = h))
  plain <- fit_at(history)
  profile <- plain$diagnostics$history_profile
  expect_equal(names(profile)[1:4], c("history", "objective", "posterior_objective", "convergence"))
  expect_equal(plain$settings$history, history[[which.min(profile$posterior_objective)]])
  for (i in seq_along(history)) {
    alone <- fit_at(history[[i]])
    expect_equal(profile$objective[[i]], alone$diagnostics$objective)
    expect_equal(profile$posterior_objective[[i]], alone$diagnostics$posterior_objective)
  }
  # The candidate nlminb's objective prefers is reported at a worse mode: the
  # other one is kept, and nlminb's objectives are unchanged.
  best <- which.min(profile$objective)
  real <- sparse_joint_inner_mode
  calls <- 0L
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    value <- real(objective, par_fixed, start)
    if (calls == best) structure(as.numeric(value) + 1000, restart = attr(value, "restart")) else value
  })
  worse <- fit_at(history)
  expect_equal(calls, length(history))
  moved <- worse$diagnostics$history_profile
  expect_equal(moved$objective, profile$objective)
  expect_equal(moved$posterior_objective[[best]], profile$posterior_objective[[best]] + 1000)
  expect_equal(worse$settings$history, history[[-best]])
})

test_that("a candidate whose scale optimisation gave no finite objective is not compared", {
  sim <- simulate_sparse_joint_counts()
  history <- c(sim$dt, 300)
  fit_at <- function(h) suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta,
                                                          major_frequency = 0.05, max_cn = 4, history = h))
  best <- which.min(fit_at(history)$diagnostics$history_profile$posterior_objective)
  # The candidate the posterior-mode objective prefers has no finite nlminb
  # objective; its posterior mode is still found.
  real <- sparse_joint_optimise_scales
  calls <- 0L
  local_mocked_bindings(sparse_joint_optimise_scales = function(...) {
    calls <<- calls + 1L
    optimum <- real(...)
    if (calls == best) optimum$objective <- NaN
    optimum
  })
  fit <- fit_at(history)
  expect_equal(calls, length(history))
  profile <- fit$diagnostics$history_profile
  expect_true(is.nan(profile$objective[[best]]) && is.finite(profile$posterior_objective[[best]]))
  expect_equal(fit$settings$history, history[[-best]])
})

test_that("the trends of a cohort are compared on the objective at the posterior mode", {
  sim <- simulate_cohort_counts()
  fit <- function() suppressWarnings(fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta,
                                                             major_frequency = 0.05, max_cn = 4, trend = "select"))
  plain <- fit()
  posterior <- plain$diagnostics$trend_posterior_objective
  expect_equal(names(posterior), c("shared", "sample"))
  best <- which.min(posterior)
  expect_equal(plain$diagnostics$trend, names(posterior)[[best]])
  # the trend preferred is reported at a worse mode: the other is kept, and
  # nlminb's objectives are unchanged (one posterior per trend, shared first)
  real_mode <- sparse_joint_inner_mode
  calls <- 0L
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    value <- real_mode(objective, par_fixed, start)
    if (calls == best) structure(as.numeric(value) + 1000, restart = attr(value, "restart")) else value
  })
  worse <- fit()
  expect_equal(calls, 2L)
  expect_equal(worse$diagnostics$trend_objective, plain$diagnostics$trend_objective)
  expect_equal(worse$diagnostics$trend_posterior_objective[[best]], posterior[[best]] + 1000)
  expect_equal(worse$diagnostics$trend, names(posterior)[[-best]])
  # the trend preferred has no finite nlminb objective: not compared
  local_mocked_bindings(sparse_joint_inner_mode = real_mode)
  real_scales <- sparse_joint_optimise_scales
  calls <- 0L
  local_mocked_bindings(sparse_joint_optimise_scales = function(...) {
    calls <<- calls + 1L
    optimum <- real_scales(...)
    if (calls == best) optimum$objective <- NaN
    optimum
  })
  unusable <- fit()
  expect_equal(calls, 2L)
  expect_true(is.nan(unusable$diagnostics$trend_objective[[best]]))
  expect_equal(unusable$diagnostics$trend, names(posterior)[[-best]])
})
