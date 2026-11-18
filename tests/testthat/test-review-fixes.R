# Regression tests for the fixes of the code review of 2026-11-17.

test_that("an explicit missegregation rate replaces the one stored in prepared data", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_karst_data(sim$counts, dt = sim$dt, beta = sim$beta)
  # without `beta`, the data's own rate is fitted, through the candidate
  # histories too, and gives the fit of the count matrix at that rate
  own <- fit_sparse_joint(data, major_frequency = 0.05, max_cn = 4)
  expect_equal(own$data$beta, sim$beta)
  direct <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4)
  expect_equal(own$summary, direct$summary)
  given <- fit_sparse_joint(data, beta = 2 * sim$beta, major_frequency = 0.05, max_cn = 4, history = sim$dt)
  expect_equal(given$data$beta, 2 * sim$beta)
  # a range of rates refits prepared data at each rate
  rates <- fit_sparse_joint_rates(data, major_frequency = 0.05, max_cn = 4, history = sim$dt,
                                  beta_range = c(sim$beta / 2, 2 * sim$beta))
  expect_equal(vapply(rates$fits, function(f) f$data$beta, numeric(1)), rates$rates)
})

test_that("the two-layer pipeline refuses more than two time points", {
  sim <- simulate_series_counts()
  data <- prepare_karst_data(sim$counts, beta = sim$beta, times = sim$times)
  expect_error(build_karyotype_graph(data, shell_depth = 1, min_cn = 1, max_cn = 4), "two time points")
  expect_error(fit_local_posterior(data, graph = NULL), "two time points")
})

test_that("every status of a local fit has a default covariance inflation", {
  status <- c("TMB_sdreport", "untrusted_gradient", "untrusted_inner_gradient", "untrusted_nonconverged",
              "untrusted_sdreport_missing", "untrusted_scale_boundary", "untrusted_hessian",
              "untrusted_sdreport_nonfinite", "untrusted_lognormal_moments", NA)
  multiplier <- anchor_covariance_multiplier(status)
  expect_length(multiplier, length(status))
  expect_true(all(is.finite(multiplier) & multiplier >= 1))
})

test_that("cross-validation of a pairwise fit centres each pair's unseen karyotypes as the fit does", {
  # thinly sampled, so that some neighbours in the pairs' domains are unseen
  sim <- simulate_series_counts(depth = 150)
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                          series = "pairwise")
  problem <- sparse_joint_problem(fit$data, fit$settings)
  y_start <- problem$tmb_data$y_start
  y_end <- problem$tmb_data$y_end
  centre <- log(colSums(y_start)) + fit$settings$initial_frequency[[1L]]
  unseen_any <- FALSE
  for (k in seq_len(ncol(y_start))) {
    unseen <- problem$tmb_data$member[, k] > 0 & y_start[, k] == 0 & y_end[, k] == 0
    if (!any(unseen)) next
    unseen_any <- TRUE
    expect_equal(unname(problem$tmb_data$zeta_prior_mean[unseen, k]), rep(unname(centre[[k]]), sum(unseen)))
  }
  expect_true(unseen_any)
  expect_true(is.finite(as_alfak(fit, minobs = 20, nboot = 5, n0 = 1e5, nb = 1e7)$xval$R2R))
})
