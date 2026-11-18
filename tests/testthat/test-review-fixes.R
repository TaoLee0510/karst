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

# The prior that the pairwise cross-validation gives a held-out fold is
# checked on the calls the cross-validation makes, in test-review-fixes-2.R
# ("cross-validation holds out each fold's prior ..."); the check here only
# compared the fit's own prior with its formula.
