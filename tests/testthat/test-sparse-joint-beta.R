# The missegregation rate must reach the likelihood through every entry point,
# from count matrices and from prepared data (`karst_data`) alike. Twice it
# was held fixed without an error:
# - alfak2 up to b9a275e ignored a rate given with prepared data (fixed in
#   e7e32f8);
# - karst e7e32f8 to cfabe12 fitted prepared data passed to fit_karst() at
#   the default rate instead of the data's own (fixed in 42d13a9).
# Each path below is compared with the fit of the count matrix at the rate it
# must use, and the fits at the two rates must differ, so a rate held fixed
# anywhere fails here. The tests in test-review-fixes*.R only read the stored
# rate (`data$beta`), not the fit.

beta_low <- 5e-4
beta_high <- 4 * beta_low

# The fits at the two rates differ in the likelihood and in fitness (about 3
# to 16 in the objective and 0.28 to 0.44 in the neighbours' fitness here).
expect_rates_matter <- function(low, high, table = "summary") {
  expect_gt(abs(low$diagnostics$objective - high$diagnostics$objective), 1)
  expect_identical(low[[table]]$karyotype, high[[table]]$karyotype)
  expect_gt(max(abs(low[[table]]$relative_fitness - high[[table]]$relative_fitness)), 0.05)
}

# The same fit as the reference: the same likelihood and the same tables.
expect_same_fit <- function(fit, reference, tables = "summary") {
  expect_equal(fit$diagnostics$objective, reference$diagnostics$objective)
  for (table in tables) expect_equal(fit[[table]], reference[[table]])
}

test_that("every two-time-point entry point fits at the rate it is given", {
  sim <- simulate_sparse_joint_counts()
  control <- list(major_frequency = 0.05, max_cn = 4, history = sim$dt)
  reference <- lapply(c(low = beta_low, high = beta_high), function(rate) {
    do.call(fit_sparse_joint, c(list(sim$counts, dt = sim$dt, beta = rate), control))
  })
  expect_rates_matter(reference$low, reference$high)
  prepared_low <- prepare_karst_data(sim$counts, dt = sim$dt, beta = beta_low)
  prepared_high <- prepare_karst_data(sim$counts, dt = sim$dt, beta = beta_high)

  # fit_sparse_joint(): the prepared data's own rate, or the one given
  own <- do.call(fit_sparse_joint, c(list(prepared_high), control))
  expect_equal(own$data$beta, beta_high)
  expect_same_fit(own, reference$high)
  expect_same_fit(do.call(fit_sparse_joint, c(list(prepared_low, beta = beta_high), control)), reference$high)
  expect_same_fit(do.call(fit_sparse_joint, c(list(prepared_high, beta = beta_low), control)), reference$low)

  # fit_karst(): from a matrix, from prepared data without a rate (its own),
  # and from prepared data with a rate (that one)
  expect_same_fit(fit_karst(sim$counts, dt = sim$dt, beta = beta_high, sparse_joint_control = control),
                  reference$high)
  expect_same_fit(fit_karst(prepared_high, sparse_joint_control = control), reference$high)
  expect_same_fit(fit_karst(prepared_low, sparse_joint_control = control), reference$low)
  expect_same_fit(fit_karst(prepared_low, beta = beta_high, sparse_joint_control = control), reference$high)

  # fit_sparse_joint_rates(): each rate of the range, from a matrix and from
  # prepared data stored at the lowest rate
  for (input in list(list(sim$counts, dt = sim$dt), list(prepared_low))) {
    rates <- do.call(fit_sparse_joint_rates, c(input, control, list(beta_range = c(beta_low, beta_high), n_rates = 2L)))
    expect_same_fit(rates$fits[[1L]], reference$low)
    expect_same_fit(rates$fits[[2L]], reference$high)
  }
})

test_that("every method for more than two time points fits at the rate it is given", {
  sim <- simulate_series_counts()
  prepared_low <- prepare_karst_data(sim$counts, beta = beta_low, times = sim$times)
  prepared_high <- prepare_karst_data(sim$counts, beta = beta_high, times = sim$times)
  for (method in c("anchored", "pairwise", "state_space")) {
    fit <- function(counts, ...) {
      suppressWarnings(fit_sparse_joint(counts, ..., max_cn = 4, history = 60, series = method))
    }
    reference <- list(low = fit(sim$counts, times = sim$times, beta = beta_low),
                      high = fit(sim$counts, times = sim$times, beta = beta_high))
    expect_rates_matter(reference$low, reference$high)
    expect_same_fit(fit(prepared_high), reference$high)
    expect_same_fit(fit(prepared_low, beta = beta_high), reference$high)
  }
})

test_that("a cohort fits at the rate it is given, or at its prepared datasets' own", {
  sim <- simulate_cohort_counts()
  tables <- c("landscape", "samples")
  fit <- function(samples, ...) fit_sparse_joint_cohort(samples, ..., major_frequency = 0.05, max_cn = 4)
  reference <- list(low = fit(sim$samples, dt = sim$dt, beta = beta_low),
                    high = fit(sim$samples, dt = sim$dt, beta = beta_high))
  expect_rates_matter(reference$low, reference$high, table = "landscape")
  prepare_at <- function(rate) Map(function(x, dt) prepare_karst_data(x, dt = dt, beta = rate), sim$samples, sim$dt)
  prepared_low <- prepare_at(beta_low)
  prepared_high <- prepare_at(beta_high)

  expect_same_fit(fit(prepared_high), reference$high, tables)
  expect_same_fit(fit(prepared_low, beta = beta_high), reference$high, tables)
  # a count matrix among prepared datasets takes their rate
  mixed <- list(A = sim$samples$A, B = prepared_high$B)
  expect_same_fit(fit(mixed, dt = sim$dt), reference$high, tables)

  rates <- fit_sparse_joint_rates(prepared_low, major_frequency = 0.05, max_cn = 4, beta_range = c(beta_low, beta_high),
                                  n_rates = 2L, fit_function = fit_sparse_joint_cohort)
  expect_same_fit(rates$fits[[1L]], reference$low, tables)
  expect_same_fit(rates$fits[[2L]], reference$high, tables)
})

test_that("alfak_compat() fits at its `pm`", {
  sim <- simulate_sparse_joint_counts()
  x <- sim$counts
  colnames(x) <- c("0", "1")
  fitted_at <- function(pm) {
    dir <- tempfile("karst_alfak_")
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    set.seed(5)
    alfak_compat(list(x = x, dt = sim$dt), dir, minobs = 20, nboot = 5, pm = pm, max_cn = 4,
                 major_frequency = 0.05, history = sim$dt)
    readRDS(file.path(dir, "karst_fit.Rds"))
  }
  for (rate in c(beta_low, beta_high)) {
    reference <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = rate, major_frequency = 0.05, max_cn = 4,
                                  history = sim$dt)
    fit <- fitted_at(rate)
    expect_equal(fit$data$beta, rate)
    expect_same_fit(fit, reference)
  }
})
