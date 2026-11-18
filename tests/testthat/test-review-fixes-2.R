# Regression tests for the fixes of the code review of 2026-11-18.

test_that("fit_karst() fits prepared data at the data's own missegregation rate", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_karst_data(sim$counts, dt = sim$dt, beta = 2 * sim$beta)
  control <- list(major_frequency = 0.05, max_cn = 4, history = sim$dt)
  fit <- fit_karst(data, sparse_joint_control = control)
  expect_equal(fit$data$beta, 2 * sim$beta)
  direct <- do.call(fit_sparse_joint, c(list(data), control))
  expect_equal(fit$summary, direct$summary)
  # a rate given to fit_karst() still replaces the data's
  given <- fit_karst(data, beta = sim$beta, sparse_joint_control = control)
  expect_equal(given$data$beta, sim$beta)
  expect_error(fit_karst(data, beta = sim$beta, sparse_joint_control = c(control, beta = sim$beta)), "not both")
})

test_that("karyotype labels are read one way only", {
  counts <- rbind("2.2.2" = c(10, 12), "2.03.2" = c(3, 1), "02.2.1" = c(0, 2))
  data <- prepare_karst_data(counts, dt = 1)
  expect_equal(rownames(data$counts), c("2.2.2", "2.3.2", "2.2.1"))
  for (bad in c("2.2.2.", " 2.2.2", "+2.2.2", "2..2.2", "2.2.-2", "2.2.2a")) {
    malformed <- counts
    rownames(malformed)[[1L]] <- bad
    expect_error(prepare_karst_data(malformed, dt = 1), "Malformed karyotype label", fixed = TRUE)
    expect_error(parse_karyotypes(bad), "Malformed karyotype label", fixed = TRUE)
  }
  same <- rbind("2.2.2" = c(10, 12), "2.02.2" = c(3, 1))
  expect_error(prepare_karst_data(same, dt = 1), "are the same karyotype")
  expect_equal(canonical_karyotype_labels(c("0.10.003", "00.0")), c("0.10.3", "0.0"))

  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                          history = sim$dt)
  major <- fit$summary$karyotype[[1L]]
  padded <- gsub("(^|\\.)([0-9])", "\\10\\2", major)
  prediction <- predict(fit, c(padded, major, major))
  expect_equal(prediction$karyotype, rep(major, 3L))
  expect_true(all(prediction$in_domain))
  expect_equal(prediction[2:3, ], prediction[c(1, 1), ], ignore_attr = TRUE)
  # a factor, and repeated karyotypes outside the domain
  expect_equal(predict(fit, factor(c(major, "1.1.1.1", "1.1.1.1"))),
               predict(fit, c(major, "1.1.1.1", "1.1.1.1")))
})

test_that("the missegregation rate of a fit must be positive", {
  sim <- simulate_sparse_joint_counts()
  expect_error(prepare_karst_data(sim$counts, dt = sim$dt, beta = 0), "positive missegregation rate")
  data <- prepare_karst_data(sim$counts, dt = sim$dt, beta = sim$beta)
  expect_error(fit_sparse_joint(data, beta = 0), "positive missegregation rate")
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = 0), "positive missegregation rate")
})

test_that("a free balance level is refused for the pairwise and state-space methods", {
  sim <- simulate_series_counts()
  for (method in c("anchored", "pairwise", "state_space")) {
    expect_error(fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                  series = method, balance_level = "free"),
                 "implemented for two time points only")
  }
})

test_that("a failed inner optimisation is an error, not a mode", {
  objective <- list(fn = function(par) NaN, env = list(last.par = c(g = 1)))
  expect_error(sparse_joint_inner_mode(objective, numeric(0)), class = "karst_inner_failure")
  expect_invisible(sparse_joint_inner_mode(list(fn = function(par) 3), numeric(0)))
})

test_that("a fit with no usable candidate history stops", {
  sim <- simulate_sparse_joint_counts()
  local_mocked_bindings(sparse_joint_choose_candidate = function(objective, convergence) NA_integer_)
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                                history = c(sim$dt, 300)),
               "No candidate history")
})

test_that("time points left out keep their names, their karyotypes and alfak's outputs", {
  sim <- simulate_series_counts()
  counts <- sim$counts
  # t1 has cells on neither the majors nor the domain: one far karyotype only.
  counts[, 2L] <- 0
  counts <- rbind(counts, "4.4.4.1" = c(0, 1, 0, 0))
  fit <- suppressWarnings(fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 5, history = 60))
  expect_equal(fit$diagnostics$dropped_samples$sample, 2L)
  expect_true(all(c("count_t0", "count_t2", "count_t3", "expected_t2") %in% names(fit$summary)))
  expect_false("count_t1" %in% names(fit$summary))
  expect_equal(colnames(fit$data$counts), c("t0", "t2", "t3"))
  # seen only at the time point left out: observed, at the observed level
  expect_true("4.4.4.1" %in% fit$dropped)
  expect_equal(fit$diagnostics$dropped_karyotypes, length(fit$dropped))
  far <- predict(fit, "4.4.4.1")
  expect_false(far$in_domain)
  expect_equal(log(far$fitness_mean / (2 * sparse_joint_self_probability(parse_karyotypes("4.4.4.1"), sim$beta) - 1)),
               unname(fit$trend$observed_level[["log_birth"]]), tolerance = 1e-10)
  out <- as_alfak(fit, minobs = 20, nboot = 3, xval = FALSE)
  expect_true(all(is.finite(out$landscape$mean)))

  # the pairwise and state-space methods report the same columns with and
  # without a time point left out
  for (method in c("pairwise", "state_space")) {
    whole <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                              series = method)
    gap <- suppressWarnings(fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 5, history = 60,
                                             series = method))
    expect_equal(nrow(whole$diagnostics$dropped_samples), 0L)
    expect_gt(nrow(gap$diagnostics$dropped_samples), 0L)
    expect_equal(vapply(whole$diagnostics$dropped_samples, class, character(1)),
                 vapply(gap$diagnostics$dropped_samples, class, character(1)))
  }
})

test_that("karyotypes seen only between the first and the last time point count as unobserved", {
  # thinly sampled, so that some neighbours in the domain are unseen
  sim <- simulate_series_counts(depth = 150)
  fit0 <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60)
  unseen <- fit0$summary$karyotype[!fit0$summary$major & sparse_joint_summary_total(fit0$summary) == 0]
  expect_gt(length(unseen), 0L)
  counts <- rbind(sim$counts, matrix(c(0, 2, 1, 0), 1L, dimnames = list(unseen[[1L]], NULL)))
  fit <- fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60)
  row <- match(unseen[[1L]], fit$summary$karyotype)
  s <- fit$summary
  calibration <- fit$trend$calibration
  expected_sd <- sqrt(s$relative_sd_model[row]^2 + (calibration[["deficit_unobserved"]] * s$relative_fitness[row])^2 +
                        (calibration[["floor"]] * fit$trend$birth_rate * (1 + s$shell[row]))^2)
  expect_equal(s$relative_sd[row], expected_sd, tolerance = 1e-12)
  used <- sparse_joint_problem(fit$data, fit$settings)$used
  expect_equal(fit$trend$observed_level, sparse_joint_observed_level(s, rowSums(used) > 0))
  expect_false(isTRUE(all.equal(fit$trend$observed_level, sparse_joint_observed_level(s))))
})

test_that("the efflux correction is alfak's: the x0-weighted mean at time 0 is g0", {
  sim <- simulate_sparse_joint_counts(beta = 2e-3)
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                          history = sim$dt)
  g0 <- log(1e7 / 1e5) / sim$dt
  plain <- as_alfak(fit, minobs = 20, nboot = 1)
  for (first_time in c(0, 40)) {
    one <- as_alfak(fit, minobs = 20, nboot = 1, correct_efflux = TRUE, first_time = first_time)
    boot <- one$bootstrap_res
    expect_equal(sum(boot$final_frequencies * boot$final_fitness), g0, tolerance = 1e-9)
    # (f - k) / v with one k for the landscape, the cross-validation and
    # landscape_data
    land <- one$landscape
    v <- 2 * sparse_joint_self_probability(parse_karyotypes(land$k), sim$beta) - 1
    relative <- predict(fit, land$k)$relative_fitness
    k <- relative - land$mean * v
    expect_lt(diff(range(k)), 1e-10)
    tested <- match(one$xval$xval_data$k, land$k)
    xval_relative <- plain$xval$xval_data$prediction - g0
    expect_equal(one$xval$xval_data$prediction, (xval_relative - k[[1L]]) / v[tested], tolerance = 1e-10)
    expect_equal(one$xval$xval_data$observation, land$mean[tested])
  }
  # landscape_data predicts the landscape, also for non-integer input (rounded with a warning)
  out <- as_alfak(fit, minobs = 20, nboot = 3, correct_efflux = TRUE, landscape_data_output = TRUE, xval = FALSE)
  model <- out$landscape_data$mean
  expect_equal(predict(model, out$landscape$k)[, 1L], out$landscape$mean, tolerance = 1e-10)
  v <- 2 * sparse_joint_self_probability(parse_karyotypes(out$landscape$k), sim$beta) - 1
  expect_equal(out$landscape$sd, predict(fit, out$landscape$k)$relative_sd / v, tolerance = 1e-10)
  cn <- parse_karyotypes(out$landscape$k[1:3])
  expect_warning(shifted <- predict(model, cn + 0.2), "rounded")
  expect_equal(shifted[, 1L], out$landscape$mean[1:3], tolerance = 1e-10)
  # without n0 and nb, the model's own net growth rates divided by v
  own <- as_alfak(fit, minobs = 20, nboot = 1, n0 = NULL, nb = NULL, correct_efflux = TRUE, xval = FALSE)
  v <- 2 * sparse_joint_self_probability(parse_karyotypes(own$landscape$k), sim$beta) - 1
  expect_equal(own$landscape$mean, predict(fit, own$landscape$k)$fitness_mean / v, tolerance = 1e-10)
})

test_that("cross-validation holds out each fold's prior and drift weight as the fit treats unseen states", {
  sim <- simulate_series_counts(depth = 150)
  calls <- list()
  rebuild <- sparse_joint_rebuild
  local_mocked_bindings(sparse_joint_rebuild = function(fit, weight = NULL, zeta_prior_mean = NULL,
                                                        zeta_prior_sd = NULL, design = NULL, posterior = TRUE,
                                                        drift_weight = NULL) {
    if (!posterior) calls[[length(calls) + 1L]] <<- list(weight = weight, mean = zeta_prior_mean,
                                                         sd = zeta_prior_sd, drift = drift_weight)
    rebuild(fit, weight, zeta_prior_mean, zeta_prior_sd, design, posterior, drift_weight)
  })
  for (method in c("pairwise", "state_space")) {
    fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            series = method)
    problem <- sparse_joint_problem(fit$data, fit$settings)
    calls <- list()
    out <- as_alfak(fit, minobs = 20, nboot = 3)
    expect_true(is.finite(out$xval$R2R))
    expect_gt(length(calls), 0L)
    depth <- colSums(problem$tmb_data[["y"]])
    for (call in calls) {
      held <- call$weight == 0
      expect_true(any(held))
      mean <- as.matrix(call$mean)
      expected <- as.matrix(problem$tmb_data$zeta_prior_mean)
      if (method == "pairwise") {
        centre <- log(colSums(problem$tmb_data$y_start)) + fit$settings$initial_frequency[[1L]]
        expected[held, ] <- matrix(centre, sum(held), ncol(expected), byrow = TRUE)
        expect_null(call$drift)
      } else {
        expected[held, ] <- log(depth[[1L]]) + fit$settings$initial_frequency[[1L]]
        drift <- problem$tmb_data$drift_weight
        drift[held, ] <- matrix(depth[-length(depth)] + depth[-1L], sum(held), ncol(drift), byrow = TRUE)
        expect_equal(call$drift, drift)
      }
      expect_equal(mean, expected)
      expect_true(all(as.matrix(call$sd)[held, ] == fit$settings$initial_frequency[[2L]]))
    }
  }
})

test_that("the simulators keep their parameters and their samples", {
  land <- simulate_grf_landscape(n_chr = 3, n_centroids = 6, seed = 1, max_cn = 4)
  k <- matrix(c(2, 2, 2, 3, 2, 1), nrow = 2, byrow = TRUE)
  first <- predict_landscape_fitness(land, k)
  changed <- land
  changed$lambda <- 2 * land$lambda
  expect_equal(predict_landscape_fitness(changed, k), predict_landscape_fitness(changed, k, use_cache = FALSE))
  expect_false(isTRUE(all.equal(predict_landscape_fitness(changed, k), first)))
  expect_equal(predict_landscape_fitness(land, k), first)
  # one copy of one chromosome: only that copy can missegregate
  set.seed(1)
  daughters <- grf_missegregate_daughters(c(0L, 0L, 3L), beta = 0.5)
  expect_equal(daughters[[1L]][1:2], c(0L, 0L))
  set.seed(2)
  daughters <- grf_missegregate_daughters(c(0L, 1L, 0L), beta = 0.9)
  expect_equal(sort(c(daughters[[1L]][[2L]], daughters[[2L]][[2L]])), c(0L, 2L))
  one <- simulate_grf_landscape(n_chr = 2, n_centroids = 3, centroid_min = 5, centroid_max = 5)
  expect_true(all(one$centroids == 5L))
  # no karyotype without a cell, also with detection_threshold = 0
  sim <- simulate_sparse_counts(land, n0 = 20, n1 = 20, seed = 3, detection_threshold = 0, beta = 0.01)
  expect_true(all(rowSums(sim$counts) > 0))
  set.seed(4)
  expect_error(grf_multinomial_draw(rep(0.1, 10), 100, overdispersion = 1e-300), "too small")
})
