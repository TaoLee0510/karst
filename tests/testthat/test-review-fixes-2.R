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

# A stand-in for a TMB objective: `values` are what successive calls of fn()
# return; each call records the inner start it would use.
mock_inner_objective <- function(values) {
  env <- new.env()
  env$par <- c(scale = 0, g = 0)
  env$random <- 2L
  env$last.par.best <- c(scale = 0, g = 5)
  env$random.start <- expression(last.par.best[random])
  env$starts <- numeric(0)
  fn <- function(par) {
    env$starts <- c(env$starts, unname(eval(env$random.start, env)))
    values[[min(length(env$starts), length(values))]]
  }
  list(fn = fn, env = env)
}

test_that("the posterior mode is searched once more from the initial values before it fails", {
  # found at once: one search, from TMB's best mode
  first <- mock_inner_objective(3)
  expect_invisible(sparse_joint_inner_mode(first, numeric(0)))
  expect_equal(first$env$starts, 5)
  # failed from TMB's best mode, found from the initial values
  second <- mock_inner_objective(c(NaN, 4))
  found <- sparse_joint_inner_mode(second, numeric(0))
  expect_equal(as.numeric(found), 4)
  expect_true(isTRUE(attr(found, "restart")))
  expect_equal(second$env$starts, c(5, 0))
  expect_identical(second$env$random.start, expression(last.par.best[random]))
  # failed from both: an error, and TMB's own start restored
  failed <- mock_inner_objective(c(NaN, Inf))
  expect_error(sparse_joint_inner_mode(failed, numeric(0)), class = "karst_inner_failure")
  expect_equal(failed$env$starts, c(5, 0))
  expect_identical(failed$env$random.start, expression(last.par.best[random]))
})

test_that("a candidate history whose posterior mode cannot be found is not chosen", {
  sim <- simulate_sparse_joint_counts()
  real <- sparse_joint_inner_mode
  calls <- 0L
  # The first candidate's posterior fails, the second's is found.
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    if (calls == 1L) stop(errorCondition("mocked", class = "karst_inner_failure"))
    real(objective, par_fixed, start)
  })
  fit <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = c(sim$dt, 300)))
  profile <- fit$diagnostics$history_profile
  expect_equal(profile$inner_failure, c(TRUE, FALSE))
  expect_true(is.na(profile$objective[[1L]]) && is.na(profile$posterior_objective[[1L]]) &&
                is.na(profile$convergence[[1L]]))
  expect_true(is.finite(profile$objective[[2L]]) && is.finite(profile$posterior_objective[[2L]]))
  expect_equal(fit$settings$history, 300)
  alone <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                             max_cn = 4, history = 300))
  expect_equal(fit$summary, alone$summary)
  # Every candidate failed: the fit stops, with the class of the failure.
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    stop(errorCondition("mocked", class = "karst_inner_failure"))
  })
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                                history = c(sim$dt, 300)),
               "posterior mode at the fitted scales could not be found at any", class = "karst_inner_failure")
  # One history: no candidate to fall back on.
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                                history = sim$dt),
               class = "karst_inner_failure")
})

test_that("the second search uses TMB's real objective and the problem's own initial values", {
  sim <- simulate_sparse_joint_counts()
  fit <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = sim$dt))
  problem <- sparse_joint_problem(fit$data, fit$settings)
  make <- function() TMB::MakeADFun(problem$tmb_data, problem$parameters, map = problem$map,
                                    random = problem$random, DLL = "karst", silent = TRUE,
                                    inner.control = fit$settings$inner_control)
  # In a fit: TMB's best mode so far far off, the search from it fails.
  objective <- make()
  env <- objective$env
  scales <- env$par[-env$random]
  reference <- make()$fn(scales)   # a fresh objective starts at the initial values
  bad <- env$last.par.best
  bad[env$random] <- 1e4
  assign("last.par.best", bad, envir = env)
  tmb_start <- env$random.start
  expect_false(is.finite(objective$fn(scales)))
  found <- sparse_joint_inner_mode(objective, scales)
  expect_true(isTRUE(attr(found, "restart")))
  expect_equal(as.numeric(found), as.numeric(reference))
  expect_identical(env$random.start, tmb_start)
  # In a rebuild TMB starts from the stored mode (the objective is taped
  # there); the rebuild passes on the problem's own initial values for the
  # second search.
  real <- sparse_joint_inner_mode
  passed <- NULL
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    passed <<- start
    real(objective, par_fixed, start)
  })
  rebuilt <- sparse_joint_rebuild(fit)
  map <- problem$map
  map$log_scale <- factor(rep(NA, length(problem$parameters$log_scale)))
  env <- rebuilt$objective$env
  expect_identical(passed, unname(sparse_joint_latent_vector(problem$parameters, map,
                                                             unique(names(env$par)[env$random]))))
  expect_identical(length(passed), length(env$random))
  expect_null(rebuilt$solver$restart)
  # ... and when the search from the stored mode fails, the second finds the
  # mode that a search from the initial values finds.
  bad <- env$last.par.best
  bad[env$random] <- 1e4
  assign("last.par.best", bad, envir = env)
  found <- real(rebuilt$objective, numeric(0), passed)
  expect_true(isTRUE(attr(found, "restart")))
  initial <- problem$parameters
  initial$log_scale <- fit$posterior$log_scale
  fresh <- TMB::MakeADFun(problem$tmb_data, initial, map = map, random = problem$random, DLL = "karst",
                          silent = TRUE, inner.control = fit$settings$inner_control)
  expect_equal(as.numeric(found), as.numeric(fresh$fn(numeric(0))))
  # posterior = FALSE passes the same start
  passed <- NULL
  sparse_joint_rebuild(fit, posterior = FALSE)
  expect_identical(passed, unname(sparse_joint_latent_vector(problem$parameters, map,
                                                             unique(names(env$par)[env$random]))))
})

test_that("the initial latent values are ordered as in TMB's parameter vector", {
  check <- function(fit) {
    problem <- sparse_joint_problem(fit$data, fit$settings)
    objective <- TMB::MakeADFun(problem$tmb_data, problem$parameters, map = problem$map,
                                random = problem$random, DLL = "karst", silent = TRUE)
    latent <- unique(names(objective$env$par)[objective$env$random])
    expect_setequal(latent, problem$random)
    expect_identical(sparse_joint_latent_vector(problem$parameters, problem$map, latent),
                     unname(objective$env$par[objective$env$random]))
  }
  two <- simulate_sparse_joint_counts()
  check(suppressWarnings(fit_sparse_joint(two$counts, dt = two$dt, beta = two$beta, major_frequency = 0.05,
                                          max_cn = 4, history = two$dt)))
  sim <- simulate_series_counts()
  for (method in c("anchored", "pairwise", "state_space")) {
    check(suppressWarnings(fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4,
                                            history = 60, series = method)))
  }
  # in the order given; a mapped parameter: one value per level, entries
  # mapped to NA left out
  expect_equal(sparse_joint_latent_vector(list(a = 1:2, z = c(1, 2, 3, 4), s = 9),
                                          list(z = factor(c(2, NA, 1, NA))), c("z", "a")),
               c(3, 1, 1, 2))
})

test_that("candidates that failed beside ones whose optimisation failed, and the restart record", {
  sim <- simulate_sparse_joint_counts()
  real <- sparse_joint_inner_mode
  calls <- 0L
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    if (calls == 1L) stop(errorCondition("mocked", class = "karst_inner_failure"))
    real(objective, par_fixed, start)
  })
  # the free level of the balance prior: NA for the failed candidate
  free <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                            max_cn = 4, history = c(sim$dt, 300), balance_level = "free"))
  expect_true(is.na(free$diagnostics$history_profile$balance_level[[1L]]))
  expect_true(is.finite(free$diagnostics$history_profile$balance_level[[2L]]))
  # one failed, the other's optimisation failed: an error, not of the class
  calls <- 0L
  local_mocked_bindings(sparse_joint_choose_candidate = function(objective, convergence) NA_integer_)
  err <- tryCatch(suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta,
                                                    major_frequency = 0.05, max_cn = 4,
                                                    history = c(sim$dt, 300))),
                  error = function(e) e)
  expect_match(conditionMessage(err), "could not be found at 1, and the scale optimisation failed at the others")
  expect_false(inherits(err, "karst_inner_failure"))
})

test_that("a posterior mode found by the second search is recorded", {
  sim <- simulate_sparse_joint_counts()
  real <- sparse_joint_inner_mode
  calls <- 0L
  # the first candidate's mode is reported as found by the second search
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    value <- real(objective, par_fixed, start)
    if (calls == 1L) structure(as.numeric(value), restart = TRUE) else value
  })
  fit <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = c(sim$dt, 300)))
  restart <- fit$diagnostics$history_profile$posterior_restart
  expect_true(is.finite(restart[[1L]]) && is.na(restart[[2L]]))
  expect_null(fit$diagnostics$history_profile$inner_failure)
  calls <- 0L
  one <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = sim$dt))
  expect_true(is.finite(one$diagnostics$posterior_restart))
  expect_null(one$diagnostics$history_profile)
})

test_that("the first fits of the intermediate weight and of initial_state = 'auto' fall back on failure", {
  real <- sparse_joint_inner_mode
  calls <- 0L
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    calls <<- calls + 1L
    if (calls == 1L) stop(errorCondition("mocked", class = "karst_inner_failure"))
    real(objective, par_fixed, start)
  })
  # the first fit (at the shortest history) sets the weight of the time
  # points in between: failed, the weight is 1
  sim <- simulate_series_counts()
  fit <- suppressWarnings(fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4,
                                           history = c(60, 300)))
  expect_equal(fit$diagnostics$intermediate_weight, 1)
  expect_null(fit$diagnostics$intermediate_dispersion)
  expect_null(fit$diagnostics$history_profile$inner_failure)
  # the check fit of initial_state = "auto": failed, the balance is kept
  calls <- 0L
  two <- simulate_sparse_joint_counts()
  auto <- suppressWarnings(fit_sparse_joint(two$counts, dt = two$dt, beta = two$beta, major_frequency = 0.05,
                                            max_cn = 4, history = two$dt, initial_state = "auto"))
  check <- auto$diagnostics$initial_state_check
  expect_equal(check$chosen, "balance")
  expect_true(is.na(check$level) && is.na(check$level_sd) && is.na(check$check_convergence))
  expect_equal(auto$diagnostics$initial_state, "balance")
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
