test_that("a state's own terms reproduce the fit's objective for every method", {
  sim <- simulate_sparse_joint_counts()
  two <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  expect_lt(sparse_joint_own_check(two), 1e-8)
  series <- simulate_series_counts()
  for (method in c("anchored", "pairwise", "state_space")) {
    fit <- fit_sparse_joint(series$counts, times = series$times, beta = series$beta, max_cn = 4, history = 60,
                            series = method)
    expect_lt(sparse_joint_own_check(fit), 1e-8)
  }
})

test_that("own-count fitness uses a karyotype's own counts only", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  labels <- fit$summary$karyotype
  outside <- "4.4.4.4"
  own <- sparse_joint_own_fitness(fit, c(labels, outside), resample = FALSE)
  cells <- attr(own, "cells")
  expect_true(is.na(own[[outside]]))
  expect_true(all(is.na(own[labels][cells[seq_along(labels)] == 0])))
  expect_true(all(is.finite(own[labels][cells[seq_along(labels)] > 0])))
  # The three clones the data are drawn from: their own counts give their
  # net rates' order, without the landscape prior.
  major <- labels[1:3]
  expect_equal(order(own[major]), order(sim$net[match(major, sim$labels)]))
  # A clone's estimate does not depend on which others are asked for, and the
  # resampling is reproducible and leaves the global random state alone.
  expect_equal(as.numeric(sparse_joint_own_fitness(fit, major[2], resample = FALSE)), as.numeric(own[major[2]]))
  set.seed(5)
  before <- stats::runif(1)
  set.seed(5)
  a <- sparse_joint_own_fitness(fit, labels, seed = 3)
  after <- stats::runif(1)
  expect_equal(before, after)
  b <- sparse_joint_own_fitness(fit, rev(labels), seed = 3)
  expect_equal(a[labels], b[labels])
  expect_false(isTRUE(all.equal(a[labels], own[labels])))
})

test_that("as_alfak compares held-out predictions with own-count targets on request", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  set.seed(3)
  model <- as_alfak(fit, minobs = 20, nboot = 5)
  set.seed(3)
  default <- as_alfak(fit, minobs = 20, nboot = 5, xval_target = "model")
  expect_identical(model, default)
  expect_null(model$xval$target)
  set.seed(3)
  own <- as_alfak(fit, minobs = 20, nboot = 5, xval_target = "own")
  expect_identical(own$landscape, model$landscape)
  expect_identical(own$bootstrap_res, model$bootstrap_res)
  expect_equal(own$xval$target, "own")
  # The same held-out predictions, compared with other targets, on the
  # karyotypes with cells.
  k <- own$xval$xval_data$k
  expect_true(all(k %in% model$xval$xval_data$k))
  m <- match(k, model$xval$xval_data$k)
  expect_equal(own$xval$xval_data$prediction, model$xval$xval_data$prediction[m])
  cells <- rowSums(fit$data$counts)[k]
  expect_true(all(cells > 0))
  targets <- sparse_joint_own_fitness(fit, k, resample = TRUE, seed = 1L)
  divisor <- rep(1, length(k))
  level <- log(1e7 / 1e5) / sim$dt
  expect_equal(own$xval$xval_data$observation,
               unname((as.numeric(targets) - fit$trend$reference$fitness + level) / divisor), tolerance = 1e-12)
  expect_true(is.finite(own$xval$R2R))
  point <- as_alfak(fit, minobs = 20, nboot = 5, xval_target = "own", xval_resample = FALSE)
  expect_false(isTRUE(all.equal(point$xval$xval_data$observation, own$xval$xval_data$observation)))
})
