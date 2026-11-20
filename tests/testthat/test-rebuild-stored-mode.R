# karst 1.1.3: a fit's own posterior is rebuilt at its stored mode, with no new
# search. 1.1.1 and 1.1.2 searched again; where TMB failed at the stored mode
# (a numerically singular inner Hessian) they searched from the initial latent
# values and could end at another, worse mode.

test_that("a fit's own posterior is rebuilt at its stored mode, with no search", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  real_mode <- sparse_joint_inner_mode
  # any search fails: the fit's own posterior does not need one
  local_mocked_bindings(sparse_joint_inner_mode = function(objective, par_fixed, start = NULL) {
    stop(errorCondition("mocked", class = "karst_inner_failure"))
  })
  rebuilt <- sparse_joint_rebuild(fit)
  solver <- rebuilt$solver
  expect_true(solver$stored)
  expect_null(solver$restart)
  expect_equal(solver$value, fit$diagnostics$posterior_objective)
  expect_lt(solver$gradient, 1e-4)
  latent <- unique(solver$names)
  expect_equal(solver$mode, sparse_joint_latent_vector(fit$posterior, rebuilt$problem$map, latent))
  g <- which(solver$names == "g")
  expect_equal(sqrt(solver$variance(g)), fit$summary$log_birth_sd, tolerance = 1e-8)
  mode <- sparse_joint_rebuild(fit, posterior = FALSE)
  expect_equal(mode$g, fit$posterior$g)
  expect_null(attr(mode, "restart"))
  # a cross-validation fold changes the problem: it searches, and here fails
  n <- nrow(fit$summary)
  expect_error(sparse_joint_rebuild(fit, weight = c(0, rep(1, n - 1L)), posterior = FALSE),
               class = "karst_inner_failure")
  # as_alfak(): the fold is left out and counted, the rebuild recorded
  out <- as_alfak(fit, minobs = 5, nboot = 4)
  expect_equal(out$rebuild$gradient, solver$gradient)
  expect_named(out$rebuild$xval_folds, c("rebuilt", "restarted", "failed"))
  expect_equal(out$rebuild$xval_folds[["rebuilt"]], 0L)
  expect_gt(out$rebuild$xval_folds[["failed"]], 0L)
  expect_null(attr(out$xval, "folds"))
  # with searches back, the same landscape, the folds rebuilt
  local_mocked_bindings(sparse_joint_inner_mode = real_mode)
  again <- as_alfak(fit, minobs = 5, nboot = 4)
  expect_equal(again$landscape, out$landscape)
  expect_gt(again$rebuild$xval_folds[["rebuilt"]], 0L)
  expect_equal(again$rebuild$xval_folds[["failed"]], 0L)
  expect_null(attr(again$xval, "folds"))
})

test_that("the posterior of a fit at a stored mode that is not exact is the fit's own", {
  # a stored mode moved off the mode: the rebuild keeps it and reports its
  # gradient; 1.1.2 searched from there and moved back
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  moved <- fit
  moved$posterior$g <- fit$posterior$g + 0.01
  solver <- sparse_joint_rebuild(moved)$solver
  g <- which(solver$names == "g")
  expect_equal(solver$mode[g], moved$posterior$g)
  expect_gt(solver$gradient, 1e-3)
})
