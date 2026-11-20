test_that("two time points give the two-time-point fit with every series method", {
  sim <- simulate_series_counts()
  pair <- sim$counts[, c(1L, 4L)]
  pair <- pair[rowSums(pair) > 0, ]
  base <- fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60)
  for (method in c("pairwise", "state_space")) {
    fit <- fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60, series = method)
    expect_identical(fit$summary, base$summary)
    expect_identical(fit$diagnostics$objective, base$diagnostics$objective)
    expect_null(fit$settings$series)
  }
})

test_that("pair domains give each pair of time points its own majors", {
  counts <- rbind("2.2.2" = c(30, 10, 0), "3.2.2" = c(0, 25, 30), "1.1.1" = c(0, 0, 2))
  domains <- sparse_joint_pair_domains(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(domains$labels[seq_len(domains$n_major)], c("2.2.2", "3.2.2"))
  # Both are majors of both pairs; 1.1.1 has too few cells to be one and is
  # three steps from them, outside every domain.
  expect_equal(domains$pair_major[1:2, ], cbind(c(1, 1), c(1, 1)))
  expect_identical(domains$dropped, "1.1.1")
  expect_equal(nrow(domains$member), length(domains$labels))
  expect_true(all(domains$member[seq_len(domains$n_major), ] == 1))
  expect_equal(domains$y[match("2.2.2", domains$labels), ], c(30, 10, 0))
})

test_that("the pairwise method fits each pair as two time points with one landscape", {
  sim <- simulate_series_counts()
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, initial_state = "free",
                          series = "pairwise")
  expect_identical(fit$settings$series, "pairwise")
  expect_identical(fit$diagnostics$series, "pairwise")
  expect_equal(fit$diagnostics$times, sim$times)
  expect_equal(ncol(as.matrix(fit$posterior$zeta)), 3L)
  # The pairs as a cohort that shares every birth rate and the trend. The
  # cohort centres the trend on the mean of the pairs' centres, the pairwise
  # method on the counts at all time points; the fits agree closely (and exactly with
  # the same centre, validation/sparse_joint/series_methods).
  pairs <- lapply(1:3, function(k) {
    y <- sim$counts[, c(k, k + 1L)]
    y[rowSums(y) > 0, , drop = FALSE]
  })
  names(pairs) <- paste0("p", 1:3)
  cohort <- fit_sparse_joint_cohort(pairs, dt = 20, beta = sim$beta, max_cn = 4, initial_state = "free",
                                    sharing = "full", trend = "shared")
  index <- match(fit$summary$karyotype, cohort$landscape$karyotype)
  expect_false(anyNA(index))
  expect_equal(fit$diagnostics$objective, cohort$diagnostics$objective, tolerance = 1e-4)
  expect_lt(max(abs(fit$summary$fitness_mean - cohort$landscape$fitness_mean[index])), 1e-3)
  # Balance at the start of each pair over the history alone (internal
  # option of the evaluation) rather than plus the population's age.
  old <- options(karst.pair_history = "same")
  same <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                           series = "pairwise")
  options(old)
  aged <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                           series = "pairwise")
  expect_identical(same$diagnostics$pair_history, "same")
  expect_identical(aged$diagnostics$pair_history, "age")
  expect_equal(sparse_joint_problem(same$data, same$settings)$tmb_data$pair_history, rep(60, 3))
  expect_equal(sparse_joint_problem(aged$data, aged$settings)$tmb_data$pair_history, 60 + c(0, 20, 40))
  # A rebuilt posterior reproduces the fit's own.
  rebuilt <- sparse_joint_rebuild(fit)
  expect_equal(rebuilt$solver$mode, unname(unlist(lapply(rebuilt$problem$random, function(name) {
    value <- fit$posterior[[name]]
    if (name == "zeta") value[rebuilt$problem$tmb_data$member > 0] else value
  }))), tolerance = 1e-6)
  # The stored settings recreate the fit's objective: its stored mode is a mode there.
  expect_lt(rebuilt$solver$gradient, 1e-6)
})

test_that("the state-space method follows the anchored fit when there is no drift", {
  sim <- simulate_series_counts()
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                          series = "state_space")
  expect_identical(fit$diagnostics$series, "state_space")
  expect_true(fit$diagnostics$drift > 1e-9 && fit$diagnostics$drift < 1)
  expect_identical(names(fit$trend$scales), c("dosage", "curvature", "interaction", "residual", "drift"))
  expect_equal(dim(fit$posterior$zeta_later), c(nrow(fit$summary), 3L))
  rebuilt <- sparse_joint_rebuild(fit)
  expect_equal(rebuilt$solver$mode, unname(unlist(fit$posterior[rebuilt$problem$random])), tolerance = 1e-6)
  # The stored settings recreate the fit's objective: its stored mode is a mode there.
  expect_lt(rebuilt$solver$gradient, 1e-6)
  # Nearly no drift: the closed form from time point to time point, as the anchored
  # model with every time point for every karyotype, up to the second order in
  # the mutational flux.
  old <- options(karst.state_space_drift = 1e-8)
  on.exit(options(old), add = TRUE)
  tight <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            series = "state_space")
  expect_equal(tight$diagnostics$drift, 1e-8)
  anchored <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                               intermediate = "all")
  index <- match(anchored$summary$karyotype, tight$summary$karyotype)
  expect_false(anyNA(index))
  expect_gt(stats::cor(anchored$summary$relative_fitness, tight$summary$relative_fitness[index]), 0.99)
  major <- anchored$summary$major
  expect_lt(max(abs(anchored$summary$relative_fitness[major] - tight$summary$relative_fitness[index][major])), 1e-3)
})

test_that("predictions and alfak outputs work for every series method", {
  sim <- simulate_series_counts()
  for (method in c("pairwise", "state_space")) {
    fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            series = method)
    prediction <- predict(fit, c(fit$summary$karyotype[1:2], "1.1.1.1"))
    expect_true(all(is.finite(prediction$relative_fitness)))
    expect_true(all(prediction$relative_sd > 0))
    out <- as_alfak(fit, minobs = 20, nboot = 5, n0 = 1e5, nb = 1e7)
    expect_true(all(out$landscape$k[out$landscape$fq] %in% rownames(sim$counts)))
    expect_true(is.finite(out$xval$R2R))
    # alfak's level: relative fitness plus the population growth over the
    # first interval
    g0 <- log(1e7 / 1e5) / diff(sim$times)[[1L]]
    expect_equal(out$landscape$mean, predict(fit, out$landscape$k)$relative_fitness + g0, tolerance = 1e-10)
    # frequencies at the posterior mode: the fitted abundances at the first
    # time point (observed counts outside the domain)
    one <- as_alfak(fit, minobs = 20, nboot = 1, xval = FALSE)$bootstrap_res$final_frequencies
    fq <- colnames(one)
    zeta <- sparse_joint_first_sample_zeta(fit)[match(fq, fit$summary$karyotype)]
    level <- ifelse(is.na(zeta), log(pmax(fit$data$counts[fq, 1L], .Machine$double.xmin)), zeta)
    expect_equal(as.numeric(one), unname(exp(level - max(level)) / sum(exp(level - max(level)))), tolerance = 1e-10)
    # and the fitness differences are the simulated ones
    truth <- sim$net[match(fit$summary$karyotype, sim$labels)]
    major <- fit$summary$major
    weight <- sparse_joint_summary_total(fit$summary)[major]
    truth_relative <- truth - sum(truth[major] * weight) / sum(weight)
    # (the majors differ by 0.03 per step; at this seed both methods are
    # within 0.003 of the truth)
    expect_lt(max(abs(fit$summary$relative_fitness[major] - truth_relative[major])), 0.006)
    seen <- !major & sparse_joint_summary_total(fit$summary) > 0 & !is.na(truth)
    expect_lt(mean(fit$summary$relative_fitness[seen]), 0)
    expect_gt(stats::cor(fit$summary$relative_fitness[seen], truth[seen]), 0.4)
  }
})

test_that("a first time point outside every pair's domain stops the fit", {
  # The first time point is kept; its neighbour is left out, and the first and
  # the last time point alone still have no cell of the first in the domain.
  counts <- rbind("2.2.2" = c(0, 20, 15), "3.3.2" = c(1, 0, 0), "1.1.2" = c(1, 0, 0), "2.4.2" = c(1, 0, 0),
                  "4.2.3" = c(1, 0, 0))
  expect_identical(sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)$left_out, 2L)
  expect_error(suppressWarnings(fit_sparse_joint(counts, dt = c(10, 10), beta = 0.001, series = "pairwise")),
               class = "karst_domain_coverage_error")
})

# The simulated series with time point `k` replaced by `cells`.
replace_series_sample <- function(sim, k, cells) {
  counts <- sim$counts
  extra <- setdiff(names(cells), rownames(counts))
  counts <- rbind(counts, matrix(0, length(extra), ncol(counts), dimnames = list(extra, colnames(counts))))
  counts[, k] <- 0
  counts[names(cells), k] <- cells
  counts[rowSums(counts) > 0, , drop = FALSE]
}

test_that("a time point in between with no cell in its pair's domain is left out", {
  sim <- simulate_series_counts()
  # Time point t2 holds single cells far from every major: in pair (t1, t2) it
  # has no cell in the domain. It is left out, and t1 and t3 form a pair.
  far <- c("4.4.4.4" = 1, "1.1.4.4" = 1, "4.1.1.4" = 1)
  counts <- replace_series_sample(sim, 3L, far)
  pairs <- sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(pairs$left_out, 3L)
  expect_identical(pairs$uncovered, "t2")
  kept <- counts[, -3L]
  kept <- kept[rowSums(kept) > 0, ]
  for (method in c("pairwise", "state_space")) {
    expect_warning(fit <- fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                           series = method),
                   class = "karst_samples_left_out_warning")
    direct <- fit_sparse_joint(kept, times = sim$times[-3L], beta = sim$beta, max_cn = 4, history = 60,
                               series = method)
    # The same fit; its columns are named after the time points of the input.
    expect_identical(names(fit$summary)[2:4], c("count_t0", "count_t1", "count_t3"))
    expect_identical(fit$summary, renamed_like(direct$summary, fit$summary))
    expect_identical(fit$diagnostics$dropped_samples$sample, 3L)
    expect_identical(fit$diagnostics$dropped_samples$uncovered, "t2")
    expect_equal(fit$input$counts, counts)
  }
})

test_that("a last time point with no cell in its pair's domain pairs with a time point further in", {
  sim <- simulate_series_counts()
  # t2 has swept to a far clone; t3 holds two single cells of the earlier
  # majors. In (t2, t3) the last time point has no cell in the domain, so t2 is
  # left out and t3 pairs with t1, where those karyotypes are majors.
  counts <- replace_series_sample(sim, 3L, c("4.4.4.4" = 300, "4.4.4.3" = 20))
  counts <- replace_series_sample(list(counts = counts), 4L, c("2.2.2.2" = 1, "2.3.2.2" = 1))
  pairs <- sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(pairs$left_out, 3L)
  expect_identical(pairs$uncovered, "t3")
  # Every pair covered: nothing left out.
  expect_identical(sparse_joint_pair_left_out(sim$counts, 0.01, 3, 5000, 1L, 4L)$left_out, integer(0))
})
