# Shared code of the beta study (dev/realdata/beta/PLAN.md). Sourced after
# dev/R/common.R with karst loaded.

beta_grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
beta_default <- 5e-5
beta_label <- function(m) paste0("b", format(m, scientific = FALSE, drop0trailing = TRUE))
beta_results <- file.path(dev_results, "realdata", "beta")

# Counts and expectations of the first two time points of a fit's summary
# (columns are named after the input's time points).
beta_two <- function(s) {
  count <- grep("^count_", names(s), value = TRUE)[1:2]
  expected <- grep("^expected_", names(s), value = TRUE)[1:2]
  list(y0 = s[[count[1]]], y1 = s[[count[2]]], e0 = s[[expected[1]]], e1 = s[[expected[2]]])
}

beta_seed <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L

# The three criteria of beta1 for one fit: the marginal objective, the data term
# (Poisson negative log-likelihood of the domain states at the posterior mode)
# and the held-out prediction of the non-major states' counts in 5 folds
# (sparse_joint_rebuild() at the fit's scales, the fold's counts at weight 0).
beta_criteria <- function(fit, seed, folds = 5L) {
  s <- fit$summary
  two <- beta_two(s)
  data_nll <- -sum(stats::dpois(two$y0, two$e0, log = TRUE)) - sum(stats::dpois(two$y1, two$e1, log = TRUE))
  minor <- which(!s$major)
  cv <- NA_real_; cv_failed <- NA_integer_
  if (length(minor) >= folds) {
    set.seed(seed)
    fold <- sample(rep(seq_len(folds), length.out = length(minor)))
    cv <- 0; cv_failed <- 0L
    for (k in seq_len(folds)) {
      out <- minor[fold == k]
      w <- rep(1, nrow(s)); w[out] <- 0
      rb <- tryCatch(suppressWarnings(karst:::sparse_joint_rebuild(fit, weight = w)), error = function(e) e)
      if (inherits(rb, "error")) { cv_failed <- cv_failed + 1L; next }
      labels <- rb$problem$domain$labels
      if (!identical(labels, s$karyotype)) stop("The rebuilt domain differs from the fit's.")
      e0 <- rb$solver$expected_t0[out]; e1 <- rb$solver$expected_t1[out]
      cv <- cv - sum(stats::dpois(two$y0[out], e0, log = TRUE)) - sum(stats::dpois(two$y1[out], e1, log = TRUE))
    }
    if (cv_failed > 0L) cv <- NA_real_
  }
  d <- fit$diagnostics
  list(objective = d$objective, posterior_objective = d$posterior_objective, data_nll = data_nll, cv_nll = cv,
       cv_failed = cv_failed, n_minor = length(minor), n_states = nrow(s), n_major = sum(s$major),
       history = if (is.null(d$history)) NA else d$history[[1L]], convergence = d$convergence)
}

beta_fit_timed <- function(arguments) {
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
  list(fit = fit, seconds = proc.time()[["elapsed"]] - started)
}

beta_failure <- function(fit) {
  list(status = "failed", condition_class = class(fit)[[1L]], message = conditionMessage(fit))
}
