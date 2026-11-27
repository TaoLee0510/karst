# Fitting one input with alfak (dev2) or karst, writing alfak's files.
#
# Both write, per minobs, bootstrap_res.Rds, landscape.Rds,
# landscape_posterior_samples.Rds and xval.Rds into
# <outdir>/minobs_<m>/, plus run.Rds: status ("ok" or "failed"), the error
# message and the time taken. karst is fitted once and its outputs are
# derived for each minobs, which gives the same files as calling
# karst::alfak_compat() once per minobs. karst's minobs folders appear only
# after its fit, so a task killed during the fit leaves none; the collect
# steps find such fits from the round's task tables (R/tasks.R).

# saveRDS() to a temporary name next to `path`, then renamed: a task killed
# while writing leaves no truncated file under the final name (run_task.R
# skips a step whose file exists). The temporary name starts with a dot, so
# that list.files() does not list it.
save_rds_atomic <- function(object, path) {
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".", Sys.getpid(), ".tmp"))
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not write ", path)
  }
  invisible(path)
}

write_run <- function(dir, status, seconds, message = NA_character_, extra = list()) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  save_rds_atomic(c(list(status = status, seconds = seconds, message = message, finished = Sys.time()), extra),
                  file.path(dir, "run.Rds"))
}

# Sampling times of an input, as alfak and karst read them (passage_times
# NULL): column names times yi$dt.
input_times <- function(yi) as.numeric(colnames(yi$x)) * yi$dt

# Peak resident memory of this R process so far, in MB (VmHWM; NA where
# /proc is not available).
peak_rss_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  s <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
  v <- grep("^VmHWM:", s, value = TRUE)
  if (!length(v)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", v)) / 1024
}

# One value of a diagnostic, NA when absent (run.Rds of earlier commits).
first_or_na <- function(x, na = NA_real_) if (is.null(x) || !length(x)) na else x[[1L]]

fit_alfak <- function(yi, outdir, minobs, n0, nb, pm, settings = bench_settings) {
  dir <- file.path(outdir, paste0("minobs_", minobs))
  start <- proc.time()[["elapsed"]]
  result <- tryCatch({
    alfakR::alfak(yi, dir, passage_times = NULL, minobs = minobs, nboot = settings$nboot, n0 = n0, nb = nb,
                  pm = pm, correct_efflux = settings$correct_efflux)
    NULL
  }, error = function(e) conditionMessage(e))
  write_run(dir, if (is.null(result)) "ok" else "failed", proc.time()[["elapsed"]] - start,
            if (is.null(result)) NA_character_ else result, list(package = find.package("alfakR")))
  invisible(is.null(result))
}

# karst's ways of fitting more than two time points (fit_sparse_joint(series =
# ...), c7f809c on). Method "karst" is the default, "anchored"; the other two
# are methods of their own, fitted to the designs with more than two time points
# only (with two time points all three give the same fit, bit for bit).
karst_series_methods <- c(karst_pairwise = "pairwise", karst_state_space = "state_space")
# methods compared with the default ("anchored") as well as with alfak. (The
# alfak2 era's "auto" arm, karst_auto, one of the three "all" fits picked per
# input by a rule, failed its criterion in round 05a8855 and is not a method
# here any more; scripts/auto_picks.R, compose_auto.R and auto_criterion.R
# are kept for the record only.)
karst_other_methods <- names(karst_series_methods)
# karst's choice of the missegregation rate (fit_sparse_joint(beta = "auto"),
# ac36fd1 on: the fit of lowest objective over 5e-5 times each multiple).
# Each candidate rate is a task of its own, written to
# BENCH_RESULTS/candidates/<source>/<dataset>/<design>/<label>/
# (scripts/make_beta_tasks.R); per input the candidate of lowest objective is
# linked in as karst_beta_auto (scripts/beta_picks.R, compose_beta.R), which
# equals beta = "auto". karst_beta_fixed is the candidate at 5e-5 where it
# is not the round's karst fit: the missegregation simulations, whose karst
# fits take the true rate. Unlike the series methods, both have fits of their
# own on every design ("2 time points" and anchored "all time points").
karst_beta_default <- 5e-5
karst_beta_multiples <- c(1 / 16, 1 / 4, 1, 4, 16, 64)
karst_beta_methods <- c(karst_beta_auto = "auto", karst_beta_fixed = "5e-5")
beta_label <- function(multiple) ifelse(multiple < 1, paste0("x1_", round(1 / multiple)), paste0("x", multiple))
# Variants of the default fitted for checks outside the paper's tests, into
# a results folder of their own (BENCH_RESULTS), so that the round's tables
# do not include them: karst_free, initial_state = "free", the Salehi fits of
# the horizon test (validation/sparse_joint/salehi_horizon).
karst_variant_methods <- list(karst_free = list(initial_state = "free"))
karst_methods <- c("karst", karst_other_methods, names(karst_beta_methods), names(karst_variant_methods))
is_karst_method <- function(method) method %in% karst_methods

# `...` goes to karst::fit_sparse_joint() (min_cn = 0 for data with
# nullisomic karyotype calls). Any number of time points: dt is the vector of
# intervals between consecutive time points (one interval for two time points;
# karst commits before the multi-timepoint version fit two time points only and
# stop otherwise, which is recorded as a failed fit).
fit_karst <- function(yi, outdir, minobs, n0, nb, pm, settings = bench_settings, ...) {
  x <- as.matrix(yi$x)
  storage.mode(x) <- "double"
  times <- input_times(yi)
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(karst::fit_sparse_joint(x, dt = diff(times), beta = pm, ...), error = function(e) e)
  fit_seconds <- proc.time()[["elapsed"]] - start
  # peak memory of the process up to the end of the fit (input and fit)
  fit_peak_mb <- peak_rss_mb()
  failed <- inherits(fit, "error")
  # karst after dea7bd8 reports bootstrap_res frequencies at time 0 of the
  # time axis, as alfak does, given the first time point's time; dea7bd8 reports
  # them at the first time point. run.Rds records which.
  time0 <- "first_time" %in% names(formals(karst::as_alfak))
  # A failed fit records the error's class (karst_domain_coverage_error: no
  # cell of a time point inside the model domain, from 2d86b39 on); a fit records
  # the share of each time point inside the domain and, with candidate histories
  # (af362af on), the chosen history and the profile. With more than two
  # time points (21aedca on), the weight of the time points in between and the majors'
  # dispersion it comes from; `optimisation` is the convergence check (21aedca
  # on; `convergence` is its verdict there, nlminb's code before). With
  # `series` (c7f809c on), the method and the state space's drift; a
  # multi-time-point fit left with two time points (05a8855 on: pairwise or state space
  # after leaving time points out) is a two-time-point fit and has no `series`
  # (`samples_fitted` 2). From the ac36fd1 round on, `objective` (the fit's
  # negative log marginal likelihood, nlminb's; it picked the history up to
  # karst 1.1.1, and picks the missegregation rate in scripts/beta_picks.R) and `beta` (the rate fitted) are recorded too. From karst
  # 1.1.2 on, `posterior_objective` (the objective at the posterior mode the
  # fit reports, which picks the history; also a column of `history_profile`);
  # from 1.1.3 on, each minobs's `rebuild` (as_alfak()'s record of how the
  # posterior was rebuilt: `gradient`, the largest gradient over the latent
  # variables at the fit's stored mode, where the posterior is taken, and
  # `xval_folds`, the cross-validation folds whose mode was found, of them
  # those found only by the second search, and those left out). From
  # 827c4b5 on, `initial_state` (the state fitted: "balance" or "free") and,
  # for two time points under initial_state = "auto", `initial_state_check`
  # (the cloud's log level over the balance, its sd, the threshold, the
  # choice and the check fit's convergence).
  info <- if (failed) list(condition_class = class(fit)[[1L]], samples = ncol(x), fit_peak_mb = fit_peak_mb,
                           beta = pm) else
    list(objective = fit$diagnostics$objective, posterior_objective = fit$diagnostics$posterior_objective,
         beta = fit$data$beta, beta_profile = fit$diagnostics$beta_profile,
         initial_state = fit$diagnostics$initial_state, initial_state_check = fit$diagnostics$initial_state_check,
         domain_coverage = fit$diagnostics$domain_coverage, history = fit$diagnostics$history,
         history_profile = fit$diagnostics$history_profile, convergence = fit$diagnostics$convergence,
         optimisation = fit$diagnostics$optimisation,
         intermediate_weight = fit$diagnostics$intermediate_weight,
         intermediate_dispersion = fit$diagnostics$intermediate_dispersion,
         dropped_samples = fit$diagnostics$dropped_samples,
         series = fit$diagnostics$series, drift = fit$diagnostics$drift,
         samples = ncol(x), samples_fitted = ncol(x) - NROW(fit$diagnostics$dropped_samples),
         fit_peak_mb = fit_peak_mb)
  for (m in minobs) {
    dir <- file.path(outdir, paste0("minobs_", m))
    if (failed) {
      write_run(dir, "failed", fit_seconds, conditionMessage(fit), c(info, list(package = find.package("karst"))))
      next
    }
    start <- proc.time()[["elapsed"]]
    # karst 1.1.3 on: how as_alfak() rebuilt the posterior (as_alfak()'s `rebuild`)
    rebuild <- NULL
    result <- tryCatch({
      args <- list(fit, minobs = m, nboot = settings$nboot, n0 = n0, nb = nb,
                   correct_efflux = settings$correct_efflux)
      if (time0) args$first_time <- times[1]
      out <- do.call(karst::as_alfak, args)
      rebuild <- out$rebuild
      karst::write_alfak_outputs(out, dir)
      NULL
    }, error = function(e) conditionMessage(e))
    write_run(dir, if (is.null(result)) "ok" else "failed", proc.time()[["elapsed"]] - start,
              if (is.null(result)) NA_character_ else result,
              c(list(fit_seconds = fit_seconds, freq_origin = if (time0) "time0" else "first_sample",
                     package = find.package("karst"), rebuild = rebuild), info))
  }
  invisible(!failed)
}
