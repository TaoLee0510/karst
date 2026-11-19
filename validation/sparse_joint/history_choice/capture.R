# History choice study (2026-11-19): every candidate history of a fit, with
# the two objectives the choice could use.
#
# A fit with several candidate histories fits each one (a nested
# fit_sparse_joint() call with a single history) and keeps the one with the
# lowest objective. Since karst 1.1.1 that objective is nlminb's: the Laplace
# objective at the inner mode the scale optimisation ended with. The posterior
# the fit reports sits at the mode the last inner search found (from TMB's
# best mode so far, or, when that search failed, from the initial latent
# values), whose objective can differ. The rules compared in analyse.R choose
# among the same candidate fits; this file records them in one real fit:
#   - fit_sparse_joint() is wrapped in the karst namespace, so that every
#     nested single-history fit is kept, with its history;
#   - sparse_joint_inner_mode() is wrapped, so that the objective at the
#     posterior mode of each of those fits (and whether the second search
#     found it) is kept.
# The package is not changed; the wrappers live in this R process only.

hc_state <- new.env()

hc_install <- function() {
  ns <- asNamespace("karst")
  real_inner <- get("sparse_joint_inner_mode", envir = ns)
  real_fit <- get("fit_sparse_joint", envir = ns)
  inner <- function(objective, par_fixed, start = NULL) {
    value <- real_inner(objective, par_fixed, start)
    hc_state$last_inner <- list(value = as.numeric(value), restart = isTRUE(attr(value, "restart")))
    value
  }
  # The real function runs as itself (its sys.function() and formals are its
  # own); its recursive calls find this wrapper by name in the namespace.
  fit <- function(...) {
    hc_state$last_inner <- NULL
    out <- real_fit(...)
    if (inherits(out, "karst_sparse_joint_fit") && is.null(out$diagnostics$history_profile) &&
        !is.null(hc_state$last_inner)) {
      hc_state$seq <- hc_state$seq + 1L
      hc_state$fits[[length(hc_state$fits) + 1L]] <- list(
        seq = hc_state$seq, history = out$settings$history, fit = out,
        posterior_value = hc_state$last_inner$value, restart = hc_state$last_inner$restart)
    }
    out
  }
  for (name in c("sparse_joint_inner_mode", "fit_sparse_joint")) {
    unlockBinding(name, ns)
  }
  assign("sparse_joint_inner_mode", inner, envir = ns)
  assign("fit_sparse_joint", fit, envir = ns)
  lockBinding("sparse_joint_inner_mode", ns)
  lockBinding("fit_sparse_joint", ns)
  invisible(NULL)
}

# Runs `call` (a function of no arguments that makes the real fit) and
# returns the fit and its candidates: one row per candidate history of the
# kept fit's history_profile, with nlminb's objective, the objective at the
# posterior mode, whether the second search found it, the convergence check,
# and the candidate's fit object (NULL for a candidate whose posterior mode
# could not be found). For one history (no profile), the fit itself.
hc_fit_with_candidates <- function(call) {
  hc_state$fits <- list()
  hc_state$seq <- 0L
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(call(), error = function(e) e)
  seconds <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    return(list(status = "failed", condition_class = class(fit)[[1L]], message = conditionMessage(fit),
                seconds = seconds))
  }
  profile <- fit$diagnostics$history_profile
  histories <- if (is.null(profile)) fit$settings$history else profile$history
  # The last record per history: the first fit of the intermediate weight (at
  # the shortest history) comes before the candidates.
  records <- hc_state$fits
  record_history <- vapply(records, function(r) r$history, numeric(1))
  candidates <- lapply(seq_along(histories), function(i) {
    at <- which(abs(record_history - histories[[i]]) <= 1e-12 * max(1, abs(histories[[i]])))
    r <- if (length(at)) records[[max(at)]] else NULL
    optimisation <- if (is.null(r)) NULL else r$fit$diagnostics$optimisation
    list(
      row = data.frame(
        history = histories[[i]],
        objective = if (is.null(r)) NA_real_ else r$fit$diagnostics$objective,
        posterior_value = if (is.null(r)) NA_real_ else r$posterior_value,
        restart = if (is.null(r)) NA else r$restart,
        convergence = if (is.null(r)) NA_integer_ else as.integer(r$fit$diagnostics$convergence),
        reasons = if (is.null(optimisation$reasons)) "" else paste(optimisation$reasons, collapse = "; "),
        profile_objective = if (is.null(profile)) NA_real_ else profile$objective[[i]],
        stringsAsFactors = FALSE),
      fit = if (is.null(r)) NULL else r$fit)
  })
  rows <- do.call(rbind, lapply(candidates, `[[`, "row"))
  rows$kept <- abs(rows$history - fit$settings$history) <= 1e-12 * max(1, abs(fit$settings$history))
  # The kept candidate's record must be the fit returned.
  kept <- candidates[[which(rows$kept)[[1L]]]]$fit
  identical_kept <- !is.null(kept) && identical(kept$summary, fit$summary)
  list(status = "ok", fit = fit, candidates = rows, fits = lapply(candidates, `[[`, "fit"),
       identical_kept = identical_kept, seconds = seconds)
}
