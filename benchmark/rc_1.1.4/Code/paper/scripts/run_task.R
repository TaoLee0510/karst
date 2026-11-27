# Run one task of a task table: fit, forecast and score.
#
#   Rscript scripts/run_task.R TASK_FILE TASK_ID
#
# Outputs go to BENCH_RESULTS/fits/<source>/<dataset>/<design>/<method>/
# minobs_<m>/: alfak's files, run.Rds (fit status and time), abm_preds/
# (forecast replicates), forecast.Rds (forecast time) and, for the ABM data,
# score.Rds (the list of S04's abm_preds_summary.Rds). Steps already done are
# skipped, so a task can be run again after a failure: a karst task refits
# only the minobs without run.Rds (the others keep their files), a forecast
# that did not finish is resumed, and a fit whose forecast is (re)run is
# scored again. These files are written whole or not at all
# (save_rds_atomic(), R/fit.R). A score that cannot be computed ends the
# task with exit status 1 (after the peak memory is recorded), so that
# Slurm shows it as FAILED.
#
# A task table with a `beta_label` column (scripts/make_beta_tasks.R) holds
# candidate missegregation rates of karst (pm = the candidate): the fit and
# its alfak files go to BENCH_RESULTS/candidates/<source>/<dataset>/
# <design>/<beta_label>/minobs_<m>/, without forecasts. Only the candidate
# picked for an input gets its forecasts, through its link in fits/
# (karst_beta_auto, a task of scripts/compose_beta.R's table, with pm = the
# picked rate). A candidate's ABM landscape scores are kept as
# score_landscape.Rds, so that score.Rds is written with the forecasts; the
# missegregation simulations have no forecasts (except with MISSEG_FORECAST=TRUE,
# below), and their score.Rds is final.

args <- commandArgs(TRUE)
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
tasks <- utils::read.delim(args[[1]], stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
task <- tasks[tasks$task == as.integer(args[[2]]), ]
stopifnot(nrow(task) == 1L)
minobs <- as.numeric(strsplit(task$minobs, ",", fixed = TRUE)[[1]])
# A candidate rate (and the rate picked from them) as karst's grid computes
# it, not as printed in the table.
if (!is.null(task$beta_multiple) && !is.na(task$beta_multiple)) task$pm <- karst_beta_default * task$beta_multiple

cat("Task", args[[2]], ":", paste(unlist(task), collapse = " | "), "\n")
for (p in c("alfakR", "karst")) {
  cat(p, as.character(utils::packageVersion(p)), find.package(p), "\n")
}
cat("SIF:", Sys.getenv("BENCH_SIF", "none"), "\n")

candidate <- !is.null(task$beta_label) && !is.na(task$beta_label) && nzchar(task$beta_label)
outdir <- if (candidate) {
  file.path(paths$results, "candidates", task$source, task$dataset, task$design, task$beta_label)
} else {
  file.path(paths$results, "fits", task$source, task$dataset, task$design, task$method)
}
minobs_dir <- function(m) file.path(outdir, paste0("minobs_", m))
fit_ok <- function(m) {
  path <- file.path(minobs_dir(m), "run.Rds")
  file.exists(path) && identical(readRDS(path)$status, "ok")
}
fit_done <- function(m) file.exists(file.path(minobs_dir(m), "run.Rds"))

## Input --------------------------------------------------------------------
if (task$source %in% c("abm", "misseg")) {
  raw <- if (task$source == "abm") paths$abm_raw else paths$misseg_raw
  sim <- read_abm_sim(task$dataset, raw)
  design <- abm_designs[abm_designs$design == task$design, ]
  yi <- abm_input(sim, design)
} else {
  yi <- salehi_input(task$dataset, ends = task$design == "ends", paths)
}
cat("Input time points:", colnames(yi$x), " dt:", yi$dt, "\n")

## Fit ----------------------------------------------------------------------
todo <- minobs[!vapply(minobs, fit_done, logical(1))]
if (length(todo)) {
  if (task$method == "alfak") {
    for (m in todo) fit_alfak(yi, outdir, m, task$n0, task$nb, task$pm)
  } else {
    extra <- if (task$source == "salehi") bench_settings$salehi$karst_args else list()
    if (task$method %in% names(karst_series_methods)) extra$series <- karst_series_methods[[task$method]]
    if (task$method %in% names(karst_variant_methods)) extra <- c(extra, karst_variant_methods[[task$method]])
    # the minobs not done only: the finished ones keep their outputs (new
    # bootstrap draws would no longer match their forecasts and scores)
    do.call(fit_karst, c(list(yi, outdir, todo, task$n0, task$nb, task$pm), extra))
  }
}
for (m in minobs) {
  r <- readRDS(file.path(minobs_dir(m), "run.Rds"))
  cat(sprintf("Fit minobs %s: %s (%.0f s) %s\n", m, r$status, r$seconds, ifelse(is.na(r$message), "", r$message)))
}

## Forecast and score -------------------------------------------------------
# A forecast that cannot run is recorded in forecast.Rds (status "failed"),
# as the paper leaves no forecast there, and the task goes on. Replicates
# written before a failure stay for a rerun but are never scored (only
# forecasts with status "ok" are: forecast_complete(), R/forecast.R). A
# score.Rds from an earlier attempt is removed, so that the fit is scored on
# the new forecast.
forecast <- function(dir, run, extra = list()) {
  unlink(file.path(dir, "score.Rds"))
  start <- proc.time()[["elapsed"]]
  msg <- tryCatch({
    run()
    NA_character_
  }, error = function(e) conditionMessage(e))
  seconds <- proc.time()[["elapsed"]] - start
  save_rds_atomic(c(list(status = if (is.na(msg)) "ok" else "failed", seconds = seconds, message = msg), extra),
                  file.path(dir, "forecast.Rds"))
  cat(sprintf("Forecast %s: %s (%.0f s) %s\n", basename(dir), if (is.na(msg)) "ok" else "failed", seconds,
              ifelse(is.na(msg), "", msg)))
}
if (task$source == "salehi") {
  lineages <- readRDS(paths$salehi_lineages)
  meta <- utils::read.csv(paths$salehi_meta)
  sim_time <- salehi_forecast_time(task$dataset, lineages, meta)
}
scoring_failed <- 0L
for (m in minobs[vapply(minobs, fit_ok, logical(1))]) {
  dir <- minobs_dir(m)
  if (candidate) {
    if (task$source %in% c("abm", "misseg")) {
      name <- if (task$source == "abm") "score_landscape.Rds" else "score.Rds"
      if (!file.exists(file.path(dir, name))) {
        info <- cbind(abm_sim_info(task$dataset), design = task$design, method = task$method, minobs = m,
                      beta_label = task$beta_label)
        score <- tryCatch(abm_fit_summary(dir, sim, info, with_forecast = FALSE), error = function(e) {
          cat("Scoring failed for", dir, ":", conditionMessage(e), "\n")
          NULL
        })
        if (is.null(score)) scoring_failed <- scoring_failed + 1L else save_rds_atomic(score, file.path(dir, name))
      }
    }
    next
  }
  if (!forecast_complete(dir)) {
    # the missegregation simulations get forecasts only for a study that asks
    # for them (MISSEG_FORECAST=TRUE, scripts/fig2d_truth.R --source=misseg)
    if (task$source == "abm" || (task$source == "misseg" && identical(Sys.getenv("MISSEG_FORECAST"), "TRUE"))) {
      forecast(dir, function() forecast_abm_dir(dir, sim$abm_output$x, task$pm))
    } else if (task$source == "salehi" && !is.null(sim_time)) {
      forecast(dir, function() forecast_salehi_dir(dir, yi, sim_time, task$pm, task$nb), list(sim_time = sim_time))
    }
  }
  # ABM scores; without forecasts (failed) only the landscape is scored.
  if (task$source %in% c("abm", "misseg") && !file.exists(file.path(dir, "score.Rds"))) {
    info <- cbind(abm_sim_info(task$dataset), design = task$design, method = task$method, minobs = m)
    with_forecast <- task$source == "abm" && forecast_complete(dir)
    score <- tryCatch(abm_fit_summary(dir, sim, info, with_forecast = with_forecast),
                      error = function(e) {
                        cat("Scoring failed for", dir, ":", conditionMessage(e), "\n")
                        NULL
                      })
    if (is.null(score)) scoring_failed <- scoring_failed + 1L else save_rds_atomic(score, file.path(dir, "score.Rds"))
  }
}

## Peak memory ---------------------------------------------------------------
# Peak resident memory of this task's process (fit, forecasts and scores),
# kept in every minobs folder of the task as the largest over the task's
# attempts (an attempt that finds the steps done uses little). Slurm's
# accounting samples memory every 30 s or so and misses tasks shorter than
# that.
peak <- peak_rss_mb()
for (m in minobs) {
  path <- file.path(minobs_dir(m), "task_peak.Rds")
  if (!dir.exists(minobs_dir(m)) || is.na(peak)) next
  old <- if (file.exists(path)) readRDS(path)$peak_mb else NA_real_
  if (is.na(old) || peak > old) save_rds_atomic(list(peak_mb = peak, recorded = Sys.time()), path)
}
cat(sprintf("Peak memory of the task: %.0f MB\n", peak))
if (scoring_failed > 0L) {
  cat("Scoring failed for", scoring_failed, "fit folder(s); rerun the task after fixing the cause.\n")
  quit(save = "no", status = 1L)
}
cat("Done.\n")
