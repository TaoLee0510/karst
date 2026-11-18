# Forecasts with alfakR::predict_evo (ABM mode), the same for both methods:
# only the landscape differs.

# ABM data, as S03_predict_alfak_abm_sims.R: start from the time point at the last
# sampling time before `last_time` (restricted to the landscape's karyotypes),
# report at the next `horizon` sampling times. Replicate r uses seed r.
forecast_abm <- function(lscape, sim_x, pm, rep, last_time = 120, settings = bench_settings$forecast_abm) {
  times <- as.numeric(colnames(sim_x))
  start <- max(times[times < last_time])
  report <- utils::head(times[times > start], settings$horizon) - start
  x0 <- sim_x[, times == start]
  names(x0) <- rownames(sim_x)
  x0 <- x0[names(x0) %in% lscape$k]
  x0 <- x0 / sum(x0)
  init <- 0 * lscape$mean
  names(init) <- lscape$k
  init[names(x0)] <- x0
  y <- alfakR::predict_evo(lscape = lscape, p = pm, times = c(0, report), x0 = init, prediction_type = "ABM",
                           abm_pop_size = settings$abm_pop_size, abm_max_pop = settings$abm_max_pop,
                           abm_delta_t = settings$abm_delta_t, abm_culling_survival = settings$abm_culling_survival,
                           abm_record_interval = settings$abm_record_interval, abm_seed = rep)
  ix <- vapply(report, function(i) which.min(abs(i - y$time)), integer(1))
  y[c(1L, ix), ]
}

# All replicates for one fit directory, saved as abm_preds/rep_XX.Rds. A
# replicate already written is kept (a task can be rerun after it was
# killed); each is written whole or not at all (save_rds_atomic(), R/fit.R).
forecast_abm_dir <- function(fit_dir, sim_x, pm, last_time = 120, settings = bench_settings$forecast_abm) {
  lscape <- readRDS(file.path(fit_dir, "landscape.Rds"))
  out <- file.path(fit_dir, "abm_preds")
  dir.create(out, showWarnings = FALSE)
  start <- proc.time()[["elapsed"]]
  for (rep in seq_len(settings$reps)) {
    path <- file.path(out, sprintf("rep_%02d.Rds", rep))
    if (file.exists(path)) next
    save_rds_atomic(forecast_abm(lscape, sim_x, pm, rep, last_time, settings), path)
  }
  proc.time()[["elapsed"]] - start
}

# Whether the forecasts of this benchmark's fit folder are complete:
# run_task.R writes forecast.Rds once the last replicate is written, with
# status "ok", or "failed" when the forecast stopped. Replicates written
# before a failure (or before the task was killed) stay for a rerun, which
# skips them, but are never scored.
forecast_complete <- function(fit_dir) {
  path <- file.path(fit_dir, "forecast.Rds")
  file.exists(path) && identical(readRDS(path)$status, "ok")
}
