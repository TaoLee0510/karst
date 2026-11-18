# Collect the ABM and missegregation scores into tables.
#
#   Rscript scripts/collect_abm.R
#
# Reads every score.Rds under BENCH_RESULTS/fits/{abm,misseg} and writes to
# BENCH_RESULTS/tables/:
#   abm_dfl.Rds, misseg_dfl.Rds  landscape scores, one row per fit (the dfl
#                                table of the paper's ABM_summaries.Rds plus
#                                design, method and fit status columns)
#   abm_dfp.Rds                  forecast scores in S05's long format: per
#                                fit, passage and metric, prediction,
#                                baseline, win and angle
#   abm_initial.Rds              Supp. Note 3.2: estimated vs sampled
#                                frequencies of the frequent karyotypes at
#                                the first time point (S04B), one row per fit
#   abm_runs.Rds                 status, time, error class and chosen history
#                                of every fit and forecast
#   abm_inputs.Rds               per simulation and design, before any fit:
#                                the share of each time point inside karst's
#                                model domain and whether the input is
#                                analysable (R/coverage.R)
# The paper's own tables are added as method "alfak_paper" (ntp 2, 4, 8 as
# designs w2_ends, w4_all, w8_all).
# First every fit of the round's task tables is checked (R/tasks.R): run.Rds,
# and for a fit that ran score.Rds and, ABM, forecast.Rds. Anything missing
# stops the step (tables/missing_abm.tsv) unless ALLOW_MISSING=1, which
# records a fit without run.Rds as status "missing".

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
out <- file.path(paths$results, "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fit_dirs <- function(source) Sys.glob(file.path(paths$results, "fits", source, "*", "*", "*", "minobs_*"))
dir_info <- function(d) {
  p <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
  data.frame(source = p[5], sim_id = p[4], design = p[3], method = p[2], minobs = sub("minobs_", "", p[1]),
             stringsAsFactors = FALSE)
}

## The round's tasks ------------------------------------------------------------
not_run <- check_round_outputs(round_expected_fits(c("abm", "misseg"), paths), "abm", function(row, run) {
  if (!identical(run$status, "ok")) return(character())
  c(if (!file.exists(file.path(row$dir, "score.Rds"))) "score.Rds",
    if (row$source == "abm" && !file.exists(file.path(row$dir, "forecast.Rds"))) "forecast.Rds")
}, paths)

## Status and time ------------------------------------------------------------
# A fit folder without run.Rds: the task was killed (memory, time); a task
# killed before karst's fit ended leaves no folder (listed from the task
# tables, ALLOW_MISSING=1).
expected_dirs <- if (is.null(not_run)) character() else not_run$dir[!dir.exists(not_run$dir)]
runs <- do.call(rbind, lapply(c(fit_dirs("abm"), fit_dirs("misseg"), expected_dirs), function(d) {
  run <- file.path(d, "run.Rds")
  r <- if (file.exists(run)) readRDS(run) else list(status = "missing", seconds = NA_real_, message = NA_character_)
  fc <- file.path(d, "forecast.Rds")
  cbind(dir_info(d), status = r$status, fit_seconds = r$seconds,
        karst_fit_seconds = first_or_na(r$fit_seconds),
        forecast_seconds = if (file.exists(fc)) readRDS(fc)$seconds else NA_real_,
        message = first_or_na(r$message, NA_character_),
        condition_class = first_or_na(r$condition_class, NA_character_),
        # the first value where a multi-time-point fit reports one per interval
        history = first_or_na(r$history), convergence = first_or_na(r$convergence),
        samples = first_or_na(r$samples, NA_integer_), karst_fit_peak_mb = first_or_na(r$fit_peak_mb),
        intermediate_weight = first_or_na(r$intermediate_weight),
        intermediate_dispersion = first_or_na(r$intermediate_dispersion),
        series = first_or_na(r$series, NA_character_), drift = first_or_na(r$drift), samples_fitted = first_or_na(r$samples_fitted, NA_integer_),
        dropped_samples = if (is.null(r$dropped_samples)) NA_integer_ else NROW(r$dropped_samples),
        # candidate histories (af362af on): the first is dt, the fit of earlier
        # commits; how many candidates' optimisations did not converge
        dt_convergence = first_or_na(r$history_profile$convergence),
        candidates = if (is.null(r$history_profile)) NA_integer_ else nrow(r$history_profile),
        candidates_unconverged = if (is.null(r$history_profile)) NA_integer_ else
          sum(r$history_profile$convergence != 0),
        coverage = if (is.null(r$domain_coverage)) NA_real_ else min(r$domain_coverage$fraction),
        stringsAsFactors = FALSE)
}))
saveRDS(runs, file.path(out, "abm_runs.Rds"))
inputs <- rbind(abm_input_coverage("abm", paths), abm_input_coverage("misseg", paths))
saveRDS(inputs, file.path(out, "abm_inputs.Rds"))
cat("ABM and missegregation inputs (simulation x design):", nrow(inputs), " not analysable:",
    sum(inputs$analysable %in% FALSE), " not judged (installed karst reads two time points only):",
    sum(is.na(inputs$analysable)), "\n")

## Landscape and forecast scores ----------------------------------------------
reformat <- function(df, metric) {
  colnames(df) <- c("prediction", "baseline")
  df$passage <- seq_len(nrow(df)) - 1L
  df$metric <- metric
  df
}
collect <- function(source) {
  scored <- Filter(function(d) file.exists(file.path(d, "score.Rds")), fit_dirs(source))
  res <- lapply(scored, function(d) list(info = dir_info(d), s = readRDS(file.path(d, "score.Rds"))))
  dfl <- do.call(rbind, lapply(res, function(x) cbind(x$info, w = x$s$info$w, abmrep = x$s$info$abmrep,
                                                       pmis = x$s$info$pmis, x$s$summary)))
  dfp <- NULL
  if (source == "abm") {
    dfp <- do.call(rbind, lapply(res, function(x) {
      s <- x$s
      if (is.null(s$df_so)) return(NULL)  # forecast failed: landscape scores only
      tmp <- rbind(cbind(reformat(s$df_so, "overlap"), s$df_a), cbind(reformat(s$df_sc, "cosine"), s$df_a),
                   cbind(reformat(s$df_de, "euclidean"), s$df_a), cbind(reformat(s$df_dw, "wasserstein"), s$df_a))
      cbind(tmp, x$info, w = s$info$w, Rxv = s$summary$Rxv, row.names = NULL)
    }))
    up <- dfp$metric %in% c("overlap", "cosine")
    dfp$win <- ifelse(up, dfp$prediction > dfp$baseline, dfp$prediction < dfp$baseline)
    dfp$pos_xv <- dfp$Rxv > 0
  }
  list(dfl = dfl, dfp = dfp)
}
abm <- collect("abm")
misseg <- collect("misseg")

# The paper's tables, as a reference.
paper <- readRDS(file.path(paths$data, "processed", "ABM_summaries.Rds"))
ntp_design <- c("2" = "w2_ends", "4" = "w4_all", "8" = "w8_all")
pl <- paper$dfl
pl <- data.frame(source = "abm", sim_id = sprintf("w_%s_m_5e-05_rep_%s", sub(".", "p", pl$w, fixed = TRUE), pl$abmrep),
                 design = ntp_design[pl$ntp], method = "alfak_paper", minobs = pl$minobs, w = pl$w,
                 abmrep = pl$abmrep, pmis = 5e-05, pl[c("r", "rfq", "rnn", "rd2", "rho", "rhofq", "rhonn", "rhod2", "R",
                                                        "Rfq", "Rnn", "Rd2", "Rxv", "nfq")], stringsAsFactors = FALSE)
abm$dfl <- if (is.null(abm$dfl)) pl else rbind(abm$dfl, pl[names(abm$dfl)])
saveRDS(abm$dfl, file.path(out, "abm_dfl.Rds"))
saveRDS(abm$dfp, file.path(out, "abm_dfp.Rds"))
saveRDS(misseg$dfl, file.path(out, "misseg_dfl.Rds"))

## Supp. Note 3.2: initial frequencies ------------------------------------------
# alfak's frequencies are at time 0 of the time axis and are projected to the
# first time point with each replicate's fitness (S04B). So are karst's after
# dea7bd8; dea7bd8 reports them at the first time point (run.Rds freq_origin).
initial <- do.call(rbind, lapply(c(fit_dirs("abm"), fit_dirs("misseg")), function(d) {
  info <- dir_info(d)
  b <- file.path(d, "bootstrap_res.Rds")
  run <- file.path(d, "run.Rds")
  # alfak writes bootstrap_res.Rds before its cross-validation; a killed
  # task leaves it without run.Rds.
  if (!file.exists(b) || !file.exists(run)) return(NULL)
  fit <- readRDS(b)
  raw <- if (info$source == "abm") paths$abm_raw else paths$misseg_raw
  sim <- read_abm_sim(info$sim_id, raw)
  yi <- abm_input(sim, abm_designs[abm_designs$design == info$design, ])
  tt <- input_times(yi)[1]
  origin <- readRDS(file.path(d, "run.Rds"))$freq_origin
  at_first <- is_karst_method(info$method) && (is.null(origin) || identical(origin, "first_sample"))
  freq <- if (at_first) fit$final_frequencies else {
    proj <- fit$final_frequencies * exp(fit$final_fitness * tt)
    proj / rowSums(proj)
  }
  est <- colMeans(freq, na.rm = TRUE)
  counts <- yi$x[, 1]
  names(counts) <- rownames(yi$x)
  true <- counts[names(est)] / sum(counts[names(est)])
  cbind(info, correlation = stats::cor(est, true), rmse = sqrt(mean((est - true)^2)), nfq = length(est),
        late_n = sum(true == 0), late_median_est = if (any(true == 0)) stats::median(est[true == 0]) else NA_real_)
}))
saveRDS(initial, file.path(out, "abm_initial.Rds"))

cat("ABM fits scored:", nrow(abm$dfl) - nrow(pl), " missegregation fits scored:", nrow(misseg$dfl),
    " initial-frequency rows:", nrow(initial), "\n")
print(table(runs$source, runs$method, runs$status))
