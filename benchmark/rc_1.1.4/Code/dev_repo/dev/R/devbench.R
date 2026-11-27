# Fitting and scoring one development-benchmark dataset (GRF landscapes, dev or
# heldout split) with a dev method, in the form calibration_fit.R reads:
# scored nodes with the fit's model sd, shell and own relative fitness.
# Task columns: method, lib, dataset_id, drift ("none"/"record"), c_drift,
# beta_multiple, extra. Settings: the dev benchmark's karst_sparse_joint_balance
# (birth_rate 1, min_cn 1, max_cn 8, initial_state "balance").

dev_devbench_input <- function() {
  local <- file.path(dev_inputs, "devbench")
  # rc freeze: the development split's inputs are data/own_benchmark_dev
  if (dir.exists(local)) local else file.path(Sys.getenv("RC_HOME"), "data", "own_benchmark_dev")
}

dev_devbench_source <- function() {
  if (exists("stt_score_predictions", mode = "function")) return(invisible())
  for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R")) {
    source(file.path(dev_repo, "benchmark", "sparse_two_timepoint", part))
  }
}

dev_devbench_fit <- function(task) {
  dev_devbench_source()
  input <- dev_devbench_input()
  config <- stt_config(file.path(dev_results, "abm_cache", "devbench"))
  d <- readRDS(file.path(input, "datasets", paste0(task$dataset_id, ".rds")))
  nodes <- readRDS(file.path(input, "nodes", paste0(task$dataset_id, ".rds")))
  arguments <- list(d$counts, dt = d$interval, beta = d$beta * task$beta_multiple, birth_rate = 1,
                    min_cn = 1L, max_cn = 8L, initial_state = "balance")
  drift_variance <- NA_real_
  if (identical(task$drift, "record")) {
    passages <- readRDS(file.path(dev_results, "phase1", "devbench_passages.rds"))
    row <- passages[passages$dataset_id == task$dataset_id, ]
    if (nrow(row) != 1L || !row$check) stop("No checked passage record for ", task$dataset_id)
    drift_variance <- dev_drift_variance(row$survivors[[1L]], task$c_drift)
    arguments$drift <- drift_variance
  }
  if (!is.null(task$extra) && !is.na(task$extra) && nzchar(task$extra)) {
    arguments <- c(arguments, eval(parse(text = task$extra)))
  }
  set.seed(d$seed)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
  key <- list(method = task$method, dataset_id = task$dataset_id, scenario = d$scenario,
              drift_variance = drift_variance, seconds = proc.time()[["elapsed"]] - started)
  if (inherits(fit, "error")) return(c(key, list(status = "failed", detail = conditionMessage(fit))))
  prediction <- stats::predict(fit, nodes$karyotype)
  pred <- data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                     fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                     shell = prediction$shell, stringsAsFactors = FALSE)
  truth <- stt_true_growth_rate(nodes$karyotype, d$truth)
  score <- stt_score_predictions(nodes, truth, pred, config$scoring$reference_min_count)
  table <- score$nodes
  at <- match(table$karyotype, pred$karyotype)
  table$model_sd <- pred$model_sd[at]
  table$shell <- pred$shell[at]
  table$relative <- pred$fitness[at]
  c(key, list(status = "ok", objective = fit$diagnostics$objective, history = fit$diagnostics$history,
              metrics = score$metrics, nodes = table))
}

dev_devbench_output <- function(task) {
  file.path(dev_results, "phase1", "devbench", task$method, paste0(task$dataset_id, ".rds"))
}
