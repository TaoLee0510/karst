# Fitting and scoring one dataset of the e sets (e1-e4) with a dev method.
# Sourced after common.R. A method is a row of a task table: columns
#   method        name of the fit (folder of the output)
#   set           e1 ... e4
#   dataset_id
#   beta_multiple multiple of the default rate 5e-5 (e4's true rate is
#                 4x or 16x; its fits start from the default, as in step 2)
#   drift         "none" or "record" (R1: F from the replayed passages)
#   c_drift       the passage constant c of R1 (record only)
#   extra         optional R expression of further fit_sparse_joint arguments
# The package comes from the library first on the path (DEV_LIBDIR on HPC).

dev_eset_window_id <- function(dataset_id) sub("_n[0-9]{4}_k[0-9]+$", "_n1000_k8", dataset_id)

# Passages between the first and the last time point of the dataset, from the
# replay of its window (phase 0).
dev_eset_passages <- function(set, dataset) {
  path <- file.path(dev_results, "phase0", "process_noise", "replay",
                    sprintf("replay_%s_%s.rds", set, dev_eset_window_id(dataset$dataset_id)))
  if (!file.exists(path)) stop("No replay of the window of ", dataset$dataset_id, call. = FALSE)
  replay <- readRDS(path)
  if (!all(replay$identical)) stop("The replay of ", dataset$dataset_id, " differs from the data.", call. = FALSE)
  last <- max(dataset$times)
  replay$passages[replay$passages$day > 0 & replay$passages$day <= last + 1e-9, , drop = FALSE]
}

dev_drift_variance <- function(survivors, c_drift) {
  if (!length(survivors)) return(0)
  1 - prod(1 - pmin(c_drift / survivors, 1))
}

dev_eset_load <- function(set, dataset_id) {
  config <- dev_set_config(set)
  dataset <- readRDS(file.path(config$input_dir, "datasets", paste0(dataset_id, ".rds")))
  nodes_path <- file.path(config$input_dir, "nodes", paste0(dataset_id, ".rds"))
  nodes <- if (file.exists(nodes_path)) readRDS(nodes_path) else mt_dataset_nodes(config, dataset)
  list(config = config, dataset = dataset, nodes = nodes)
}

dev_eset_fit <- function(task) {
  dev_source_validation()
  loaded <- dev_eset_load(task$set, task$dataset_id)
  config <- loaded$config; dataset <- loaded$dataset; nodes <- loaded$nodes
  counts <- mt_counts(dataset, "ends")
  times <- dataset$times[c(1L, length(dataset$times))]
  fit_beta <- 5e-5 * task$beta_multiple
  arguments <- list(counts, times = times, beta = fit_beta, birth_rate = 1, min_cn = 1L, max_cn = 8L)
  drift_variance <- NA_real_
  if (identical(task$drift, "record")) {
    passages <- dev_eset_passages(task$set, dataset)
    drift_variance <- dev_drift_variance(passages$after, task$c_drift)
    arguments$drift <- drift_variance
  }
  if (!is.null(task$extra) && !is.na(task$extra) && nzchar(task$extra)) {
    arguments <- c(arguments, eval(parse(text = task$extra)))
  }
  set.seed(dataset$seed)
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(do.call(karst::fit_sparse_joint, arguments),
                                      warning = function(w) {
                                        warnings <<- c(warnings, conditionMessage(w))
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) e)
  seconds <- proc.time()[["elapsed"]] - started
  key <- list(method = task$method, set = task$set, dataset_id = task$dataset_id,
              beta_multiple = task$beta_multiple, drift_variance = drift_variance, seconds = seconds,
              karst = as.character(utils::packageVersion("karst")), library = dirname(find.package("karst")))
  if (inherits(fit, "error")) return(c(key, list(status = "failed", detail = conditionMessage(fit),
                                                 error_class = class(fit)[[1L]])))
  prediction <- stats::predict(fit, nodes$karyotype)
  s <- fit$summary
  data_nll <- -sum(stats::dpois(s$count_t0, s$expected_t0, log = TRUE)) -
    sum(stats::dpois(s$count_t1, s$expected_t1, log = TRUE))
  pred <- data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                     fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                     shell = prediction$shell, stringsAsFactors = FALSE)
  truth <- sm_true_net_rate(nodes$karyotype, dataset$truth, dataset$beta)
  score <- stt_score_predictions(nodes, truth, pred, config$scoring$reference_min_count)
  at <- match(score$nodes$karyotype, pred$karyotype)
  model_sd <- pred$model_sd[at]
  d <- fit$diagnostics
  c(key, list(
    status = "ok", objective = d$objective, data_nll = data_nll, history = d$history,
    convergence = d$convergence, n_states = d$n_states, n_major = d$n_major, warnings = warnings,
    metrics = score$metrics,
    nodes = cbind(score$nodes[c("karyotype", "distance", "truth_relative", "pred_relative", "pred_sd",
                                "available")], model_sd = model_sd,
                  # the fit's own relative fitness and shell (for refitting calibration constants)
                  relative = pred$fitness[at], shell = pred$shell[at]),
    power_counts = sm_power_counts(score$nodes, model_sd, c(0.02, 0.05, 0.1), 0.005),
    summary = s[c("karyotype", "major", "shell", "count_t0", "count_t1", "expected_t0", "expected_t1",
                  "relative_fitness", "relative_sd", "relative_sd_model")]))
}

dev_eset_output <- function(task) {
  file.path(dev_results, "phase1", "eset", task$method, task$set, paste0(task$dataset_id, ".rds"))
}
