#!/usr/bin/env Rscript
# beta0/beta1 of dev/realdata/beta/PLAN.md on the simulations: one dataset,
# fitted at every rate of the grid (and, for C2, at its beta_sim), scored
# against the truth, with the three criteria of beta1.
#   Rscript dev/realdata/beta/beta_sim_fit.R TASK_FILE ID
# Task columns: kind (eset | c2), set (e1..e4, or a C2 grid setting), dataset_id.
# Output: DEV_RESULTS/realdata/beta/fits/<kind>/<set>/<dataset_id>/<b...>.rds
#   (status, criteria, metrics: stt_score_predictions() metrics, seconds).
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library_dir <- Sys.getenv("DEV_LIBDIR", file.path(dev_results, "lib", "karst_cd342db"))
.libPaths(c(library_dir, .libPaths()))
loadNamespace("karst")
source(file.path(dev_repo, "dev", "realdata", "beta", "R", "beta.R"))
Sys.setenv(OMP_NUM_THREADS = "1")
out_dir <- file.path(beta_results, "fits", task$kind, task$set, task$dataset_id)

if (task$kind == "eset") {
  source(file.path(dev_repo, "dev", "R", "eset.R"))
  dev_source_validation()
  loaded <- dev_eset_load(task$set, task$dataset_id)
  dataset <- loaded$dataset; nodes <- loaded$nodes
  counts <- mt_counts(dataset, "ends")
  times <- dataset$times[c(1L, length(dataset$times))]
  base_args <- list(counts, times = times, birth_rate = 1, min_cn = 1L, max_cn = 8L)
  truth <- sm_true_net_rate(nodes$karyotype, dataset$truth, dataset$beta)
  reference_min <- loaded$config$scoring$reference_min_count
  true_beta <- dataset$beta
  rates <- beta_default * beta_grid
  labels <- beta_label(beta_grid)
} else {
  source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
  c2_bench_source()
  d <- c2_dataset(task$set, task$dataset_id)
  nodes <- c2_nodes(task$set, task$dataset_id)
  base_args <- list(d$counts, dt = d$interval, birth_rate = 1, min_cn = 1L, max_cn = as.integer(d$max_cn))
  truth <- c2_true_rate(nodes$karyotype, d)
  reference_min <- 5
  true_beta <- d$truth$beta_sim
  rates <- c(beta_default * beta_grid, true_beta)
  labels <- c(beta_label(beta_grid), "bsim")
}
seed <- beta_seed(task$dataset_id)
for (j in seq_along(rates)) {
  path <- file.path(out_dir, paste0(labels[j], ".rds"))
  if (file.exists(path)) next
  set.seed(seed)
  run <- beta_fit_timed(c(base_args, list(beta = rates[j])))
  key <- list(kind = task$kind, set = task$set, dataset_id = task$dataset_id, beta = rates[j], label = labels[j],
              multiple = rates[j] / beta_default, true_beta = true_beta, seconds = run$seconds,
              karst = as.character(utils::packageVersion("karst")))
  if (inherits(run$fit, "error")) {
    record <- c(key, beta_failure(run$fit))
  } else {
    fit <- run$fit
    p <- stats::predict(fit, nodes$karyotype)
    pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
                       model_sd = p$relative_sd_model, shell = p$shell, stringsAsFactors = FALSE)
    score <- stt_score_predictions(nodes, truth, pred, reference_min)
    started <- proc.time()[["elapsed"]]
    criteria <- beta_criteria(fit, seed)
    record <- c(key, list(status = "ok", criteria = criteria, cv_seconds = proc.time()[["elapsed"]] - started,
                          metrics = score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"),
                                                  c("stratum", "n_nodes", "rmse", "rmse_zero", "bias", "spearman")]))
  }
  dev_save_rds(record, path)
}
