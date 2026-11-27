#!/usr/bin/env Rscript
# C2 fits (dev/realdata/PLAN.md, phase C): one fit of a C2 dataset, scored
# against the simulation truth on the development benchmark's nodes
# (observed karyotypes, one- and two-step neighbours).
#   Rscript dev/realdata/CE/c2_fit.R TASK_FILE ID      (ID may be a range "a-b")
# Task columns: method, lib, regime, dataset_id, extra. `regime` is H, M or W
# (accepted datasets, c2/regimes/<regime>) or a tuning setting (c2/grid/<setting>).
# `extra`: optional R expression of a list of further fit_sparse_joint
# arguments. Fit settings as the Salehi benchmark: beta 5e-5, birth rate 1;
# copy numbers 1 to the simulation's maximum.
# Output: DEV_RESULTS/realdata/CE/c2/fits/<method>/<regime>/<dataset_id>.rds with
#   status, seconds, (on failure) condition_class, message, coverage;
#   near_nodes (karyotype, distance, count_t0, count_t1, truth_relative,
#   pred_relative, pred_sd, available, model_sd, shell, relative) as far_fit.R,
#   near_metrics, summary (fit$summary columns), diagnostics (domain coverage,
#   history, ...), dispersion (second time point, Pearson), meta.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
args <- commandArgs(TRUE)
tasks <- utils::read.delim(args[[1L]], stringsAsFactors = FALSE)
ids <- if (grepl("-", args[[2L]], fixed = TRUE)) {
  r <- as.integer(strsplit(args[[2L]], "-", fixed = TRUE)[[1L]]); seq(r[[1L]], r[[2L]])
} else as.integer(args[[2L]])
library_dir <- Sys.getenv("DEV_LIBDIR", "")
if (!nzchar(library_dir)) library_dir <- file.path(dev_results, "lib", tasks$lib[[ids[[1L]]]])
.libPaths(c(library_dir, .libPaths()))
loadNamespace("karst")
Sys.setenv(OMP_NUM_THREADS = "1")
c2_bench_source()

for (i in ids[ids <= nrow(tasks)]) {
  task <- tasks[i, , drop = FALSE]
  if (basename(library_dir) != task$lib) stop("Task ", i, " needs library ", task$lib, ", loaded ", library_dir)
  path <- file.path(c2_root, "fits", task$method, task$regime, paste0(task$dataset_id, ".rds"))
  if (file.exists(path)) next
  d <- c2_dataset(task$regime, task$dataset_id)
  nodes <- c2_nodes(task$regime, task$dataset_id)
  extra <- if (!is.na(task$extra) && nzchar(task$extra)) eval(parse(text = task$extra)) else list()
  arguments <- utils::modifyList(list(d$counts, dt = d$interval, beta = d$beta, birth_rate = 1, min_cn = 1L,
                                      max_cn = as.integer(d$max_cn)), extra)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
  seconds <- proc.time()[["elapsed"]] - started
  key <- list(method = task$method, regime = task$regime, dataset_id = task$dataset_id, seconds = seconds,
              library = library_dir,
              meta = list(setting = d$setting, landscape_id = d$landscape_id, wavelength = d$wavelength,
                          interval = d$interval, depth = c(d$depth_t0, d$depth_t1)))
  if (inherits(fit, "error")) {
    record <- c(key, list(status = "failed", condition_class = class(fit)[[1L]], message = conditionMessage(fit),
                          coverage = fit$coverage))
  } else {
    p <- stats::predict(fit, nodes$karyotype)
    pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
                       model_sd = p$relative_sd_model, shell = p$shell, stringsAsFactors = FALSE)
    truth <- c2_true_rate(nodes$karyotype, d)
    score <- stt_score_predictions(nodes, truth, pred, 5)
    table <- score$nodes
    at <- match(table$karyotype, pred$karyotype)
    table$model_sd <- pred$model_sd[at]
    table$shell <- pred$shell[at]
    table$relative <- pred$fitness[at]
    s <- fit$summary
    keep <- intersect(c("karyotype", "major", "shell", "count_t0", "count_t1", "expected_t0", "expected_t1",
                        "fitness_mean", "relative_fitness", "relative_sd", "relative_sd_model"), names(s))
    dg <- fit$diagnostics
    record <- c(key, list(status = "ok", near_metrics = score$metrics, near_nodes = table, summary = s[keep],
                          diagnostics = dg[intersect(c("domain_coverage", "history", "convergence", "objective",
                                                       "n_states", "n_major", "initial_state"), names(dg))],
                          dispersion = c2_dispersion(s)))
  }
  dev_save_rds(record, path)
}
