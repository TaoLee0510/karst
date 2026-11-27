#!/usr/bin/env Rscript
# Step 2-3 of dev/realdata/noise/PLAN_v4.md: one fit of one screened v3
# simulation (first and last time point; the screen's counts, i.e. after its
# injection), scored against the truth.
#   T0: beta 5e-5; T1: beta = beta_hat; T2: beta 5e-5 and noise_phi = phi_hat
#   (library karst_v4); T2r: beta 5e-5 and noise_phi = phi_r, the robust
#   variance inflation at beta_hat (PLAN_v4.md correction 16:24);
#   T2_phi1: noise_phi = 1 (identity check only).
# Settings as the sets' earlier fits: e sets birth rate 1, copy numbers 1-8,
# times of the first and last time point; C2 dt = interval, copy numbers 1 to
# the simulation's maximum. Seed: the dataset's seed.
# Scoring: stt_score_predictions() on the v3 nodes (e sets) or
# stt_scoring_nodes() of the original counts (C2), truth sm_true_net_rate()
# or c2_true_rate(). Also the share of nodes (d0-d2) whose 95% interval
# excludes 0 ("false selection" on e3).
#   Rscript dev/realdata/noise/v4/fit_task.R TASK_FILE ID
# Task columns: kind, set, variant, dataset_id, method, keep (1: store the
#   predictions, for the identity check).
# Output: DEV_RESULTS/realdata/noise/v4/fits/<method>[__<library>]/<kind>_<set>_<variant>/<dataset_id>.rds
#   (the library suffix only when it is not karst_v4)
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library_dir <- Sys.getenv("DEV_LIBDIR", "")
if (!nzchar(library_dir) || !dir.exists(file.path(library_dir, "karst"))) stop("No karst in ", library_dir)
.libPaths(c(library_dir, .libPaths()))
loadNamespace("karst")
lib <- basename(library_dir)
Sys.setenv(OMP_NUM_THREADS = "1")
dev_source_validation()
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
c2_bench_source()
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
tag <- sprintf("%s_%s_%s", task$kind, task$set, task$variant)
method_dir <- if (lib == "karst_v4") task$method else paste0(task$method, "__", lib)
out_file <- file.path(dev_results, "realdata", "noise", "v4", "fits", method_dir, tag, paste0(task$dataset_id, ".rds"))
if (file.exists(out_file)) quit(save = "no")
scr <- readRDS(file.path(dev_results, "realdata", "noise", "v4", "screen", tag, paste0(task$dataset_id, ".rds")))
scr$source <- rc_as_run_path(scr$source)  # rc freeze
d <- readRDS(scr$source)
key <- list(kind = task$kind, set = task$set, variant = task$variant, dataset_id = task$dataset_id, method = task$method,
            library = lib, half = scr$half, landscape_id = scr$landscape_id, class = scr$class, depth = scr$depth,
            bottleneck = scr$bottleneck, truth_beta = scr$truth, beta_hat = scr$beta_hat, phi_hat = scr$phi_hat)
save <- function(rec) dev_save_rds(c(key, rec), out_file)
if (!identical(scr$status, "ok")) { save(list(status = "no screen")); quit(save = "no") }

## counts: the original first and last time point, or the screen's injected ones
original <- if (task$kind == "eset") mt_counts(d, "ends") else d$counts
counts <- if (task$variant == "none") original else {
  u <- union(names(scr$x0), names(scr$x1))
  m <- matrix(0, length(u), 2L, dimnames = list(u, colnames(original)[c(1L, ncol(original))]))
  m[names(scr$x0), 1L] <- scr$x0; m[names(scr$x1), 2L] <- scr$x1
  storage.mode(m) <- storage.mode(original)
  m
}
beta <- switch(task$method, T0 = 5e-5, T1 = scr$beta_hat, T2 = 5e-5, T2r = 5e-5, T2_phi1 = 5e-5)
if (!is.finite(beta)) { save(list(status = "no beta_hat")); quit(save = "no") }
arguments <- if (task$kind == "eset") {
  list(counts, times = d$times[c(1L, length(d$times))], beta = beta, birth_rate = 1, min_cn = 1L, max_cn = 8L)
} else list(counts, dt = d$interval, beta = beta, birth_rate = 1, min_cn = 1L, max_cn = as.integer(d$max_cn))
if (task$method == "T2") arguments$noise_phi <- scr$phi_hat
if (task$method == "T2_phi1") arguments$noise_phi <- 1
if (task$method == "T2r") {
  robust <- env_phi_robust(scr$x0, scr$x1, scr$dt, 1, scr$beta_hat)
  key$phi_r <- robust$phi_r; key$phi_r_testable <- robust$n_testable
  arguments$noise_phi <- robust$phi_r
}

## nodes and truth
if (task$kind == "eset") {
  nodes <- readRDS(file.path(dirname(dirname(scr$source)), "nodes", paste0(task$dataset_id, ".rds")))
  truth <- sm_true_net_rate(nodes$karyotype, d$truth, d$beta)
} else {
  nodes <- stt_scoring_nodes(original, 1L, as.integer(d$max_cn), 2000L, 2000L, d$seed)
  d_truth <- d
  d_truth$truth <- readRDS(file.path(dirname(dirname(scr$source)), "truth", paste0(d$landscape_id, ".rds")))
  truth <- c2_true_rate(nodes$karyotype, d_truth)
}

set.seed(d$seed)
started <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
seconds <- proc.time()[["elapsed"]] - started
if (inherits(fit, "error")) {
  save(list(status = "failed", seconds = seconds, condition_class = class(fit)[[1L]], message = conditionMessage(fit)))
  quit(save = "no")
}
p <- stats::predict(fit, nodes$karyotype)
pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
                   stringsAsFactors = FALSE)
score <- stt_score_predictions(nodes, truth, pred, 5)
nt <- score$nodes
use <- nt$distance <= 2L & nt$available & is.finite(nt$pred_sd)
false_sel <- if (any(use)) mean(abs(nt$pred_relative[use]) > stats::qnorm(0.975) * nt$pred_sd[use]) else NA_real_
rec <- list(status = "ok", seconds = seconds, objective = fit$diagnostics$objective,
            metrics = score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), ],
            false_selection = false_sel, n_scored = sum(use),
            domain_coverage = fit$diagnostics$domain_coverage)
if (isTRUE(task$keep == 1)) rec$predictions <- p
save(rec)
