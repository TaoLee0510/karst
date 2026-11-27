#!/usr/bin/env Rscript
# Step 4 of dev/realdata/noise/PLAN_v7.md (only if step 3 passes): one fit of
# one v7 dataset, T0 (beta 5e-5) or T6 (beta = beta_used), library 1.1.3,
# scored against the truth.
#   e sets: times 0 and interval, birth rate 1, copy numbers 1-8 (as v4);
#   C2: dt = interval, birth rate 1, copy numbers 1 to the simulation's maximum.
#   Nodes stt_scoring_nodes() (2000 / 2000); truth sm_true_net_rate() (e) or
#   c2_true_rate() (C2).
#   Rscript dev/realdata/noise/v7/fit_task.R TASK_FILE ID   (columns: family, dataset_id, method, beta)
# Output: DEV_RESULTS/realdata/noise/v7/fits/<method>/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library_dir <- Sys.getenv("DEV_LIBDIR", "")
if (!nzchar(library_dir) || !dir.exists(file.path(library_dir, "karst"))) stop("No karst in ", library_dir)
.libPaths(c(library_dir, .libPaths()))
loadNamespace("karst")
Sys.setenv(OMP_NUM_THREADS = "1")
dev_source_validation()
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
c2_bench_source()
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
out_file <- v7_path("fits", task$method, paste0(task$dataset_id, ".rds"))
if (file.exists(out_file)) quit(save = "no")
d <- readRDS(v7_check_read(v7_path("sims", task$family, "datasets", paste0(task$dataset_id, ".rds"))))
beta <- as.numeric(task$beta)
key <- list(family = d$family, design = d$design, dataset_id = d$dataset_id, method = task$method, beta = beta,
            multiple = d$multiple, half = d$half, landscape_id = d$landscape_id, karst = as.character(utils::packageVersion("karst")))
counts <- d$counts
arguments <- if (d$family == "e") {
  list(counts, times = c(0, d$interval), beta = beta, birth_rate = 1, min_cn = 1L, max_cn = 8L)
} else list(counts, dt = d$interval, beta = beta, birth_rate = 1, min_cn = 1L, max_cn = as.integer(d$max_cn))
max_cn <- if (d$family == "e") 8L else as.integer(d$max_cn)
nodes <- stt_scoring_nodes(counts, 1L, max_cn, 2000L, 2000L, d$seed)
truth <- if (d$family == "e") sm_true_net_rate(nodes$karyotype, d$truth, d$beta) else
  c2_true_rate(nodes$karyotype, list(truth = readRDS(v7_check_read(d$truth_file))))
set.seed(d$seed)
started <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
seconds <- proc.time()[["elapsed"]] - started
if (inherits(fit, "error")) {
  dev_save_rds(c(key, list(status = "failed", seconds = seconds, message = conditionMessage(fit))), out_file)
  quit(save = "no")
}
p <- stats::predict(fit, nodes$karyotype)
pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd, stringsAsFactors = FALSE)
score <- stt_score_predictions(nodes, truth, pred, 5)
dev_save_rds(c(key, list(status = "ok", seconds = seconds, objective = fit$diagnostics$objective,
                         metrics = score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), ],
                         domain_coverage = fit$diagnostics$domain_coverage)), out_file)
