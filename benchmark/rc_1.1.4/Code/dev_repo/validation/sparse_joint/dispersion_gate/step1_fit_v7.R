#!/usr/bin/env Rscript
# Step 1 of PLAN.md: one KARST fit (library 1.1.3) of one v7 dataset
# (dev/realdata/noise/v7 simulations, read only) at one multiple of 5e-5,
# scored on stt_scoring_nodes() against the net rate at the true rate. The
# fit and the scoring are those of dev/improve/fit_task.R (kind "sim"),
# copied here; the objective and posterior objective are kept as well, for
# the pick step.
#   Rscript validation/sparse_joint/dispersion_gate/step1_fit_v7.R TASK_FILE ID STAGE V7_RESULTS OUT
# TASK_FILE columns: family (e | c2), dataset_id, multiple. STAGE holds the
# copied dev/ and validation code (as a repository tree), V7_RESULTS the dev
# results folder whose realdata/noise/v7/sims hold the datasets, OUT the
# results folder: OUT/x<m>/<dataset_id>.rds.
args <- commandArgs(TRUE)
task_file <- args[[1]]; index <- args[[2]]; stage <- args[[3]]; v7_results <- args[[4]]; out <- args[[5]]
Sys.setenv(DEV_REPO = stage, DEV_RESULTS = v7_results, OMP_NUM_THREADS = "1")
source(file.path(stage, "dev", "R", "common.R"))
task <- dev_task_row(task_file, index)
loadNamespace("karst")
m <- as.numeric(task$multiple); beta <- 5e-5 * m
out_file <- file.path(out, sprintf("x%02d", as.integer(m)), paste0(gsub("[^A-Za-z0-9_.-]", "_", task$dataset_id), ".rds"))
if (file.exists(out_file)) quit(save = "no")
key <- list(family = task$family, dataset_id = task$dataset_id, multiple = m, beta = beta,
            karst = as.character(utils::packageVersion("karst")), karst_path = find.package("karst"))
save <- function(rec) dev_save_rds(c(key, rec), out_file)

dev_source_validation()
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
c2_bench_source()
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
d <- readRDS(v7_path("sims", task$family, "datasets", paste0(task$dataset_id, ".rds")))
key <- c(key, list(design = d$design, half = d$half, landscape_id = d$landscape_id, wavelength = d$wavelength,
                   truth_multiple = d$multiple, bottleneck = d$bottleneck, interval = d$interval,
                   depth_t0 = d$depth_t0, depth_t1 = d$depth_t1))
counts <- d$counts
max_cn <- if (d$family == "e") 8L else as.integer(d$max_cn)
arguments <- if (d$family == "e") list(counts, times = c(0, d$interval), beta = beta, birth_rate = 1, min_cn = 1L, max_cn = 8L) else
  list(counts, dt = d$interval, beta = beta, birth_rate = 1, min_cn = 1L, max_cn = max_cn)
nodes <- stt_scoring_nodes(counts, 1L, max_cn, 2000L, 2000L, d$seed)
truth <- if (d$family == "e") sm_true_net_rate(nodes$karyotype, d$truth, d$beta) else
  c2_true_rate(nodes$karyotype, list(truth = readRDS(v7_check_read(d$truth_file))))
set.seed(d$seed)
started <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
seconds <- proc.time()[["elapsed"]] - started
if (inherits(fit, "error")) {
  save(list(status = if (inherits(fit, "karst_domain_coverage_error")) "domain error" else "failed",
            message = conditionMessage(fit), seconds = seconds))
  quit(save = "no")
}
p <- stats::predict(fit, nodes$karyotype)
pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd, stringsAsFactors = FALSE)
score <- stt_score_predictions(nodes, truth, pred, 5)
nt <- score$nodes
keep <- nt$distance <= 2L
compact <- data.frame(distance = as.integer(nt$distance[keep]), truth = nt$truth_relative[keep],
                      pred = nt$pred_relative[keep], sd = nt$pred_sd[keep], available = nt$available[keep])
met <- score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse", "rmse_zero")]
r_of <- function(dd) { k <- compact$distance %in% dd & is.finite(compact$pred) & is.finite(compact$truth)
  if (sum(k) < 3) NA_real_ else suppressWarnings(stats::cor(compact$pred[k], compact$truth[k])) }
save(list(status = "ok", seconds = seconds, metrics = met, r = c(all = r_of(0:2), d0 = r_of(0), d1 = r_of(1), d2 = r_of(2)),
          nodes = compact, objective = fit$diagnostics$objective, posterior_objective = fit$diagnostics$posterior_objective,
          n_major = fit$diagnostics$n_major))
