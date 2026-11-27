#!/usr/bin/env Rscript
# Analysis of dev/realdata/noise/PLAN_v3.md on one chunk of the new
# simulations: R/envelope.R unchanged, on the first and last time point, after
# an optional injection (variant none | dm5 | dm20 | reg05 | reg10), as v2_task.R.
#   Rscript dev/realdata/noise/v3/v3_task.R TASK_FILE ID
# Task columns: kind (eset | c2), set (e1-e4 or a C2 setting), variant, chunk, chunks.
# Output: DEV_RESULTS/realdata/noise/v3/analysis/<kind>_<set>_<variant>/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
dev_source_validation()
sims <- file.path(dev_results, "realdata", "noise", "v3", "sims")
out_dir <- file.path(dev_results, "realdata", "noise", "v3", "analysis", sprintf("%s_%s_%s", task$kind, task$set, task$variant))
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
files <- if (task$kind == "eset") {
  sort(Sys.glob(file.path(sims, task$set, "*", "datasets", "*.rds")))
} else sort(Sys.glob(file.path(sims, "c2", "grid", task$set, "datasets", "*.rds")))
files <- files[(seq_along(files) - 1L) %% task$chunks + 1L == task$chunk]
for (f in files) {
  d <- readRDS(f)
  out_file <- file.path(out_dir, paste0(d$dataset_id, ".rds"))
  if (file.exists(out_file)) next
  if (task$kind == "eset") {
    x <- mt_counts(d, "ends"); dt <- diff(d$times[c(1L, length(d$times))]); truth <- d$beta
    depth <- d$depth; bottleneck <- d$bottleneck
  } else {
    x <- d$counts; dt <- d$interval; truth <- d$truth$beta_sim; depth <- d$depth; bottleneck <- d$truth$bottleneck
  }
  x0 <- col(x, 1); x1 <- col(x, ncol(x))
  if (task$variant != "none") {
    set.seed(seed_of(paste(d$dataset_id, task$variant)))
    kind <- if (startsWith(task$variant, "dm")) "dm" else "region"
    level <- switch(task$variant, dm5 = 5, dm20 = 20, reg05 = 0.5, reg10 = 1.0)
    fam <- if (kind == "region") env_families(x0[x0 > 0], x1[x1 > 0]) else NULL
    x0 <- env_inject(x0[x0 > 0], kind, level, fam); x1 <- env_inject(x1[x1 > 0], kind, level, fam)
  }
  started <- proc.time()[["elapsed"]]
  res <- tryCatch(env_analyse(x0, x1, dt, 1, seed = seed_of(d$dataset_id)), error = function(e) e)
  rec <- if (inherits(res, "error")) list(status = "failed", message = conditionMessage(res)) else c(list(status = "ok"), res)
  dev_save_rds(c(list(kind = task$kind, set = task$set, dataset_id = d$dataset_id, variant = task$variant, truth = truth,
                      depth = depth, bottleneck = bottleneck, dt = dt, seconds = proc.time()[["elapsed"]] - started), rec),
               out_file)
}
cat("done", length(files), "\n")
