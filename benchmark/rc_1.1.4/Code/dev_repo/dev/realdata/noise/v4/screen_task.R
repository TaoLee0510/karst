#!/usr/bin/env Rscript
# Step 1 of dev/realdata/noise/PLAN_v4.md on one chunk of the v3 simulations
# (300 and 1000 cells per time point): optional injection as v3_task.R (same
# seeds), env_analyse() (R/envelope.R, unchanged), env_by_beta() and the class
# with beta_hat and phi_hat (R/envelope_v4.R). The counts the fits use (after
# the injection) are stored with the record.
#   Rscript dev/realdata/noise/v4/screen_task.R TASK_FILE ID
# Task columns: kind (eset | c2), set (e1-e4 or a C2 setting), variant (none | dm20 | reg10), chunk, chunks.
# Output: DEV_RESULTS/realdata/noise/v4/screen/<kind>_<set>_<variant>/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
dev_source_validation()
sims <- file.path(dev_results, "realdata", "noise", "v3", "sims")
out_dir <- file.path(dev_results, "realdata", "noise", "v4", "screen", sprintf("%s_%s_%s", task$kind, task$set, task$variant))
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
files <- if (task$kind == "eset") {
  sort(Sys.glob(file.path(sims, task$set, "*", "datasets", "*.rds")))
} else sort(Sys.glob(file.path(sims, "c2", "grid", task$set, "datasets", "*.rds")))
files <- files[grepl("_n(0300|1000)", basename(files))]
files <- files[(seq_along(files) - 1L) %% task$chunks + 1L == task$chunk]
for (f in files) {
  d <- readRDS(f)
  out_file <- file.path(out_dir, paste0(d$dataset_id, ".rds"))
  if (file.exists(out_file)) next
  if (task$kind == "eset") {
    x <- mt_counts(d, "ends"); dt <- diff(d$times[c(1L, length(d$times))]); truth <- d$beta
    bottleneck <- d$bottleneck
  } else {
    x <- d$counts; dt <- d$interval; truth <- d$truth$beta_sim; bottleneck <- d$truth$bottleneck
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
  rec <- tryCatch({
    env <- env_analyse(x0, x1, dt, 1, seed = seed_of(d$dataset_id))
    bb <- env_by_beta(x0, x1, dt, 1)
    cl <- env_v4_class(env, bb)
    list(status = "ok", summary = env$summary, by_beta = bb, class = cl$class, beta_hat = cl$beta_hat,
         phi_hat = cl$phi_hat)
  }, error = function(e) list(status = "failed", message = conditionMessage(e)))
  dev_save_rds(c(list(kind = task$kind, set = task$set, dataset_id = d$dataset_id, landscape_id = d$landscape_id,
                      half = v4_half(d$landscape_id), variant = task$variant, truth = truth, depth = d$depth,
                      bottleneck = bottleneck, dt = dt, source = f, x0 = x0, x1 = x1,
                      seconds = proc.time()[["elapsed"]] - started), rec),
               out_file)
}
cat("done", length(files), "\n")
