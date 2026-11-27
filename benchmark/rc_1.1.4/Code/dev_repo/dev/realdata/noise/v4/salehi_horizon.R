#!/usr/bin/env Rscript
# Step 4 of dev/realdata/noise/PLAN_v4.md: one forecast of one T0/T1/T2 fit of
# a Salehi trajectory with validation/sparse_joint/salehi_horizon/forecast_score.R,
# unchanged (forecast rate pm and nb from the task; common support = the
# trajectory's three fits under the v4 root).
#   Rscript dev/realdata/noise/v4/salehi_horizon.R TASK_FILE TASK_ID
# Task columns: id, method, minobs, design, root, support_roots, pm, nb.
# Output: DEV_RESULTS/realdata/noise/v4/salehi/horizon/
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
rd_load_karst()
.horizon_args <- c(normalizePath(commandArgs(TRUE)[[1]]), commandArgs(TRUE)[[2]],
                   file.path(dev_results, "realdata", "noise", "v4", "salehi", "horizon"))
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "v4", "salehi", "bench_unused"))
commandArgs <- function(trailingOnly = FALSE) .horizon_args
source(file.path(dev_repo, "validation", "sparse_joint", "salehi_horizon", "forecast_score.R"))
