#!/usr/bin/env Rscript
# D3 W3: validation/sparse_joint/salehi_horizon/forecast_score.R, unchanged,
# on the D3 fits of dev/realdata/CE/d3_salehi.R (as dev/realdata/B/horizon_task.R).
#   Rscript dev/realdata/CE/d3_horizon.R TASK_FILE TASK_ID
# Rows go to DEV_RESULTS/realdata/CE/d3/horizon/.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
rd_load_karst()
.horizon_args <- c(normalizePath(commandArgs(TRUE)[[1]]), commandArgs(TRUE)[[2]], file.path(rd_results, "CE", "d3", "horizon"))
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(rd_results, "CE", "d3", "bench_unused"))
# forecast_score.R reads its arguments with commandArgs(TRUE)
commandArgs <- function(trailingOnly = FALSE) .horizon_args
source(file.path(dev_repo, "validation", "sparse_joint", "salehi_horizon", "forecast_score.R"))
