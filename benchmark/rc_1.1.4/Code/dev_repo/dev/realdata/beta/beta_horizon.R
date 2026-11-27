#!/usr/bin/env Rscript
# E4 and S1 of dev/realdata/beta/PLAN.md: validation/sparse_joint/salehi_horizon/forecast_score.R,
# unchanged, on the beta-grid Salehi fits of beta_salehi.R (minobs 10; the
# forecast's missegregation rate is the fit's beta, column pm). Run after every
# fit of the grid exists, so that the common support is the same for all rates.
#   Rscript dev/realdata/beta/beta_horizon.R TASK_FILE TASK_ID
# Rows go to DEV_RESULTS/realdata/beta/horizon/.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
rd_load_karst()
.horizon_args <- c(normalizePath(commandArgs(TRUE)[[1]]), commandArgs(TRUE)[[2]],
                   file.path(dev_results, "realdata", "beta", "horizon"))
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "beta", "salehi", "bench_unused"))
commandArgs <- function(trailingOnly = FALSE) .horizon_args
source(file.path(dev_repo, "validation", "sparse_joint", "salehi_horizon", "forecast_score.R"))
