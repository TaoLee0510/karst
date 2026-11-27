#!/usr/bin/env Rscript
# decomp/PLAN.md: one (L, F) forecast of one Salehi trajectory with
# validation/sparse_joint/salehi_horizon/forecast_score.R, unchanged. The fit is
# the alfak output of rate L (column method), the forecast rate is column pm;
# rows go to DEV_RESULTS/realdata/beta/decomp/<out>/ (out = F's label).
#   Rscript dev/realdata/beta/decomp/decomp_horizon.R TASK_FILE TASK_ID
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
rd_load_karst()
.task_file <- normalizePath(commandArgs(TRUE)[[1]])
.task_id <- commandArgs(TRUE)[[2]]
.task_out <- utils::read.delim(.task_file, stringsAsFactors = FALSE)$out[as.integer(.task_id)]
.horizon_args <- c(.task_file, .task_id, file.path(dev_results, "realdata", "beta", "decomp", .task_out))
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "beta", "salehi", "bench_unused"))
commandArgs <- function(trailingOnly = FALSE) .horizon_args
source(file.path(dev_repo, "validation", "sparse_joint", "salehi_horizon", "forecast_score.R"))
