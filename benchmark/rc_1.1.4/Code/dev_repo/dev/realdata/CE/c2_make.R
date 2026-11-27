#!/usr/bin/env Rscript
# C2 trajectories and datasets (dev/realdata/PLAN.md, phase C): one task per
# setting and landscape.
#   Rscript dev/realdata/CE/c2_make.R TASK_FILE ID     (columns: setting, landscape_id)
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1L]], args[[2L]])
setting <- c2_setting(task$setting)
ids <- c2_make_datasets(setting, task$landscape_id)
truth <- c2_simulate(setting, task$landscape_id)
cat(task$setting, task$landscape_id, "datasets", paste(ids, collapse = " "), "| seconds", round(truth$seconds), "\n")
