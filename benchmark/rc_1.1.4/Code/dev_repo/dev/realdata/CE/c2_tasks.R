#!/usr/bin/env Rscript
# Task tables of a C2 round: c2_make_<round>.tsv (setting x landscape) and
# c2_fit_<round>.tsv (ce_base fits of every dataset).
#   Rscript dev/realdata/CE/c2_tasks.R ROUND
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
round <- commandArgs(TRUE)[[1L]]
s <- c2_settings()
s <- s[as.character(s$round) == round, ]
make <- expand.grid(landscape_id = c2_landscapes$landscape_id, setting = s$setting, stringsAsFactors = FALSE)[c("setting", "landscape_id")]
fit <- expand.grid(interval = c2_intervals, landscape_id = c2_landscapes$landscape_id, setting = s$setting,
                   stringsAsFactors = FALSE)
fit <- data.frame(method = "ce_base", lib = "karst_cd342db", regime = fit$setting,
                  dataset_id = c2_dataset_id(fit$landscape_id, fit$interval), extra = "", stringsAsFactors = FALSE)
dir <- file.path(dev_repo, "dev", "realdata", "CE", "tasks")
utils::write.table(make, file.path(dir, paste0("c2_make_", round, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(fit, file.path(dir, paste0("c2_fit_", round, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(make), "trajectories,", nrow(fit), "fits\n")
