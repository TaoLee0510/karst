#!/usr/bin/env Rscript
# ce_base fit tasks of the accepted regimes: tasks/c2_fit_regimes.tsv.
#   Rscript dev/realdata/CE/c2_regime_tasks.R REGIME...
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
regimes <- commandArgs(TRUE)
rows <- do.call(rbind, lapply(regimes, function(r) {
  ids <- sub("\\.rds$", "", list.files(file.path(c2_root, "regimes", r, "datasets")))
  data.frame(method = "ce_base", lib = "karst_cd342db", regime = r, dataset_id = ids, extra = "", stringsAsFactors = FALSE)
}))
path <- file.path(dev_repo, "dev", "realdata", "CE", "tasks", paste0("c2_fit_regimes_", paste(regimes, collapse = ""), ".tsv"))
utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(rows), "tasks ->", path, "\n")
