#!/usr/bin/env Rscript
# Copy the accepted setting of a regime (dev/realdata/PLAN.md, phase C) to
# c2/regimes/<regime>/ (datasets with `regime` set, nodes) and append its
# ce_base fits to tasks/c2_fit_regimes.tsv.
#   Rscript dev/realdata/CE/c2_promote.R REGIME SETTING
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
a <- commandArgs(TRUE)
regime <- a[[1L]]; setting <- a[[2L]]
from <- file.path(c2_root, "grid", setting)
to <- file.path(c2_root, "regimes", regime)
ids <- sub("\\.rds$", "", list.files(file.path(from, "datasets")))
for (id in ids) {
  d <- readRDS(file.path(from, "datasets", paste0(id, ".rds")))
  d$regime <- regime
  dev_save_rds(d, file.path(to, "datasets", paste0(id, ".rds")))
  dir.create(file.path(to, "nodes"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(from, "nodes", paste0(id, ".rds")), file.path(to, "nodes", paste0(id, ".rds")), overwrite = TRUE)
}
writeLines(setting, file.path(to, "SETTING"))
cat(regime, "<-", setting, ":", length(ids), "datasets\n")
