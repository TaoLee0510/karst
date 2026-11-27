#!/usr/bin/env Rscript
# C2 H and M: base fits at beta 5e-5 (ce_base) against fits at the simulated
# rate beta_sim (ce_bsim; H 5e-4, M 1.5e-4). Medians over datasets of the
# RMSE ratio to predicting zero, bias and Spearman, by distance.
#   Rscript dev/realdata/CE/c2_bsim_compare.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
root <- file.path(dev_results, "realdata", "CE", "c2", "fits")
rows <- do.call(rbind, lapply(c("ce_base", "ce_bsim"), function(m) do.call(rbind, lapply(c("H", "M"), function(g) {
  files <- list.files(file.path(root, m, g), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)
    if (!identical(r$status, "ok")) return(data.frame(method = m, regime = g, dataset_id = r$dataset_id, stratum = NA,
                                                      ratio = NA, bias = NA, spearman = NA, coverage = NA, ok = FALSE))
    x <- r$near_metrics[r$near_metrics$stratum %in% c("d0", "d1", "d2"), ]
    data.frame(method = m, regime = g, dataset_id = r$dataset_id, stratum = x$stratum, ratio = x$rmse / x$rmse_zero,
               bias = x$bias, spearman = x$spearman, coverage = x$coverage95, ok = TRUE)
  }))
}))))
cat("fits:\n"); print(stats::aggregate(ok ~ method + regime, unique(rows[c("method", "regime", "dataset_id", "ok")]), sum))
ok <- rows[rows$ok, ]
med <- stats::aggregate(cbind(ratio, bias, spearman, coverage) ~ regime + stratum + method, ok, stats::median)
options(width = 200)
print(med[order(med$regime, med$stratum, med$method), ], digits = 3, row.names = FALSE)
# paired: share of datasets where beta_sim is closer to the truth
w <- stats::reshape(ok[c("method", "regime", "dataset_id", "stratum", "ratio")], idvar = c("regime", "dataset_id", "stratum"),
                    timevar = "method", direction = "wide")
w$better <- w$ratio.ce_bsim < w$ratio.ce_base
print(stats::aggregate(better ~ regime + stratum, w, mean), digits = 3)
dev_save_rds(list(rows = rows, medians = med), file.path(dev_results, "realdata", "CE", "c2_bsim_compare.rds"))
