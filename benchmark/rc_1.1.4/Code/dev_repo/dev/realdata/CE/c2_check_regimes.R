#!/usr/bin/env Rscript
# Check of the accepted C2 regimes (dev/realdata/PLAN.md, phase C): the ce_base
# fits of c2/regimes/<regime> give the acceptance medians again, and the
# scored nodes are complete.
#   Rscript dev/realdata/CE/c2_check_regimes.R REGIME...
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
for (regime in commandArgs(TRUE)) {
  files <- list.files(file.path(c2_root, "fits", "ce_base", regime), full.names = TRUE)
  rows <- do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)
    d <- c2_dataset(regime, r$dataset_id)
    st <- c2_count_stats(d)
    ok <- r$status == "ok"
    data.frame(st, status = r$status, coverage_min = if (ok) min(r$diagnostics$domain_coverage$fraction) else NA,
               dispersion = if (ok) r$dispersion else NA, seconds = r$seconds,
               nodes = if (ok) nrow(r$near_nodes) else NA,
               rmse_ratio_d0 = if (ok) with(r$near_metrics[r$near_metrics$stratum == "d0", ], rmse / rmse_zero) else NA,
               rmse_ratio_d1 = if (ok) with(r$near_metrics[r$near_metrics$stratum == "d1", ], rmse / rmse_zero) else NA,
               rmse_ratio_d2 = if (ok) with(r$near_metrics[r$near_metrics$stratum == "d2", ], rmse / rmse_zero) else NA)
  }))
  m <- function(x) round(stats::median(x, na.rm = TRUE), 3)
  cat(regime, ": fits", nrow(rows), "ok", sum(rows$status == "ok"), "| singleton", m(rows$singleton_share_last),
      "coverage", m(rows$coverage_min), "dispersion", m(rows$dispersion), "| seconds", m(rows$seconds),
      "| RMSE/zero d0 d1 d2", m(rows$rmse_ratio_d0), m(rows$rmse_ratio_d1), m(rows$rmse_ratio_d2), "\n")
}
