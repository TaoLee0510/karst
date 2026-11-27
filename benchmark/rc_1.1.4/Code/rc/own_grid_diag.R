# Part 2 of PLAN.md, follow-up to own_grid_summary.R: why the error ratio moves with the
# multiple. Core designs, distances 0-2, same aggregation (designs within landscape, then
# landscapes): bias, error sd, truth sd, slope (cov(pred, truth) / var(truth); 1 = right
# scale), share of |pred| > 1, robust error ratio (median abs error / null median abs
# error). Also the failure messages of KARST fits at every multiple.
#   Rscript own_grid_diag.R TASK_FILE ID OUT_DIR
args <- commandArgs(TRUE)
out <- args[[3]]
sdir <- file.path(out, "summary")
m <- utils::read.csv(file.path(out, "metrics_heldout_grid.csv"), stringsAsFactors = FALSE)
m <- m[m$role == "core" & m$stratum %in% c("d0", "d1", "d2"), ]
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
cols <- c("rmse", "rmse_zero", "bias", "error_sd", "truth_sd", "slope", "extreme_fraction", "median_abs_error", "median_abs_error_zero")
cell <- stats::aggregate(m[cols], m[c("method", "multiple", "wavelength_class", "landscape_id", "scenario", "stratum")], mean_na)
land <- stats::aggregate(cell[cols], cell[c("method", "multiple", "wavelength_class", "landscape_id", "stratum")], mean_na)
pool <- stats::aggregate(land[cols], land[c("method", "multiple", "stratum")], mean_na)
pool$ratio <- pool$rmse / pool$rmse_zero
pool$robust_ratio <- pool$median_abs_error / pool$median_abs_error_zero
pool <- pool[order(pool$method, pool$stratum, pool$multiple), ]
utils::write.table(pool, file.path(sdir, "core_scale_diag.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(pool[c("method", "stratum", "multiple", "ratio", "robust_ratio", "bias", "error_sd", "truth_sd", "slope", "extreme_fraction")], digits = 3, row.names = FALSE)

# per-landscape robust ratio: share of landscapes where KARST beats ALFA-K at the same multiple
lr <- land; lr$robust <- lr$median_abs_error / lr$median_abs_error_zero; lr$ratio <- lr$rmse / lr$rmse_zero
k <- merge(lr[lr$method == "karst_sparse_joint_balance", c("multiple", "landscape_id", "wavelength_class", "stratum", "ratio", "robust")],
           lr[lr$method == "alfak_minobs10", c("multiple", "landscape_id", "stratum", "ratio", "robust")],
           by = c("multiple", "landscape_id", "stratum"), suffixes = c(".karst", ".alfak"))
share <- stats::aggregate(cbind(ratio = ratio.karst < ratio.alfak, robust = robust.karst < robust.alfak) ~ multiple + stratum + wavelength_class, k, mean)
utils::write.table(share, file.path(sdir, "core_karst_better_share_by_wave.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats::reshape(share[share$stratum == "d1", c("multiple", "wavelength_class", "ratio")], idvar = "wavelength_class", timevar = "multiple", direction = "wide"), digits = 2, row.names = FALSE)

# failures of KARST fits
fit_dir <- function(k) if (k == 1) Sys.getenv("OWN_ROUND_FITS") else file.path(out, sprintf("x%02d", k), "fits")  # rc freeze
f <- unique(utils::read.csv(file.path(out, "metrics_heldout_grid.csv"), stringsAsFactors = FALSE)[c("dataset_id", "method", "multiple", "status", "scenario", "wavelength_class")])
f <- f[f$status != "ok" & f$method == "karst_sparse_joint_balance", ]
f$detail <- vapply(seq_len(nrow(f)), function(i) {
  d <- readRDS(file.path(fit_dir(f$multiple[i]), sprintf("%s__%s.rds", f$dataset_id[i], f$method[i])))$detail
  substr(gsub("[[:space:]]+", " ", paste(d, collapse = " ")), 1, 160)
}, "")
utils::write.table(f, file.path(sdir, "karst_failures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(f$multiple, f$scenario))
print(table(f$multiple, f$wavelength_class))
print(sort(table(sub("[0-9.e+-]+", "#", f$detail)), decreasing = TRUE)[1:10])
