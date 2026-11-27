# How much information do the Salehi inputs carry? (descriptive, 2026-11-26)
# Per window: per-karyotype information (karyotypes seen at both time points),
# per-chromosome mean copy-number shifts t0 -> t1 (z), and whether continuing
# the training-period centroid shift predicts the daughter's direction (the
# Fig. 2c angle of a "trend continuation" forecast).
here <- "/Users/4482173/Documents/GitHub/karst/patch/outside_rank"
source(file.path(here, "R", "common.R"))
root <- "/Users/4482173/Documents/GitHub/karst/benchmark/paper"; data_root <- "/Users/4482173/Documents/GitHub/ALFA-K/data"
pairs <- utils::read.delim(file.path(root, "results_6ae2a24/analysis/fig2c_two_timepoints/shared_full_parent/daughter_angles.tsv"), stringsAsFactors = FALSE)
meta <- utils::read.csv(file.path(data_root, "raw/salehi/metadata.csv"), stringsAsFactors = FALSE)
inputs_dir <- file.path(data_root, "processed/salehi/alfak_inputs"); inputs <- list.files(inputs_dir)
s3 <- utils::read.delim(file.path(here, "results/step3/salehi_rows.tsv"), stringsAsFactors = FALSE)
cmean <- function(v) od_centre(v)
csd <- function(v) { cn <- od_parse(names(v)); p <- v / sum(v); sqrt(pmax(colSums(p * cn^2) - colSums(p * cn)^2, 0)) }
rows <- list()
for (w in unique(s3$window)) {
  x <- readRDS(file.path(inputs_dir, paste0(w, ".Rds")))$x
  x <- as.matrix(x[, c(1L, ncol(x)), drop = FALSE]); x <- x[rowSums(x) > 0, , drop = FALSE]
  v0 <- x[x[, 1] > 0, 1]; v1 <- x[x[, 2] > 0, 2]; n0 <- sum(v0); n1 <- sum(v1)
  both3 <- rownames(x)[x[, 1] >= 3 & x[, 2] >= 3]
  fit <- readRDS(file.path(here, "data", "salehi_fits", paste0(w, ".Rds")))
  maj <- fit$summary$karyotype[fit$summary$major]
  shift <- cmean(v1) - cmean(v0)
  z <- shift / sqrt(csd(v0)^2 / n0 + csd(v1)^2 / n1)
  for (j in which(pairs$window == w)) {
    fj <- paste(meta[meta$uid == pairs$daughter_uid[j], c("datasetname", "timepoint")], collapse = "_")
    tru <- readRDS(file.path(inputs_dir, inputs[grepl(fj, inputs)][1]))$x
    d <- stats::setNames(as.numeric(tru[, ncol(tru)]), rownames(tru)); d <- d[d > 0]
    par <- d; par <- v1  # parent = last training sample (as the shared full parent)
    obs <- cmean(d) - cmean(v1)
    sel <- abs(z) > 3
    rows[[length(rows) + 1L]] <- data.frame(window = w, daughter = pairs$daughter[j],
      family = if (grepl("^SA906", w)) "SA906" else if (grepl("^SA609|^SA000", w)) "SA609" else "other",
      cells_t0 = n0, cells_t1 = n1, karyotypes_t1 = length(v1), singleton_share_t1 = sum(v1[v1 == 1]) / n1,
      both3 = length(both3), both3_share_t1 = sum(x[both3, 2]) / n1, majors = length(maj), major_share_t1 = sum(v1[names(v1) %in% maj]) / n1,
      chr_shift_z3 = sum(sel), angle_trend = od_angle(shift, obs),
      angle_trend_sig = if (any(sel)) od_angle(ifelse(sel, shift, 0), obs) else NA_real_,
      angle_karst = s3$angle[s3$window == w & s3$daughter == pairs$daughter[j] & s3$arm == "karst_current"][1],
      angle_alfak = s3$angle[s3$window == w & s3$daughter == pairs$daughter[j] & s3$arm == "alfak"][1], stringsAsFactors = FALSE)
  }
}
r <- do.call(rbind, rows)
options(width = 200)
cat("daughters", nrow(r), " windows", length(unique(r$window)), "\n\n== per window (medians by family)\n")
num <- c("cells_t1", "karyotypes_t1", "singleton_share_t1", "majors", "major_share_t1", "both3", "both3_share_t1", "chr_shift_z3")
print(stats::aggregate(r[!duplicated(r$window), num], list(family = r$family[!duplicated(r$window)]), stats::median), digits = 3, row.names = FALSE)
cat("\n== Fig. 2c angle (median) of: continuing the training shift / only chromosomes with |z| > 3 / KARST / ALFA-K\n")
print(stats::aggregate(cbind(angle_trend, angle_trend_sig, angle_karst, angle_alfak) ~ family, r, stats::median, na.action = na.pass), digits = 3, row.names = FALSE)
cat("all:", sapply(r[, c("angle_trend", "angle_karst", "angle_alfak")], stats::median), "\n")
cat("trend continuation smaller angle than KARST in", sum(r$angle_trend < r$angle_karst), "of", nrow(r), "; than ALFA-K in", sum(r$angle_trend < r$angle_alfak), "\n")
cat("angle < 90 deg (direction informative): trend", sum(r$angle_trend < 90), " KARST", sum(r$angle_karst < 90), " ALFA-K", sum(r$angle_alfak < 90), "\n")
utils::write.table(r, "/Users/4482173/Documents/GitHub/karst/patch/forecast_rate/diag/info_check_rows.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
