# The paper's tests in version A (the round as it is: KARST aligned to alfak,
# on alfak's karyotypes of each minobs) and version B (alfak aligned to KARST,
# ALIGN=karst, R/align.R: every landscape on the input's extended list E)
# side by side (docs/align_karst_6ae2a24.md).
#
#   Rscript scripts/align_report.R
#
# Reads BENCH_RESULTS/analysis/ (A) and BENCH_RESULTS/align_karst/analysis/
# (B) and B's task records (align_karst/align_tasks.Rds,
# scripts/align_summary.R). Writes BENCH_RESULTS/analysis_align_compare/:
# compare_<test>.tsv with a column `version` ("A", "B"), and lists.tsv (B's
# list sizes, the share of each alfak landscape that Kriging fills, the
# Kriging's agreement with alfak's mean on its own list, times). In B the
# landscape accuracy columns rfq, rnn and rd2 hold the distance classes 0, 1
# and 2 to the observed karyotypes, and the frequent karyotypes of every
# landscape are KARST's majors.

for (f in list.files("R", full.names = TRUE)) source(f)
Sys.setenv(ALIGN = "")
paths <- bench_paths()
A <- file.path(paths$results, "analysis")
B <- file.path(paths$results, "align_karst", "analysis")
out <- file.path(paths$results, "analysis_align_compare")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tab <- function(d, name) utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE,
                                                  row.names = FALSE)
both <- function(rel, name) {
  read <- function(root) {
    p <- file.path(root, rel)
    if (!file.exists(p)) return(NULL)
    utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  }
  a <- read(A)
  b <- read(B)
  if (is.null(a) || is.null(b)) {
    cat("missing:", rel, if (is.null(a)) "(A)", if (is.null(b)) "(B)", "\n")
    return(invisible(NULL))
  }
  d <- bind_rows_fill(list(cbind(version = "A", a), cbind(version = "B", b)))
  write_tab(d, name)
  invisible(d)
}
both("abm/supp2a_landscape_accuracy.tsv", "supp2a")
both("abm/paired_landscape.tsv", "paired_landscape")
both("abm/paired_forecast.tsv", "paired_forecast_abm")
both("abm/fig1g_beat_baseline.tsv", "fig1g")
both("abm/supp3_misseg_pearson.tsv", "supp3")
both("abm/supp4b_cv_flags_accuracy.tsv", "supp4b")
both("abm/supp4c_selected_by_cv.tsv", "supp4c")
both("salehi/paired_forecast.tsv", "paired_forecast_salehi")
both("salehi/fig2c_angle_test.tsv", "fig2c")
both("salehi/fig2d_wasserstein.tsv", "fig2d")
both("salehi/fig2f_beat_baseline_eval_pooled.tsv", "fig2f")
both("salehi/fig3ef_novel_karyotypes.tsv", "fig3ef")
both("salehi/fig3_paired.tsv", "fig3_paired")
both("salehi/fit_failures.tsv", "salehi_fit_failures")
both("fig456/conclusions.tsv", "fig456")
rec <- file.path(paths$results, "align_karst", "align_tasks.Rds")
if (file.exists(rec)) {
  r <- readRDS(rec)
  f <- r$fits
  f$kind <- ifelse(f$method == "alfak", "alfak", "KARST")
  q <- function(x, p) round(stats::quantile(x, p, na.rm = TRUE, names = FALSE), 3)
  lists <- do.call(rbind, lapply(split(f, f[c("source", "kind")], drop = TRUE), function(x) data.frame(
    source = x$source[1], fit = x$kind[1], landscapes = nrow(x), E_median = stats::median(x$E), E_max = max(x$E),
    kriged_share_median = if (x$kind[1] == "alfak") q(x$kriged / x$E, 0.5) else NA,
    kriging_agreement_median = if (x$kind[1] == "alfak") q(x$kriging_agreement, 0.5) else NA,
    kriging_agreement_q10 = if (x$kind[1] == "alfak") q(x$kriging_agreement, 0.1) else NA,
    landscape_seconds_median = q(x$landscape_seconds, 0.5), forecast_seconds_median = q(x$forecast_seconds, 0.5),
    forecast_seconds_max = q(x$forecast_seconds, 1), stringsAsFactors = FALSE)))
  write_tab(lists, "lists")
}
cat("written to", out, "\n")
