# The multi-timepoint session's criterion for karst's "auto" arm, registered
# before any E4 result, checked on BENCH_RESULTS/analysis/summary/
# summary_accuracy.tsv (scripts/summary_comparison.R: every method on the same
# inputs of a group, so auto and anchored are paired). Auto passes if:
#   ABM and missegregation, all time points, windows 4 and 8: the median landscape
#     Pearson r >= anchored - 0.01, and (ABM) the share of forecasts beating
#     the baseline at passage 10 >= anchored - 0.02
#   Salehi, all time points: the share of forecasts beating the parent and the
#     share of novel-karyotype regressions where fitness predicts emergence
#     >= anchored - 0.02
#   Salehi sweeps (anchored left out time points): the share of forecasts beating
#     the parent >= anchored
# Forecast lines are checked with the Euclidean and the Wasserstein distance,
# each as a line of its own.
#
#   Rscript scripts/auto_criterion.R
#
# Writes BENCH_RESULTS/analysis/summary/auto_criterion.tsv: line, panel, group,
# pairs, auto, anchored, difference, allowed margin, pass.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
dir <- file.path(paths$results, "analysis", "summary")
acc <- utils::read.delim(file.path(dir, "summary_accuracy.tsv"), stringsAsFactors = FALSE)
picks_path <- file.path(paths$results, "tables", "auto_picks.tsv")
picks <- if (file.exists(picks_path)) utils::read.delim(picks_path) else NULL

line <- function(label, panel, x, margin) {
  a <- acc[acc$panel == panel & acc$x == x & acc$method == "karst_auto", ]
  b <- acc[acc$panel == panel & acc$x == x & acc$method == "karst", ]
  if (nrow(a) != 1L || nrow(b) != 1L) {
    return(data.frame(line = label, panel = panel, group = x, pairs = NA, auto = NA, anchored = NA, difference = NA,
                      margin = margin, pass = NA))
  }
  d <- a$value - b$value
  data.frame(line = label, panel = panel, group = x, pairs = a$pairs, auto = a$value, anchored = b$value,
             difference = round(d, 3), margin = margin, pass = d >= -margin - 1e-9)
}
rows <- list()
for (src in c("ABM", "Missegregation")) for (w in c(4, 8)) {
  rows[[length(rows) + 1]] <- line(paste(src, "landscape r, window", w), paste0(src, ": landscape Pearson r"),
                                   paste0("all time points window ", w), 0.01)
}
for (m in c("euclidean", "wasserstein")) for (w in c(4, 8)) {
  rows[[length(rows) + 1]] <- line(paste("ABM forecast beats baseline (", m, "), window", w),
                                   paste0("ABM: forecast beats baseline, passage 10 (", m, ")"),
                                   paste0("all time points window ", w), 0.02)
}
for (m in c("euclidean", "wasserstein")) {
  panel <- paste0("Salehi: share of daughters where the forecast beats the parent (", m, ")")
  rows[[length(rows) + 1]] <- line(paste("Salehi forecast beats parent (", m, ")"), panel,
                                   "all time points (3-10 passages)", 0.02)
  rows[[length(rows) + 1]] <- line(paste("Salehi sweeps, forecast beats parent (", m, ")"), panel,
                                   "all time points (sweeps)", 0)
}
rows[[length(rows) + 1]] <- line("Salehi novel karyotypes", "Salehi: novel karyotypes, fitness predicts emergence",
                                 "all time points (3-10 passages)", 0.02)
crit <- do.call(rbind, rows)
if (!is.null(picks)) {
  crit$share_cut <- picks$share_cut[1]
  crit$weight_cut <- picks$weight_cut[1]
}
utils::write.table(crit, file.path(dir, "auto_criterion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(crit, row.names = FALSE)
cat("\nAll lines pass:", all(crit$pass), " (lines not evaluated:", sum(is.na(crit$pass)), ")\n")
