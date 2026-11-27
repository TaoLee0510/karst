# Criterion for karst's missegregation rate picked per input (round ac36fd1,
# beta = "auto"), fixed before any scoring: the 13 lines of
# scripts/auto_criterion.R with the same margins, for each new arm against
# its fixed-rate counterpart on the same inputs (summary_accuracy.tsv of
# scripts/summary_comparison.R: every method on the same inputs of a group,
# so the arms are paired):
#   all time points  karst_beta_auto against karst (anchored, 5e-5)
#   2 time points    karst_beta_auto against karst (first and last, 5e-5),
#                    the same tests on the two-time-point inputs: ABM and
#                    missegregation windows 4 and 8 (w*_ends), every Salehi
#                    trajectory, and the sweep trajectories' two-time-point
#                    fits
# The counterpart on the missegregation simulations is karst_beta_fixed
# (5e-5), not karst, whose fits there take the true rate. The lines:
#   ABM and missegregation, windows 4 and 8: the median landscape Pearson r
#     >= counterpart - 0.01; ABM: the share of forecasts beating the baseline
#     at passage 10 >= counterpart - 0.02 (Euclidean, Wasserstein)
#   Salehi: the share of forecasts beating the parent (Euclidean,
#     Wasserstein) and the share of novel-karyotype regressions where fitness
#     predicts emergence >= counterpart - 0.02
#   Salehi sweeps: the share of forecasts beating the parent >= counterpart
# Reported besides, not part of the criterion: window 2 (two time points
# only), and on the missegregation simulations karst_beta_auto and 5e-5
# against the true rate (karst). Fig. 4d is reported by fig456_summary.R and
# is not a criterion.
#
#   Rscript scripts/beta_criterion.R
#
# Writes BENCH_RESULTS/analysis/summary/beta_criterion.tsv: arm, line, panel,
# group, method, reference, pairs, value, reference value, difference,
# allowed margin, pass, criterion (FALSE: reported only); the shares of
# paired inputs where the arm is better and worse, and the paired test.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
dir <- file.path(paths$results, "analysis", "summary")
acc <- utils::read.delim(file.path(dir, "summary_accuracy.tsv"), stringsAsFactors = FALSE)

line <- function(arm, label, panel, x, margin, method = "karst_beta_auto", reference = NULL, criterion = TRUE) {
  a <- acc[acc$panel == panel & acc$x == x & acc$method == method, ]
  if (is.null(reference)) reference <- if (startsWith(panel, "Missegregation")) "karst_beta_fixed" else "karst"
  b <- acc[acc$panel == panel & acc$x == x & acc$method == reference, ]
  empty <- data.frame(arm = arm, line = label, panel = panel, group = x, method = method, reference = reference,
                      pairs = NA, value = NA, reference_value = NA, difference = NA, margin = margin, pass = NA,
                      criterion = criterion, better = NA, worse = NA, p = NA)
  if (nrow(a) != 1L || nrow(b) != 1L) return(empty)
  d <- a$value - b$value
  # the paired shares and test against the arm's own reference (summary_comparison.R)
  same <- identical(a$reference, reference)
  data.frame(arm = arm, line = label, panel = panel, group = x, method = method, reference = reference,
             pairs = a$pairs, value = a$value, reference_value = b$value, difference = round(d, 3), margin = margin,
             pass = d >= -margin - 1e-9, criterion = criterion,
             better = if (same) a$better_than_anchored else NA, worse = if (same) a$worse_than_anchored else NA,
             p = if (same) a$p_vs_anchored else NA)
}
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- line(...)
for (g in c("all time points", "2 time points")) {
  arm <- paste("KARST beta auto,", g)
  for (src in c("ABM", "Missegregation")) for (w in c(4, 8)) {
    add(arm, paste(src, "landscape r, window", w), paste0(src, ": landscape Pearson r"), paste(g, "window", w), 0.01)
  }
  for (m in c("euclidean", "wasserstein")) for (w in c(4, 8)) {
    add(arm, paste0("ABM forecast beats baseline (", m, "), window ", w),
        paste0("ABM: forecast beats baseline, passage 10 (", m, ")"), paste(g, "window", w), 0.02)
  }
  sal_x <- if (g == "all time points") "all time points (3-10 passages)" else "2 time points"
  for (m in c("euclidean", "wasserstein")) {
    panel <- paste0("Salehi: share of daughters where the forecast beats the parent (", m, ")")
    add(arm, paste0("Salehi forecast beats parent (", m, ")"), panel, sal_x, 0.02)
    add(arm, paste0("Salehi sweeps, forecast beats parent (", m, ")"), panel, paste(g, "(sweeps)"), 0)
  }
  add(arm, "Salehi novel karyotypes", "Salehi: novel karyotypes, fitness predicts emergence", sal_x, 0.02)
}
crit <- do.call(rbind, rows)
# reported only
info <- list()
addi <- function(...) info[[length(info) + 1L]] <<- line(..., criterion = FALSE)
addi("KARST beta auto, 2 time points", "ABM landscape r, window 2", "ABM: landscape Pearson r", "2 time points window 2",
     0.01)
addi("KARST beta auto, 2 time points", "Missegregation landscape r, window 2", "Missegregation: landscape Pearson r",
     "2 time points window 2", 0.01)
for (m in c("euclidean", "wasserstein")) {
  addi("KARST beta auto, 2 time points", paste0("ABM forecast beats baseline (", m, "), window 2"),
       paste0("ABM: forecast beats baseline, passage 10 (", m, ")"), "2 time points window 2", 0.02)
}
for (g in c("all time points", "2 time points")) for (w in if (g == "2 time points") c(2, 4, 8) else c(4, 8)) {
  for (mt in c("karst_beta_auto", "karst_beta_fixed")) {
    addi(paste(method_label[[mt]], "vs true rate,", g), paste("Missegregation landscape r, window", w),
         "Missegregation: landscape Pearson r", paste(g, "window", w), 0.01, method = mt, reference = "karst")
  }
}
out <- rbind(crit, do.call(rbind, info))
utils::write.table(out, file.path(dir, "beta_criterion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 220)
print(out[c("arm", "line", "pairs", "value", "reference", "reference_value", "difference", "margin", "pass",
            "criterion", "better", "worse", "p")], row.names = FALSE)
for (a in unique(crit$arm)) {
  x <- crit[crit$arm == a, ]
  cat(sprintf("\n%s: %d of %d lines pass, %d not evaluated -> %s", a, sum(x$pass %in% TRUE), nrow(x), sum(is.na(x$pass)),
              if (all(x$pass %in% TRUE)) "PASS" else "FAIL"))
}
cat("\n")
