# Part 1 of PLAN.md: the paper benchmark's Salehi Fig. 2c / 2d for KARST and
# ALFA-K at every missegregation multiple, from the beta-grid study's tables
# (benchmark/paper/results_6ae2a24/beta_grid/analysis; 1.1.3 = 1.1.4 by
# default). No new fits.
#   Rscript benchmark/rc_1.1.4/paper_grid_summary.R
a <- "/Users/4482173/Documents/GitHub/karst/benchmark/paper/results_6ae2a24/beta_grid/analysis"
out <- "/Users/4482173/Documents/GitHub/karst/benchmark/rc_1.1.4/results/paper_grid"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) utils::read.delim(file.path(a, f), stringsAsFactors = FALSE)
f2c <- rd("fig2c_paired.tsv"); f2r <- rd("fig2c_rule.tsv"); f2d <- rd("fig2d_rule.tsv"); fn <- rd("fig2d_noise_differences.tsv"); fl <- rd("failures.tsv")
fmt <- function(x, d = 1) formatC(x, format = "f", digits = d)
tab <- function(x, cols) { x <- x[order(x$design, x$multiple), cols]; x }
t1 <- tab(transform(f2c, KARST = fmt(angle_karst), ALFA_K = fmt(angle_alfak), diff_ALFA_K_minus_KARST = sprintf("%s (%s, %s)", fmt(diff), fmt(lo), fmt(hi))),
          c("design", "multiple", "pairs", "KARST", "ALFA_K", "diff_ALFA_K_minus_KARST", "significant"))
t2 <- tab(transform(f2r, KARST = fmt(angle_karst), ALFA_K = fmt(angle_alfak), diff_ALFA_K_minus_KARST = sprintf("%s (%s, %s)", fmt(diff), fmt(lo), fmt(hi))),
          c("design", "multiple", "n_karst", "n_alfak", "KARST", "ALFA_K", "diff_ALFA_K_minus_KARST", "significant"))
t3 <- tab(transform(f2d, KARST = fmt(rho_karst, 2), ALFA_K = fmt(rho_alfak, 2), diff_KARST_minus_ALFA_K = sprintf("%s (%s, %s)", fmt(diff, 2), fmt(lo, 2), fmt(hi, 2))),
          c("design", "multiple", "n_karst", "n_alfak", "KARST", "ALFA_K", "diff_KARST_minus_ALFA_K", "karst_ahead"))
fn <- fn[fn$subset == "every trajectory" & fn$arm == "karst", ]
t4 <- fn[order(fn$group, fn$measure, fn$multiple), c("group", "measure", "multiple", "difference_from_alfak", "lo", "hi")]
t5 <- stats::reshape(fl, idvar = c("multiple", "design"), timevar = "method", direction = "wide")
t5 <- t5[order(t5$design, t5$multiple), ]
for (nm in c("fig2c_no_selection", "fig2c_paper_rule", "fig2d_paper_rule", "fig2d_no_selection", "failures"))
  utils::write.table(get(c(fig2c_no_selection = "t1", fig2c_paper_rule = "t2", fig2d_paper_rule = "t3", fig2d_no_selection = "t4", failures = "t5")[[nm]]),
                     file.path(out, paste0(nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
md <- function(x) c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
                    apply(x, 1, function(r) paste("|", paste(trimws(r), collapse = " | "), "|")))
writeLines(c("# Paper benchmark: Salehi Fig. 2c / 2d by missegregation multiple (KARST 1.1.4 rc = 1.1.3 fits)", "",
             "Multiples of 5e-5, the same rate for fitting and forecasting, for both methods. Source: `benchmark/paper/results_6ae2a24/beta_grid/analysis` (beta-grid study of 2026-11-23/24).", "",
             "## Fig. 2c without selection (paired median angle; ALFA-K minus KARST > 0: KARST closer)", "", md(t1), "",
             "## Fig. 2c by the paper rule (median of lineage-type medians)", "", md(t2), "",
             "## Fig. 2d by the paper rule (Spearman)", "", md(t3), "",
             "## Fig. 2d without selection, every trajectory (KARST minus ALFA-K)", "", md(t4), "",
             "## Failure rates", "", md(t5)), file.path(out, "paper_grid.md"))
cat(readLines(file.path(out, "paper_grid.md")), sep = "\n")
