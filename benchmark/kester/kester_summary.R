# Kester summary (benchmark/rc_1.1.4/PLAN.md, part 3): tables from scores.tsv.
#   Rscript kester_summary.R SCORES_TSV MANIFEST_TSV OUT_DIR
args <- commandArgs(TRUE)
s <- utils::read.delim(args[[1]], stringsAsFactors = FALSE); man <- utils::read.delim(args[[2]], stringsAsFactors = FALSE)
out <- args[[3]]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 220)
s <- merge(s, man[, c("input", "kind", "culture", "design", "target_cells")], by = "input")
w <- function(x, d = 3) formatC(x, format = "f", digits = d)
md <- function(x) c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
                    apply(x, 1, function(r) paste("|", paste(trimws(r), collapse = " | "), "|")))
lines <- c("# Kester: KARST 1.1.4 rc against ALFA-K by missegregation multiple", "",
           "Fits per fold on training passages only; 1 passage = 7 days; multiples of 5e-5 for fitting and forecasting.", "")

## Salehi-like, engine on list E ----------------------------------------------------------------
e <- s[s$kind == "salehi" & s$type == "engine" & s$status == "ok", ]
arms <- c("karst", "alfak_m5", "alfak_m10", "alfak_m20", "neutral", "persistence")
agg <- function(x, v) stats::aggregate(x[[v]], list(design = x$design, multiple = x$multiple, arm = x$arm), mean)
for (v in c("log_score", "angle", "w1")) {
  a <- agg(e[e$arm %in% arms, ], v)
  t <- stats::reshape(a, idvar = c("design", "multiple"), timevar = "arm", direction = "wide")
  names(t) <- sub("^x\\.", "", names(t)); t <- t[order(t$design, t$multiple), ]
  utils::write.table(t, file.path(out, paste0("salehi_engine_", v, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  t[-(1:2)] <- lapply(t[-(1:2)], w, d = if (v == "log_score") 3 else if (v == "angle") 1 else 3)
  lines <- c(lines, sprintf("## Salehi-like folds, engine on list E: mean %s over folds (6 all / 3 ends)", v), "", md(t), "")
}
# paired KARST - ALFA-K (minobs 10) and KARST - neutral, per fold, all-design folds
pk <- function(a, b, v, design = "all") {
  x <- e[e$arm == a & e$design == design, c("input", "multiple", v)]; y <- e[e$arm == b & e$design == design, c("input", "multiple", v)]
  m <- merge(x, y, by = c("input", "multiple")); m$d <- m[[paste0(v, ".x")]] - m[[paste0(v, ".y")]]
  stats::aggregate(d ~ multiple, m, function(z) sprintf("%s (%d/%d)", w(mean(z)), sum(z > 0), length(z)))
}
pt <- Reduce(function(a, b) merge(a, b, by = "multiple"), list(
  stats::setNames(pk("karst", "alfak_m10", "log_score"), c("multiple", "KARST - ALFA-K m10")),
  stats::setNames(pk("karst", "alfak_m5", "log_score"), c("multiple", "KARST - ALFA-K m5")),
  stats::setNames(pk("karst", "neutral", "log_score"), c("multiple", "KARST - neutral")),
  stats::setNames(pk("karst", "persistence", "log_score"), c("multiple", "KARST - persistence")),
  stats::setNames(pk("alfak_m10", "neutral", "log_score"), c("multiple", "ALFA-K m10 - neutral"))))
lines <- c(lines, "## Salehi-like, all-passage folds: paired log-score differences per cell (mean; folds where the first is better / folds)", "", md(pt), "")

## Salehi-like, paper-style ABM ----------------------------------------------------------------
p <- s[s$kind == "salehi" & s$type == "paper_abm" & s$status == "ok", ]
if (nrow(p)) {
  a <- stats::aggregate(cbind(angle, wasserstein) ~ design + multiple + arm, p, mean)
  t <- stats::reshape(a, idvar = c("design", "multiple"), timevar = "arm", direction = "wide"); t <- t[order(t$design, t$multiple), ]
  utils::write.table(t, file.path(out, "salehi_paper_abm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  keep <- c("design", "multiple", grep("^angle", names(t), value = TRUE)); t2 <- t[, keep]; t2[-(1:2)] <- lapply(t2[-(1:2)], w, d = 1)
  lines <- c(lines, "## Salehi-like, paper-style ABM forecast (each method's landscape): mean Fig. 2c angle over folds", "", md(t2), "")
}

## all passages against the endpoints (three folds) ----------------------------------------------
ae <- e[e$arm %in% c("karst", "alfak_m10") & sub("_(all|ends)$", "", sub("_P[0-9-]+_to", "_to", e$input)) %in% c("APKS2_to_P18", "APKS3_to_P19", "APKS3_to_P23"), ]
if (nrow(ae)) {
  ae$fold <- sub("_(all|ends)$", "", sub("_P[0-9-]+_to", "_to", ae$input))
  ae <- ae[ae$input %in% c("APKS2_P5-8-14_to_P18_all", "APKS3_P3-7-15_to_P19_all", "APKS3_P3-7-15-19_to_P23_all",
                           "APKS2_P5-14_to_P18_ends", "APKS3_P3-15_to_P19_ends", "APKS3_P3-19_to_P23_ends"), ]
  a <- stats::aggregate(log_score ~ arm + design + multiple, ae, mean)
  t <- stats::reshape(a, idvar = c("arm", "multiple"), timevar = "design", direction = "wide"); t <- t[order(t$arm, t$multiple), ]
  t$all_minus_ends <- t$log_score.all - t$log_score.ends
  utils::write.table(t, file.path(out, "salehi_all_vs_ends.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  t[-(1:2)] <- lapply(t[-(1:2)], w)
  lines <- c(lines, "## All training passages against the endpoints (3 folds, same targets): mean log score per cell", "", md(t), "")
}

## GBM-like: within-harvest cell holdout ---------------------------------------------------------
g <- s[s$kind == "gbm" & s$status == "ok" & s$design != "full", ]
if (nrow(g)) {
  g$score <- g$nll_per_cell
  a <- stats::aggregate(cbind(nll_per_cell, baseline_nll_per_cell, scored_coverage) ~ culture + type + multiple + arm, g, mean)
  utils::write.table(a, file.path(out, "gbm_cell_holdout.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (cu in unique(a$culture)) {
    x <- a[a$culture == cu & a$type %in% c("engine", "passage_2"), c("multiple", "arm", "nll_per_cell")]
    t <- stats::reshape(x, idvar = "multiple", timevar = "arm", direction = "wide"); names(t) <- sub("^nll_per_cell\\.", "", names(t))
    t <- t[order(t$multiple), ]; t[-1] <- lapply(t[-1], w)
    lines <- c(lines, sprintf("## %s: within-harvest cell holdout, later passage, negative log likelihood per held-out cell (lower is better; mean of 3 folds)", cu), "", md(t), "")
  }
}
fails <- s[s$status != "ok", ]
lines <- c(lines, "## Failures", "", if (nrow(fails)) md(stats::aggregate(multiple ~ input + arm + status, fails, function(z) paste(sort(unique(z)), collapse = ","))) else "none", "")
writeLines(lines, file.path(out, "kester_summary.md"))
cat(lines, sep = "\n")
