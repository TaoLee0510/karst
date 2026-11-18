#!/usr/bin/env Rscript
# Analysis of scripts/fig2d_check.R: why karst with two time points ranks the
# size of change badly in Fig. 2d (Salehi, round 05a8855).
#
#   BENCH_RESULTS=results_05a8855 SALEHI_META=.../metadata.csv Rscript scripts/fig2d_check_analyse.R
#
# Reads tables/salehi_scores.Rds, tables/salehi_inputs.Rds and the chunks in
# analysis/fig2d_check/{salehi,abm}/, and writes the tables (tsv) and
# figures of analysis/fig2d_check/.
#
# Hypotheses (plan of 2026-11-16):
#   H1  the paper's selection by CV score picks different trajectories;
#   H2  the karyotypes the axes are computed on (cells outside the landscape or
#       outside karst's model domain);
#   H3  karst's predicted change is set by the amount of data or by
#       missegregation, not by the strength of selection;
#   H4  a property of the model, seen on the ABM simulations with the truth.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
res <- Sys.getenv("BENCH_RESULTS")
tab <- file.path(res, "tables")
out <- file.path(res, "analysis", "fig2d_check")
fig <- file.path(out, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
theme_set(theme_classic(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top"))
write_tab <- function(d, name) {
  utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(d)
}
sp <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman", use = "complete.obs"))
arms <- c("KARST 2 time points", "alfak dev2@f93b2d2 2 time points")
arm_of <- c(karst = arms[[1]], alfak = arms[[2]])
cols <- c(karst = "#D55E00", alfak = "#0072B2")

# Spearman with a percentile interval from resampling clusters (lineages).
boot_sp <- function(d, x, y, cluster, B = 4000) {
  idx <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  bs <- replicate(B, {
    s <- unlist(idx[sample(length(idx), replace = TRUE)], use.names = FALSE)
    sp(d[[x]][s], d[[y]][s])
  })
  data.frame(n = nrow(d), clusters = length(idx), spearman = round(sp(d[[x]], d[[y]]), 3),
             low = round(stats::quantile(bs, 0.025, na.rm = TRUE), 3),
             high = round(stats::quantile(bs, 0.975, na.rm = TRUE), 3), row.names = NULL)
}
# Paired difference of two Spearman correlations with the same x, by cluster.
boot_diff <- function(d, x, y1, y2, cluster, B = 4000) {
  idx <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  bs <- replicate(B, {
    s <- unlist(idx[sample(length(idx), replace = TRUE)], use.names = FALSE)
    sp(d[[x]][s], d[[y1]][s]) - sp(d[[x]][s], d[[y2]][s])
  })
  data.frame(n = nrow(d), difference = round(sp(d[[x]], d[[y1]]) - sp(d[[x]], d[[y2]]), 3),
             low = round(stats::quantile(bs, 0.025, na.rm = TRUE), 3),
             high = round(stats::quantile(bs, 0.975, na.rm = TRUE), 3),
             share_above_0 = round(mean(bs > 0, na.rm = TRUE), 3), row.names = NULL)
}
# Partial Spearman of x and y given z (ranks, linear residuals).
partial_sp <- function(d, x, y, z) {
  r <- function(v) rank(v)
  ok <- stats::complete.cases(d[c(x, y, z)])
  d <- d[ok, ]
  zz <- as.data.frame(lapply(d[z], r))
  rx <- stats::resid(stats::lm(r(d[[x]]) ~ ., data = zz))
  ry <- stats::resid(stats::lm(r(d[[y]]) ~ ., data = zz))
  round(stats::cor(rx, ry), 3)
}

## Step 1: reproduce the paper's selection (analyse_salehi.R, verbatim) -------------
meta <- utils::read.csv(Sys.getenv("SALEHI_META"))
scores <- label_salehi_fits(readRDS(file.path(tab, "salehi_scores.Rds")))
scores <- scores[scores$fit %in% arms, ]
scores$type <- salehi_lineage_fig2(scores$fi, meta)
inputs_all <- readRDS(file.path(tab, "salehi_inputs.Rds"))
if (is.null(inputs_all$design)) inputs_all$design <- "ends"
excluded_keys <- not_analysable_keys(inputs_all, c("fi", "design"))
scores <- scores[!input_key(scores$fi, scores$design) %in% excluded_keys, ]
select_paper <- function(x) {
  x <- x[!is.na(x$xv) & x$xv > 0 & !is.na(x$test_treat) & x$test_treat != "NaN", ]
  if (!nrow(x)) return(NULL)
  x <- select_by_cv(x, "fi", "min_obs", "xv")
  x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
  x <- select_by_cv(x, "clean_id", "fi", "xv")
  rownames(x) <- NULL
  x
}
eval_rows <- function(x) {
  x <- split(x, x$fi)
  x <- x[vapply(x, function(d) 1 %in% d$tt, logical(1))]
  x <- do.call(rbind, lapply(x, function(d) {
    maxp <- if (isTRUE(d$type[1] %in% c("p53 k.o", "p53 w.t"))) 25 else 15
    d <- d[d$tt %in% 1:maxp, ]
    d$eval_tp <- d$tt == maxp
    d
  }))
  rownames(x) <- NULL
  x
}
sel <- do.call(rbind, lapply(arms, function(a) eval_rows(select_paper(scores[scores$fit == a, ]))))
sel <- sel[sel$eval_tp & sel$metric == "wasserstein", ]
repro <- do.call(rbind, lapply(arms, function(a) {
  x <- sel[sel$fit == a, ]
  data.frame(fit = a, n = nrow(x), spearman = round(sp(x$base, x$pred0), 3),
             pearson = round(stats::cor(x$base, x$pred0), 3))
}))
write_tab(repro, "s1_reproduction")

## The recomputed rows ----------------------------------------------------------------
read_chunks <- function(stage) {
  files <- list.files(file.path(out, stage), pattern = "^chunk_.*\\.Rds$", full.names = TRUE)
  do.call(rbind, lapply(files, readRDS))
}
rows <- read_chunks("salehi")
rows <- rows[!input_key(rows$fi, "ends") %in% excluded_keys, ]
rows$coverage <- inputs_all$coverage[match(input_key(rows$fi, "ends"),
                                           input_key(inputs_all$fi, inputs_all$design))]
rows$eval_ok <- !is.na(rows$test_treat) & rows$test_treat != "NaN"
# Each selected row of the paper's selection, found among the recomputed rows.
key3 <- function(d) paste(d$fi, d$method, d$min_obs)
match_sel <- vapply(seq_len(nrow(sel)), function(i) {
  cand <- which(key3(rows) == key3(sel[i, ]) & abs(rows$base - sel$base[i]) < 1e-6 &
                  abs(rows$pred0 - sel$pred0[i]) < 1e-6)
  if (length(cand)) cand[[1]] else NA_integer_
}, integer(1))
check <- data.frame(selected_rows = nrow(sel), found_with_same_axes = sum(!is.na(match_sel)))
write_tab(check, "s1_recompute_check")
sel_rows <- rows[match_sel[!is.na(match_sel)], ]
sel_rows$fit <- arm_of[sel_rows$method]

# The rows both methods have (same trajectory, minobs, daughter passage), with
# a daughter passage the paper evaluates.
both <- merge(rows[rows$method == "karst" & rows$eval_ok, ], rows[rows$method == "alfak" & rows$eval_ok, ],
              by = c("fi", "min_obs", "dj", "type", "tt"), suffixes = c("_a2", "_aR"))
same_base <- data.frame(common_rows = nrow(both),
                        same_observed_axis = sum(abs(both$base_a2 - both$base_aR) < 1e-9, na.rm = TRUE),
                        same_landscape_size = sum(both$landscape_size_a2 == both$landscape_size_aR))
write_tab(same_base, "s1_common_rows")
both$base <- both$base_a2
both$coverage <- both$coverage_a2

## Step 1: selection (H1) ----------------------------------------------------------
h1 <- rbind(
  cbind(set = "paper's selection, each its own", method = "karst",
        boot_sp(sel_rows[sel_rows$method == "karst", ], "base", "pred0", "type")),
  cbind(set = "paper's selection, each its own", method = "alfak",
        boot_sp(sel_rows[sel_rows$method == "alfak", ], "base", "pred0", "type")))
for (mo in c("5", "10", "20", "all")) {
  d <- if (mo == "all") both else both[both$min_obs == mo, ]
  h1 <- rbind(h1, cbind(set = paste("common rows, minobs", mo), method = "karst", boot_sp(d, "base", "pred0_a2", "type")),
              cbind(set = paste("common rows, minobs", mo), method = "alfak", boot_sp(d, "base", "pred0_aR", "type")))
}
# each method on the other's selection (same trajectory, minobs, daughter)
sel_key <- function(m) {
  x <- sel_rows[sel_rows$method == m, ]
  paste(x$fi, x$min_obs, x$dj)
}
bk <- paste(both$fi, both$min_obs, both$dj)
for (m in c("karst", "alfak")) {
  d <- both[bk %in% sel_key(m), ]
  h1 <- rbind(h1, cbind(set = paste0(m, "'s selection"), method = "karst", boot_sp(d, "base", "pred0_a2", "type")),
              cbind(set = paste0(m, "'s selection"), method = "alfak", boot_sp(d, "base", "pred0_aR", "type")))
}
write_tab(h1, "s1_selection")
both$sample <- sub("_X[0-9]+_l_.*$", "", both$fi)
h1_diff <- do.call(rbind, lapply(c("5", "10", "20", "all"), function(mo) {
  d <- if (mo == "all") both else both[both$min_obs == mo, ]
  cbind(set = paste("common rows, minobs", mo), cluster = "lineage (6)",
        boot_diff(d, "base", "pred0_a2", "pred0_aR", "type"))
}))
h1_diff <- rbind(h1_diff,
  cbind(set = "common rows, minobs all", cluster = "sample (12)", boot_diff(both, "base", "pred0_a2", "pred0_aR", "sample")),
  do.call(rbind, lapply(c("karst", "alfak"), function(m) {
    d <- both[bk %in% sel_key(m), ]
    rbind(cbind(set = paste0(m, "'s selection"), cluster = "lineage (6)",
                boot_diff(d, "base", "pred0_a2", "pred0_aR", "type")),
          cbind(set = paste0(m, "'s selection"), cluster = "sample (12)",
                boot_diff(d, "base", "pred0_a2", "pred0_aR", "sample")))
  })))
write_tab(h1_diff, "s1_paired_difference")
# by trajectory length (passages fitted, first and last only): 2 = one interval
both$len <- as.integer(sub(".*_l_([0-9]+)_d1_.*", "\\1", both$fi))
both$len_group <- ifelse(both$len >= 5, "5+", as.character(both$len))
h1_len <- do.call(rbind, lapply(split(both, both$len_group), function(d) {
  cbind(length = d$len_group[1], karst = round(sp(d$base, d$pred0_a2), 3), alfak = round(sp(d$base, d$pred0_aR), 3),
        boot_diff(d, "base", "pred0_a2", "pred0_aR", "sample"))
}))
write_tab(h1_len, "s1_by_length")

## Step 2: karyotypes the axes see (H2) ------------------------------------------------
cut_cov <- function(v) cut(v, c(0, 0.1, 0.2, 0.5, 1), right = FALSE, include.lowest = TRUE,
                           labels = c("(0, 0.1)", "[0.1, 0.2)", "[0.2, 0.5)", "[0.5, 1]"))
both$cov_bin <- cut_cov(both$coverage)
h2 <- do.call(rbind, lapply(split(both, both$cov_bin, drop = TRUE), function(d) {
  rbind(cbind(cov_bin = as.character(d$cov_bin[1]), method = "karst", boot_sp(d, "base", "pred0_a2", "type", 1000)),
        cbind(cov_bin = as.character(d$cov_bin[1]), method = "alfak", boot_sp(d, "base", "pred0_aR", "type", 1000)))
}))
write_tab(h2, "s2_by_coverage")
h2_axes <- rbind(
  data.frame(axes = "paper (landscape karyotypes)", karst = round(sp(both$base, both$pred0_a2), 3),
             alfak = round(sp(both$base, both$pred0_aR), 3)),
  data.frame(axes = "every observed karyotype", karst = round(sp(both$base_full_a2, both$pred0_full_a2), 3),
             alfak = round(sp(both$base_full_aR, both$pred0_full_aR), 3)),
  data.frame(axes = "Euclidean (mean karyotype)", karst = round(sp(both$base_euc_a2, both$pred0_euc_a2), 3),
             alfak = round(sp(both$base_euc_aR, both$pred0_euc_aR), 3)))
cov_tab <- data.frame(
  quantity = c("parent cells in the landscape", "daughter cells in the landscape",
               "parent cells in karst's model domain", "input coverage (karst's domain)"),
  median = round(c(stats::median(both$cov_parent_a2), stats::median(both$cov_daughter_a2),
                   stats::median(both$domain_parent_a2, na.rm = TRUE), stats::median(both$coverage)), 3))
write_tab(h2_axes, "s2_axes")
write_tab(cov_tab, "s2_coverage_medians")

## Step 3: what sets the predicted change (H3) -------------------------------------------
feat <- c("base", "sel0", "var_f", "price", "entropy_parent", "kary_parent_in", "cells_parent", "cov_parent",
          "domain_parent", "coverage", "n_fq", "n_major", "history", "cells_t0", "fit_dt")
long <- rbind(
  data.frame(both[c("fi", "min_obs", "dj", "type", "base", "coverage")], method = "karst",
             pred0 = both$pred0_a2, sel0 = both$sel0_a2, var_f = both$var_f_a2, price = both$price_a2,
             entropy_parent = both$entropy_parent_a2, kary_parent_in = both$kary_parent_in_a2,
             cells_parent = both$cells_parent_a2, cov_parent = both$cov_parent_a2,
             domain_parent = both$domain_parent_a2, n_fq = both$n_fq_a2, n_major = both$n_major_a2,
             history = both$history_a2, cells_t0 = both$cells_t0_a2, fit_dt = both$fit_dt_a2),
  data.frame(both[c("fi", "min_obs", "dj", "type", "base", "coverage")], method = "alfak",
             pred0 = both$pred0_aR, sel0 = both$sel0_aR, var_f = both$var_f_aR, price = both$price_aR,
             entropy_parent = both$entropy_parent_aR, kary_parent_in = both$kary_parent_in_aR,
             cells_parent = both$cells_parent_aR, cov_parent = both$cov_parent_aR,
             domain_parent = both$domain_parent_a2, n_fq = both$n_fq_aR, n_major = both$n_major_a2,
             history = both$history_a2, cells_t0 = both$cells_t0_a2, fit_dt = both$fit_dt_a2))
h3 <- do.call(rbind, lapply(c("karst", "alfak"), function(m) {
  d <- long[long$method == m, ]
  data.frame(method = m, feature = feat,
             with_pred0 = vapply(feat, function(f) round(sp(d$pred0, d[[f]]), 3), numeric(1)),
             with_base = vapply(feat, function(f) round(sp(d$base, d[[f]]), 3), numeric(1)),
             with_var_f = vapply(feat, function(f) round(sp(d$var_f, d[[f]]), 3), numeric(1)), row.names = NULL)
}))
write_tab(h3, "s3_spearman_features")
h3_partial <- do.call(rbind, lapply(c("karst", "alfak"), function(m) {
  d <- long[long$method == m, ]
  data.frame(method = m,
             pred0_base = round(sp(d$pred0, d$base), 3),
             given_entropy = partial_sp(d, "pred0", "base", "entropy_parent"),
             given_entropy_karyotypes = partial_sp(d, "pred0", "base", c("entropy_parent", "kary_parent_in")),
             given_entropy_karyotypes_cells = partial_sp(d, "pred0", "base",
                                                         c("entropy_parent", "kary_parent_in", "cells_parent")),
             sel0_base = round(sp(d$sel0, d$base), 3),
             pred0_over_sel0 = round(stats::median(d$pred0 / d$sel0, na.rm = TRUE), 3),
             pred0_sel0 = round(sp(d$pred0, d$sel0), 3))
}))
write_tab(h3_partial, "s3_partial")
# Scale of the fitness differences among the parent's cells, paired
scale_tab <- data.frame(
  median_var_ratio_a2_over_aR = round(stats::median(both$var_f_a2 / both$var_f_aR, na.rm = TRUE), 3),
  spearman_var_a2_var_aR = round(sp(both$var_f_a2, both$var_f_aR), 3),
  spearman_pred0_a2_pred0_aR = round(sp(both$pred0_a2, both$pred0_aR), 3),
  median_pred0_ratio_a2_over_aR = round(stats::median(both$pred0_a2 / both$pred0_aR, na.rm = TRUE), 3),
  median_pred0_over_base_a2 = round(stats::median(both$pred0_a2 / both$base, na.rm = TRUE), 3),
  median_pred0_over_base_aR = round(stats::median(both$pred0_aR / both$base, na.rm = TRUE), 3),
  spearman_ratio_var_with_cells_t0 = round(sp(both$var_f_a2 / both$var_f_aR, both$cells_t0_a2), 3),
  spearman_ratio_var_with_n_major = round(sp(both$var_f_a2 / both$var_f_aR, both$n_major_a2), 3),
  spearman_ratio_var_with_kary = round(sp(both$var_f_a2 / both$var_f_aR, both$kary_parent_in_a2), 3))
write_tab(scale_tab, "s3_scale")

# Exploratory (not in the plan): karst's chosen history, relative to the
# interval it fitted, by trajectory length. With one interval it scales the
# predicted change without following the observed change.
both$hist_rel <- both$history_a2 / both$fit_dt_a2
h3_hist <- do.call(rbind, lapply(split(both, both$len_group), function(d) {
  m <- stats::lm(log(pred0_a2) ~ log(hist_rel), d)
  data.frame(length = d$len_group[1], n = nrow(d),
             history_vs_pred0_a2 = round(sp(d$hist_rel, d$pred0_a2), 3),
             history_vs_base = round(sp(d$hist_rel, d$base), 3),
             history_vs_pred0_aR = round(sp(d$hist_rel, d$pred0_aR), 3),
             slope_log_pred0_a2_on_log_history = round(unname(stats::coef(m)[2]), 3),
             karst = round(sp(d$base, d$pred0_a2), 3),
             karst_given_history = partial_sp(d, "pred0_a2", "base", "hist_rel"),
             alfak = round(sp(d$base, d$pred0_aR), 3),
             alfak_given_history = partial_sp(d, "pred0_aR", "base", "hist_rel"))
}))
write_tab(h3_hist, "s3_history_by_length")
d2 <- both[both$len_group == "2", ]
d2$history_group <- ifelse(d2$hist_rel <= 1.01, "1 interval", ifelse(d2$hist_rel <= 2.01, "2 intervals", "> 2 intervals"))
h3_hist_groups <- do.call(rbind, lapply(split(d2, d2$history_group), function(x) data.frame(
  history = x$history_group[1], n = nrow(x), karst = round(sp(x$base, x$pred0_a2), 3),
  alfak = round(sp(x$base, x$pred0_aR), 3), median_pred0_a2 = round(stats::median(x$pred0_a2), 3),
  median_pred0_aR = round(stats::median(x$pred0_aR), 3),
  median_var_ratio_a2_over_aR = round(stats::median(x$var_f_a2 / x$var_f_aR), 3))))
write_tab(h3_hist_groups, "s3_history_groups_length2")

## Step 4: ABM with the truth (H4) ------------------------------------------------------
abm <- read_chunks("abm")
if (!is.null(abm) && nrow(abm)) {
  h4 <- do.call(rbind, lapply(split(abm, abm[c("design", "w", "method", "min_obs")], drop = TRUE), function(d) {
    data.frame(design = d$design[1], w = d$w[1], method = d$method[1], min_obs = d$min_obs[1], n = nrow(d),
               fig2d = round(sp(d$base, d$pred0), 3),
               truth_ceiling = round(sp(d$base, d$true_sel0), 3),
               sel_vs_true_sel = round(sp(d$sel0, d$true_sel0), 3),
               var_vs_true_var = round(sp(d$var_f, d$true_var_f), 3),
               price_vs_true_price = round(sp(d$price, d$true_price), 3),
               median_var_ratio = round(stats::median(d$var_f / d$true_var_f, na.rm = TRUE), 3),
               median_price_angle = round(stats::median(d$angle_price_truth, na.rm = TRUE), 1),
               median_spearman_f_parent = round(stats::median(d$spearman_f_parent, na.rm = TRUE), 3),
               median_pred0_over_sel0 = round(stats::median(d$pred0 / d$sel0, na.rm = TRUE), 3))
  }))
  h4 <- h4[order(h4$design, h4$w, h4$min_obs, h4$method), ]
  write_tab(h4, "s4_abm")
  # pooled over wavelengths, minobs 5
  h4p <- do.call(rbind, lapply(split(abm[abm$min_obs == "5", ], abm[abm$min_obs == "5", c("design", "method")],
                                     drop = TRUE), function(d) {
    data.frame(design = d$design[1], method = d$method[1], n = nrow(d), fig2d = round(sp(d$base, d$pred0), 3),
               truth_ceiling = round(sp(d$base, d$true_sel0), 3), var_vs_true_var = round(sp(d$var_f, d$true_var_f), 3),
               sel_vs_true_sel = round(sp(d$sel0, d$true_sel0), 3),
               median_var_ratio = round(stats::median(d$var_f / d$true_var_f, na.rm = TRUE), 3))
  }))
  write_tab(h4p, "s4_abm_pooled_minobs5")
  a <- abm[abm$method == "karst" & abm$min_obs == "5", ]
  a$hist_rel <- a$history / a$fit_dt
  h4h <- do.call(rbind, lapply(split(a, a$design), function(d) data.frame(
    design = d$design[1], n = nrow(d), history_equals_interval = round(mean(abs(d$hist_rel - 1) < 1e-6), 3),
    history_vs_pred0 = round(sp(d$hist_rel, d$pred0), 3), history_vs_base = round(sp(d$hist_rel, d$base), 3),
    history_vs_true_var = round(sp(d$hist_rel, d$true_var_f), 3))))
  write_tab(h4h, "s4_abm_history")
}

# The length-2 Salehi rows, coloured by karst's history (exploratory)
pl2 <- rbind(data.frame(d2[c("base", "history_group")], pred0 = d2$pred0_a2, method = "karst"),
             data.frame(d2[c("base", "history_group")], pred0 = d2$pred0_aR, method = "alfak"))
pl2$history_group <- factor(paste("karst history:", pl2$history_group),
                            paste("karst history:", c("1 interval", "2 intervals", "> 2 intervals")))
p2 <- ggplot(pl2, aes(base, pred0, colour = method)) +
  geom_point(size = 0.7, alpha = 0.7) +
  facet_wrap(~ history_group) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(NULL, values = cols) +
  labs(x = "observed change (Wasserstein)",
       y = "predicted change (Wasserstein)",
       subtitle = "trajectories of two passages (one interval), common rows")
ggsave(file.path(fig, "salehi_fig2d_length2_by_history.png"), p2, width = 180, height = 70, units = "mm", dpi = 200,
       bg = "white")
ggsave(file.path(fig, "salehi_fig2d_length2_by_history.pdf"), p2, width = 180, height = 70, units = "mm")

## Figures ----------------------------------------------------------------------------
pl <- rbind(data.frame(both[c("base", "min_obs")], pred0 = both$pred0_a2, method = "karst"),
            data.frame(both[c("base", "min_obs")], pred0 = both$pred0_aR, method = "alfak"))
p <- ggplot(pl, aes(base, pred0, colour = method)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ paste("minobs", min_obs)) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(NULL, values = cols) +
  labs(x = "observed change (Wasserstein)",
       y = "predicted change (Wasserstein)")
ggsave(file.path(fig, "salehi_fig2d_common_rows.png"), p, width = 180, height = 70, units = "mm", dpi = 200, bg = "white")
ggsave(file.path(fig, "salehi_fig2d_common_rows.pdf"), p, width = 180, height = 70, units = "mm")

cat("== reproduction\n"); print(repro, row.names = FALSE)
cat("== recompute check\n"); print(check, row.names = FALSE); print(same_base, row.names = FALSE)
cat("== H1\n"); print(h1, row.names = FALSE); print(h1_diff, row.names = FALSE)
cat("== H2\n"); print(h2, row.names = FALSE); print(h2_axes, row.names = FALSE); print(cov_tab, row.names = FALSE)
cat("== H3\n"); print(h3, row.names = FALSE); print(h3_partial, row.names = FALSE); print(scale_tab, row.names = FALSE)
if (exists("h4p")) { cat("== H4\n"); print(h4p, row.names = FALSE) }
