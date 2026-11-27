# Tests of the ALFA-K paper on the Salehi data, for every fit (04, 05).
#
#   Rscript scripts/analyse_salehi.R
#
# Reads BENCH_RESULTS/tables/salehi_{scores,reference,novel}.Rds
# (scripts/collect_salehi.R). Writes tables to BENCH_RESULTS/analysis/salehi/
# and figures to BENCH_RESULTS/figures/.
#
# Fits ("fit" column, R/labels.R), each trajectory being 2 to 10 consecutive
# passages:
#   karst 2 time points    karst, first and last passage of the trajectory
#   karst all time points  karst, every passage (multi-timepoint karst; for
#                       2-passage trajectories the same fit as above)
#   alfak 2 time points    alfak, the same two passages
#   alfak all time points  alfak, every passage (the paper's design; for
#                       2-passage trajectories the same fit as above)
#   paper alfak        the paper's fits (all passages), reference only
# The paper's rule selects, for each fit, the minobs and the trajectory length
# with the best CV score, so each fit is tested on different trajectories.
# The paired comparison takes every trajectory and minobs on which both
# fits have a forecast, without selection (pairs: `comparisons`). A fit
# failed when its run did (run.Rds, R/metrics_salehi.R salehi_fit_ok()), the
# definition of summary_comparison.R too; a fit that ran but has no CV score
# is left out only where the paper selects by CV score.
# Inputs that karst cannot fit, decided from the data before any fit for
# each design on the passages it fits (tables/salehi_inputs.Rds,
# R/coverage.R), are left out of every method's analysis on that design;
# excluded_inputs.tsv shows what the other fits give there. A sensitivity
# analysis splits the comparison by the share of the input inside karst's
# model domain (coverage_*.tsv, figures/coverage_salehi).

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
tab <- file.path(paths$results, "tables")
out <- file.path(paths$results, "analysis", "salehi")
fig <- file.path(paths$results, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
set.seed(1)

label_fits <- label_salehi_fits
# The paired comparisons: fit, reference.
comparisons <- list(c("KARST 2 time points", "alfak dev2@f93b2d2 2 time points"), c("KARST 2 time points", "alfak dev2@f93b2d2 all time points"),
                    c("KARST all time points", "alfak dev2@f93b2d2 all time points"), c("KARST all time points", "KARST 2 time points"))
# karst's other multi-time-point methods (c7f809c on; no rows in rounds without
# them): against alfak and against the default, "anchored".
for (f in c("KARST pairwise all time points", "KARST state space all time points")) {
  comparisons <- c(comparisons, list(c(f, "alfak dev2@f93b2d2 all time points"), c(f, "KARST all time points")))
}
# karst with its missegregation rate picked per input (ac36fd1 round; no rows
# in rounds without it), on each design: against alfak, against karst at
# 5e-5 (on the missegregation simulations "KARST beta 5e-5"; their "karst"
# takes the true rate) and against karst as fitted in the round.
for (g in c("2 time points", "all time points")) {
  f <- paste("KARST beta auto", g)
  comparisons <- c(comparisons, list(c(f, paste(method_label[["alfak"]], g)), c(f, paste("KARST", g)),
                                     c(f, paste("KARST beta 5e-5", g)), c(paste("KARST beta 5e-5", g), paste("KARST", g))))
}
write_tab <- function(d, name) {
  utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(d)
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = width, height = height, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = width, height = height, units = "mm", dpi = 200, bg = "white")
}
theme_set(theme_classic(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top"))
by_fit <- function(d, fun) {
  do.call(rbind, lapply(levels(d$fit), function(fl) {
    x <- d[d$fit == fl, ]
    if (!nrow(x)) return(NULL)
    r <- fun(x)
    if (is.null(r)) NULL else cbind(fit = fl, r)
  }))
}

meta <- utils::read.csv(paths$salehi_meta)
runs <- readRDS(file.path(tab, "salehi_runs.Rds"))
scores <- readRDS(file.path(tab, "salehi_scores.Rds"))
scores$fit_ok <- salehi_fit_ok(scores, runs)  # before label_fits(), which copies rows under other methods
# A fit of the run table without scores (a task killed before its fit ended,
# collected with ALLOW_MISSING=1, R/tasks.R) counts as failed.
run_key <- function(d) paste(d$fi, d$design, d$method, d$min_obs)
gone <- runs[!run_key(runs) %in% run_key(scores), c("fi", "design", "method", "min_obs")]
if (nrow(gone)) {
  scores <- bind_rows_fill(list(scores, data.frame(gone, xv = NaN, fit_ok = FALSE, metric = "NaN", test_treat = "NaN",
                                                   stringsAsFactors = FALSE)))
}
scores <- label_fits(scores)
ref <- readRDS(file.path(tab, "salehi_reference.Rds"))
novel <- label_fits(readRDS(file.path(tab, "salehi_novel.Rds")))
scores$type <- salehi_lineage_fig2(scores$fi, meta)
scores$gap <- scores$ntrain - 1L  # passages between the first and last sample
# Scores computed before the denoised baseline and entropies were added
for (n in c("denoised", "entropy_pred", "entropy_base", "entropy_true")) {
  if (is.null(scores[[n]])) scores[[n]] <- NA_real_
}

## Inputs that cannot be analysed -------------------------------------------------
# Decided from the data before any fit (R/coverage.R), for each design on the
# passages it fits: "ends" the first and last passage, "all" every passage
# (for 2-passage trajectories the same two). An input with a passage that has
# no cell inside karst's model domain is not linked at the karyotype level
# (no cell of that passage lies within one step of a major karyotype), and
# karst cannot fit it. A fit is analysed only when its own input is
# analysable, for every method: two-time-point fits by the "ends" rule (24
# trajectories with 8327282 and ff13cf4), all-time-point fits, the paper's line
# included, by the "all" rule. exclusion_check confirms, per design, that the
# set equals karst's failures for this reason; inputs the installed karst
# cannot read (more than two passages before the multi-timepoint version) are
# not judged and are kept.
inputs_all <- readRDS(file.path(tab, "salehi_inputs.Rds"))
if (is.null(inputs_all$design)) inputs_all$design <- "ends"  # tables of earlier rounds
inputs <- inputs_all[inputs_all$design == "ends", ]  # one row per trajectory
excluded_keys <- not_analysable_keys(inputs_all, c("fi", "design"))
excluded_row <- function(d) input_key(d$fi, d$design) %in% excluded_keys
excluded_ids <- inputs$fi[inputs$analysable %in% FALSE]  # by the "ends" rule
a2_runs <- runs[runs$method == "karst", ]
domain_class <- a2_runs$condition_class %in% "karst_domain_coverage_error"
domain_text <- grepl("no cell inside the model domain", a2_runs$message, fixed = TRUE)
domain_fail <- a2_runs$status != "ok" & (domain_class | domain_text)
exclusion_check <- do.call(rbind, lapply(c("ends", "all"), function(dg) {
  x <- inputs_all[inputs_all$design == dg, ]
  # 2-passage trajectories have no "all" fit of their own
  if (dg == "all") x <- x[x$gap > 1, ]
  if (!nrow(x)) return(NULL)
  bad <- x$fi[x$analysable %in% FALSE]
  failed <- unique(a2_runs$fi[domain_fail & a2_runs$design == dg])
  data.frame(design = dg, trajectories = nrow(x), judged = sum(!is.na(x$analysable)),
             not_analysable = length(bad), karst_domain_failures = length(failed), same_set = setequal(bad, failed),
             failed_fits_by_class = sum(a2_runs$status != "ok" & a2_runs$design == dg & domain_class),
             failed_fits_by_text_only = sum(a2_runs$status != "ok" & a2_runs$design == dg & domain_text & !domain_class),
             not_analysable_also_by_ends_rule = length(intersect(bad, excluded_ids)))
}))
write_tab(exclusion_check, "exclusion_check")
write_tab(inputs_all[inputs_all$analysable %in% FALSE, ], "excluded_trajectories")
# Candidate history chosen by karst (af362af on) and convergence of the
# chosen fit, by design and gap: one row per trajectory (minobs 5, the fit is
# shared). For all time points the chosen value (its first value where there is
# one per interval) is one of the fixed times or "other".
h <- a2_runs[a2_runs$status == "ok" & a2_runs$min_obs == "5" & !is.na(a2_runs$history), ]
hist_tab <- NULL
if (nrow(h)) {
  ends <- h$design == "ends"
  h$chosen <- karst_history_class(h$history)
  h$chosen[ends] <- karst_chosen_history(h$history[ends], h$candidates[ends])
  h$chosen <- factor(h$chosen, c("dt", "30", "100", "300", "1000", "3000", "other"))
  h$gap_group <- cut(inputs$gap[match(h$fi, inputs$fi)], c(0, 1, 2, 4, Inf), labels = c("1", "2", "3-4", "5+"))
  hist_tab <- bind_rows_fill(lapply(split(h, h[c("design", "gap_group")], drop = TRUE), function(x) data.frame(
    design = x$design[1], gap_group = x$gap_group[1], trajectories = nrow(x),
    t(as.matrix(table(x$chosen))),
    chosen_unconverged = round(mean(x$convergence != 0), 3), dt_unconverged = round(mean(x$dt_convergence != 0), 3),
    any_candidate_unconverged = round(mean(x$candidates_unconverged > 0), 3), check.names = FALSE)))
  counts <- intersect(names(hist_tab), levels(h$chosen))
  hist_tab[counts][is.na(hist_tab[counts])] <- 0L
  write_tab(hist_tab, "karst_history")
}
# For the coverage analysis: coverage of the input and its bin, gap groups
# and the sample (lineage branch) the trajectory comes from.
cut_coverage <- function(v) cut(v, c(0, 0.1, 0.2, 0.5, 1), right = FALSE, include.lowest = TRUE,
                                labels = c("(0, 0.1)", "[0.1, 0.2)", "[0.2, 0.5)", "[0.5, 1]"))
scores$coverage <- inputs$coverage[match(scores$fi, inputs$fi)]
scores$cov_bin <- cut_coverage(scores$coverage)
scores$gap_group <- cut(scores$gap, c(0, 1, 2, 4, Inf), labels = c("1", "2", "3-4", "5+"))
scores$sample <- sub("_X[0-9]+_l_.*$", "", scores$fi)

fail_table <- function(s) {
  f <- s[!duplicated(s[c("fit", "fi", "min_obs")]), ]
  by_fit(f, function(x) do.call(rbind, lapply(split(x, x$min_obs), function(y) data.frame(
    min_obs = y$min_obs[1], fits = nrow(y), failed = sum(!y$fit_ok), share_failed = round(mean(!y$fit_ok), 3),
    no_cv_score = sum(y$fit_ok & is.na(y$xv))))))
}
# Failures over all inputs, as before the exclusion, and by cause.
fail_all <- fail_table(scores)
write_tab(fail_all, "fit_failures_all_inputs")
fail_cause <- runs[runs$status != "ok", ]
# (karst's pairwise and state-space methods stop on an analysable input when
# a time point of a pair has no cell inside that pair's domain: "pair check".)
fail_cause$cause <- ifelse(excluded_row(fail_cause), "input not analysable",
                           ifelse(fail_cause$condition_class %in% "karst_domain_coverage_error", "pair check",
                                  ifelse(grepl("no frequent karyotypes", fail_cause$message, fixed = TRUE),
                                         "no frequent karyotype", "other")))
fail_cause <- as.data.frame(table(method = fail_cause$method, design = fail_cause$design,
                                  min_obs = fail_cause$min_obs, cause = fail_cause$cause))
write_tab(fail_cause[fail_cause$Freq > 0, ], "fit_failures_by_cause")

# What every fit gives on its excluded inputs (called after eval_rows is
# defined below).
excluded_summary <- function(s) {
  s <- s[excluded_row(s), ]
  f <- s[!duplicated(s[c("fit", "fi", "min_obs")]), ]
  head <- by_fit(f, function(x) {
    ok <- x[x$fit_ok, ]
    xv <- ok$xv[!is.na(ok$xv)]
    data.frame(trajectories = length(unique(x$fi)), fits = nrow(x), failed = sum(!x$fit_ok),
               fitted_trajectories = length(unique(ok$fi)),
               xv_median = if (length(xv)) round(stats::median(xv), 3) else NA_real_,
               xv_share_pos = if (length(xv)) round(mean(xv > 0), 3) else NA_real_)
  })
  ev <- s[s$fit_ok & !is.na(s$test_treat) & s$test_treat != "NaN", ]
  wins <- if (nrow(ev)) by_fit(ev, function(x) {
    e <- eval_rows(x)
    e <- e[e$eval_tp, ]
    if (!nrow(e)) return(NULL)
    w <- stats::aggregate(win ~ metric, e, mean)
    data.frame(fits_with_forecast = length(unique(paste(e$fi, e$min_obs))),
               t(stats::setNames(round(w$win, 3), paste0("win_", w$metric))))
  }) else NULL
  if (is.null(wins)) head else merge(head, wins, all.x = TRUE)
}

scores_all <- scores
scores <- scores[!excluded_row(scores), ]
novel <- novel[!excluded_row(novel), ]

## Fit failures (analysable inputs) ------------------------------------------------
fits <- scores[!duplicated(scores[c("fit", "fi", "min_obs")]), ]
fail <- fail_table(scores)
write_tab(fail, "fit_failures")

## Fig. 2b: CV score against training length ----------------------------------------
# 04: one row per fit and minobs that did not fail and has a CV score, best
# minobs per trajectory
x0 <- fits[fits$fit_ok & !is.na(fits$xv), ]
x0 <- select_by_cv(x0, c("fit", "fi"), "min_obs", "xv")
f2b <- by_fit(x0, function(x) do.call(rbind, lapply(split(x, x$ntrain), function(y) data.frame(
  ntrain = y$ntrain[1], gap = y$gap[1], trajectories = nrow(y), xv_median = round(stats::median(y$xv), 3),
  xv_share_pos = round(mean(y$xv > 0), 3)))))
write_tab(f2b, "fig2b_cv_by_length")
# 04: median CV score of trajectories longer than 3 passages, lineage bootstrap
f2b_ci <- by_fit(x0, function(x) {
  if (!any(x$ntrain > 3)) return(NULL)
  s <- lapply(split(x$xv[x$ntrain > 3], x$type[x$ntrain > 3]), identity)
  set.seed(1)  # per fit, so that a fit's interval does not depend on the other fits
  ci <- replicate(10000, stats::median(unlist(s[sample(length(s), length(s), replace = TRUE)])))
  data.frame(t(round(stats::quantile(ci, c(0.025, 0.5, 0.975)), 3)), check.names = FALSE)
})
write_tab(f2b_ci, "fig2b_cv_median_ci")
p <- ggplot(x0, aes(pmax(-1, xv), factor(ntrain), colour = fit)) +
  facet_grid(. ~ fit) +
  geom_violin(linewidth = 0.3) + geom_jitter(width = 0, height = 0.2, size = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey50") +
  scale_colour_manual(NULL, values = fit_colours, guide = "none") +
  labs(x = "CV score (best minobs; not comparable between methods)", y = "trajectory length (passages)")
save_plot(p, "fig2b_cv_by_length", 200, 70)

## Paper's selection and evaluation passage (04) ------------------------------------
# The best minobs of each trajectory, then the best trajectory of each start;
# ties resolved (R/select.R).
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
  x$win <- x$base > x$pred
  up <- x$metric %in% c("overlap", "cosine")
  x$win[up] <- (x$base < x$pred)[up]
  rownames(x) <- NULL
  x
}
sel <- by_fit(scores, function(x) {
  s <- select_paper(x)
  if (is.null(s)) NULL else eval_rows(s)
})
sel$fit <- factor(sel$fit, fit_levels)
sel_eval <- sel[sel$eval_tp, ]

# Fig. 2e-f: fraction beating the parent-passage baseline
f2e <- stats::aggregate(win ~ fit + type + metric + tt, data = sel, FUN = mean)
write_tab(f2e, "fig2e_beat_baseline_by_day")
f2f <- merge(stats::aggregate(win ~ fit + type + metric, data = sel_eval, FUN = mean),
             stats::aggregate(cbind(n = win) ~ fit + type + metric, data = sel_eval, FUN = length))
write_tab(f2f, "fig2f_beat_baseline_eval")
f2f_all <- merge(stats::aggregate(win ~ fit + metric, data = sel_eval, FUN = mean),
                 stats::aggregate(cbind(n = win) ~ fit + metric, data = sel_eval, FUN = length))
write_tab(f2f_all, "fig2f_beat_baseline_eval_pooled")
p <- ggplot(f2f, aes(type, win, fill = fit)) +
  facet_grid(. ~ metric) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(NULL, values = fit_colours) +
  coord_flip() + labs(x = NULL, y = "fraction beating the baseline at the evaluation passage",
                      subtitle = "paper's rule: minobs and trajectory length with the best CV score (> 0) per fit")
save_plot(p, "fig2f_beat_baseline", 200, 80)

# Fig. 2c: angle between predicted and observed change, cluster test
euc <- sel_eval[sel_eval$metric == "euclidean", ]
null_test <- function(angles, cluster) {
  set.seed(1)  # per test, independent of the other fits
  ok <- !is.na(cluster) & is.finite(angles)
  r <- sphere_null_test(angles[ok], cluster[ok])
  data.frame(n = sum(ok), clusters = length(unique(cluster[ok])), median_of_medians = round(r$stat, 1), p = r$p)
}
f2c <- rbind(
  by_fit(euc, function(x) cbind(group = "predicted", null_test(x$angle, x$type))),
  cbind(fit = "none", group = "forked sisters", null_test(ref$sisters$value[ref$sisters$metric == "angle"],
                                                        ref$sisters$lineage[ref$sisters$metric == "angle"])),
  cbind(fit = "none", group = "unrelated", null_test(ref$unrelated$angle,
                                                   paste(pmin(ref$unrelated$PDX1, ref$unrelated$PDX2),
                                                         pmax(ref$unrelated$PDX1, ref$unrelated$PDX2), sep = "--"))))
write_tab(f2c, "fig2c_angle_test")
ang <- rbind(data.frame(group = as.character(euc$fit), angle = euc$angle),
             data.frame(group = "forked sisters", angle = ref$sisters$value[ref$sisters$metric == "angle"]),
             data.frame(group = "unrelated", angle = ref$unrelated$angle))
p <- ggplot(ang, aes(angle, colour = group)) +
  stat_ecdf(linewidth = 0.4) +
  scale_colour_manual(NULL, values = c(fit_colours, "forked sisters" = "#009E73", "unrelated" = "#000000")) +
  labs(x = "angle (degrees)", y = "cumulative fraction")
save_plot(p, "fig2c_angle", 120, 80)

# Fig. 2d: predicted vs observed Wasserstein distance from the parent passage
f2d <- by_fit(sel_eval[sel_eval$metric == "wasserstein", ], function(x) data.frame(
  n = nrow(x), spearman = round(stats::cor(x$base, x$pred0, method = "spearman"), 3),
  pearson = round(stats::cor(x$base, x$pred0), 3)))
write_tab(f2d, "fig2d_wasserstein")

# Fig. 2g-h: forked sister passages as a reference
sis <- ref$sisters[ref$sisters$metric != "angle", ]
sis$dir <- ifelse(sis$metric %in% c("cosine", "overlap"), -1, 1)
sis_win <- c(sis$dir * (sis$value - sis$baseline1), sis$dir * (sis$value - sis$baseline2))
f2g <- stats::aggregate(list(win = sis_win < 0), by = list(lineage = c(sis$lineage, sis$lineage),
                                                           metric = c(sis$metric, sis$metric)), FUN = mean)
write_tab(f2g, "fig2g_sister_reference")

## Paired comparison: same trajectory and minobs, no selection ----------------------
# Per fit: the share of daughter passages beaten at the evaluation passage,
# averaged over the daughters. Every fit that ran (run status) and has a
# forecast, whether or not it has a CV score (no selection here).
all_eval <- by_fit(scores[scores$fit_ok & !is.na(scores$test_treat) & scores$test_treat != "NaN", ], eval_rows)
all_eval <- all_eval[all_eval$eval_tp, ]
# Trajectories of time points outside the first 12 lineages of Fig. 2 (04)
# have no lineage type; aggregate() would drop them from the pairs.
all_eval$type[is.na(all_eval$type)] <- "none (not in Fig. 2)"
# The parent restricted to its frequent karyotypes (R/metrics_salehi.R):
# does it beat the parent, and does the forecast beat it?
up <- all_eval$metric %in% c("overlap", "cosine")
all_eval$dwin <- ifelse(up, all_eval$denoised > all_eval$base, all_eval$denoised < all_eval$base)
all_eval$beats_dn <- ifelse(up, all_eval$pred > all_eval$denoised, all_eval$pred < all_eval$denoised)
fit_keys <- c("fi", "min_obs", "metric", "gap", "type", "coverage", "cov_bin", "gap_group", "sample")
per_fit <- stats::aggregate(cbind(win, angle, dwin, beats_dn) ~ fit + fi + min_obs + metric + gap + type + coverage +
                              cov_bin + gap_group + sample, data = all_eval, FUN = mean, na.action = stats::na.pass)
saveRDS(per_fit, file.path(out, "per_fit.Rds"))  # for scripts/summary_comparison.R
# The columns *.a2 / *_a2 hold the first fit of a comparison, *.ref / *_ref
# its reference.
paired_fits <- function(fit, reference) {
  a <- per_fit[per_fit$fit == fit, ]
  b <- per_fit[per_fit$fit == reference, ]
  merge(a, b, by = fit_keys, suffixes = c(".a2", ".ref"))
}
paired_fc <- function(fit, reference, by) {
  m <- paired_fits(fit, reference)
  do.call(rbind, lapply(split(m, m[by], drop = TRUE), function(x) cbind(
    fit = fit, reference = reference, x[1, by, drop = FALSE], fits = nrow(x),
    fwin_a2 = round(mean(x$win.a2), 3), fwin_ref = round(mean(x$win.ref), 3),
    a2_better = round(mean(x$win.a2 > x$win.ref), 3), ref_better = round(mean(x$win.a2 < x$win.ref), 3),
    angle_a2 = round(stats::median(x$angle.a2, na.rm = TRUE), 1),
    angle_ref = round(stats::median(x$angle.ref, na.rm = TRUE), 1),
    p = signif(tryCatch(stats::wilcox.test(x$win.a2, x$win.ref, paired = TRUE)$p.value,
                        error = function(e) NA_real_), 2))))
}
over <- function(cmps, by) do.call(rbind, lapply(cmps, function(cmp) paired_fc(cmp[1], cmp[2], by)))
pair_all <- over(comparisons, "metric")
write_tab(pair_all, "paired_forecast")
pair_gap <- over(comparisons, c("metric", "gap"))
write_tab(pair_gap, "paired_forecast_by_gap")
pg <- pair_gap[pair_gap$fit == "KARST 2 time points" & pair_gap$reference == "alfak dev2@f93b2d2 2 time points", ]
pg <- rbind(data.frame(pg[c("metric", "gap", "fits")], fit = "KARST 2 time points", fwin = pg$fwin_a2),
            data.frame(pg[c("metric", "gap", "fits")], fit = "alfak dev2@f93b2d2 2 time points", fwin = pg$fwin_ref))
p <- ggplot(pg, aes(gap, fwin, colour = fit)) +
  facet_grid(. ~ metric) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  geom_line() + geom_point(aes(size = fits)) +
  scale_size_area(max_size = 2.5) +
  scale_colour_manual(NULL, values = fit_colours) +
  scale_x_continuous("passages between the two time points", breaks = 1:9) +
  labs(y = "fraction beating baseline\n(paired fits, all minobs)")
save_plot(p, "stability_salehi_by_gap", 200, 70)

## Excluded trajectories: what the other fits give there ---------------------------
excluded_tab <- excluded_summary(scores_all)
write_tab(excluded_tab, "excluded_inputs")

## Sensitivity: how much of the input karst's model domain holds -------------------
# Coverage = the lower of the shares of the first and last passage inside
# karst's model domain (R/coverage.R); a description of the input, not an
# exclusion criterion. Coverage falls with the gap, and karst's advantage
# grows with the gap, so the paired comparison is also split by gap and a
# regression takes both, with a bootstrap over the samples the trajectories
# come from.
cov_inputs <- inputs[inputs$analysable %in% TRUE, ]
cov_inputs$cov_bin <- cut_coverage(cov_inputs$coverage)
cov_desc <- merge(stats::aggregate(cbind(trajectories = coverage) ~ cov_bin, cov_inputs, length),
                  stats::aggregate(cbind(gap_median = gap) ~ cov_bin, cov_inputs, stats::median))
write_tab(cov_desc, "coverage_bins")
cov_by_gap <- stats::aggregate(cbind(coverage_median = coverage) ~ gap, cov_inputs, function(v) round(stats::median(v), 3))
write_tab(cov_by_gap, "coverage_by_gap")
cov_rho <- stats::cor(cov_inputs$coverage, cov_inputs$gap, method = "spearman")

cov_unpaired <- merge(stats::aggregate(win ~ fit + cov_bin + metric, all_eval, function(v) round(mean(v), 3)),
                      stats::aggregate(cbind(rows = win) ~ fit + cov_bin + metric, all_eval, length))
write_tab(cov_unpaired, "coverage_unpaired")
# karst against alfak only
cov_cmps <- comparisons[vapply(comparisons, function(cmp) fit_method(cmp[2]) == "alfak", logical(1))]
pair_cov <- over(cov_cmps, c("metric", "cov_bin"))
write_tab(pair_cov, "paired_forecast_by_coverage")
pair_cov_gap <- over(cov_cmps, c("metric", "gap_group", "cov_bin"))
write_tab(pair_cov_gap, "paired_forecast_by_coverage_gap")

# Paired difference (karst minus reference, share of daughters beaten) ~
# log10(coverage) + gap. Predicted difference at coverage 0.05 and 1 and gap
# 1 and 3; 95% intervals from a bootstrap over samples (clusters).
cov_regression <- do.call(rbind, lapply(cov_cmps, function(cmp) {
  reference <- cmp[2]
  m <- paired_fits(cmp[1], reference)
  if (nrow(m) < 10) return(NULL)
  m$d <- m$win.a2 - m$win.ref
  do.call(rbind, lapply(split(m, m$metric), function(x) {
    grid <- data.frame(coverage = c(0.05, 1, 0.05, 1), gap = c(1, 1, 3, 3))
    estimate <- function(d) {
      f <- stats::lm(d ~ log10(coverage) + gap, data = d)
      c(slope_log10_coverage = unname(stats::coef(f)[2]), slope_gap = unname(stats::coef(f)[3]),
        stats::setNames(stats::predict(f, grid), sprintf("diff_cov%s_gap%d", grid$coverage, grid$gap)))
    }
    est <- estimate(x)
    clusters <- split(seq_len(nrow(x)), x$sample)
    set.seed(1)  # per comparison and metric
    boot <- replicate(2000, {
      take <- unlist(clusters[sample(length(clusters), replace = TRUE)], use.names = FALSE)
      estimate(x[take, ])
    })
    ci <- apply(boot, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE)
    data.frame(fit = cmp[1], reference = reference, metric = x$metric[1], fits = nrow(x), clusters = length(clusters),
               term = names(est), estimate = round(est, 3), lower = round(ci[1, ], 3), upper = round(ci[2, ], 3),
               row.names = NULL)
  }))
}))
write_tab(cov_regression, "coverage_regression")

# Denoising: does a forecast win because it concentrates on the parent's
# frequent karyotypes? Per fit and coverage bin: the share of daughters each
# forecast beats, the share the parent restricted to its frequent karyotypes
# beats (no model), the share where the forecast beats that restricted
# parent, and the median entropy change from parent to forecast and from
# parent to daughter.
denoise <- all_eval[!is.na(all_eval$denoised), ]
denoise_tab <- if (nrow(denoise)) {
  merge(stats::aggregate(cbind(win, dwin, beats_dn) ~ fit + cov_bin + metric, denoise,
                         function(v) round(mean(v), 3)),
        merge(stats::aggregate(cbind(rows = win) ~ fit + cov_bin + metric, denoise, length),
              stats::aggregate(cbind(entropy_change_forecast = entropy_pred - entropy_base,
                                     entropy_change_daughter = entropy_true - entropy_base) ~ fit + cov_bin + metric,
                               denoise, function(v) round(stats::median(v), 3))))
} else NULL
if (!is.null(denoise_tab)) write_tab(denoise_tab, "coverage_denoise")

pc <- pair_cov[pair_cov$fit == "KARST 2 time points" & pair_cov$reference == "alfak dev2@f93b2d2 2 time points", ]
pc <- rbind(data.frame(pc[c("metric", "cov_bin", "fits")], fit = "KARST 2 time points", fwin = pc$fwin_a2),
            data.frame(pc[c("metric", "cov_bin", "fits")], fit = "alfak dev2@f93b2d2 2 time points", fwin = pc$fwin_ref))
if (!is.null(denoise_tab)) {
  dn <- denoise_tab[denoise_tab$fit == "KARST 2 time points", ]
  pc <- rbind(pc, data.frame(metric = dn$metric, cov_bin = dn$cov_bin, fits = NA, fit = "parent, frequent only",
                             fwin = dn$dwin))
}
p <- ggplot(pc, aes(cov_bin, fwin, colour = fit, group = fit)) +
  facet_grid(. ~ metric) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  geom_line() + geom_point() +
  scale_colour_manual(NULL, values = c(fit_colours, "parent, frequent only" = "#009E73")) +
  labs(x = "share of the input inside karst's model domain (lower of the two time points)",
       y = "fraction beating baseline\n(paired fits, all minobs)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, "coverage_salehi", 200, 75)

## Fig. 3e-f: novel karyotype emergence (05) -------------------------------------------
novel_summary <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  z <- stats::aggregate(list(frac = df$Estimate > 0 & df$p < 0.05), by = list(var = df$ids), FUN = mean)
  # 05: per trajectory, the most significant positive variable; its first row
  # (the intercept of the first daughter's regression) is dropped first.
  w <- do.call(rbind, lapply(split(df, df$fi), function(wi) {
    wi <- wi[-1, ]
    wi <- wi[wi$Estimate > 0, ]
    wi[which.min(wi$p), ]
  }))
  data.frame(regressions = sum(df$ids == "(Intercept)"), trajectories = length(unique(df$fi)),
             f_significant = round(z$frac[z$var == "f"], 3), d1_significant = round(z$frac[z$var == "d1"], 3),
             f_most_significant = round(mean(w$ids == "f"), 3),
             d1_most_significant = round(mean(w$ids == "d1"), 3))
}
paper_rule <- function(x) {
  m <- scores[scores$fit == x$fit[1] & !is.na(scores$xv) & scores$tt %in% 1 & scores$metric %in% "cosine" &
                scores$xv > 0, ]
  if (!nrow(m)) return(NULL)
  m$tail_id <- vapply(strsplit(m$fi, "_l_", fixed = TRUE), `[`, character(1), 1)
  # the best fit and minobs of each start, ties resolved (R/select.R): the
  # lowest minobs, then the first trajectory
  m$choice <- paste(sprintf("%04d", as.integer(m$min_obs)), m$fi)
  m <- select_by_cv(m, "tail_id", "choice", "xv")
  # keep the row order: 05 drops the first row of each trajectory (an intercept)
  x[paste(x$fi, x$min_obs) %in% paste(m$fi, m$min_obs), ]
}
with_rule <- function(rule, s) if (is.null(s)) NULL else cbind(rule = rule, s)
f3 <- rbind(by_fit(novel, function(x) with_rule("paper", novel_summary(paper_rule(x)))),
            by_fit(novel, function(x) with_rule("all fits", novel_summary(x))))
write_tab(f3, "fig3ef_novel_karyotypes")
# paired: the same trajectory, minobs and daughter fitted by both
nf <- novel[novel$ids == "f", c("fit", "fi", "min_obs", "dec_id", "Estimate", "p")]
nf$sig <- nf$Estimate > 0 & nf$p < 0.05
pair_novel <- do.call(rbind, lapply(comparisons, function(cmp) {
  reference <- cmp[2]
  m <- merge(nf[nf$fit == cmp[1], ], nf[nf$fit == reference, ], by = c("fi", "min_obs", "dec_id"),
             suffixes = c(".a2", ".ref"))
  if (!nrow(m)) return(NULL)
  data.frame(fit = cmp[1], reference = reference, regressions = nrow(m), f_sig_a2 = round(mean(m$sig.a2), 3),
             f_sig_ref = round(mean(m$sig.ref), 3), only_a2 = sum(m$sig.a2 & !m$sig.ref),
             only_ref = sum(!m$sig.a2 & m$sig.ref),
             mcnemar_p = signif(stats::mcnemar.test(table(factor(m$sig.a2, c(FALSE, TRUE)),
                                                          factor(m$sig.ref, c(FALSE, TRUE))))$p.value, 2))
}))
write_tab(pair_novel, "fig3_paired")

## Console summary -----------------------------------------------------------------------
cat("\n== Inputs not analysable (decided before fitting) and karst's domain failures\n")
print(exclusion_check, row.names = FALSE)
cat("\n== Fits on the excluded trajectories\n"); print(excluded_tab, row.names = FALSE)
if (!is.null(hist_tab)) {
  cat("\n== karst: chosen history and convergence, by design and gap\n"); print(hist_tab, row.names = FALSE)
}
cat("\n== Fit failures, all inputs (before the exclusion)\n"); print(fail_all, row.names = FALSE)
cat("\n== Fit failures by cause (this benchmark's fits)\n")
print(fail_cause[fail_cause$Freq > 0, ], row.names = FALSE)
cat("\n== Fit failures (analysable inputs)\n"); print(fail, row.names = FALSE)
cat("\n== Fig. 2b CV score by trajectory length (best minobs)\n"); print(f2b, row.names = FALSE)
cat("\n== Fig. 2b median CV score (ntrain > 3), lineage bootstrap\n"); print(f2b_ci, row.names = FALSE)
cat("\n== Fig. 2f fraction beating baseline at the evaluation passage (paper's rule, pooled)\n")
print(f2f_all, row.names = FALSE)
cat("\n== Fig. 2c angle test\n"); print(f2c, row.names = FALSE)
cat("\n== Fig. 2d Wasserstein\n"); print(f2d, row.names = FALSE)
cat("\n== Fig. 2g sister reference (pooled)\n")
print(stats::aggregate(win ~ metric, f2g, mean), row.names = FALSE)
cat("\n== Paired forecasts (no selection)\n"); print(pair_all, row.names = FALSE)
cat("\n== Paired forecasts by gap, euclidean\n")
print(pair_gap[pair_gap$metric == "euclidean", ], row.names = FALSE)
cat("\n== Coverage: trajectories per bin; median by gap; Spearman with gap", round(cov_rho, 3), "\n")
print(cov_desc, row.names = FALSE); print(cov_by_gap, row.names = FALSE)
cat("\n== Unpaired, by coverage (wasserstein)\n")
print(cov_unpaired[cov_unpaired$metric == "wasserstein", ], row.names = FALSE)
cat("\n== Paired forecasts by coverage\n"); print(pair_cov, row.names = FALSE)
cat("\n== Paired forecasts by gap and coverage, euclidean and wasserstein, karst vs alfak 2 time points\n")
print(pair_cov_gap[pair_cov_gap$fit == "KARST 2 time points" & pair_cov_gap$reference == "alfak dev2@f93b2d2 2 time points" &
                     pair_cov_gap$metric %in% c("euclidean", "wasserstein"), ], row.names = FALSE)
cat("\n== Paired difference ~ log10(coverage) + gap, sample bootstrap\n"); print(cov_regression, row.names = FALSE)
if (!is.null(denoise_tab)) {
  cat("\n== Denoising: forecast win, parent-frequent-only win, forecast beats it, entropy changes\n")
  print(denoise_tab[denoise_tab$fit != "paper (alfak main@4ad82b0)", ], row.names = FALSE)
}
cat("\n== Fig. 3e-f\n"); print(f3, row.names = FALSE)
cat("\n== Fig. 3 paired\n"); print(pair_novel, row.names = FALSE)
