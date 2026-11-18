# Tests of the ALFA-K paper on the ABM data, for every fit (S05, X01).
#
#   Rscript scripts/analyse_abm.R
#
# Reads BENCH_RESULTS/tables/abm_{dfl,dfp,initial,runs}.Rds and
# misseg_dfl.Rds (scripts/collect_abm.R) and the paper's ABM_summaries.Rds.
# Writes tables to BENCH_RESULTS/analysis/abm/ and figures to
# BENCH_RESULTS/figures/.
#
# Fits ("fit" column, R/labels.R):
#   karst 2 time points    karst, first and last time point of the window
#   karst all time points  karst, every time point of the window (multi-timepoint
#                       karst; for the 2-passage window the same fit as
#                       "KARST 2 time points")
#   alfak 2 time points    alfak, the same two time points
#   alfak all time points  alfak, every time point of the window (for the 2-passage
#                       window the same fit as "alfak dev2@f93b2d2 2 time points")
#   paper alfak        the paper's fits (all time points), reference only
# Two ways of reporting, as in the README: the paper's rule (failed fits
# dropped; where the paper selects, the minobs with the best CV score of each
# fit) and a paired comparison of two fits on the same input: karst against
# alfak on the same time points (2 time points; all time points), karst on two time points
# against alfak on all, and karst on all time points against karst on two.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
tab <- file.path(paths$results, "tables")
out <- file.path(paths$results, "analysis", "abm")
fig <- file.path(paths$results, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
options(width = 200)

label_fits <- label_abm_fits
# The paired comparisons: fit, reference.
comparisons <- list(c("KARST 2 time points", "alfak dev2@f93b2d2 2 time points"), c("KARST 2 time points", "alfak dev2@f93b2d2 all time points"),
                    c("KARST all time points", "alfak dev2@f93b2d2 all time points"), c("KARST all time points", "KARST 2 time points"))
# karst's other multi-time-point methods (c7f809c on; no rows in rounds without
# them): against alfak and against the default, "anchored".
for (f in c("KARST pairwise all time points", "KARST state space all time points", "KARST auto all time points")) {
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
q3 <- function(x) stats::setNames(stats::quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE),
                                  c("q10", "median", "q90"))
summarise_by <- function(d, by, vars, fun = q3) {
  g <- split(d, d[by], drop = TRUE)
  do.call(rbind, lapply(g, function(x) {
    row <- x[1, by, drop = FALSE]
    for (v in vars) {
      s <- fun(x[[v]])
      for (n in names(s)) row[[paste(v, n, sep = "_")]] <- round(s[[n]], 4)
    }
    row$n <- nrow(x)
    row
  }))
}
write_tab <- function(d, name) {
  by <- intersect(c("fit", "window", "w", "minobs", "passage", "metric"), names(d))
  # (order() of no column gives no rows)
  if (length(by)) d <- d[do.call(order, unname(as.list(d[by]))), ]
  utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(d)
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = width, height = height, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = width, height = height, units = "mm", dpi = 200, bg = "white")
}
theme_set(theme_classic(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top"))
# Best CV score among the minobs of one fit of one input (S05: ties kept).
select_best_xv <- function(d, by = c("fit", "sim_id", "window")) {
  d <- d[is.finite(d$Rxv), ]
  do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(x) x[x$Rxv == max(x$Rxv), , drop = FALSE]))
}

## Data -----------------------------------------------------------------------
dfl <- label_fits(readRDS(file.path(tab, "abm_dfl.Rds")))
lscape_cols <- c("r", "rfq", "rnn", "rd2", "rho", "rhofq", "rhonn", "rhod2", "R", "Rfq", "Rnn", "Rd2")
for (v in c(lscape_cols, "Rxv")) dfl[[v]][!is.finite(dfl[[v]])] <- -Inf
clip <- function(x) pmax(-1, x)

dfp <- readRDS(file.path(tab, "abm_dfp.Rds"))
paper_dfp <- readRDS(file.path(paths$data, "processed", "ABM_summaries.Rds"))$dfp
paper_dfp <- data.frame(prediction = paper_dfp$prediction, baseline = paper_dfp$baseline,
                        passage = paper_dfp$passage - 1L, metric = paper_dfp$metric, angle = paper_dfp$angle,
                        source = "abm",
                        sim_id = sprintf("w_%s_m_5e-05_rep_%s", sub(".", "p", paper_dfp$w, fixed = TRUE),
                                         paper_dfp$abmrep),
                        design = c("2" = "w2_ends", "4" = "w4_all", "8" = "w8_all")[as.character(paper_dfp$ntp)],
                        method = "alfak_paper", minobs = paper_dfp$minobs, w = paper_dfp$w, Rxv = paper_dfp$Rxv,
                        win = paper_dfp$win, pos_xv = paper_dfp$pos_xv, stringsAsFactors = FALSE)
dfp <- label_fits(rbind(dfp[names(paper_dfp)], paper_dfp))

runs <- readRDS(file.path(tab, "abm_runs.Rds"))

## Fit failures ----------------------------------------------------------------
fail <- label_fits(runs[runs$source == "abm", ])
fail <- do.call(rbind, lapply(split(fail, fail[c("fit", "window", "minobs")], drop = TRUE), function(x)
  data.frame(fit = x$fit[1], window = x$window[1], minobs = x$minobs[1], fits = nrow(x),
             failed = sum(x$status != "ok"), share_failed = round(mean(x$status != "ok"), 3))))
write_tab(fail, "fit_failures")

# Inputs karst cannot fit (a time point with no cell inside its model domain),
# decided before any fit for each design on the time points it fits
# (R/coverage.R): first and last time point, or every time point. A fit is analysed
# only when its own input is analysable, for every method. Per data set and
# kind of design, the set must equal the karst fits that failed for this
# reason; inputs the installed karst cannot read (more than two time points
# before the multi-timepoint version) are not judged and are kept.
inputs <- readRDS(file.path(tab, "abm_inputs.Rds"))
key_cols <- c("source", "sim_id", "design")
excluded_keys <- not_analysable_keys(inputs, key_cols)
domain_fail <- runs$method == "karst" & runs$status != "ok" &
  (runs$condition_class %in% "karst_domain_coverage_error" |
     grepl("no cell inside the model domain", runs$message, fixed = TRUE))
failed_keys <- unique(input_key(runs$source, runs$sim_id, runs$design)[domain_fail])
kind <- function(design) ifelse(grepl("_all$", design), "all time points", "2 time points")
exclusion_check <- do.call(rbind, lapply(c("abm", "misseg"), function(s) do.call(rbind, lapply(
  c("2 time points", "all time points"), function(g) {
    x <- inputs[inputs$source == s & kind(inputs$design) == g, ]
    if (!nrow(x)) return(NULL)
    keys <- input_key(x$source, x$sim_id, x$design)
    bad <- keys[x$analysable %in% FALSE]
    failed <- intersect(failed_keys, keys)
    data.frame(source = s, inputs = g, judged = sum(!is.na(x$analysable)), not_judged = sum(is.na(x$analysable)),
               not_analysable = length(bad), karst_domain_failures = length(failed), same_set = setequal(bad, failed),
               coverage_median = round(stats::median(x$coverage, na.rm = TRUE), 3),
               coverage_min = round(suppressWarnings(min(x$coverage, na.rm = TRUE)), 3))
  }))))
write_tab(exclusion_check, "exclusion_check")
drop_excluded <- function(d) d[!input_key(d$source, d$sim_id, d$design) %in% excluded_keys, ]
dfl <- drop_excluded(dfl)
dfp <- drop_excluded(dfp)
# Candidate history chosen by karst (af362af on) and whether the chosen fit's
# optimisation converged, by design: one row per input (minobs 10). For two
# time points the candidates are dt (the interval, the fit of earlier commits)
# and 30, 100, 300, 1000, 3000 above it; for all time points the chosen value
# (its first value where there is one per interval) is one of the fixed
# times or "other".
h <- runs[runs$method == "karst" & runs$status == "ok" & runs$minobs == "10" & !is.na(runs$history), ]
if (nrow(h)) {
  ends <- grepl("_ends$", h$design)
  h$chosen <- karst_history_class(h$history)
  h$chosen[ends] <- karst_chosen_history(h$history[ends], h$candidates[ends])
  h$chosen <- factor(h$chosen, c("dt", "30", "100", "300", "1000", "3000", "other"))
  hist_tab <- bind_rows_fill(lapply(split(h, h[c("source", "design")], drop = TRUE), function(x) data.frame(
    source = x$source[1], design = x$design[1], inputs = nrow(x),
    t(as.matrix(table(x$chosen))),
    chosen_unconverged = round(mean(x$convergence != 0), 3),
    dt_unconverged = round(mean(x$dt_convergence != 0), 3),
    any_candidate_unconverged = round(mean(x$candidates_unconverged > 0), 3), check.names = FALSE)))
  counts <- intersect(names(hist_tab), levels(h$chosen))
  hist_tab[counts][is.na(hist_tab[counts])] <- 0L
  write_tab(hist_tab, "karst_history")
}

## Supp. Fig. 2a: landscape accuracy by wavelength, minobs and window ----------
s2a <- dfl
for (v in lscape_cols) s2a[[v]] <- clip(s2a[[v]])
s2a <- summarise_by(s2a, c("fit", "window", "w", "minobs"), lscape_cols)
write_tab(s2a, "supp2a_landscape_accuracy")
long <- do.call(rbind, lapply(c("r", "R", "rho"), function(v) data.frame(
  s2a[c("fit", "window", "w", "minobs")], metric = v, lo = s2a[[paste0(v, "_q10")]],
  med = s2a[[paste0(v, "_median")]], hi = s2a[[paste0(v, "_q90")]])))
long$metric <- factor(long$metric, c("r", "R", "rho"), c("Pearson r", "rescaled R^2", "Spearman rho"))
long$minobs <- factor(long$minobs, c("5", "10", "20"))
p <- ggplot(long, aes(x = factor(window), colour = fit)) +
  facet_grid(metric ~ paste("lambda =", w), scales = "free_y") +
  geom_errorbar(aes(ymin = lo, ymax = hi, group = interaction(minobs, fit)), width = 0,
                position = position_dodge(width = 0.85), alpha = 0.6) +
  geom_point(aes(y = med, shape = minobs, group = interaction(minobs, fit)), position = position_dodge(width = 0.85),
             size = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  scale_colour_manual(NULL, values = fit_colours) +
  scale_shape_manual("minobs", values = c("5" = 16, "10" = 17, "20" = 15)) +
  labs(x = "window (passages)", y = "median and 10-90% range over simulations")
save_plot(p, "supp2a_landscape_accuracy", 200, 130)

## Supp. Fig. 2b, 2f: frequent-karyotype R^2, number of frequent karyotypes ---
s2b <- dfl
s2b$nfq_gt8 <- s2b$nfq > 8
s2b$Rfq_pos <- s2b$Rfq > 0
s2b$R <- clip(s2b$R)
s2b_tab <- summarise_by(s2b, c("fit", "window", "Rfq_pos", "nfq_gt8"), "R")
write_tab(s2b_tab, "supp2b_R_by_nfq_and_Rfq")
s2f <- do.call(rbind, lapply(split(dfl, dfl[c("fit", "window", "w")], drop = TRUE), function(x) data.frame(
  fit = x$fit[1], window = x$window[1], w = x$w[1], fits = nrow(x),
  share_Rfq_pos = round(mean(x$Rfq > 0), 3), share_nfq_gt8 = round(mean(x$nfq > 8), 3),
  share_R_pos = round(mean(x$R > 0), 3),
  share_R_pos_given_Rfq_pos_nfq_gt8 = round(mean(x$R[x$Rfq > 0 & x$nfq > 8] > 0), 3),
  share_R_pos_otherwise = round(mean(x$R[!(x$Rfq > 0 & x$nfq > 8)] > 0), 3))))
write_tab(s2f, "supp2f_flows")

## Supp. Fig. 3: missegregation rate -------------------------------------------
mis <- drop_excluded(label_fits(readRDS(file.path(tab, "misseg_dfl.Rds"))))
mis$r[!is.finite(mis$r)] <- NA
s3 <- summarise_by(mis, c("fit", "window", "pmis"), "r")
write_tab(s3, "supp3_misseg_pearson")
p <- ggplot(mis, aes(x = factor(pmis), y = r, colour = fit)) +
  facet_grid(. ~ paste("window", window)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  scale_colour_manual(NULL, values = fit_colours) +
  labs(x = "missegregation rate (per chromosome per division)", y = "Pearson r, true vs inferred fitness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "supp3_misseg_pearson", 200, 80)

## Supp. Fig. 4: the CV score as a diagnostic -----------------------------------
s4 <- dfl
s4$Rxv <- clip(s4$Rxv)
s4a <- summarise_by(s4, c("fit", "window", "w", "minobs"), "Rxv", function(x) c(q3(x), share_pos = mean(x > 0)))
write_tab(s4a, "supp4a_cv_score")
# 4b: does a positive CV score flag a landscape with R^2 > 0?
s4b <- do.call(rbind, lapply(split(dfl, dfl[c("fit", "window", "w")], drop = TRUE), function(x) {
  pos <- x$Rxv > 0
  data.frame(fit = x$fit[1], window = x$window[1], w = x$w[1], fits = nrow(x), share_xv_pos = round(mean(pos), 3),
             P_R_pos_if_xv_pos = round(mean(x$R[pos] > 0), 3), P_R_pos_if_xv_neg = round(mean(x$R[!pos] > 0), 3),
             spearman_xv_r = round(suppressWarnings(stats::cor(clip(x$Rxv), x$r, method = "spearman",
                                                                 use = "complete")), 3))
}))
write_tab(s4b, "supp4b_cv_flags_accuracy")
# 4c: best CV score per input; accuracy split by its sign
best <- select_best_xv(dfl)
best$xv_pos <- best$Rxv > 0
for (v in c("r", "R", "rho")) best[[v]] <- clip(best[[v]])
s4c <- summarise_by(best, c("fit", "window", "xv_pos"), c("r", "R", "rho"))
write_tab(s4c, "supp4c_selected_by_cv")

## Supp. Fig. 5, Fig. 1g: forecasts ----------------------------------------------
sel <- select_best_xv(dfp)
s5b <- do.call(rbind, lapply(split(sel, sel[c("fit", "window", "metric", "passage", "pos_xv")], drop = TRUE),
                             function(x) data.frame(fit = x$fit[1], window = x$window[1], metric = x$metric[1],
                                                    passage = x$passage[1], xv_pos = x$pos_xv[1], n = nrow(x),
                                                    fwin = round(mean(x$win, na.rm = TRUE), 3))))
write_tab(s5b, "supp5b_beat_baseline_selected")
ang <- sel[sel$metric == "euclidean" & sel$passage %in% c(1, 10), ]
s5a <- summarise_by(ang, c("fit", "window", "passage", "pos_xv"), "angle")
write_tab(s5a, "supp5a_angle_selected")
# (X01 counts passages from 1 at the last time point and plots passage - 1; here
# they count from 0 already.)
f1g <- sel[sel$w %in% c(0.4, 0.8) & sel$pos_xv & sel$metric == "euclidean" & sel$passage >= 1, ]
f1g <- do.call(rbind, lapply(split(f1g, f1g[c("fit", "window", "w", "passage")], drop = TRUE), function(x)
  data.frame(fit = x$fit[1], window = x$window[1], w = x$w[1], passage = x$passage[1], n = nrow(x),
             fwin = round(mean(x$win, na.rm = TRUE), 3))))
write_tab(f1g, "fig1g_beat_baseline")
p <- ggplot(s5b[s5b$xv_pos, ], aes(passage, fwin, colour = fit)) +
  facet_grid(paste("window", window) ~ metric) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line() +
  scale_colour_manual(NULL, values = fit_colours) +
  scale_y_continuous("fraction beating the no-evolution baseline", limits = c(0, 1)) +
  labs(x = "passages after the last time point", subtitle = "minobs with the best CV score of each fit; CV score > 0")
save_plot(p, "supp5b_beat_baseline", 200, 120)
p <- ggplot(ang, aes(angle, colour = fit)) +
  facet_grid(paste("passage", passage) ~ paste("window", window)) +
  stat_ecdf(linewidth = 0.4) +
  geom_line(data = local({
    N <- 22
    a <- 0:180
    dens <- function(t) sin(t)^(N - 2) / stats::integrate(function(t) sin(t)^(N - 2), 0, pi)$value
    data.frame(angle = a, CDF = vapply(a * pi / 180, function(t) stats::integrate(dens, 0, t)$value, numeric(1)))
  }), aes(angle, CDF), inherit.aes = FALSE, linetype = 2, colour = "black", linewidth = 0.3) +
  scale_colour_manual(NULL, values = fit_colours) +
  labs(x = "angle between predicted and observed change (degrees)", y = "cumulative fraction",
       subtitle = "minobs with the best CV score of each fit; dashed: random directions in 22 dimensions")
save_plot(p, "supp5a_angle", 200, 100)

## Supp. Note 3.2: frequencies at the first time point --------------------------------
ini <- drop_excluded(label_fits(readRDS(file.path(tab, "abm_initial.Rds"))))
ini <- ini[ini$source == "abm", ]
ini$late_share <- ini$late_n / ini$nfq
write_tab(summarise_by(ini, c("fit", "window"), c("correlation", "rmse", "late_share", "late_median_est")),
          "note32_initial_frequencies")

## Paired comparison: same input, no selection ---------------------------------
# One input = simulation, window and minobs; each pair of `comparisons`. The
# columns *_a2 hold the first fit (`fit`), *_ref the reference.
paired <- function(fit, ref) {
  key <- c("sim_id", "window", "minobs")
  a <- dfl[dfl$fit == fit, c(key, "w", lscape_cols, "nfq")]
  b <- dfl[dfl$fit == ref, c(key, lscape_cols)]
  m <- merge(a, b, by = key, suffixes = c(".a2", ".ref"))
  do.call(rbind, lapply(split(m, m[c("window", "w")], drop = TRUE), function(x) {
    row <- data.frame(fit = fit, reference = ref, window = x$window[1], w = x$w[1], inputs = nrow(x))
    for (v in c("r", "R", "rho", "rfq", "rnn", "rd2", "Rfq", "Rnn", "Rd2")) {
      a2 <- clip(x[[paste0(v, ".a2")]])
      rf <- clip(x[[paste0(v, ".ref")]])
      ok <- is.finite(a2) & is.finite(rf)
      row[[paste0(v, "_a2")]] <- round(stats::median(a2[ok]), 3)
      row[[paste0(v, "_ref")]] <- round(stats::median(rf[ok]), 3)
      row[[paste0(v, "_win")]] <- round(mean(a2[ok] > rf[ok]), 3)
      row[[paste0(v, "_p")]] <- if (sum(ok) > 5) signif(suppressWarnings(stats::wilcox.test(a2[ok], rf[ok],
                                                                                    paired = TRUE)$p.value), 2) else NA
    }
    row
  }))
}
pair_tab <- do.call(rbind, lapply(comparisons, function(cmp) paired(cmp[1], cmp[2])))
utils::write.table(pair_tab, file.path(out, "paired_landscape.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Forecasts, paired: same input and minobs, both forecasts exist.
fk <- c("sim_id", "window", "minobs", "metric", "passage")
pf <- function(fit, ref) {
  a <- dfp[dfp$fit == fit, c(fk, "w", "win", "angle")]
  b <- dfp[dfp$fit == ref, c(fk, "win", "angle")]
  m <- merge(a, b, by = fk, suffixes = c(".a2", ".ref"))
  m <- m[m$passage %in% c(1, 5, 10), ]
  do.call(rbind, lapply(split(m, m[c("window", "metric", "passage")], drop = TRUE), function(x) data.frame(
    fit = fit, reference = ref, window = x$window[1], metric = x$metric[1], passage = x$passage[1], inputs = nrow(x),
    fwin_a2 = round(mean(x$win.a2), 3), fwin_ref = round(mean(x$win.ref), 3),
    angle_a2 = round(stats::median(x$angle.a2, na.rm = TRUE), 1),
    angle_ref = round(stats::median(x$angle.ref, na.rm = TRUE), 1),
    mcnemar_p = signif(tryCatch(stats::mcnemar.test(table(factor(x$win.a2, c(FALSE, TRUE)),
                                                          factor(x$win.ref, c(FALSE, TRUE))))$p.value,
                                error = function(e) NA_real_), 2))))
}
pair_fc <- do.call(rbind, lapply(comparisons, function(cmp) pf(cmp[1], cmp[2])))
utils::write.table(pair_fc, file.path(out, "paired_forecast.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Stability with the time between the two time points: paired medians by window.
st <- pair_tab[pair_tab$fit == "KARST 2 time points" & pair_tab$reference == "alfak dev2@f93b2d2 2 time points", ]
st <- do.call(rbind, lapply(c("r", "R"), function(v) rbind(
  data.frame(metric = v, window = st$window, w = st$w, fit = "KARST 2 time points", value = st[[paste0(v, "_a2")]]),
  data.frame(metric = v, window = st$window, w = st$w, fit = "alfak dev2@f93b2d2 2 time points", value = st[[paste0(v, "_ref")]]))))
p <- ggplot(st, aes(window, value, colour = fit)) +
  facet_grid(metric ~ paste("lambda =", w), scales = "free_y") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_line() + geom_point(size = 1) +
  scale_colour_manual(NULL, values = fit_colours) +
  scale_x_continuous("window (passages; time points 1, 3 and 7 passages apart)", breaks = c(2, 4, 8)) +
  labs(y = "median over paired inputs (all minobs)")
save_plot(p, "stability_by_window", 200, 90)

## Console summary ---------------------------------------------------------------
cat("\n== Inputs karst cannot fit (decided before fitting) and karst's domain failures\n")
print(exclusion_check, row.names = FALSE)
if (exists("hist_tab")) {
  cat("\n== karst: chosen history and convergence (minobs 10)\n"); print(hist_tab, row.names = FALSE)
}
cat("\n== Fit failures (share)\n")
print(stats::xtabs(share_failed ~ fit + window, aggregate(share_failed ~ fit + window, fail, mean)), digits = 3)
cat("\n== Paired landscape accuracy: median fit / median reference / share fit better\n")
print(pair_tab[c("fit", "reference", "window", "w", "inputs", "r_a2", "r_ref", "r_win", "R_a2", "R_ref", "R_win",
                 "rho_win", "rfq_win", "rnn_win", "rd2_win")], row.names = FALSE)
cat("\n== Paired forecasts: share beating baseline (fit / reference)\n")
print(pair_fc[pair_fc$metric %in% c("euclidean", "wasserstein"), ], row.names = FALSE)
cat("\n== Fig. 1g (selected by CV score, CV > 0, euclidean): fraction beating baseline\n")
print(f1g[f1g$passage %in% c(1, 5, 10), ], row.names = FALSE)
cat("\n== Supp. Fig. 3 (median Pearson r by missegregation rate)\n")
print(stats::reshape(s3[c("fit", "window", "pmis", "r_median")], idvar = c("fit", "window"), timevar = "pmis",
                     direction = "wide"), row.names = FALSE)
cat("\n== Supp. Note 3.2 (initial frequencies)\n")
print(read.delim(file.path(out, "note32_initial_frequencies.tsv"))[, 1:8], row.names = FALSE)
