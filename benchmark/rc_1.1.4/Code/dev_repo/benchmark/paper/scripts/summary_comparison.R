# karst against alfak at a glance: two time points and all time points.
#
#   Rscript scripts/summary_comparison.R
#
# After collect_abm.R, collect_salehi.R, analyse_salehi.R (per_fit.Rds) and
# collect_resources.R. Writes to BENCH_RESULTS/analysis/summary/ and the
# figure BENCH_RESULTS/figures/summary_comparison.{pdf,png}.
#
# Two groups of inputs, each method fitted to the same time points:
#   2 time points    first and last time point: ABM and missegregation windows of 2,
#                4 and 8 passages (w*_ends), every Salehi trajectory
#   all time points  every time point, inputs with more than two time points only: ABM
#                and missegregation windows of 4 and 8 passages (w4_all,
#                w8_all), Salehi trajectories of 3 to 10 passages
# Accuracy is compared on the inputs where every method of the group has a
# result (paired); failures, time and memory over every input. With karst's
# other multi-time-point methods in the round (karst_pairwise,
# karst_state_space; "all time points" only, R/labels.R) the methods are karst
# (2 time points: first and last; all time points: the default, "anchored"), karst
# pairwise, karst state space and alfak; Salehi "all time points" is also shown
# on the trajectories where anchored left out time points in between (sweeps). Inputs that karst
# cannot fit (R/coverage.R; per design, on the time points it fits) are left out
# for both methods, as in analyse_abm.R and analyse_salehi.R. With karst's
# missegregation-rate arms (ac36fd1 round: karst_beta_auto on both groups,
# karst_beta_fixed = 5e-5 on the missegregation simulations, whose karst
# fits take the true rate) each is also compared with its fixed-rate
# counterpart (`reference`: karst_beta_fixed where present, else karst), the
# Salehi sweeps are shown for both groups, and the time of beta = "auto" per
# input is that of all its candidates (tables/beta_picks.tsv).
#   summary_accuracy.tsv     per metric, group (and ABM window) and method:
#                            pairs, median (or mean share) and its range, the
#                            shares of pairs where the method is better and
#                            worse than alfak (the rest: ties) and a paired
#                            test (better_than_alfak, worse_than_alfak,
#                            p_vs_alfak), and for pairwise and state space
#                            the same against anchored, for the
#                            missegregation-rate arms against `reference`
#                            (columns *_anchored)
#   summary_performance.tsv  per data set, group and method: fits, failed
#                            fits, fit time per input, peak memory per input
#                            (of the task: fit, forecasts and scores; alfak
#                            one task per minobs, the largest is taken;
#                            karst one task for all minobs)
#   summary_dropped.tsv      karst's multi-time-point fits per data set, design
#                            and method: inputs fitted, inputs with time points
#                            left out, time points left out, inputs left with
#                            the first and last time point only (fitted as two
#                            time points)
#   summary_cv.tsv           CV score (R2R) per data set, group and method:
#                            median and share > 0. R2R is measured against
#                            each method's own full-data estimates, so it is
#                            not comparable between the two methods.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
tab <- file.path(paths$results, "tables")
out <- file.path(paths$results, "analysis", "summary")
fig <- file.path(paths$results, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
set.seed(1)
methods <- c("karst", karst_other_methods, names(karst_beta_methods), "alfak")
method_colours <- c(karst = "#D55E00", karst_pairwise = "#CC79A7", karst_state_space = "#7B3294",
                    karst_beta_auto = "#1B7837", karst_beta_fixed = "#762A83", alfak = "#0072B2")
# the fixed-rate fit a method is also compared with, NA for none
reference_of <- function(mt, arms) {
  if (mt %in% karst_other_methods || mt == "karst_beta_fixed") return("karst")
  if (mt == "karst_beta_auto") return(if ("karst_beta_fixed" %in% arms) "karst_beta_fixed" else "karst")
  NA_character_
}
groups <- c("2 time points", "all time points")
write_tab <- function(d, name) {
  # labels with line breaks (figure axes) are written on one line
  for (v in names(d)[vapply(d, is.character, logical(1))]) d[[v]] <- gsub("\n", " ", d[[v]], fixed = TRUE)
  utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(d)
}
clip <- function(x) pmax(-1, x)
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- 1.959964
  p <- k / n
  c(p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)),
    p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)
}
paired_p <- function(a, b) signif(tryCatch(stats::wilcox.test(a, b, paired = TRUE)$p.value,
                                            error = function(e) NA_real_), 2)

## Paired accuracy rows ----------------------------------------------------------------
# d: one row per fit with columns fit, the pairing keys and `value`; `higher`:
# a larger value is better. Continuous values: median and 25-75% range;
# shares (`share` TRUE, value 0/1 or a per-fit share): mean and a 95% interval.
arm_fit <- function(mt, g) paste(method_label[[mt]], g)
accuracy_rows <- function(d, keys, panel, x, share = FALSE, higher = TRUE) {
  do.call(rbind, lapply(groups, function(g) {
    arms <- methods[vapply(methods, function(mt) any(d$fit == arm_fit(mt, g)), logical(1))]
    if (!all(c("karst", "alfak") %in% arms)) return(NULL)
    vals <- lapply(arms, function(mt) {
      v <- d[d$fit == arm_fit(mt, g), c(keys, "value")]
      names(v)[names(v) == "value"] <- mt
      v
    })
    m <- Reduce(function(a, b) merge(a, b, by = keys), vals)
    m <- m[Reduce(`&`, lapply(arms, function(mt) is.finite(m[[mt]]))), , drop = FALSE]
    do.call(rbind, lapply(arms, function(mt) {
      v <- m[[mt]]
      n <- length(v)
      est <- if (!n) NA_real_ else if (share) mean(v) else stats::median(v)
      rng <- if (!n) c(NA_real_, NA_real_) else if (share && all(v %in% c(0, 1))) wilson(sum(v), n) else if (share)
        mean(v) + c(-1, 1) * 1.959964 * stats::sd(v) / sqrt(n) else stats::quantile(v, c(0.25, 0.75), names = FALSE)
      # shares of pairs where this method is better and worse than `ref`
      # (the rest are ties), paired test
      vs <- function(ref) {
        if (mt == ref || !n || !ref %in% arms) return(rep(NA_real_, 3))
        d <- if (higher) v - m[[ref]] else m[[ref]] - v
        c(round(mean(d > 0), 3), round(mean(d < 0), 3), if (n > 5) paired_p(v, m[[ref]]) else NA_real_)
      }
      a <- vs("alfak")
      ref <- reference_of(mt, arms)
      b <- if (!is.na(ref)) vs(ref) else rep(NA_real_, 3)
      data.frame(panel = panel, x = x(g), group = g, method = mt, pairs = n, value = est, lo = rng[1], hi = rng[2],
                 better_than_alfak = a[1], worse_than_alfak = a[2], p_vs_alfak = a[3],
                 better_than_anchored = b[1], worse_than_anchored = b[2], p_vs_anchored = b[3],
                 reference = if (!is.na(ref) && ref %in% arms) ref else NA_character_, stringsAsFactors = FALSE)
    }))
  }))
}

## Inputs karst cannot fit ---------------------------------------------------------------
# Per input and design (R/coverage.R): a fit is left out when its own input
# is not analysable, for both methods (as analyse_abm.R, analyse_salehi.R).
abm_inputs <- readRDS(file.path(tab, "abm_inputs.Rds"))
sal_inputs <- readRDS(file.path(tab, "salehi_inputs.Rds"))
if (is.null(sal_inputs$design)) sal_inputs$design <- "ends"  # tables of earlier rounds
sal_inputs$source <- "salehi"
excluded_keys <- c(not_analysable_keys(abm_inputs, c("source", "sim_id", "design")),
                   not_analysable_keys(sal_inputs, c("source", "fi", "design")))
excluded <- function(source, dataset, design) input_key(source, dataset, design) %in% excluded_keys

## ABM and missegregation landscapes --------------------------------------------------
abm_x <- function(window) function(g) paste0(g, "\nwindow ", window)
landscape <- function(file, source_name, metrics) {
  d <- label_abm_fits(readRDS(file.path(tab, file)))
  d <- d[d$method %in% methods & !excluded(d$source, d$sim_id, d$design), ]
  do.call(rbind, lapply(names(metrics), function(v) {
    do.call(rbind, lapply(c(2L, 4L, 8L), function(w) {
      x <- d[d$window == w, ]
      x$value <- clip(x[[v]])
      r <- accuracy_rows(x, c("sim_id", "minobs"), paste0(source_name, ": ", metrics[[v]]), abm_x(w))
      # windows of two time points have no separate "all time points" fit
      if (w == 2L) r <- r[r$group == "2 time points", ]
      if (!is.null(r)) {
        r$window <- w
        r$subset <- "all"
      }
      r
    }))
  }))
}
acc <- rbind(
  landscape("abm_dfl.Rds", "ABM", c(r = "landscape Pearson r", R = "landscape rescaled R^2",
                                    rho = "landscape Spearman rho")),
  landscape("misseg_dfl.Rds", "Missegregation", c(r = "landscape Pearson r")))

# ABM forecasts: does the forecast beat the no-evolution baseline 10 passages
# after the last time point (Euclidean and Wasserstein distance)?
dfp <- label_abm_fits(readRDS(file.path(tab, "abm_dfp.Rds")))
dfp <- dfp[dfp$method %in% methods & dfp$passage == 10 & !excluded(dfp$source, dfp$sim_id, dfp$design), ]
for (mt in c("euclidean", "wasserstein")) {
  for (w in c(2L, 4L, 8L)) {
    x <- dfp[dfp$window == w & dfp$metric == mt, ]
    x$value <- as.numeric(x$win)
    r <- accuracy_rows(x, c("sim_id", "minobs"), paste0("ABM: forecast beats baseline, passage 10 (", mt, ")"),
                       abm_x(w), share = TRUE)
    if (w == 2L) r <- r[r$group == "2 time points", ]
    r$window <- w
    r$subset <- "all"
    acc <- rbind(acc, r)
  }
}

## Salehi ------------------------------------------------------------------------------
# per_fit.Rds: analysable inputs only (analyse_salehi.R)
multi <- function(fi) !grepl("_l_2_", fi)  # trajectories of more than 2 passages
salehi_x <- function(g) if (g == "2 time points") "2 time points" else "all time points\n(3-10 passages)"
# Sweeps: trajectories where anchored left out time points in between (no cell on
# a major; diagnostics$dropped_samples of the "all" fit).
sweep_x <- function(g) paste0(g, "\n(sweeps)")
beta_round <- any(names(karst_beta_methods) %in% unique(readRDS(file.path(tab, "salehi_runs.Rds"))$method))
sal_runs_all <- readRDS(file.path(tab, "salehi_runs.Rds"))
sweeps <- unique(sal_runs_all$fi[sal_runs_all$method == "karst" & sal_runs_all$design == "all" &
                                   sal_runs_all$dropped_samples %in% seq_len(100)])
salehi_rows <- function(x, keys, panel, share = TRUE) {
  r <- accuracy_rows(x, keys, panel, salehi_x, share = share)
  if (length(sweeps)) {
    rs <- accuracy_rows(x[x$fi %in% sweeps, ], keys, panel, sweep_x, share = share)
    # the two-time-point fits of the same trajectories too, for the
    # missegregation-rate arms, which have their own
    if (!is.null(rs)) r <- rbind(r, rs[rs$group == "all time points" | beta_round, ])
  }
  r$subset <- ifelse(grepl("sweeps", r$x), "sweeps", "all")
  r
}
per_fit <- readRDS(file.path(paths$results, "analysis", "salehi", "per_fit.Rds"))
per_fit <- per_fit[!(fit_group(per_fit$fit) == "all time points" & !multi(per_fit$fi)), ]
for (mt in c("wasserstein", "euclidean")) {
  x <- per_fit[per_fit$metric == mt, ]
  x$value <- x$win
  r <- salehi_rows(x, c("fi", "min_obs"), paste0("Salehi: share of daughters where the forecast beats the parent (",
                                                 mt, ")"))
  r$window <- NA_integer_
  acc <- rbind(acc, r)
}
# Both Salehi panels pair the fits that ran (run status; per_fit.Rds of
# analyse_salehi.R likewise), whatever their CV score.
novel <- readRDS(file.path(tab, "salehi_novel.Rds"))
run_ok <- with(sal_runs_all[sal_runs_all$status == "ok", ], paste(fi, design, method, min_obs))
novel <- novel[paste(novel$fi, novel$design, novel$method, novel$min_obs) %in% run_ok, ]
novel <- label_salehi_fits(novel)
novel <- novel[novel$ids == "f" & novel$method %in% methods & !excluded("salehi", novel$fi, novel$design), ]
novel <- novel[!(fit_group(novel$fit) == "all time points" & !multi(novel$fi)), ]
novel$value <- as.numeric(novel$Estimate > 0 & novel$p < 0.05)
r <- salehi_rows(novel, c("fi", "min_obs", "dec_id"), "Salehi: novel karyotypes, fitness predicts emergence")
r$window <- NA_integer_
acc <- rbind(acc, r)
acc$value <- round(acc$value, 3)
if (is.null(acc$subset)) acc$subset <- "all"
acc$lo <- round(acc$lo, 3)
acc$hi <- round(acc$hi, 3)
write_tab(acc, "summary_accuracy")

## Failures, time and memory ---------------------------------------------------------
design_group <- function(design) ifelse(grepl("all$", design), "all time points", "2 time points")
abm_runs <- readRDS(file.path(tab, "abm_runs.Rds"))
abm_runs <- abm_runs[!excluded(abm_runs$source, abm_runs$sim_id, abm_runs$design), ]
salehi_runs <- readRDS(file.path(tab, "salehi_runs.Rds"))
salehi_runs <- salehi_runs[!excluded("salehi", salehi_runs$fi, salehi_runs$design), ]
fits <- rbind(
  data.frame(source = abm_runs$source, dataset = abm_runs$sim_id, design = abm_runs$design, method = abm_runs$method,
             status = abm_runs$status, stringsAsFactors = FALSE),
  data.frame(source = "salehi", dataset = salehi_runs$fi, design = salehi_runs$design, method = salehi_runs$method,
             status = salehi_runs$status, stringsAsFactors = FALSE))
fits <- fits[fits$method %in% methods, ]
fits$group <- design_group(fits$design)
source_name <- c(abm = "ABM", misseg = "Missegregation", salehi = "Salehi")
resources_path <- file.path(tab, "resources.Rds")
res <- if (file.exists(resources_path)) readRDS(resources_path) else NULL
if (!is.null(res)) {
  res <- res[res$method %in% methods & !excluded(res$source, res$dataset, res$design), ]
  res$group <- design_group(res$design)
  # beta = "auto" fits every candidate rate: its time is theirs together
  # (each candidate's fit plus the picked one's outputs)
  picks_path <- file.path(tab, "beta_picks.tsv")
  if (file.exists(picks_path)) {
    bp <- utils::read.delim(picks_path, stringsAsFactors = FALSE)
    i <- which(res$method == "karst_beta_auto")
    j <- match(paste(res$source[i], res$dataset[i], res$design[i]), paste(bp$source, bp$dataset, bp$design))
    res$fit_seconds[i] <- bp$auto_seconds[j]
  }
}
q <- function(v, p) if (any(is.finite(v))) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
# Peak memory per input (collect_resources.R); below 100 MB the task was not
# measured (an R process alone takes more): Slurm's sampling missed it.
mem_of <- function(r) {
  m <- if (is.null(r$peak_mb)) r$task_peak_mb else r$peak_mb
  m[!is.na(m) & m < 100] <- NA
  m
}
perf <- do.call(rbind, lapply(names(source_name), function(s) do.call(rbind, lapply(groups, function(g) {
  do.call(rbind, lapply(methods, function(mt) {
    f <- fits[fits$source == s & fits$group == g & fits$method == mt, ]
    r <- if (is.null(res)) NULL else res[res$source == s & res$group == g & res$method == mt, ]
    ci <- wilson(sum(f$status != "ok"), nrow(f))
    data.frame(source = s, group = g, method = mt, inputs = length(unique(paste(f$dataset, f$design))),
               fits = nrow(f), failed = sum(f$status != "ok"),
               share_failed = if (nrow(f)) round(mean(f$status != "ok"), 4) else NA_real_,
               share_failed_lo = round(ci[1], 4), share_failed_hi = round(ci[2], 4),
               time_median_s = round(q(r$fit_seconds, 0.5), 1), time_q10_s = round(q(r$fit_seconds, 0.1), 1),
               time_q90_s = round(q(r$fit_seconds, 0.9), 1), time_max_s = round(q(r$fit_seconds, 1), 1),
               memory_inputs = sum(is.finite(mem_of(r))),
               memory_median_mb = round(q(mem_of(r), 0.5)), memory_q10_mb = round(q(mem_of(r), 0.1)),
               memory_q90_mb = round(q(mem_of(r), 0.9)), memory_max_mb = round(q(mem_of(r), 1)),
               karst_fit_memory_median_mb = round(q(r$karst_fit_peak_mb, 0.5)),
               karst_fit_memory_max_mb = round(q(r$karst_fit_peak_mb, 1)), stringsAsFactors = FALSE)
  }))
}))))
perf <- perf[perf$fits > 0, ]
write_tab(perf, "summary_performance")

## Time points left out by karst's multi-time-point fits --------------------------------
dropped_rows <- function(r, src, ds) {
  r <- r[is_karst_method(r$method) & r$status == "ok" & !is.na(r$samples) & r$samples > 2, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  r$dataset <- ds(r)
  r$source <- src(r)
  r <- r[!duplicated(r[c("source", "dataset", "design", "method")]), , drop = FALSE]
  n_left <- ifelse(is.na(r$dropped_samples), 0L, r$dropped_samples)
  fitted <- if (is.null(r$samples_fitted)) r$samples - n_left else ifelse(is.na(r$samples_fitted), r$samples - n_left,
                                                                          r$samples_fitted)
  do.call(rbind, lapply(split(seq_len(nrow(r)), r[c("source", "design", "method")], drop = TRUE), function(i) data.frame(
    source = r$source[i[1]], design = r$design[i[1]], method = r$method[i[1]], inputs = length(i),
    inputs_left_out = sum(n_left[i] > 0), samples_left_out = sum(n_left[i]), two_samples_left = sum(fitted[i] == 2),
    stringsAsFactors = FALSE)))
}
dropped <- rbind(dropped_rows(abm_runs, function(r) r$source, function(r) r$sim_id),
                 dropped_rows(salehi_runs, function(r) rep("salehi", nrow(r)), function(r) r$fi))
if (!is.null(dropped)) write_tab(dropped, "summary_dropped")

## CV score (not comparable between the methods) ------------------------------------------
cv_rows <- function(d, s, value) {
  d$value <- value
  d <- d[d$method %in% methods & is.finite(d$value), ]
  if (!nrow(d)) return(NULL)
  d$group <- fit_group(d$fit)
  do.call(rbind, lapply(split(d, d[c("method", "group")], drop = TRUE), function(x) data.frame(
    source = s, method = x$method[1], group = x$group[1], fits = nrow(x),
    cv_median = round(stats::median(x$value), 3), cv_share_pos = round(mean(x$value > 0), 3),
    stringsAsFactors = FALSE)))
}
abm_dfl <- label_abm_fits(readRDS(file.path(tab, "abm_dfl.Rds")))
abm_dfl <- abm_dfl[!(fit_group(abm_dfl$fit) == "all time points" & abm_dfl$window == 2L) &
                     !excluded(abm_dfl$source, abm_dfl$sim_id, abm_dfl$design), ]
mis_dfl <- label_abm_fits(readRDS(file.path(tab, "misseg_dfl.Rds")))
mis_dfl <- mis_dfl[!(fit_group(mis_dfl$fit) == "all time points" & mis_dfl$window == 2L) &
                     !excluded(mis_dfl$source, mis_dfl$sim_id, mis_dfl$design), ]
sal <- label_salehi_fits(readRDS(file.path(tab, "salehi_scores.Rds")))
sal <- sal[!excluded("salehi", sal$fi, sal$design) & !(fit_group(sal$fit) == "all time points" & !multi(sal$fi)), ]
sal <- sal[!duplicated(sal[c("fit", "fi", "min_obs")]), ]
cv <- rbind(cv_rows(abm_dfl, "abm", abm_dfl$Rxv), cv_rows(mis_dfl, "misseg", mis_dfl$Rxv),
            cv_rows(sal, "salehi", sal$xv))
cv$note <- "R2R against each method's own full-data estimates; not comparable between methods"
write_tab(cv, "summary_cv")

## Figure ------------------------------------------------------------------------------
# Time and memory differ by orders of magnitude between the methods: log10.
short_name <- c(abm = "ABM", misseg = "Misseg.", salehi = "Salehi")
perf_x <- function(s, g) paste0(short_name[s], "\n", g)
time_panel <- "Fit time per input, log10 minutes (median, 10-90%)"
memory_panel <- "Peak memory per input, log10 GB (median, 10-90%; x = max)"
perf_long <- rbind(
  data.frame(panel = "Failed fits (share of input x minobs)", x = perf_x(perf$source, perf$group),
             method = perf$method, value = perf$share_failed, lo = perf$share_failed_lo, hi = perf$share_failed_hi),
  data.frame(panel = time_panel, x = perf_x(perf$source, perf$group), method = perf$method,
             value = log10(perf$time_median_s / 60), lo = log10(perf$time_q10_s / 60),
             hi = log10(perf$time_q90_s / 60)),
  data.frame(panel = memory_panel, x = perf_x(perf$source, perf$group), method = perf$method,
             value = log10(perf$memory_median_mb / 1024), lo = log10(perf$memory_q10_mb / 1024),
             hi = log10(perf$memory_q90_mb / 1024)))
peak <- data.frame(panel = memory_panel, x = perf_x(perf$source, perf$group), method = perf$method,
                   value = log10(perf$memory_max_mb / 1024))
show <- acc[acc$panel %in% c("ABM: landscape Pearson r", "ABM: landscape rescaled R^2",
                             "ABM: forecast beats baseline, passage 10 (euclidean)",
                             "Missegregation: landscape Pearson r",
                             "Salehi: share of daughters where the forecast beats the parent (wasserstein)",
                             "Salehi: share of daughters where the forecast beats the parent (euclidean)",
                             "Salehi: novel karyotypes, fitness predicts emergence"), ]
plot_d <- rbind(show[c("panel", "x", "method", "value", "lo", "hi")], perf_long)
plot_d <- plot_d[is.finite(plot_d$value), ]
plot_d$lo[!is.finite(plot_d$lo)] <- NA
plot_d$hi[!is.finite(plot_d$hi)] <- NA
panel_levels <- c(unique(show$panel), unique(perf_long$panel))
wrap <- function(s) vapply(strwrap(s, 48, simplify = FALSE), paste, character(1), collapse = "\n")
plot_d$panel <- factor(wrap(plot_d$panel), wrap(panel_levels))
peak <- peak[is.finite(peak$value), ]
peak$panel <- factor(wrap(peak$panel), wrap(panel_levels))
x_levels <- unique(c(paste0("2 time points\nwindow ", c(2, 4, 8)), paste0("all time points\nwindow ", c(4, 8)),
                     "2 time points", "2 time points\n(sweeps)", "all time points\n(3-10 passages)",
                     "all time points\n(sweeps)",
                     as.vector(t(outer(short_name, groups, function(s, g) paste0(s, "\n", g))))))
plot_d$x <- factor(plot_d$x, x_levels)
peak$x <- factor(peak$x, x_levels)
p <- ggplot(plot_d, aes(x, value, colour = method)) +
  facet_wrap(~panel, scales = "free", ncol = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = position_dodge(width = 0.5), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 1.6) +
  geom_point(data = peak, shape = 4, size = 1.4, position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_colour_manual(NULL, values = method_colours, breaks = methods,
                      labels = unname(c(if (beta_round) "KARST (beta 5e-5; misseg.: true rate)" else
                        "KARST (all time points: default, anchored)", method_label[setdiff(methods, "karst")]))) +
  scale_y_continuous(NULL) + scale_x_discrete(NULL) +
  labs(caption = paste("Accuracy: paired inputs; points median (bars 25-75%) or mean share (bars 95%).",
                       "Each method fitted to the same time points; 'all time points' = inputs with more than two time points;",
                       "accuracy on the inputs where every method shown has a result. Sweeps: Salehi trajectories where",
                       "anchored (karst's default) left out time points in between.",
                       if (beta_round) paste("karst beta auto: missegregation rate picked per input from 5e-5 x",
                                             "1/16 ... 64 (time: every candidate); all time points: anchored.") else NULL,
                       "CV scores: summary_cv.tsv (not comparable between methods).", sep = "\n")) +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), legend.position = "top", plot.caption = element_text(hjust = 0))
n_panels <- length(unique(plot_d$panel))
save_plot_dims <- c(200, 55 * ceiling(n_panels / 2) + 25)
ggsave(file.path(fig, "summary_comparison.pdf"), p, width = save_plot_dims[1], height = save_plot_dims[2], units = "mm")
ggsave(file.path(fig, "summary_comparison.png"), p, width = save_plot_dims[1], height = save_plot_dims[2], units = "mm",
       dpi = 200, bg = "white")

## Console -----------------------------------------------------------------------------
cat("\n== Accuracy (paired; value = median or mean share)\n")
print(acc[c("panel", "group", "subset", "window", "method", "pairs", "value", "lo", "hi", "better_than_alfak",
           "worse_than_alfak", "p_vs_alfak", "better_than_anchored", "worse_than_anchored", "p_vs_anchored",
           "reference")],
      row.names = FALSE)
cat("\n== Failures, time, memory\n"); print(perf, row.names = FALSE)
if (!is.null(dropped)) { cat("\n== Time points left out by karst's multi-time-point fits\n")
  print(stats::setNames(dropped, sub("samples", "time_points", names(dropped))), row.names = FALSE) }
cat("\n== CV score (not comparable between methods)\n"); print(cv[setdiff(names(cv), "note")], row.names = FALSE)
