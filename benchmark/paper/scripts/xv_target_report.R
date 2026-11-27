# Tables of the cross-validation study (docs/xv_target_6ae2a24.md): KARST's
# cross-validation scored against its default targets ("own": its estimate
# from all the data), against alfak's ("alfak_target",
# scripts/xv_alfak_target.R) and against targets from each karyotype's own
# counts ("own_target", resampled once; "own_point", not resampled;
# scripts/xv_own_target.R), and the paper's tests that select fits by the
# score, from the default analyses (analysis/) and those run with
# XV_SOURCE=<source> (analysis_xv_<source>/).
#
#   Rscript scripts/xv_target_report.R
#
# Uses every source whose table exists (and, for the paper's tests, whose
# analyses exist). Writes BENCH_RESULTS/analysis_xv_targets/, each table with
# a column `xv` naming the score:
#   scores.tsv            the scores by data, method and design: fits by
#                         status, medians and shares above 0, Pearson's r of
#                         target and prediction
#   truth_scores.tsv      simulations: Spearman's rho of each score with the
#                         landscape's true accuracy (Pearson r), within each
#                         wavelength and window (median over them), and pooled
#   truth_targets.tsv     simulations: Spearman's rho of each kind of target
#                         with the true fitness, per fit (median, quartiles):
#                         KARST's own-count target, without resampling, its
#                         default target, its held-out prediction; alfak's
#                         target (all it tested; those KARST tested with cells)
#   forecast_salehi.tsv   Salehi: Spearman's rho of each score with the share
#                         of daughter passages a fit's forecast beats
#                         (Wasserstein, the paper's evaluation day)
#   fig2d_swap.tsv        Fig. 2d (Spearman of observed and predicted change)
#                         on each method's selection, on the other's, and on
#                         every row both have
#   compare_<test>.tsv    the paper's tests that select by the score, from
#                         every run
# and stops if a row of alfak or of the paper's fits differs between runs:
# only KARST's scores change.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
tabdir <- file.path(paths$results, "tables")
root_of <- function(src) {
  Sys.setenv(XV_SOURCE = if (src == "own") "" else src)
  on.exit(Sys.setenv(XV_SOURCE = ""))
  analysis_root(paths)
}
spec <- xv_source_tables
sources <- names(spec)[vapply(names(spec), function(s) file.exists(file.path(tabdir, paste0(spec[[s]][["table"]], ".Rds"))),
                              logical(1))]
with_analysis <- sources[vapply(sources, function(s) dir.exists(root_of(s)), logical(1))]
cat("score sources:", paste(sources, collapse = ", "), " with analyses:", paste(with_analysis, collapse = ", "), "\n")
out <- file.path(paths$results, "analysis_xv_targets")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tab <- function(d, name) {
  utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(d)
}
med <- function(x) round(stats::median(x, na.rm = TRUE), 3)
q <- function(x, p) round(stats::quantile(x, p, na.rm = TRUE, names = FALSE), 3)
pos <- function(x) round(mean(x > 0, na.rm = TRUE), 3)
rho <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5) return(NA_real_)
  round(suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman")), 3)
}
group_of <- function(design) ifelse(grepl("ends$", design), "2 time points", "all time points")
excluded_keys <- c(not_analysable_keys(readRDS(file.path(tabdir, "abm_inputs.Rds")), c("source", "sim_id", "design")),
                   not_analysable_keys(transform(readRDS(file.path(tabdir, "salehi_inputs.Rds")), source = "salehi"),
                                       c("source", "fi", "design")))
keep <- function(d) d[!input_key(d$source, d$dataset, d$design) %in% excluded_keys, , drop = FALSE]

# One row per KARST fit and minobs: the new score of each source, with
# KARST's default score (xv_own) and, from the alfak-target table, alfak's.
alfak_tab <- if ("alfak_target" %in% sources) keep(readRDS(file.path(tabdir, "xv_alfak_target.Rds"))) else NULL
score_rows <- lapply(sources, function(src) {
  d <- keep(readRDS(file.path(tabdir, paste0(spec[[src]][["table"]], ".Rds"))))
  col <- spec[[src]][["column"]]
  pear <- switch(src, alfak_target = "pearson_alfak_target", own_target = "pearson_own_target",
                 own_point = "pearson_own_point")
  data.frame(xv = src, source = d$source, dataset = d$dataset, design = d$design, method = d$method,
             minobs = d$minobs, group = group_of(d$design), status = d$status, score = d[[col]], pearson = d[[pear]],
             xv_own = d$xv_own, stringsAsFactors = FALSE)
})
names(score_rows) <- sources
key5 <- function(d) paste(d$source, d$dataset, d$design, d$method, d$minobs)

## Scores -----------------------------------------------------------------------------
scores <- do.call(rbind, lapply(score_rows, function(d) {
  do.call(rbind, lapply(split(d, d[c("source", "method", "group")], drop = TRUE), function(x) {
    ok <- x[x$status == "ok", ]
    a <- if (is.null(alfak_tab)) NULL else alfak_tab[match(key5(ok), key5(alfak_tab)), ]
    data.frame(xv = x$xv[1], source = x$source[1], method = x$method[1], group = x$group[1], fits = nrow(x),
               ok = nrow(ok), without_score = nrow(x) - nrow(ok), own_median = med(ok$xv_own), own_pos = pos(ok$xv_own),
               score_median = med(ok$score), score_pos = pos(ok$score), pearson_median = med(ok$pearson),
               alfak_own_median = if (is.null(a)) NA else med(a$alfak_xv_own),
               alfak_own_pos = if (is.null(a)) NA else pos(a$alfak_xv_own),
               rho_own_score = rho(ok$xv_own, ok$score), stringsAsFactors = FALSE)
  }))
}))
write_tab(scores, "scores")

## Do the scores follow the true accuracy? (simulations) -----------------------------------
truth <- rbind(readRDS(file.path(tabdir, "abm_dfl.Rds")), readRDS(file.path(tabdir, "misseg_dfl.Rds")))
truth <- truth[c("source", "sim_id", "design", "method", "minobs", "w", "r", "Rxv")]
truth_key <- paste(truth$source, truth$sim_id, truth$design, truth$method, truth$minobs)
truth_scores <- do.call(rbind, lapply(score_rows, function(d) {
  d <- d[d$source %in% c("abm", "misseg") & d$status == "ok", ]
  m <- match(key5(d), truth_key)
  d$r <- truth$r[m]
  d$w <- truth$w[m]
  am <- match(paste(d$source, d$dataset, d$design, "alfak", d$minobs), truth_key)
  d$r_alfak <- truth$r[am]
  d$xv_alfak <- truth$Rxv[am]
  d$window <- as.integer(sub("^w(\\d+)_.*$", "\\1", d$design))
  do.call(rbind, lapply(split(d, d[c("source", "method", "group")], drop = TRUE), function(x) {
    cells <- split(x, x[c("w", "window")], drop = TRUE)
    within <- function(sc, r) round(stats::median(vapply(cells, function(y) rho(y[[sc]], y[[r]]), numeric(1)),
                                                   na.rm = TRUE), 3)
    data.frame(xv = x$xv[1], source = x$source[1], method = x$method[1], group = x$group[1], fits = nrow(x),
               cells = length(cells), score_within = within("score", "r"), own_within = within("xv_own", "r"),
               alfak_within = within("xv_alfak", "r_alfak"), score_pooled = rho(x$score, x$r),
               own_pooled = rho(x$xv_own, x$r), alfak_pooled = rho(x$xv_alfak, x$r_alfak), stringsAsFactors = FALSE)
  }))
}))
write_tab(truth_scores, "truth_scores")

## Do the targets follow the true fitness? (simulations, own-count table) --------------------
if ("own_target" %in% sources) {
  o <- keep(readRDS(file.path(tabdir, "xv_own_target.Rds")))
  o <- o[o$source %in% c("abm", "misseg") & o$status == "ok", ]
  o$group <- group_of(o$design)
  o$w <- vapply(o$dataset, function(s) abm_sim_info(s)$w, numeric(1))
  cols <- c(own = "rho_truth_own", own_point = "rho_truth_point", model = "rho_truth_model",
            prediction = "rho_truth_prediction", alfak = "rho_truth_alfak", alfak_common = "rho_truth_alfak_common")
  tt <- do.call(rbind, lapply(split(o, o[c("source", "method", "group", "w")], drop = TRUE), function(x) {
    do.call(rbind, lapply(names(cols), function(nm) data.frame(
      source = x$source[1], method = x$method[1], group = x$group[1], w = x$w[1], target = nm, fits = sum(is.finite(x[[cols[[nm]]]])),
      q25 = q(x[[cols[[nm]]]], 0.25), median = med(x[[cols[[nm]]]]), q75 = q(x[[cols[[nm]]]], 0.75),
      stringsAsFactors = FALSE)))
  }))
  write_tab(tt, "truth_targets")
}

## Do the scores follow the forecasts? (Salehi) ----------------------------------------------
sal <- readRDS(file.path(tabdir, "salehi_scores.Rds"))
meta <- utils::read.csv(paths$salehi_meta)
sal$type <- salehi_lineage_fig2(sal$fi, meta)
eval_day <- ifelse(sal$type %in% c("p53 k.o", "p53 w.t"), 25, 15)
ev <- sal[sal$metric %in% "wasserstein" & sal$tt == eval_day & !is.na(sal$test_treat) & sal$test_treat != "NaN", ]
ev$win <- ev$base > ev$pred
per_fit <- stats::aggregate(win ~ design + method + fi + min_obs, ev, mean)
per_fit <- per_fit[!input_key("salehi", per_fit$fi, per_fit$design) %in% excluded_keys, ]
per_fit$group <- group_of(per_fit$design)
alfak_xv <- sal[sal$method == "alfak" & !duplicated(sal[c("design", "fi", "min_obs")]), ]
al <- per_fit[per_fit$method == "alfak", ]
al$score <- alfak_xv$xv[match(paste(al$design, al$fi, al$min_obs), paste(alfak_xv$design, alfak_xv$fi, alfak_xv$min_obs))]
fc_row <- function(xv, x) data.frame(xv = xv, method = x$method[1], group = x$group[1], fits = sum(is.finite(x$score)),
                                     rho_score_win = rho(x$score, x$win),
                                     win_if_pos = round(mean(x$win[x$score > 0], na.rm = TRUE), 3),
                                     win_if_not = round(mean(x$win[x$score <= 0], na.rm = TRUE), 3), stringsAsFactors = FALSE)
forecast <- list(do.call(rbind, lapply(split(al, al$group), function(x) fc_row("alfak (own)", x))))
ka <- per_fit[per_fit$method %in% karst_xv_methods, ]
pk <- function(d) paste("salehi", d$fi, d$design, d$method, d$min_obs)
first <- score_rows[[1]]
x0 <- ka
x0$score <- first$xv_own[match(pk(x0), key5(first))]
forecast[[2]] <- do.call(rbind, lapply(split(x0, x0[c("method", "group")], drop = TRUE), function(x) fc_row("own", x)))
for (src in sources) {
  d <- score_rows[[src]]
  x <- ka
  x$score <- d$score[match(pk(x), key5(d))]
  forecast[[length(forecast) + 1L]] <- do.call(rbind, lapply(split(x, x[c("method", "group")], drop = TRUE),
                                                             function(y) fc_row(src, y)))
}
write_tab(do.call(rbind, forecast), "forecast_salehi")

## Fig. 2d on swapped selections ----------------------------------------------------------------
# The paper's rule (analyse_salehi.R): fits with a score > 0 and a test
# passage, the best minobs of each trajectory, then the best trajectory of
# each start; Fig. 2d at the evaluation day (Wasserstein).
select_paper <- function(x) {
  x <- x[!is.na(x$xv) & x$xv > 0 & !is.na(x$test_treat) & x$test_treat != "NaN", ]
  if (!nrow(x)) return(NULL)
  x <- select_by_cv(x, "fi", "min_obs", "xv")
  x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
  select_by_cv(x, "clean_id", "fi", "xv")
}
fig2d_rows <- function(d) {
  d <- d[d$metric %in% "wasserstein" & !is.na(d$test_treat) & d$test_treat != "NaN", ]
  d[d$tt == ifelse(d$type %in% c("p53 k.o", "p53 w.t"), 25, 15), ]
}
runs <- readRDS(file.path(tabdir, "salehi_runs.Rds"))
swap <- list()
for (src in c("own", sources)) {
  Sys.setenv(XV_SOURCE = if (src == "own") "" else src)
  s <- read_bench_table("salehi_scores", paths)
  Sys.setenv(XV_SOURCE = "")
  s$fit_ok <- salehi_fit_ok(s, runs)
  s <- label_salehi_fits(s)
  s$type <- salehi_lineage_fig2(s$fi, meta)
  s <- s[!input_key("salehi", s$fi, s$design) %in% excluded_keys, ]
  for (g in c("2 time points", "all time points")) {
    fits <- list(KARST = s[s$fit == paste("KARST", g), ], alfak = s[s$fit == paste(method_label[["alfak"]], g), ])
    sel <- lapply(fits, select_paper)
    pair <- function(rows_of) {
      chosen <- unique(rows_of[c("fi", "min_obs")])
      one <- function(d) {
        d <- fig2d_rows(d[d$fit_ok & paste(d$fi, d$min_obs) %in% paste(chosen$fi, chosen$min_obs), ])
        d$bk <- round(d$base, 8)
        d
      }
      m <- merge(one(fits$KARST), one(fits$alfak), by = c("fi", "min_obs", "tt", "test_treat", "type", "bk"),
                 suffixes = c(".k", ".a"))
      if (anyDuplicated(m[c("fi", "min_obs", "bk")])) stop("Fig. 2d rows of a daughter passage are not unique")
      list(m = m, windows = nrow(chosen), one_interval = sum(grepl("_l_2_", chosen$fi)))
    }
    sp <- function(m, side) rho(m[[paste0("base.", side)]], m[[paste0("pred0.", side)]])
    boot <- function(m) {
      set.seed(1)
      g2 <- split(seq_len(nrow(m)), m$type)
      d <- replicate(2000, {
        i <- unlist(g2[sample(length(g2), length(g2), TRUE)])
        sp(m[i, ], "k") - sp(m[i, ], "a")
      })
      paste(round(stats::quantile(d, c(0.025, 0.975), na.rm = TRUE), 2), collapse = " to ")
    }
    e <- lapply(sel, function(x) if (is.null(x)) NULL else fig2d_rows(x))
    swap[[length(swap) + 1L]] <- data.frame(xv = src, group = g, rows_from = "each its own selection", windows = NA,
      one_interval = NA, rows = paste(nrow(e$KARST), nrow(e$alfak), sep = " / "),
      KARST = rho(e$KARST$base, e$KARST$pred0), alfak = rho(e$alfak$base, e$alfak$pred0), difference = NA,
      interval = NA, stringsAsFactors = FALSE)
    for (who in c("KARST", "alfak", "every row")) {
      rows_of <- if (who == "every row") {
        merge(unique(fits$KARST[fits$KARST$fit_ok, c("fi", "min_obs")]), unique(fits$alfak[fits$alfak$fit_ok, c("fi", "min_obs")]))
      } else sel[[who]]
      p <- pair(rows_of)
      swap[[length(swap) + 1L]] <- data.frame(xv = src, group = g,
        rows_from = if (who == "every row") "every row both have" else paste(who, "selection"),
        windows = p$windows, one_interval = p$one_interval, rows = as.character(nrow(p$m)),
        KARST = sp(p$m, "k"), alfak = sp(p$m, "a"), difference = round(sp(p$m, "k") - sp(p$m, "a"), 3),
        interval = boot(p$m), stringsAsFactors = FALSE)
    }
  }
}
write_tab(do.call(rbind, swap), "fig2d_swap")

## The paper's tests in every run ------------------------------------------------------------------
# alfak's and the paper's rows must be the same in every run.
fixed_fit <- function(x) grepl("^alfak|^paper", x)
compare <- function(rel, name, fit_col = "fit") {
  read <- function(src) utils::read.delim(file.path(root_of(src), rel), check.names = FALSE, stringsAsFactors = FALSE)
  base <- read("own")
  ord <- function(d) {
    d <- d[fixed_fit(d[[fit_col]]), , drop = FALSE]
    d <- d[do.call(order, unname(as.list(d))), , drop = FALSE]
    rownames(d) <- NULL
    d
  }
  all_runs <- list(cbind(xv = "own", base))
  for (src in with_analysis) {
    b <- read(src)
    if (!isTRUE(all.equal(ord(base), ord(b), check.attributes = FALSE))) {
      stop(rel, ": alfak's or the paper's rows differ between the default run and ", src)
    }
    all_runs[[length(all_runs) + 1L]] <- cbind(xv = src, b)
  }
  write_tab(do.call(rbind, all_runs), paste0("compare_", name))
}
if (length(with_analysis)) {
  compare("salehi/fig2b_cv_median_ci.tsv", "fig2b")
  compare("salehi/fig2c_angle_test.tsv", "fig2c")
  compare("salehi/fig2d_wasserstein.tsv", "fig2d")
  compare("salehi/fig2f_beat_baseline_eval_pooled.tsv", "fig2f")
  compare("salehi/fig3ef_novel_karyotypes.tsv", "fig3ef")
  compare("fig456/conclusions.tsv", "fig456", fit_col = "arm_label")
  compare("abm/supp4a_cv_score.tsv", "supp4a")
  compare("abm/supp4b_cv_flags_accuracy.tsv", "supp4b")
  compare("abm/supp4c_selected_by_cv.tsv", "supp4c")
  compare("abm/supp5b_beat_baseline_selected.tsv", "supp5b")
  compare("abm/fig1g_beat_baseline.tsv", "fig1g")
}
cat("written to", out, "\n")
