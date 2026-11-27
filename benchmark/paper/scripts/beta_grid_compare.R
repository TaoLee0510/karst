# KARST against alfak on the Salehi tests of Fig. 2c and 2d across a grid of
# missegregation rates (docs/beta_grid_salehi_6ae2a24.md). Both methods fit
# and forecast at the same rate; x01 is the round's own fits (linked).
#
#   Rscript scripts/beta_grid_compare.R --root=DIR [--multiples=1,2,4,8,16,32,64] [--dispersion=FILE] [--boot=2000]
#
# DIR holds one folder per rate, x<mm>/, each a results folder after
# scripts/score_salehi.R, collect_salehi.R, analyse_salehi.R and
# fig2d_noise.R (hpc/bg_job.sh): analysis/salehi/{fig2c_angle_test.tsv,
# fig2d_wasserstein.tsv, per_fit.Rds, sel_eval.Rds}, analysis/fig2d_noise/
# {fig2d_noise.tsv, fig2d_noise_differences.tsv}, tables/salehi_runs.Rds.
# --dispersion: the dispersion table of scripts/fig2d_truth.R
# --stage=dispersion (Salehi rows: line, once).
# Writes DIR/analysis/: fig2c_paired.tsv (no selection: the angle of the same
# trajectory, minobs and evaluation passage, alfak minus KARST, median and
# bootstrap interval over Salehi lines), fig2c_rule.tsv (the paper's rule:
# each fit's median angle, the difference and its interval), fig2d.tsv (the
# paper's rule: Spearman per fit, the difference and its interval; and
# fig2d_noise's raw / noise-held-fixed rho with its intervals), failures.tsv,
# line_best.tsv (per Salehi line, the rate with KARST's lowest median angle,
# against the line's dispersion), summary.tsv, figures/*.pdf.
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
root <- arg("root")
multiples <- as.numeric(strsplit(arg("multiples", "1,2,4,8,16,32,64"), ",")[[1]])
boot <- as.integer(arg("boot", "2000"))
out <- file.path(root, "analysis")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(ggplot2))
alfak <- "alfak dev2@f93b2d2"
designs <- c("2 time points", "all time points")
line_of <- function(fi) sub("_.*$", "", fi)
folder <- function(m) file.path(root, sprintf("x%02d", m))
have <- multiples[vapply(multiples, function(m) file.exists(file.path(folder(m), "analysis", "salehi", "per_fit.Rds")), logical(1))]
if (!length(have)) stop("no rate folder has its analysis yet")
cat("rates with results:", paste0("x", have, collapse = " "), "\n")

# bootstrap over Salehi lines: `stat(rows)` on lines drawn with replacement
line_boot <- function(d, stat) {
  lines <- unique(d$line)
  set.seed(1)
  v <- replicate(boot, {
    pick <- sample(lines, length(lines), replace = TRUE)
    stat(do.call(rbind, lapply(pick, function(l) d[d$line == l, , drop = FALSE])))
  })
  c(est = stat(d), lo = unname(stats::quantile(v, 0.025, na.rm = TRUE)), hi = unname(stats::quantile(v, 0.975, na.rm = TRUE)))
}

paired <- list(); rule <- list(); f2d <- list(); fails <- list(); per_line <- list(); noise <- list()
for (m in have) {
  a <- file.path(folder(m), "analysis", "salehi")
  pf <- readRDS(file.path(a, "per_fit.Rds"))
  se <- readRDS(file.path(a, "sel_eval.Rds"))
  pf <- pf[pf$metric == "euclidean", ]
  for (g in designs) {
    k <- paste("KARST", g); r <- paste(alfak, g)
    keys <- c("fi", "min_obs", "metric", "gap", "type", "coverage", "cov_bin", "gap_group", "sample")
    x <- merge(pf[pf$fit == k, c(keys, "angle")], pf[pf$fit == r, c(keys, "angle")], by = keys, suffixes = c("_karst", "_alfak"))
    x <- x[is.finite(x$angle_karst) & is.finite(x$angle_alfak), ]
    x$line <- line_of(x$fi); x$diff <- x$angle_alfak - x$angle_karst
    b <- line_boot(x, function(d) stats::median(d$diff))
    paired[[length(paired) + 1L]] <- data.frame(multiple = m, design = g, pairs = nrow(x), lines = length(unique(x$line)),
      angle_karst = stats::median(x$angle_karst), angle_alfak = stats::median(x$angle_alfak),
      diff = b[["est"]], lo = b[["lo"]], hi = b[["hi"]], karst_smaller = mean(x$diff > 0), significant = b[["lo"]] > 0)
    # per line: KARST's and alfak's median angle (no selection)
    pl <- stats::aggregate(cbind(angle_karst, angle_alfak) ~ line, x, stats::median)
    per_line[[length(per_line) + 1L]] <- data.frame(multiple = m, design = g, pl, n = as.vector(table(x$line)[pl$line]))
    # the paper's rule (each fit on its own selected rows)
    # Fig. 2c's statistic: the median over lineage types of each type's median angle
    s <- se[se$fit %in% c(k, r) & se$metric == "euclidean" & is.finite(se$angle), c("fit", "fi", "type", "angle")]
    s$line <- line_of(s$fi)
    mom <- function(d) { d <- d[!is.na(d$type), ]; if (!nrow(d)) NA_real_ else stats::median(tapply(d$angle, d$type, stats::median)) }
    if (nrow(s) && all(c(k, r) %in% s$fit)) {
      b <- line_boot(s, function(d) mom(d[d$fit == r, ]) - mom(d[d$fit == k, ]))
      rule[[length(rule) + 1L]] <- data.frame(multiple = m, design = g, n_karst = sum(s$fit == k), n_alfak = sum(s$fit == r),
        angle_karst = mom(s[s$fit == k, ]), angle_alfak = mom(s[s$fit == r, ]),
        diff = b[["est"]], lo = b[["lo"]], hi = b[["hi"]], significant = b[["lo"]] > 0)
    }
    w <- se[se$fit %in% c(k, r) & se$metric == "wasserstein" & is.finite(se$base) & is.finite(se$pred0), c("fit", "fi", "base", "pred0")]
    w$line <- line_of(w$fi)
    if (nrow(w) && all(c(k, r) %in% w$fit)) {
      sp <- function(d, f) { z <- d[d$fit == f, ]; if (nrow(z) < 3) NA_real_ else stats::cor(z$base, z$pred0, method = "spearman") }
      b <- line_boot(w, function(d) sp(d, k) - sp(d, r))
      f2d[[length(f2d) + 1L]] <- data.frame(multiple = m, design = g, selection = "paper rule", n_karst = sum(w$fit == k),
        n_alfak = sum(w$fit == r), rho_karst = sp(w, k), rho_alfak = sp(w, r), diff = b[["est"]], lo = b[["lo"]], hi = b[["hi"]],
        karst_ahead = b[["lo"]] > 0)
    }
  }
  nf <- file.path(folder(m), "analysis", "fig2d_noise")
  if (file.exists(file.path(nf, "fig2d_noise_differences.tsv"))) {
    nd <- utils::read.delim(file.path(nf, "fig2d_noise_differences.tsv"), stringsAsFactors = FALSE)
    nr <- utils::read.delim(file.path(nf, "fig2d_noise.tsv"), stringsAsFactors = FALSE)
    noise[[length(noise) + 1L]] <- list(m = m, diff = nd, rho = nr)
  }
  runs <- readRDS(file.path(folder(m), "tables", "salehi_runs.Rds"))
  runs <- runs[runs$method %in% c("karst", "alfak"), ]
  fails[[length(fails) + 1L]] <- data.frame(multiple = m, stats::aggregate(list(failed = runs$status != "ok"),
                                                                            list(method = runs$method, design = runs$design), mean))
}
paired <- do.call(rbind, paired); rule <- do.call(rbind, rule); f2d <- do.call(rbind, f2d); per_line <- do.call(rbind, per_line)
fails <- do.call(rbind, lapply(fails, function(x) x[, c("multiple", "method", "design", "failed")]))
wt <- function(d, name) utils::write.table(d, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
wt(paired, "fig2c_paired"); wt(rule, "fig2c_rule"); wt(fails, "failures"); wt(per_line, "per_line_angles")
nz <- do.call(rbind, lapply(noise, function(z) {
  d <- z$diff[z$diff$arm == "karst", ]
  r <- z$rho
  cbind(multiple = z$m, d)
}))
if (!is.null(nz)) wt(nz, "fig2d_noise_differences")
nzr <- do.call(rbind, lapply(noise, function(z) cbind(multiple = z$m, z$rho)))
if (!is.null(nzr)) wt(nzr, "fig2d_noise_rho")
wt(f2d, "fig2d_rule")
options(width = 220)
cat("\n== Fig. 2c, no selection: paired angle (alfak - KARST, degrees), line bootstrap\n"); print(paired, row.names = FALSE, digits = 3)
cat("\n== Fig. 2c, paper rule (median over lineage types of the type medians; alfak - KARST)\n"); print(rule, row.names = FALSE, digits = 3)
cat("\n== Fig. 2d, paper rule: Spearman (KARST - alfak), line bootstrap\n"); print(f2d, row.names = FALSE, digits = 3)
if (!is.null(nz)) { cat("\n== Fig. 2d (fig2d_noise): KARST - alfak, raw and noise held fixed\n"); print(nz, row.names = FALSE, digits = 3) }
cat("\n== failures (share of fits)\n"); print(reshape(fails, idvar = c("method", "design"), timevar = "multiple", direction = "wide"), row.names = FALSE, digits = 3)

# best rates: same rate for both (curve), each method's own best, per line
best <- do.call(rbind, lapply(designs, function(g) {
  p <- paired[paired$design == g, ]
  data.frame(design = g, default_diff = p$diff[p$multiple == 1], default_sig = p$significant[p$multiple == 1],
             best_same_rate = p$multiple[which.max(p$diff)], best_same_diff = max(p$diff),
             karst_best_rate = p$multiple[which.min(p$angle_karst)], karst_best_angle = min(p$angle_karst),
             alfak_best_rate = p$multiple[which.min(p$angle_alfak)], alfak_best_angle = min(p$angle_alfak),
             karst_default_angle = p$angle_karst[p$multiple == 1], alfak_default_angle = p$angle_alfak[p$multiple == 1])
}))
wt(best, "summary")
cat("\n== best rates (Fig. 2c median angle, no selection)\n"); print(best, row.names = FALSE, digits = 3)

lb <- do.call(rbind, lapply(split(per_line[per_line$design == "2 time points", ], per_line$line[per_line$design == "2 time points"]),
                            function(d) data.frame(line = d$line[1], n = d$n[d$multiple == min(d$multiple)][1],
                                                   karst_best_rate = d$multiple[which.min(d$angle_karst)],
                                                   karst_default = d$angle_karst[d$multiple == 1][1], karst_best = min(d$angle_karst),
                                                   alfak_best_rate = d$multiple[which.min(d$angle_alfak)],
                                                   alfak_default = d$angle_alfak[d$multiple == 1][1], alfak_best = min(d$angle_alfak))))
disp_file <- arg("dispersion")
if (!is.null(disp_file) && file.exists(disp_file)) {
  ds <- utils::read.delim(disp_file, stringsAsFactors = FALSE)
  ds <- ds[ds$set == "Salehi", ]
  ds$line <- line_of(ds$line)
  dl <- stats::aggregate(cbind(once, simpson) ~ line, ds, stats::median)
  lb <- merge(lb, dl, by = "line", all.x = TRUE)
}
if ("once" %in% names(lb)) lb$karst_gain <- lb$karst_default - lb$karst_best
wt(lb, "line_best")
cat("\n== per Salehi line (2 time points): KARST's best rate against dispersion\n"); print(lb[order(lb$once), ], row.names = FALSE, digits = 3)
if (all(c("once") %in% names(lb))) {
  ct <- suppressWarnings(stats::cor.test(log2(lb$karst_best_rate), lb$once, method = "spearman"))
  cat("Spearman(log2 best rate, once):", round(ct$estimate, 3), " p =", signif(ct$p.value, 2), " lines:", sum(is.finite(lb$once)), "\n")
}

if ("once" %in% names(lb)) {
  ok <- is.finite(lb$once)
  cg <- suppressWarnings(stats::cor.test(lb$karst_gain[ok], lb$once[ok], method = "spearman"))
  cat("Spearman(KARST's gain from its best rate, once):", round(cg$estimate, 3), " p =", signif(cg$p.value, 2), "\n")
}

## each method at its own best rate (in-sample) ---------------------------------
own <- list()
keys <- c("fi", "min_obs", "metric", "gap", "type", "coverage", "cov_bin", "gap_group", "sample")
pf_of <- function(m, f) { d <- readRDS(file.path(folder(m), "analysis", "salehi", "per_fit.Rds")); d[d$fit == f & d$metric == "euclidean", c(keys, "angle")] }
se_of <- function(m, f, metric) { d <- readRDS(file.path(folder(m), "analysis", "salehi", "sel_eval.Rds")); d[d$fit == f & d$metric == metric, ] }
for (g in designs) {
  k <- paste("KARST", g); r <- paste(alfak, g)
  p <- paired[paired$design == g, ]
  for (cmp in list(c("own best", p$multiple[which.min(p$angle_karst)], p$multiple[which.min(p$angle_alfak)]),
                   c("KARST best, alfak default", p$multiple[which.min(p$angle_karst)], 1))) {
    mk <- as.numeric(cmp[2]); ma <- as.numeric(cmp[3])
    x <- merge(pf_of(mk, k), pf_of(ma, r), by = keys, suffixes = c("_karst", "_alfak"))
    x <- x[is.finite(x$angle_karst) & is.finite(x$angle_alfak), ]
    x$line <- line_of(x$fi); x$diff <- x$angle_alfak - x$angle_karst
    b <- line_boot(x, function(d) stats::median(d$diff))
    own[[length(own) + 1L]] <- data.frame(test = "Fig. 2c, no selection", comparison = cmp[1], design = g, karst_rate = mk, alfak_rate = ma,
      n = nrow(x), karst = stats::median(x$angle_karst), alfak = stats::median(x$angle_alfak), diff = b[["est"]], lo = b[["lo"]], hi = b[["hi"]])
  }
  q <- f2d[f2d$design == g, ]
  for (cmp in list(c("own best", q$multiple[which.max(q$rho_karst)], q$multiple[which.max(q$rho_alfak)]),
                   c("KARST best, alfak default", q$multiple[which.max(q$rho_karst)], 1))) {
    mk <- as.numeric(cmp[2]); ma <- as.numeric(cmp[3])
    w <- rbind(se_of(mk, k, "wasserstein")[, c("fit", "fi", "base", "pred0")], se_of(ma, r, "wasserstein")[, c("fit", "fi", "base", "pred0")])
    w <- w[is.finite(w$base) & is.finite(w$pred0), ]; w$line <- line_of(w$fi)
    sp <- function(d, f) { z <- d[d$fit == f, ]; if (nrow(z) < 3) NA_real_ else stats::cor(z$base, z$pred0, method = "spearman") }
    b <- line_boot(w, function(d) sp(d, k) - sp(d, r))
    own[[length(own) + 1L]] <- data.frame(test = "Fig. 2d, paper rule", comparison = cmp[1], design = g, karst_rate = mk, alfak_rate = ma,
      n = nrow(w), karst = sp(w, k), alfak = sp(w, r), diff = b[["est"]], lo = b[["lo"]], hi = b[["hi"]])
  }
}
own <- do.call(rbind, own)
wt(own, "own_best")
cat("\n== each method at its own best rate (in-sample, optimistic); Fig. 2c diff = alfak - KARST, Fig. 2d diff = KARST - alfak\n")
print(own, row.names = FALSE, digits = 3)

## figures
theme_set(theme_classic(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top"))
cols <- c(KARST = "#D55E00", alfak = "#0072B2")
pa <- rbind(data.frame(paired[, c("multiple", "design")], method = "KARST", angle = paired$angle_karst),
            data.frame(paired[, c("multiple", "design")], method = "alfak", angle = paired$angle_alfak))
g1 <- ggplot(pa, aes(multiple, angle, colour = method)) + geom_line() + geom_point(size = 1.4) + facet_wrap(~design) +
  scale_x_log10(breaks = multiples) + scale_colour_manual(NULL, values = cols) +
  labs(x = "missegregation rate (multiple of 5e-5), both methods", y = "median angle (degrees, lower is better)",
       caption = "Fig. 2c, every trajectory and minobs both methods forecast (no selection).")
g2 <- ggplot(paired, aes(multiple, diff)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) + facet_wrap(~design) + scale_x_log10(breaks = multiples) +
  labs(x = "missegregation rate (multiple of 5e-5)", y = "alfak - KARST angle (degrees)",
       caption = "Median paired difference; 95% bootstrap over Salehi lines. Above 0: KARST closer to the observed change.")
ggsave(file.path(out, "figures", "fig2c_grid.pdf"), gridExtra::arrangeGrob(g1, g2, ncol = 1), width = 160, height = 140, units = "mm")
r2 <- rbind(data.frame(rule[, c("multiple", "design")], method = "KARST", angle = rule$angle_karst),
            data.frame(rule[, c("multiple", "design")], method = "alfak", angle = rule$angle_alfak))
g5 <- ggplot(r2, aes(multiple, angle, colour = method)) + geom_line() + geom_point(size = 1.4) + facet_wrap(~design) +
  scale_x_log10(breaks = multiples) + scale_colour_manual(NULL, values = cols) +
  labs(x = "missegregation rate (multiple of 5e-5), both methods", y = "median angle (degrees)",
       caption = "Fig. 2c, the paper's rule (each fit's minobs and trajectory length with the best CV score).")
ggsave(file.path(out, "figures", "fig2c_rule_grid.pdf"), g5, width = 160, height = 70, units = "mm")
if (!is.null(nzr)) {
  r <- nzr[nzr$arm %in% c("karst", "alfak") & nzr$measure %in% c("raw", "noise held fixed"), ]
  r$method <- ifelse(r$arm == "karst", "KARST", "alfak")
  g3 <- ggplot(r, aes(multiple, rho, colour = method)) + geom_line() + geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.05, linewidth = 0.3) +
    facet_grid(measure ~ paste(group, subset, sep = ", ")) + scale_x_log10(breaks = multiples) +
    scale_colour_manual(NULL, values = cols) +
    labs(x = "missegregation rate (multiple of 5e-5), both methods", y = "Spearman, observed vs predicted change",
         caption = "Fig. 2d (scripts/fig2d_noise.R): rows every arm has; intervals by lineage-type bootstrap.")
  ggsave(file.path(out, "figures", "fig2d_grid.pdf"), g3, width = 200, height = 110, units = "mm")
}
if ("once" %in% names(lb)) {
  g4 <- ggplot(lb, aes(once, karst_best_rate, label = line)) + geom_point() + geom_text(size = 2, vjust = -0.6) +
    scale_y_log10(breaks = multiples) + labs(x = "dispersion of the line (share of cells whose karyotype occurs once, 300 cells)",
                                             y = "KARST's best rate (lowest median angle)", caption = "2 time points, no selection.")
  ggsave(file.path(out, "figures", "line_best_rate.pdf"), g4, width = 120, height = 90, units = "mm")
}
cat("\nwritten:", out, "\n")
