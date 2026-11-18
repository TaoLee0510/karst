# The paper's Fig. 4-6 conclusions for every arm (after fig4_cna_effects.R,
# fig5_wgd.R and fig6_steady_states.R): do the fits of each method lead to
# the conclusions the paper drew from its own fits?
#
#   Rscript scripts/fig456_summary.R
#
# Reference: Fig. 4, the paper's stored bootstrap ("paper, published",
# data/processed/salehi/landscape_analysis_bootstrap.Rds, 200 draws); Fig. 5
# (none stored) and Fig. 6, the paper's own fits on every input through the
# same code ("paper (alfak main@4ad82b0), every input"). For Fig. 6 the stored results
# (steadyStatePredictions.Rds) are not the reference: the same landscapes give
# other steady states with alfakR_dev2's find_steady_state() than with the
# version the paper used (fig6_reproduction.tsv), so every arm is compared
# with the paper's fits through the same function, and the stored results are
# listed as "paper, published". The arms use the analysable
# trajectories only (the paper line too, "paper (alfak main@4ad82b0)"), except the
# reproduction rows (*_all_inputs).
#
# Writes BENCH_RESULTS/analysis/fig456/:
#   fig4_estimates.tsv   per arm and parameter: draws, 2.5/50/97.5% quantiles
#                        (4d: PDX_v, cisplatin_v, PDX_kv, cisplatin_kv, and
#                        the variance models with na.rm, PDX_kv_narm,
#                        cisplatin_kv_narm; 4f: the similarity model's
#                        coefficients)
#   fig4_kvar_dropped.tsv  per arm: frequent karyotypes per draw and the
#                        share the variance models of 4d lose (a neighbour
#                        missing from the landscape; 06 takes the variance
#                        without na.rm)
#   fig4_deltaf.tsv      example draw (seed 42): delta f by context and
#                        treatment (4c)
#   fig5_tests.tsv       (from fig5_wgd.R) and fig5_model.tsv (5b terms)
#   fig6_dominance.tsv   per arm: landscapes, landscapes whose dominant
#                        karyotype changes with the rate, median distance of
#                        the dominant karyotype from the lowest rate's
#   fig6_agreement.tsv   per arm and landscape: same dominant karyotype as the
#                        paper's at the lowest and the highest rate
#   fig6_reproduction.tsv  the paper line rerun against the stored results
#   conclusions.tsv      per statement and arm: the value, its significance
#                        and whether it agrees with the reference: the same
#                        significance and, where both are significant, the
#                        same sign; intercepts are not statements, Fig. 6 is
#                        descriptive (fig456_conclusions(), R/paper_fig456.R)
# and figures/fig4_effects.png, fig5_wgd.png, fig6_dominance.png.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
dir <- file.path(paths$results, "analysis", "fig456")
fig <- file.path(paths$results, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
write_tab <- function(d, name) utils::write.table(d, file.path(dir, paste0(name, ".tsv")), sep = "\t", quote = FALSE,
                                                  row.names = FALSE)
arm_levels <- c(fig456_arms$key, "paper_all_inputs", "paper_published")
arm_label <- c(stats::setNames(fig456_arms$fit, fig456_arms$key), paper_all_inputs = "paper (alfak main@4ad82b0), every input",
               paper_published = "paper, published")
q3 <- function(v) stats::quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)

## Fig. 4 ---------------------------------------------------------------------------
read_fig4 <- function(sub, key_suffix = "") {
  files <- list.files(file.path(dir, sub), pattern = "\\.Rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  rs <- lapply(files, readRDS)
  keys <- paste0(vapply(rs, function(r) r$arm$key, character(1)), key_suffix)
  list(df1 = do.call(rbind, Map(function(r, k) cbind(arm = k, r$df1), rs, keys)),
       df2 = do.call(rbind, Map(function(r, k) if (is.null(r$df2)) NULL else cbind(arm = k, r$df2[c("param", "Estimate")]),
                                rs, keys)),
       example = do.call(rbind, Map(function(r, k) if (is.null(r$example)) NULL else
         cbind(arm = k, r$example[c("deltaf", "context", "treat")]), rs, keys)))
}
f4 <- read_fig4("fig4")
f4_all <- read_fig4("fig4_all_inputs", "_all_inputs")
f4_all <- if (is.null(f4_all)) NULL else lapply(f4_all, function(d) if (is.null(d)) NULL else d[d$arm == "paper_all_inputs", ])
pub <- readRDS(file.path(paths$data, "processed", "salehi", "landscape_analysis_bootstrap.Rds"))
pub1 <- cbind(arm = "paper_published", do.call(rbind, lapply(pub, `[[`, "df1")))
pub2 <- cbind(arm = "paper_published", do.call(rbind, lapply(pub, function(x) data.frame(param = rownames(x$df2),
                                                                                          Estimate = x$df2$Estimate))))
df1 <- bind_rows_fill(list(f4$df1, f4_all$df1, pub1))
df2 <- bind_rows_fill(list(f4$df2, f4_all$df2, pub2))
est <- NULL
if (!is.null(df1)) {
  # the na.rm sensitivity of the variance models (fig4_cna_effects.R; rounds
  # from 2026-11-18 on)
  params1 <- intersect(c("PDX_v", "cisplatin_v", "PDX_kv", "cisplatin_kv", "PDX_kv_narm", "cisplatin_kv_narm"),
                       names(df1))
  long1 <- do.call(rbind, lapply(params1, function(v)
    data.frame(arm = df1$arm, param = v, value = df1[[v]])))
  long2 <- data.frame(arm = df2$arm, param = df2$param, value = df2$Estimate)
  long <- rbind(long1, long2)
  est <- do.call(rbind, lapply(split(long, long[c("arm", "param")], drop = TRUE), function(x) {
    q <- q3(x$value)
    data.frame(arm = x$arm[1], param = x$param[1], draws = sum(is.finite(x$value)), lo = q[1], mid = q[2], hi = q[3])
  }))
  est <- est[order(match(est$arm, arm_levels), est$param), ]
  write_tab(est, "fig4_estimates")
}
# Frequent karyotypes the variance models of 4d lose: 06 takes the variance
# of delta f without na.rm, so one neighbour missing from the landscape
# drops the karyotype (fig4_cna_effects.R; the published bootstrap does not
# record it).
kvar_tab <- NULL
if (!is.null(df1$kvar_frequent)) {
  kv <- df1[!is.na(df1$kvar_frequent) & df1$kvar_frequent > 0, ]
  kvar_tab <- do.call(rbind, lapply(split(kv, kv$arm), function(x) {
    share <- x$kvar_na / x$kvar_frequent
    q <- q3(share)
    data.frame(arm = x$arm[1], draws = nrow(x), frequent_median = stats::median(x$kvar_frequent),
               dropped_median = stats::median(x$kvar_na), share_dropped_lo = round(q[1], 3),
               share_dropped_median = round(q[2], 3), share_dropped_hi = round(q[3], 3))
  }))
  kvar_tab <- kvar_tab[order(match(kvar_tab$arm, arm_levels)), ]
  write_tab(kvar_tab, "fig4_kvar_dropped")
}
if (!is.null(f4$example)) {
  ex <- f4$example
  dtab <- do.call(rbind, lapply(split(ex, ex[c("arm", "context", "treat")], drop = TRUE), function(x) data.frame(
    arm = x$arm[1], context = x$context[1], treat = x$treat[1], n = nrow(x),
    median_deltaf = stats::median(x$deltaf, na.rm = TRUE), median_abs_deltaf = stats::median(abs(x$deltaf), na.rm = TRUE))))
  write_tab(dtab, "fig4_deltaf")
}

## Fig. 5 ---------------------------------------------------------------------------
f5 <- if (file.exists(file.path(dir, "fig5.Rds"))) readRDS(file.path(dir, "fig5.Rds")) else NULL
f5_all <- if (file.exists(file.path(dir, "fig5_all_inputs.Rds"))) readRDS(file.path(dir, "fig5_all_inputs.Rds")) else NULL
f5_tests <- bind_rows_fill(c(lapply(f5, `[[`, "tests"),
                             list(if (!is.null(f5_all$paper)) transform(f5_all$paper$tests, arm = "paper_all_inputs"))))
f5_model <- bind_rows_fill(c(lapply(f5, `[[`, "coefs"),
                             list(if (!is.null(f5_all$paper$coefs)) transform(f5_all$paper$coefs, arm = "paper_all_inputs"))))
if (!is.null(f5_model)) {
  names(f5_model)[names(f5_model) == "Pr(>|t|)"] <- "p"
  write_tab(f5_model, "fig5_model")
}

## Fig. 6 ---------------------------------------------------------------------------
dominants <- function(dfi) rownames(dfi)[apply(dfi, 2, which.max)]
read_fig6 <- function(sub) {
  files <- list.files(file.path(dir, sub), pattern = "\\.Rds$", full.names = TRUE)
  stats::setNames(lapply(files, readRDS), sub("\\.Rds$", "", basename(files)))
}
f6 <- read_fig6("fig6")
f6_all <- read_fig6("fig6_all_inputs")
pub6 <- readRDS(file.path(paths$data, "processed", "salehi", "steadyStatePredictions.Rds"))
f6_sets <- c(f6, if (!is.null(f6_all$paper)) list(paper_all_inputs = f6_all$paper), list(paper_published = pub6))
dom_tab <- do.call(rbind, lapply(names(f6_sets), function(a) {
  s <- f6_sets[[a]]
  if (!length(s)) return(NULL)
  ndom <- vapply(s, function(d) length(unique(apply(d, 2, which.max))), integer(1))
  ddom <- unlist(lapply(s, function(d) {
    k <- do.call(rbind, lapply(dominants(d), paper_s2v))
    sum(abs(k[nrow(k), ] - k[1, ]))
  }))
  data.frame(arm = a, landscapes = length(s), switching = sum(ndom > 1), share_switching = round(mean(ndom > 1), 3),
             median_distance_highest_rate = stats::median(ddom), max_distance = max(ddom))
}))
write_tab(dom_tab, "fig6_dominance")
# Same dominant karyotype as the reference's (the paper's fits through the
# same function; the stored results without that rerun) for the same landscape.
clean <- function(fi) sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", fi)
ref6 <- if (length(f6_all$paper)) f6_all$paper else pub6
pub_dom <- lapply(ref6, dominants)
names(pub_dom) <- clean(names(ref6))
agree <- do.call(rbind, lapply(names(f6_sets), function(a) {
  s <- f6_sets[[a]]
  if (!length(s) || a %in% c("paper_published", if (length(f6_all$paper)) "paper_all_inputs")) return(NULL)
  do.call(rbind, lapply(names(s), function(fi) {
    ref <- pub_dom[[clean(fi)]]
    if (is.null(ref)) return(NULL)
    d <- dominants(s[[fi]])
    data.frame(arm = a, sample = clean(fi), fi = fi, same_lowest_rate = d[1] == ref[1],
               same_highest_rate = d[length(d)] == ref[length(ref)], same_all_rates = all(d == ref))
  }))
}))
if (!is.null(agree)) write_tab(agree, "fig6_agreement")
if (length(f6_all$paper)) {
  rep6 <- do.call(rbind, lapply(names(pub6), function(fi) {
    a <- f6_all$paper[[fi]]
    b <- pub6[[fi]]
    if (is.null(a)) return(data.frame(fi = fi, rerun = FALSE, same_fit = NA, same_karyotypes = NA, max_diff = NA_real_,
                                      same_dominants = NA))
    k <- intersect(rownames(a), rownames(b))
    data.frame(fi = fi, rerun = TRUE, same_fit = identical(attr(a, "min_obs"), attr(b, "min_obs")),
               same_karyotypes = setequal(rownames(a), rownames(b)), max_diff = max(abs(a[k, ] - b[k, ])),
               same_dominants = identical(dominants(a), dominants(b)))
  }))
  write_tab(rep6, "fig6_reproduction")
}

## Conclusions ------------------------------------------------------------------------
# fig456_conclusions() (R/paper_fig456.R): Fig. 4 against the published
# bootstrap, Fig. 5 against the paper's fits on every input; agreement is the
# same significance and, where both are significant, the same sign;
# intercepts left out; Fig. 6 descriptive.
conclusions <- fig456_conclusions(est, f5_tests, f5_model, dom_tab)
if (!is.null(conclusions)) {
  conclusions$arm_label <- arm_label[conclusions$arm]
  write_tab(conclusions, "conclusions")
}

## Figures ----------------------------------------------------------------------------
fit_col <- c(fit_colours[fig456_arms$fit], "paper (alfak main@4ad82b0), every input" = "#666666", "paper, published" = "#000000")
if (!is.null(est)) {
  e <- est
  e$arm_label <- factor(arm_label[e$arm], rev(arm_label[arm_levels]))
  p <- ggplot(e, aes(mid, arm_label, colour = arm_label)) +
    facet_wrap(~param, scales = "free_x", ncol = 4) +
    geom_vline(xintercept = 0, linewidth = 0.2, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) + geom_point(size = 1) +
    scale_colour_manual(NULL, values = fit_col, guide = "none") +
    labs(x = "bootstrap median and 95% interval", y = NULL,
         subtitle = "Fig. 4d (PDX/cisplatin effects on |delta f| and its variance) and 4f (similarity model)") +
    theme_classic(base_size = 7) + theme(strip.background = element_blank())
  ggsave(file.path(fig, "fig4_effects.png"), p, width = 200, height = 120, units = "mm", dpi = 200, bg = "white")
}
if (length(f5)) {
  # every method in the same panel; panels by design (the paper's fits, all
  # time points, are drawn in both as the reference) and by WGD status
  dfl <- bind_rows_fill(lapply(f5, `[[`, "dfl"))
  design_of <- c(stats::setNames(ifelse(fig456_arms$all, "all time points", "2 time points"), fig456_arms$key))
  dfl$design <- design_of[dfl$arm]
  dfl <- rbind(dfl[dfl$arm != "paper", ], transform(dfl[dfl$arm == "paper", ], design = "2 time points"),
               dfl[dfl$arm == "paper", ])
  dfl$arm_label <- factor(arm_label[dfl$arm], arm_label[arm_levels])
  p <- ggplot(dfl, aes(df, colour = arm_label)) +
    facet_grid(wgd_str ~ design) + stat_ecdf(linewidth = 0.4) +
    scale_x_continuous("fitness effect (delta f)", limits = c(-0.05, 0.05)) +
    scale_colour_manual(NULL, values = fit_col) +
    labs(y = "cumulative fraction", subtitle = "Fig. 5c: SA906 fits (KARST one minobs per trajectory)") +
    theme_classic(base_size = 7) + theme(strip.background = element_blank(), legend.position = "top")
  ggsave(file.path(fig, "fig5_wgd.png"), p, width = 200, height = 120, units = "mm", dpi = 200, bg = "white")
}
if (!is.null(dom_tab)) {
  d <- dom_tab
  d$arm_label <- factor(arm_label[d$arm], rev(arm_label[arm_levels]))
  p <- ggplot(d, aes(share_switching, arm_label, colour = arm_label)) + geom_point(size = 2) +
    geom_text(aes(label = paste(switching, "/", landscapes)), hjust = -0.3, size = 2.2, colour = "black") +
    scale_colour_manual(NULL, values = fit_col, guide = "none") +
    scale_x_continuous("share of landscapes whose dominant karyotype changes with the rate", limits = c(0, 1)) +
    labs(y = NULL, subtitle = "Fig. 6: steady states at missegregation rates 4e-5 to 5.1e-3") +
    theme_classic(base_size = 7)
  ggsave(file.path(fig, "fig6_dominance.png"), p, width = 150, height = 60, units = "mm", dpi = 200, bg = "white")
}
options(width = 220)
if (!is.null(est)) { cat("== Fig. 4 estimates\n"); print(est, row.names = FALSE) }
if (!is.null(kvar_tab)) {
  cat("\n== Fig. 4d, frequent karyotypes the variance models lose (no na.rm, as 06)\n")
  print(kvar_tab, row.names = FALSE)
}
if (!is.null(f5_tests)) { cat("\n== Fig. 5 tests\n"); print(f5_tests, row.names = FALSE) }
cat("\n== Fig. 6 dominance\n"); print(dom_tab, row.names = FALSE)
if (exists("rep6")) {
  cat("\n== Fig. 6, the paper's fits rerun against the stored results: same fit", sum(rep6$same_fit, na.rm = TRUE),
      "of", nrow(rep6), " same karyotypes", sum(rep6$same_karyotypes, na.rm = TRUE), " same dominant karyotypes",
      sum(rep6$same_dominants, na.rm = TRUE), " largest difference", signif(max(rep6$max_diff, na.rm = TRUE), 2), "\n")
}
if (!is.null(conclusions)) {
  cat("\n== Agreement with the reference, per arm (tested statements; Fig. 6 descriptive)\n")
  print(stats::aggregate(agrees ~ arm, conclusions, function(v) sprintf("%d/%d", sum(v, na.rm = TRUE), sum(!is.na(v)))))
}
