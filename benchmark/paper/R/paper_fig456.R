# The analyses of the paper's Fig. 4-6 on the Salehi fits of every method
# (scripts/fig4_cna_effects.R, fig5_wgd.R, fig6_steady_states.R,
# fig456_summary.R), ported from the paper's scripts 06_cna_fitness_effects.R,
# 07_wgd.R and 08_misseg_rate_analysis.R. They read each fit's landscape.Rds
# (alfak's format; karst's from as_alfak()) and the fit's row of 03 (xv,
# ntrain, train_treat, test_treat), and select fits as the paper does.

# Karyotype helpers of the paper (R/utils_karyo.R), unchanged.
paper_s2v <- function(charvec) as.numeric(unlist(strsplit(charvec, split = "[.]")))
paper_gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) ids <- lapply(ids, paper_s2v)
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n)
  n <- unique(n)
  n <- n[-(1:nids), ]
  if (remove_nullisomes) n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), ]
  n
}

# The fits compared: each method's fits of the first and last passage ("2
# time points") and of every passage ("all time points"; a trajectory of two passages
# has only the former, which stands for both), and the paper's own fits
# (alfak main@4ad82b0). The fig4 array has 10 tasks per arm (tasks 1-70), the
# fig6 array one task per arm (tasks 1-7); the paper's arm is the seventh.
fig456_arms <- data.frame(
  key = c("karst_2", "karst_all", "karst_pairwise_all", "karst_state_space_all", "alfak_2", "alfak_all", "paper"),
  fit = c("KARST 2 time points", "KARST all time points", "KARST pairwise all time points",
          "KARST state space all time points", "alfak dev2@f93b2d2 2 time points", "alfak dev2@f93b2d2 all time points", "paper (alfak main@4ad82b0)"),
  method = c("karst", "karst", "karst_pairwise", "karst_state_space", "alfak", "alfak", "alfak_paper"),
  all = c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE)

# Folder of one fit of an arm.
fig456_fit_dir <- function(arm, fi, min_obs, paths = bench_paths()) {
  if (arm$method == "alfak_paper") return(file.path(paths$paper_salehi, fi, paste0("minobs_", min_obs)))
  design <- if (arm$all && !grepl("_l_2_", fi)) "all" else "ends"
  method <- if (design == "ends" && arm$method %in% karst_other_methods) "karst" else arm$method
  file.path(paths$results, "fits", "salehi", fi, design, method, paste0("minobs_", min_obs))
}

# The rows of 03 of every fit of an arm, as a list with one data frame per
# fit (as the paper's alfak_outputs_proc/ files): the paper's own files for
# the paper's fits, this round's scores (tables/salehi_scores.Rds) otherwise.
# Trajectories that are not analysable (R/coverage.R) are left out unless
# `exclude` is FALSE, each fit by the rule of the time points it fits: "all"
# for an arm of every time point (the paper's fits included) on a trajectory
# of more than two passages, "ends" otherwise (as analyse_salehi.R).
fig456_fit_rows <- function(arm, paths = bench_paths(), exclude = TRUE) {
  if (arm$method == "alfak_paper") {
    files <- list.files(file.path(paths$data, "processed", "salehi", "alfak_outputs_proc"), full.names = TRUE)
    rows <- lapply(files, readRDS)
  } else {
    s <- readRDS(file.path(paths$results, "tables", "salehi_scores.Rds"))
    two <- grepl("_l_2_", s$fi)
    want_design <- ifelse(arm$all & !two, "all", "ends")
    want_method <- ifelse(want_design == "ends" & arm$method %in% karst_other_methods, "karst", arm$method)
    s <- s[s$design == want_design & s$method == want_method, setdiff(names(s), c("design", "method"))]
    rows <- split(s, paste(s$fi, s$min_obs), drop = TRUE)
    names(rows) <- NULL
  }
  if (exclude) {
    inputs <- readRDS(file.path(paths$results, "tables", "salehi_inputs.Rds"))
    if (is.null(inputs$design)) inputs$design <- "ends"  # tables of earlier rounds: "ends" only
    bad <- not_analysable_keys(inputs, c("fi", "design"))
    has_all <- any(inputs$design == "all")
    design_of <- function(fi) if (arm$all && has_all && !grepl("_l_2_", fi)) "all" else "ends"
    rows <- Filter(function(r) !input_key(r$fi[1], design_of(r$fi[1])) %in% bad, rows)
  }
  rows
}

## Fig. 5 (07_wgd.R) --------------------------------------------------------------

# The fits of an arm for Fig. 5 (fits_raw: fi, min_obs, ntrain, xv, one row
# per fit with xv > 0). 07 keeps one fit per trajectory and minobs (the
# longest training, then the highest xv), so every minobs of a trajectory.
# alfak fits each minobs on its own; KARST derives every minobs from one fit
# (`shared`), so that its minobs repeat the same karyotypes with the same
# fitness and 07's rule would count a trajectory's karyotypes up to three
# times. For KARST one minobs per trajectory is kept, as in Fig. 4: the
# highest xv, of a tie the lowest minobs (R/select.R).
fig5_select_fits <- function(fits_raw, shared) {
  fits <- fits_raw[with(fits_raw, order(fi, min_obs, -ntrain, -xv)), , drop = FALSE]
  fits <- fits[!duplicated(fits[c("fi", "min_obs")]), , drop = FALSE]
  if (shared) {
    fits <- select_by_cv(fits, "fi", "min_obs", "xv")
    fits <- fits[!duplicated(fits$fi), , drop = FALSE]
  }
  rownames(fits) <- NULL
  fits
}

# Fig. 5b-c of an arm from its tables: `dffit` the fitness of the frequent
# karyotypes and their neighbours (fitness, dk, wgd_str, id = trajectory),
# `dfl` the fitness effects (df) of each frequent karyotype's one-chromosome
# changes (df, wgd_str, id).
#   5c  Kolmogorov-Smirnov statistic of the per-fit mean df, WGD+ against
#       WGD-, and its permutation p value (1e4 permutations, seed 1, as 07)
#   5b  piecewise mixed model fitness ~ dk * WGD + dk:(dk > 3) * WGD + (1|id),
#       p values from lmerTest
# Returns the per-fit means (cl), the KS statistic and p (ks), the medians of
# the per-fit means by WGD status and the model's coefficients (NULL when it
# cannot be fitted).
fig5_stats <- function(dffit, dfl) {
  set.seed(1)
  cl <- stats::aggregate(df ~ id + wgd_str, dfl, mean)
  groups <- split(cl$df, cl$wgd_str)
  ks <- if (length(groups) == 2L) {
    obs <- suppressWarnings(stats::ks.test(groups[[1]], groups[[2]])$statistic)
    perm <- replicate(1e4, {
      lab <- sample(cl$wgd_str)
      g <- split(cl$df, lab)
      suppressWarnings(stats::ks.test(g[[1]], g[[2]])$statistic)
    })
    c(D = unname(obs), p = mean(perm >= obs))
  } else c(D = NA_real_, p = NA_real_)
  medians <- tapply(cl$df, cl$wgd_str, stats::median)
  dffit$cut3 <- dffit$dk > 3
  mod <- tryCatch(lmerTest::lmer(fitness ~ dk * wgd_str + dk:cut3 * wgd_str + (1 | id), data = dffit),
                  error = function(e) e)
  coefs <- if (inherits(mod, "error")) NULL else {
    cf <- as.data.frame(stats::coef(summary(mod)))
    cf$term <- rownames(cf)
    cf
  }
  list(cl = cl, ks = ks, medians = medians, coefs = coefs)
}

## Conclusions of Fig. 4-6 ------------------------------------------------------------

# Does a statement agree with the reference? The same significance and,
# where both are significant, the same sign. Two results that are not
# significant agree whatever the signs of their estimates: a sign that is
# not significant is no conclusion. NA where either significance is unknown.
fig456_agrees <- function(value, significant, ref_value, ref_significant) {
  out <- significant == ref_significant & (!significant | sign(value) == sign(ref_value))
  out[is.na(significant) | is.na(ref_significant)] <- NA
  out
}

# The statements of the paper's Fig. 4-6 for every arm and whether each
# agrees with the reference (fig456_agrees()):
#   4d, 4f  bootstrap median and whether the 95% interval excludes 0
#           (est: arm, param, lo, mid, hi); reference the published
#           bootstrap (arm "paper_published")
#   5c      KS permutation p < 0.05 (a statistic without sign: significance
#           only); 5b the model terms at p < 0.05 (f5_tests, f5_model: arm,
#           term, Estimate, p); reference the paper's fits on every input
#           (arm "paper_all_inputs")
#   6       the share of landscapes whose dominant karyotype changes with the
#           missegregation rate (dom_tab), descriptive (significant and agrees
#           NA): the paper reports the share without a test, and whether any
#           landscape switches is the same for every arm by construction
# Intercepts are not statements: the paper draws no conclusion from them,
# and 5b's depends on each method's absolute fitness level. Parameters the
# reference does not have (the na.rm sensitivity of 4d) are left out.
fig456_conclusions <- function(est = NULL, f5_tests = NULL, f5_model = NULL, dom_tab = NULL) {
  concl <- list()
  if (!is.null(est)) {
    ref <- est[est$arm == "paper_published", ]
    x <- est[est$param != "(Intercept)" & est$param %in% ref$param, ]
    r <- ref[match(x$param, ref$param), ]
    sig <- function(d) ifelse(is.finite(d$lo) & is.finite(d$hi), d$lo > 0 | d$hi < 0, NA)
    concl[[length(concl) + 1L]] <- data.frame(
      figure = ifelse(x$param %in% c("PDX_v", "cisplatin_v", "PDX_kv", "cisplatin_kv"), "4d", "4f"),
      statement = x$param, arm = x$arm, value = round(x$mid, 4), interval = sprintf("[%.3g, %.3g]", x$lo, x$hi),
      significant = sig(x), reference_value = round(r$mid, 4), reference_significant = sig(r),
      agrees = fig456_agrees(x$mid, sig(x), r$mid, sig(r)), stringsAsFactors = FALSE)
  }
  if (!is.null(f5_tests) && nrow(f5_tests)) {
    r5 <- f5_tests[f5_tests$arm == "paper_all_inputs", ]
    ref_sig <- if (nrow(r5)) r5$ks_perm_p[1] < 0.05 else NA
    concl[[length(concl) + 1L]] <- data.frame(
      figure = "5c", statement = "WGD+ vs WGD- delta f (KS, permutation)", arm = f5_tests$arm,
      value = round(f5_tests$ks_D, 3), interval = sprintf("p = %.3g", f5_tests$ks_perm_p),
      significant = f5_tests$ks_perm_p < 0.05, reference_value = if (nrow(r5)) round(r5$ks_D[1], 3) else NA,
      reference_significant = ref_sig, agrees = (f5_tests$ks_perm_p < 0.05) == ref_sig, stringsAsFactors = FALSE)
  }
  if (!is.null(f5_model) && nrow(f5_model)) {
    x <- f5_model[f5_model$term != "(Intercept)", ]
    r <- f5_model[f5_model$arm == "paper_all_inputs", ]
    r <- r[match(x$term, r$term), ]
    concl[[length(concl) + 1L]] <- data.frame(
      figure = "5b", statement = x$term, arm = x$arm, value = signif(x$Estimate, 3), interval = sprintf("p = %.3g", x$p),
      significant = x$p < 0.05, reference_value = signif(r$Estimate, 3), reference_significant = r$p < 0.05,
      agrees = fig456_agrees(x$Estimate, x$p < 0.05, r$Estimate, r$p < 0.05), stringsAsFactors = FALSE)
  }
  if (!is.null(dom_tab) && nrow(dom_tab)) {
    r6 <- dom_tab[dom_tab$arm == if (any(dom_tab$arm == "paper_all_inputs")) "paper_all_inputs" else "paper_published", ]
    concl[[length(concl) + 1L]] <- data.frame(
      figure = "6", statement = "landscapes whose dominant karyotype changes with the missegregation rate (descriptive)",
      arm = dom_tab$arm, value = dom_tab$share_switching,
      interval = sprintf("%d of %d", dom_tab$switching, dom_tab$landscapes), significant = NA,
      reference_value = if (nrow(r6)) r6$share_switching[1] else NA, reference_significant = NA, agrees = NA,
      stringsAsFactors = FALSE)
  }
  if (!length(concl)) return(NULL)
  do.call(rbind, concl)
}

# Reading a landscape once per fit folder.
fig456_landscape <- local({
  cache <- new.env()
  function(dir) {
    if (!exists(dir, envir = cache, inherits = FALSE)) {
      path <- file.path(dir, "landscape.Rds")
      assign(dir, if (file.exists(path)) readRDS(path) else NULL, envir = cache)
    }
    get(dir, envir = cache, inherits = FALSE)
  }
})
