# The paper's Fig. 5 (07_wgd.R) on the Salehi fits of every arm
# (R/paper_fig456.R): the two p53 knock-out lines (SA906a, SA906b), fits with
# xv > 0, one per trajectory and minobs (the longest training, then the
# highest xv, as 07).
#   5b  fitness of the frequent karyotypes and their neighbours against the
#       number of chromosomes altered from the modal copy number (dk):
#       piecewise mixed model fitness ~ dk * WGD + dk:(dk > 3) * WGD + (1|id)
#   5c  fitness effects (delta f) of the one-chromosome changes of each
#       frequent karyotype, WGD+ (mean copy number > 3) against WGD-:
#       Kolmogorov-Smirnov statistic of the per-fit means and its permutation
#       p value (1e4 permutations, seed 1, as 07)
# 5a uses the counts only and does not depend on the fits.
#
#   Rscript scripts/fig5_wgd.R
#
# FIG456_EXCLUDE=0 keeps the trajectories that are not analysable. Writes
# BENCH_RESULTS/analysis/fig456/fig5{,_all_inputs}.Rds (per arm: tests,
# model, the delta f and fitness tables) and fig5_tests.tsv.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})
paths <- bench_paths()
exclude <- Sys.getenv("FIG456_EXCLUDE", "1") != "0"
out <- file.path(paths$results, "analysis", "fig456")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
is_wgd <- function(k, cut = 3) mean(k) > cut
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

one_arm <- function(arm) {
  rows <- fig456_fit_rows(arm, paths, exclude)
  fits_raw <- do.call(rbind, lapply(rows, function(r) utils::head(r[c("fi", "min_obs", "ntrain", "xv")], 1)))
  fits_raw <- subset(fits_raw, !is.na(xv) & xv > 0 & grepl("SA906", fi))
  if (!nrow(fits_raw)) return(NULL)
  fits_raw$traj <- ifelse(grepl("906a", fits_raw$fi, TRUE), "trajA", "trajB")
  ord <- with(fits_raw, order(fi, min_obs, -ntrain, -xv))
  fits_sorted <- fits_raw[ord, ]
  fits <- fits_sorted[!duplicated(fits_sorted[c("fi", "min_obs")]), ]
  rownames(fits) <- NULL
  dfl_list <- list()
  dffit_list <- list()
  for (i in seq_len(nrow(fits))) {
    land <- fig456_landscape(fig456_fit_dir(arm, fits$fi[i], fits$min_obs[i], paths))
    if (is.null(land)) next
    k_land <- do.call(rbind, lapply(land$k, paper_s2v))
    dk_land <- apply(k_land, 1, function(ki) sum(ki != get_mode(ki)))
    lut <- land$mean
    names(lut) <- land$k
    keep <- land$fq | land$nn
    dffit_list[[length(dffit_list) + 1]] <- data.frame(
      fitness = lut[land$k[keep]], dk = dk_land[keep],
      wgd = sapply(land$k[keep], function(j) is_wgd(paper_s2v(j))), id = fits$fi[i], k_id = land$k[keep],
      traj = fits$traj[i])
    fq <- land$k[land$fq]
    dfl_list[[length(dfl_list) + 1]] <- do.call(rbind, lapply(fq, function(fqi) {
      nn <- apply(paper_gen_all_neighbours(fqi, TRUE), 1, paste, collapse = ".")
      data.frame(df = lut[fqi] - lut[nn], wgd = is_wgd(paper_s2v(fqi)), id = fits$fi[i], traj = fits$traj[i],
                 k_id = fqi)
    }))
  }
  dfl <- do.call(rbind, dfl_list)
  dffit <- do.call(rbind, dffit_list)
  dfl$wgd_str <- ifelse(dfl$wgd, "WGD+", "WGD-")
  dffit$wgd_str <- ifelse(dffit$wgd, "WGD+", "WGD-")
  # 5c: permutation KS on the per-fit means
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
  means <- tapply(cl$df, cl$wgd_str, stats::median)
  # 5b: piecewise mixed model
  dffit$cut3 <- dffit$dk > 3
  mod <- tryCatch(lmer(fitness ~ dk * wgd_str + dk:cut3 * wgd_str + (1 | id), data = dffit), error = function(e) e)
  coefs <- if (inherits(mod, "error")) NULL else {
    cf <- as.data.frame(stats::coef(summary(mod)))
    cf$term <- rownames(cf)
    cf
  }
  tests <- data.frame(arm = arm$key, fit = arm$fit, fits = nrow(fits), trajectories = length(unique(fits$fi)),
                      landscapes = length(dffit_list), fits_wgd_plus = sum(cl$wgd_str == "WGD+"),
                      fits_wgd_minus = sum(cl$wgd_str == "WGD-"), ks_D = ks[["D"]], ks_perm_p = ks[["p"]],
                      median_df_wgd_minus = unname(means["WGD-"]), median_df_wgd_plus = unname(means["WGD+"]),
                      stringsAsFactors = FALSE)
  list(tests = tests, coefs = if (is.null(coefs)) NULL else cbind(arm = arm$key, coefs), dfl = cbind(arm = arm$key, dfl),
       dffit = cbind(arm = arm$key, dffit), fits = fits)
}

res <- lapply(seq_len(nrow(fig456_arms)), function(i) {
  arm <- fig456_arms[i, ]
  r <- tryCatch(one_arm(arm), error = function(e) {
    cat("Arm", arm$key, "failed:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) print(r$tests, row.names = FALSE)
  r
})
names(res) <- fig456_arms$key
saveRDS(res, file.path(out, if (exclude) "fig5.Rds" else "fig5_all_inputs.Rds"))
tests <- do.call(rbind, lapply(res, `[[`, "tests"))
utils::write.table(tests, file.path(out, if (exclude) "fig5_tests.tsv" else "fig5_tests_all_inputs.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
print(tests, row.names = FALSE)
