#!/usr/bin/env Rscript
# Steps 2-3 of dev/realdata/noise/PLAN_v6.md, local, from the records of
# stats_task.R (DEV_RESULTS/realdata/noise/v6/stats/<tag>.rds).
#   Rscript dev/realdata/noise/v6/calib_eval.R
# Writes DEV_RESULTS/realdata/noise/v6/calibration.rds and prints the tables.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "cloud_level.R"))
root <- file.path(dev_results, "realdata", "noise", "v6", "stats")
recs <- unlist(lapply(list.files(root, pattern = "[.]rds$", full.names = TRUE), function(f)
  lapply(readRDS(f), function(r) c(r, tag = sub("[.]rds$", "", basename(f))))), recursive = FALSE)
cat("records:", length(recs), " not ok:", sum(!vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "\n")
recs <- Filter(function(r) identical(r$status, "ok"), recs)
set_of <- function(r) {
  if (r$kind == "c2") return(sub("^c2_r1_(.*)_none$", "C2 \\1", r$tag))
  if (r$kind == "paper") return(sprintf("paper x%g", r$truth / 5e-5))
  if (r$set == "e4") return(if (abs(r$truth - 2e-4) < 1e-12) "e4 x4" else "e4 x16")
  r$set
}
D <- do.call(rbind, lapply(Filter(function(r) r$kind %in% c("eset", "c2", "paper"), recs), function(r)
  data.frame(kind = r$kind, set = set_of(r), dataset_id = r$dataset_id,
             half = if (r$kind == "paper") "paper" else r$half, landscape_id = r$landscape_id, depth = r$depth,
             dt = r$dt, window = if (is.null(r$window)) NA else r$window, truth = r$truth, m = r$truth / 5e-5,
             r$stats, stringsAsFactors = FALSE)))
D$family <- ifelse(D$kind == "c2", "C2", ifelse(D$kind == "eset", "e", "paper"))
auc <- function(s, y) {   # y TRUE = high; NA statistic counts as the lowest
  s[is.na(s)] <- -Inf
  r <- rank(s); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n0 * n1)
}
dev <- D[D$half == "dev", ]
lab <- dev$m == 1 | dev$m >= 4
aucs <- c(A = auc(dev$S_A[lab], dev$m[lab] >= 4), B = auc(dev$S_B[lab], dev$m[lab] >= 4))
cat("\n== Step 2: development half, AUC (x1 against >= x4) ==\n"); print(round(aucs, 3))
cat(sprintf("datasets: x1 %d, >= x4 %d; S_A not computable: %d (x1 %d, >= x4 %d)\n", sum(dev$m == 1), sum(dev$m >= 4),
            sum(is.na(dev$S_A)), sum(is.na(dev$S_A) & dev$m == 1), sum(is.na(dev$S_A) & dev$m >= 4)))
pass <- round(aucs, 3) >= 0.80
chosen <- if (!any(pass)) NA_character_ else {
  a <- round(aucs[pass], 3); best <- names(a)[a == max(a)]
  if ("A" %in% best) "A" else best[1]
}
verdict <- NULL; cal <- NULL; family_slopes <- NULL; tests <- NULL; real <- NULL
stat_of <- function(x, s) if (s == "A") x$S_A else x$S_B
if (is.na(chosen)) {
  verdict <- "step 2 failed: AUC below 0.80 for both statistics"
} else {
  cat("chosen statistic:", chosen, "\n")
  s <- stat_of(dev, chosen); use <- is.finite(log2(s))
  fit <- stats::lm(log2(s[use]) ~ log2(dev$m[use]))
  ci <- stats::confint(fit)[2, ]
  cal <- list(statistic = chosen, c = unname(stats::coef(fit)[1]), d = unname(stats::coef(fit)[2]),
              sigma = summary(fit)$sigma, d_ci = unname(ci), n = sum(use), r2 = summary(fit)$r.squared)
  cat(sprintf("calibration: log2 S = %.3f + %.3f log2(multiple); d 95%% CI [%.3f, %.3f]; sigma %.3f; R2 %.3f; n %d\n",
              cal$c, cal$d, ci[1], ci[2], cal$sigma, cal$r2, cal$n))
  family_slopes <- do.call(rbind, lapply(c("e", "C2"), function(fm) {
    x <- dev[dev$family == fm, ]; sx <- stat_of(x, chosen); u <- is.finite(log2(sx))
    f <- stats::lm(log2(sx[u]) ~ log2(x$m[u]))
    data.frame(family = fm, n = sum(u), c = stats::coef(f)[1], d = stats::coef(f)[2], sigma = summary(f)$sigma)
  }))
  cat("report only, per family:\n"); print(family_slopes, digits = 3, row.names = FALSE)
  if (!(ci[1] > 0)) {
    verdict <- "step 2 failed: slope d not above 0"
  } else {
    D$beta_used <- cl_beta_used(stat_of(D, chosen), cal)
    D$raised2 <- D$beta_used >= 2 * 5e-5 * (1 - 1e-9)
    crit <- function(x) {
      one <- x[x$m == 1, ]; hi <- x[x$m >= 4, ]
      data.frame(datasets = nrow(x), false_raise = mean(one$raised2), n_x1 = nrow(one), useful = mean(hi$raised2),
                 n_hi = nrow(hi), overshoot = mean(x$beta_used > 2 * x$truth * (1 + 1e-9)),
                 pass = mean(one$raised2) <= 0.10 && mean(hi$raised2) >= 0.50 &&
                   mean(x$beta_used > 2 * x$truth * (1 + 1e-9)) <= 0.05)
    }
    tests <- rbind(cbind(bed = "confirm", crit(D[D$half == "confirm", ])), cbind(bed = "paper", crit(D[D$half == "paper", ])))
    cat("\n== Step 3 ==\n"); print(tests, digits = 3, row.names = FALSE)
    verdict <- if (all(tests$pass)) "step 3 passed" else paste("step 3 failed:", paste(tests$bed[!tests$pass], collapse = ", "))
    cat("\nreport only: per set (all halves)\n")
    per <- do.call(rbind, lapply(split(D, list(D$set, D$half), drop = TRUE), function(x)
      data.frame(set = x$set[1], half = x$half[1], n = nrow(x), raised = mean(x$beta_used > 5e-5 * (1 + 1e-9)),
                 raised2 = mean(x$raised2), median_used_x = stats::median(x$beta_used / 5e-5),
                 auc_note = NA)))
    print(per[order(per$half, per$set), -7], digits = 3, row.names = FALSE)
    cat("\nreport only: raised >= 2x by depth and family\n"); print(round(tapply(D$raised2, list(D$family, D$depth), mean), 3))
    cat("report only: paper raised >= 2x by window\n"); pp <- D[D$kind == "paper", ]; print(round(tapply(pp$raised2, list(pp$m, pp$window), mean), 2))
    cat("report only: AUC (x1 vs >= x4) in the test beds:", sprintf("confirm %.3f, paper %.3f\n",
        { x <- D[D$half == "confirm" & (D$m == 1 | D$m >= 4), ]; auc(stat_of(x, chosen), x$m >= 4) },
        { x <- D[D$half == "paper" & (D$m == 1 | D$m >= 4), ]; auc(stat_of(x, chosen), x$m >= 4) }))
    R <- do.call(rbind, lapply(Filter(function(r) r$kind %in% c("salehi", "numbat"), recs), function(r)
      data.frame(kind = r$kind, line = if (r$kind == "salehi") r$line else r$cancer,
                 group = if (r$kind == "salehi") r$group else NA, dataset_id = r$dataset_id, dt = r$dt, r$stats,
                 stringsAsFactors = FALSE)))
    R$beta_used <- cl_beta_used(stat_of(R, chosen), cal)
    real <- R
    cat("\nreport only: real data, beta_used (multiple of 5e-5)\n")
    print(do.call(rbind, lapply(split(R, list(R$kind, R$line), drop = TRUE), function(x)
      data.frame(kind = x$kind[1], line = x$line[1], n = nrow(x), raised = sum(x$beta_used > 5e-5 * (1 + 1e-9)),
                 raised2 = sum(x$beta_used >= 1e-4 * (1 - 1e-9)), median_x = stats::median(x$beta_used / 5e-5),
                 max_x = max(x$beta_used / 5e-5)))), digits = 3, row.names = FALSE)
  }
}
cat("\nVERDICT:", verdict, "\n")
dev_save_rds(list(data = D, aucs = aucs, chosen = chosen, calibration = cal, family_slopes = family_slopes,
                  tests = tests, real = real, verdict = verdict),
             file.path(dev_results, "realdata", "noise", "v6", "calibration.rds"))
