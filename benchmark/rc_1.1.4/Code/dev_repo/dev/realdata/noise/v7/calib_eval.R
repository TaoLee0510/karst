#!/usr/bin/env Rscript
# Steps 2-3 of dev/realdata/noise/PLAN_v7.md, local, from the records of
# v7/stats_task.R (DEV_RESULTS/realdata/noise/v7/stats/<tag>.rds).
#   Rscript dev/realdata/noise/v7/calib_eval.R
# Writes DEV_RESULTS/realdata/noise/v7/calibration.rds and prints the tables.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "cloud_level.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
files <- list.files(v7_path("stats"), pattern = "[.]rds$", full.names = TRUE)
recs <- unlist(lapply(files, readRDS), recursive = FALSE)
cat("records:", length(recs), " not ok:", sum(!vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "\n")
recs <- Filter(function(r) identical(r$status, "ok"), recs)
sims <- Filter(function(r) r$family %in% c("e", "c2", "paper"), recs)
D <- do.call(rbind, lapply(sims, function(r) data.frame(
  family = r$family, design = r$design, dataset_id = r$dataset_id, half = r$half, landscape_id = r$landscape_id,
  wavelength = if (is.null(r$wavelength)) 0.8 else r$wavelength, m = r$multiple, truth = r$beta, dt = r$dt,
  bottleneck = if (is.null(r$bottleneck)) NA else r$bottleneck, depth_t0 = r$depth_t0, depth_t1 = r$depth_t1,
  window = if (is.null(r$window)) NA else r$window, r$stats, stringsAsFactors = FALSE)))
D$pair <- sprintf("%g/%g", D$depth_t0, D$depth_t1)
cat("datasets:\n"); print(table(D$family, D$half))
auc <- function(s, y) {
  s[is.na(s)] <- -Inf
  r <- rank(s); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n0 * n1)
}
stat_of <- function(x, s) if (s == "A") x$S_A else x$S_B
calibrate <- function(x, s) {
  v <- stat_of(x, s); u <- is.finite(log2(v))
  f <- stats::lm(log2(v[u]) ~ log2(x$m[u]))
  list(statistic = s, c = unname(stats::coef(f)[1]), d = unname(stats::coef(f)[2]), sigma = summary(f)$sigma,
       d_ci = unname(stats::confint(f)[2, ]), n = sum(u), r2 = summary(f)$r.squared)
}
crit <- function(x, used) {
  r2 <- used >= 2 * 5e-5 * (1 - 1e-9); over <- used > 2 * x$truth * (1 + 1e-9)
  one <- x$m == 1; hi <- x$m >= 4
  data.frame(datasets = nrow(x), false_raise = mean(r2[one]), n_x1 = sum(one), useful = mean(r2[hi]), n_hi = sum(hi),
             overshoot = mean(over), pass = mean(r2[one]) <= 0.10 && mean(r2[hi]) >= 0.50 && mean(over) <= 0.05)
}
dev <- D[D$half == "dev" & D$family %in% c("e", "c2"), ]
lab <- dev$m == 1 | dev$m >= 4
aucs <- c(A = auc(dev$S_A[lab], dev$m[lab] >= 4), B = auc(dev$S_B[lab], dev$m[lab] >= 4))
cat("\n== Step 2: development half (e and C2 together), AUC x1 against >= x4 ==\n"); print(round(aucs, 3))
cat(sprintf("datasets: x1 %d, >= x4 %d; S_A not computable: %d (x1 %d, >= x4 %d)\n", sum(dev$m == 1), sum(dev$m >= 4),
            sum(is.na(dev$S_A)), sum(is.na(dev$S_A) & dev$m == 1), sum(is.na(dev$S_A) & dev$m >= 4)))
pass <- round(aucs, 3) >= 0.80
chosen <- if (!any(pass)) NA_character_ else {
  a <- round(aucs[pass], 3); best <- names(a)[a == max(a)]
  if ("A" %in% best) "A" else best[1]
}
out <- list(data = D, aucs = aucs, chosen = chosen)
if (is.na(chosen)) {
  verdict <- "step 2 failed: AUC below 0.80 for both statistics"
} else {
  cal <- calibrate(dev, chosen); out$calibration <- cal
  cat(sprintf("chosen %s; calibration log2 S = %.3f + %.3f log2(multiple); d 95%% CI [%.3f, %.3f]; sigma %.3f; R2 %.3f; n %d\n",
              chosen, cal$c, cal$d, cal$d_ci[1], cal$d_ci[2], cal$sigma, cal$r2, cal$n))
  if (!(cal$d_ci[1] > 0)) {
    verdict <- "step 2 failed: slope d not above 0"
  } else {
    D$beta_used <- cl_beta_used(stat_of(D, chosen), cal)
    tests <- rbind(cbind(bed = "e confirm", crit(D[D$family == "e" & D$half == "confirm", ],
                                                 D$beta_used[D$family == "e" & D$half == "confirm"])),
                   cbind(bed = "C2 confirm", crit(D[D$family == "c2" & D$half == "confirm", ],
                                                  D$beta_used[D$family == "c2" & D$half == "confirm"])))
    cat("\n== Step 3 ==\n"); print(tests, digits = 3, row.names = FALSE)
    verdict <- if (all(tests$pass)) "step 3 passed" else paste("step 3 failed:", paste(tests$bed[!tests$pass], collapse = ", "))
    out$tests <- tests

    cat("\n== Report only: AUC per bed ==\n")
    for (b in list(c("e", "confirm"), c("c2", "confirm"), c("paper", "paper"))) {
      x <- D[D$family == b[1] & D$half == b[2] & (D$m == 1 | D$m >= 4), ]
      cat(sprintf("%s %s: %.3f\n", b[1], b[2], auc(stat_of(x, chosen), x$m >= 4)))
    }
    cat("\n== Report only: median log2 S by family and multiple (all halves) ==\n")
    print(round(tapply(log2(stat_of(D, chosen)), list(D$family, D$m), stats::median, na.rm = TRUE), 2))
    cat("\n== Report only: calibrations per family (e, C2), and their difference ==\n")
    fam <- function(x) { f <- calibrate(x, chosen); c(c = f$c, d = f$d, sigma = f$sigma, r2 = f$r2,
                                                      level_x4 = f$c + 2 * f$d) }
    both <- D[D$family %in% c("e", "c2"), ]
    per_fam <- rbind(e_dev = fam(dev[dev$family == "e", ]), c2_dev = fam(dev[dev$family == "c2", ]),
                     e_all = fam(both[both$family == "e", ]), c2_all = fam(both[both$family == "c2", ]))
    print(round(per_fam, 3))
    set.seed(7001L)
    lands <- unique(both$landscape_id)
    bs <- t(replicate(2000, {
      pick <- sample(lands, replace = TRUE)
      x <- do.call(rbind, lapply(pick, function(l) both[both$landscape_id == l, ]))
      fam(x[x$family == "e", ])[c("c", "d", "level_x4")] - fam(x[x$family == "c2", ])[c("c", "d", "level_x4")]
    }))
    diff_ci <- rbind(estimate = per_fam["e_all", c("c", "d", "level_x4")] - per_fam["c2_all", c("c", "d", "level_x4")],
                     lo = apply(bs, 2, stats::quantile, 0.025), hi = apply(bs, 2, stats::quantile, 0.975))
    cat("e minus C2 (all 20 landscapes, landscape bootstrap 2000):\n"); print(round(diff_ci, 3))
    out$per_family <- per_fam; out$family_difference <- diff_ci
    cat("\n== Report only: cross-family transfer (calibrate on one family's dev half, test the other's confirm half) ==\n")
    cross <- do.call(rbind, lapply(list(c("e", "c2"), c("c2", "e")), function(p) {
      cf <- calibrate(dev[dev$family == p[1], ], chosen)
      x <- D[D$family == p[2] & D$half == "confirm", ]
      cbind(calibrated_on = p[1], tested_on = p[2], crit(x, cl_beta_used(stat_of(x, chosen), cf)))
    }))
    print(cross, digits = 3, row.names = FALSE); out$cross <- cross
    cat("\n== Report only: raised >= 2x by family, design and multiple (confirm half) ==\n")
    cf <- D[D$half == "confirm", ]; cf$r2 <- cf$beta_used >= 1e-4 * (1 - 1e-9)
    print(round(tapply(cf$r2, list(paste(cf$family, cf$design), cf$m), mean), 2))
    cat("by depth pair:\n"); print(round(tapply(cf$r2, list(cf$pair, cf$m), mean), 2))
    cat("by wavelength:\n"); print(round(tapply(cf$r2, list(cf$wavelength, cf$m), mean), 2))
    cat("e by bottleneck (x1 only):\n"); e1 <- cf[cf$family == "e" & cf$m == 1, ]; print(round(tapply(e1$r2, e1$bottleneck, mean), 2))
    pp <- D[D$family == "paper", ]
    cat("\n== Report only: paper simulations ==\n"); print(crit(pp, pp$beta_used), digits = 3, row.names = FALSE)
    R <- do.call(rbind, lapply(Filter(function(r) r$family %in% c("salehi", "numbat"), recs), function(r)
      data.frame(family = r$family, line = if (r$family == "salehi") r$line else r$cancer,
                 group = if (r$family == "salehi") r$group else NA, dataset_id = r$dataset_id, dt = r$dt, r$stats,
                 stringsAsFactors = FALSE)))
    R$beta_used <- cl_beta_used(stat_of(R, chosen), cal); out$real <- R
    cat("\n== Report only: real data, beta_used (multiple of 5e-5) ==\n")
    print(do.call(rbind, lapply(split(R, list(R$family, R$line), drop = TRUE), function(x)
      data.frame(family = x$family[1], line = x$line[1], n = nrow(x), raised = sum(x$beta_used > 5e-5 * (1 + 1e-9)),
                 raised2 = sum(x$beta_used >= 1e-4 * (1 - 1e-9)), median_x = stats::median(x$beta_used / 5e-5),
                 max_x = max(x$beta_used / 5e-5)))), digits = 3, row.names = FALSE)
    out$data <- D
  }
}
cat("\nVERDICT:", verdict, "\n")
out$verdict <- verdict
dev_save_rds(out, v7_path("calibration.rds"))
