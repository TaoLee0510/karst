# Analysis of the Fig. 2d truth study (scripts/fig2d_truth.R): Spearman of
# observed and predicted change for KARST, alfak and the true landscape,
# raw and with the sampling noise held fixed (partial Spearman, as
# scripts/fig2d_noise.R), with intervals from resampling simulations within
# wavelengths; the landscape accuracy of the same fits.
#
#   F2T_ROOT=<root> Rscript scripts/fig2d_truth_analyse.R STEP [SIMS]
#
# STEP 1: observed change on the 1000 cells stored and on subsamples of
# n = 400, 200, 100 cells (the mean over 20 draws). STEP 2: on the
# Salehi-sized samples the fits saw and on the 1000 cells. SIMS: "all"
# (default) or "pilot" (replicates 01-05 only). Rows: every simulation,
# design (w2/w4/w8, two time points), minobs; the daughter is the sampling
# time nearest 15 days after the parent (the paper's Salehi evaluation day),
# and the first one after it as a second set. Selection: none (every minobs)
# and the paper's rule (each method's minobs of highest cross-validation
# score per simulation and design; the true landscape on that method's
# rows).
# Output: F2T_ROOT/analysis/step<k>/{fig2d.tsv, differences.tsv,
# accuracy.tsv, follow_noise.tsv, figures/*.pdf}.
#
# F2T_SOURCE=misseg (scripts/fig2d_truth.R --source=misseg): the first
# argument is the arm, beta_fixed or beta_true (the rate 5e-5's rows are
# beta_fixed's in both), SIMS "pilot" means replicates 01-02; groups are the
# missegregation rates, resampled within themselves; a fourth forecast, the
# true landscape at 5e-5 (truth_fixed). Output F2T_ROOT/analysis/<arm>/ and,
# once both arms are done, F2T_ROOT/analysis/figures/ (Fig. 2d, dispersion,
# share of cells on the list and fits that failed, against the rate;
# F2T_SALEHI_TSV, the round's analysis/fig2d_noise/fig2d_noise.tsv, gives
# the Salehi values drawn as reference lines).

for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(TRUE)
misseg <- identical(Sys.getenv("F2T_SOURCE"), "misseg")
which_sims <- if (length(args) > 1) args[[2]] else "all"
root <- Sys.getenv("F2T_ROOT", "")
if (!nzchar(root)) stop("F2T_ROOT is not set")
if (misseg) {
  arm_dir <- args[[1]]
  stopifnot(arm_dir %in% c("beta_fixed", "beta_true"))
  step <- 2L
  out <- file.path(root, "analysis", paste0(arm_dir, if (which_sims != "all") paste0("_", which_sims)))
  files <- list.files(file.path(root, arm_dir, "rows"), pattern = "[.]Rds$", full.names = TRUE)
  # the rate 5e-5: the two arms' fits are the same (beta_fixed's)
  if (arm_dir == "beta_true") {
    files <- c(files, list.files(file.path(root, "beta_fixed", "rows"), pattern = "_m_5e-05_rep_.*[.]Rds$",
                                 full.names = TRUE))
  }
} else {
  step <- as.integer(args[[1]])
  out <- file.path(root, "analysis", paste0("step", step, if (which_sims != "all") paste0("_", which_sims)))
  files <- list.files(file.path(root, paste0("step", step), "rows"), pattern = "[.]Rds$", full.names = TRUE)
}
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("F2T_BOOT", "500"))
arms <- c("karst", "alfak", "truth", if (misseg) "truth_fixed")
arm_label <- c(karst = "KARST", alfak = "alfak", truth = "true landscape",
               truth_fixed = "true landscape, forecast at 5e-5", noise = "sampling noise alone")
res <- lapply(files, readRDS)
grab <- function(name) do.call(rbind, lapply(res, `[[`, name))
preds <- grab("preds")
obs <- grab("obs")
fits <- grab("fits")
draws <- if (step == 1L) grab("draws") else NULL
if (which_sims == "pilot") {
  keep <- as.integer(obs$abmrep) <= (if (misseg) 2L else 5L)
  obs <- obs[keep, ]
}
# the groups simulations are resampled within: wavelengths, or rates
obs$stratum <- if (misseg) obs$pmis else obs$w
cat("simulations:", length(unique(obs$sim_id)), " observation rows:", nrow(obs), "\n")

key_of <- function(d) paste(d$sim_id, d$design, d$min_obs, d$j)
obs$key <- key_of(obs)
preds$key <- key_of(preds)
for (a in arms) obs[[paste0("pred_", a)]] <- preds$pred[preds$arm == a][match(obs$key, preds$key[preds$arm == a])]
# the daughter nearest 15 days after the parent, per simulation
near <- tapply(seq_len(nrow(obs)), obs$sim_id, function(i) {
  dt <- unique(obs[i, c("j", "dt")])
  dt$j[which.min(abs(dt$dt - 15))]
})
obs$j15 <- obs$j == near[obs$sim_id]
# cross-validation scores for the paper's rule
fits$fkey <- paste(fits$sim_id, fits$design, fits$min_obs)
obs$fkey <- paste(obs$sim_id, obs$design, obs$min_obs)
for (m in c("karst", "alfak")) {
  obs[[paste0("rxv_", m)]] <- fits$Rxv[fits$method == m][match(obs$fkey, fits$fkey[fits$method == m])]
}

sp <- function(a, b) {
  k <- is.finite(a) & is.finite(b)
  if (sum(k) < 5) return(NA_real_)
  stats::cor(a[k], b[k], method = "spearman")
}
# partial Spearman with z held fixed: the correlation of the ranks of y and a
# after each is regressed on the ranks of z (scripts/fig2d_noise.R), in its
# closed form from the three rank correlations
partial <- function(y, a, z) {
  ok <- is.finite(y) & is.finite(a) & is.finite(z)
  if (sum(ok) < 5) return(NA_real_)
  r <- stats::cor(cbind(y[ok], a[ok], z[ok]), method = "spearman")
  (r[1, 2] - r[1, 3] * r[2, 3]) / sqrt((1 - r[1, 3]^2) * (1 - r[2, 3]^2))
}

# observed axes: name -> list(obs = matrix rows x draws, noise = vector)
axes_for <- function(d) {
  ax <- list()
  if (step == 1L) {
    ax[["1000 cells"]] <- list(obs = matrix(d$obs_1000), noise = d$noise_1000)
    for (n in c(400L, 200L, 100L)) {
      dn <- draws[draws$n == n, ]
      dn$key <- key_of(dn)
      m <- matrix(NA_real_, nrow(d), max(dn$draw))
      for (r in seq_len(ncol(m))) {
        x <- dn[dn$draw == r, ]
        m[, r] <- x$obs[match(d$key, x$key)]
      }
      nz <- dn$noise[dn$draw == 1L][match(d$key, dn$key[dn$draw == 1L])]
      ax[[paste(n, "cells")]] <- list(obs = m, noise = nz)
    }
  } else {
    ax[["Salehi-sized samples"]] <- list(obs = matrix(d$obs_fit), noise = d$noise_fit)
    ax[["1000 cells"]] <- list(obs = matrix(d$obs_1000), noise = d$noise_1000)
  }
  ax
}

# statistics of one row set (indices `ix` into d) on one axis: per arm raw and
# partial Spearman (mean over draws), Spearman of the observed change with the
# noise alone, and of each arm's predicted change with the noise
stats_of <- function(d, ax, ix, arms_here) {
  o <- ax$obs[ix, , drop = FALSE]
  z <- ax$noise[ix]
  v <- c()
  for (a in arms_here) {
    p <- d[[paste0("pred_", a)]][ix]
    v[paste(a, "raw", sep = "|")] <- mean(apply(o, 2, sp, b = p))
    v[paste(a, "noise held fixed", sep = "|")] <- mean(apply(o, 2, function(y) partial(y, p, z)))
    v[paste(a, "step follows noise", sep = "|")] <- sp(p, z)
  }
  v["noise|raw"] <- mean(apply(o, 2, sp, b = z))
  v
}

# resample simulations within wavelengths (rates)
boot_sets <- function(d, B) {
  u <- unique(d[c("sim_id", "stratum")])
  sims <- split(u$sim_id, u$stratum)
  rows_of <- split(seq_len(nrow(d)), d$sim_id)
  set.seed(1)
  lapply(seq_len(B), function(b) unlist(lapply(sims, function(s) rows_of[sample(s, length(s), replace = TRUE)]),
                                       use.names = FALSE))
}

run_set <- function(d, label, arms_here, groups) {
  est <- list()
  dif <- list()
  for (g in names(groups)) {
    dg <- d[groups[[g]](d), , drop = FALSE]
    if (nrow(dg) < 20) next
    axes <- axes_for(dg)
    bs <- boot_sets(dg, B)
    for (an in names(axes)) {
      ax <- axes[[an]]
      full <- stats_of(dg, ax, seq_len(nrow(dg)), arms_here)
      boot <- do.call(cbind, parallel::mclapply(bs, function(ix) stats_of(dg, ax, ix, arms_here),
                                                mc.cores = as.integer(Sys.getenv("F2T_CORES", "3"))))
      ci <- apply(boot, 1L, stats::quantile, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      for (i in seq_along(full)) {
        parts <- strsplit(names(full)[i], "|", fixed = TRUE)[[1]]
        est[[length(est) + 1L]] <- data.frame(set = label, group = g, observed = an, rows = nrow(dg),
                                              simulations = length(unique(dg$sim_id)), arm = parts[1],
                                              measure = parts[2], rho = round(full[i], 3),
                                              lo = round(ci[1, i], 3), hi = round(ci[2, i], 3))
      }
      pairs <- list(c("karst", "alfak"), c("truth", "alfak"), c("truth", "karst"), c("truth", "truth_fixed"))
      for (pr in pairs) {
        if (!all(pr %in% arms_here)) next
        for (meas in c("raw", "noise held fixed")) {
          i <- paste(pr[1], meas, sep = "|")
          j <- paste(pr[2], meas, sep = "|")
          dv <- boot[i, ] - boot[j, ]
          dif[[length(dif) + 1L]] <- data.frame(set = label, group = g, observed = an, comparison =
                                                  paste(arm_label[pr[1]], "-", arm_label[pr[2]]), measure = meas,
                                                difference = round(full[[i]] - full[[j]], 3),
                                                lo = round(stats::quantile(dv, 0.025, na.rm = TRUE, names = FALSE), 3),
                                                hi = round(stats::quantile(dv, 0.975, na.rm = TRUE, names = FALSE), 3))
        }
      }
    }
  }
  list(est = do.call(rbind, est), dif = do.call(rbind, dif))
}

groups <- if (misseg) {
  rates <- sort(unique(obs$pmis))
  stats::setNames(lapply(rates, function(r) function(d) d$pmis == r), paste("rate", format(rates, scientific = TRUE)))
} else {
  c(list(`all wavelengths` = function(d) rep(TRUE, nrow(d))),
    stats::setNames(lapply(c(0.2, 0.4, 0.8, 1.6), function(w) function(d) d$w == w),
                    paste("wavelength", c(0.2, 0.4, 0.8, 1.6))))
}
# the paper's rule: pooled over wavelengths; per rate
rule_groups <- if (misseg) groups else groups[1]
complete <- Reduce(`&`, lapply(arms, function(a) is.finite(obs[[paste0("pred_", a)]])))
results <- list()
for (dset in c("15 days", "next sampling time")) {
  base <- obs[complete & (if (dset == "15 days") obs$j15 else obs$j == 1L), ]
  # no selection: every minobs, rows all three arms have
  results[[length(results) + 1L]] <- run_set(base, paste(dset, "| no selection"), arms, groups)
  # the paper's rule, per method; the true landscape on that method's rows
  for (m in c("karst", "alfak")) {
    sel <- select_by_cv(base, by = c("sim_id", "design"), choice = "min_obs", score = paste0("rxv_", m))
    results[[length(results) + 1L]] <- run_set(sel, paste(dset, "| paper rule,", arm_label[m], "rows"),
                                               c(m, intersect(c("truth", "truth_fixed"), arms)), rule_groups)
  }
}
est <- do.call(rbind, lapply(results, `[[`, "est"))
dif <- do.call(rbind, lapply(results, `[[`, "dif"))
utils::write.table(est, file.path(out, "fig2d.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(dif, file.path(out, "differences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# landscape accuracy of the fits behind the rows (no selection), paired
acc_rows <- unique(obs[complete & obs$j15, c("sim_id", "w", "pmis", "design", "min_obs")])
acc_rows$fkey <- paste(acc_rows$sim_id, acc_rows$design, acc_rows$min_obs)
acc <- list()
for (g in names(groups)) {
  a <- acc_rows[groups[[g]](acc_rows), ]
  for (s in c("r", "rfq", "rnn", "rd2")) {
    vk <- fits[[s]][fits$method == "karst"][match(a$fkey, fits$fkey[fits$method == "karst"])]
    va <- fits[[s]][fits$method == "alfak"][match(a$fkey, fits$fkey[fits$method == "alfak"])]
    ok <- is.finite(vk) & is.finite(va)
    acc[[length(acc) + 1L]] <- data.frame(group = g, score = s, fits = sum(ok),
                                          karst_median = round(stats::median(vk[ok]), 3),
                                          alfak_median = round(stats::median(va[ok]), 3),
                                          karst_higher = round(mean(vk[ok] > va[ok]), 3))
  }
}
acc <- do.call(rbind, acc)
utils::write.table(acc, file.path(out, "accuracy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# fits per input (simulation, design, minobs) and the share that worked, and
# the share of the parent's cells on the list (missegregation simulations)
if (misseg) {
  sims_per <- tapply(obs$sim_id, obs$pmis, function(x) length(unique(x)))
  fit_ok <- do.call(rbind, lapply(names(sims_per), function(r) {
    f <- fits[fits$pmis == as.numeric(r), ]
    n <- sims_per[[r]] * 3L * length(unique(obs$min_obs))
    data.frame(rate = r, inputs = n, karst_ok = round(sum(f$method == "karst") / n, 3),
               alfak_ok = round(sum(f$method == "alfak") / n, 3),
               rows_complete = round(mean(complete[obs$pmis == as.numeric(r) & obs$j15]), 3),
               on_list_parent = round(stats::median(obs$in_list_parent_fit[obs$pmis == as.numeric(r) & obs$j15]), 3),
               on_list_daughter = round(stats::median(obs$in_list_daughter_fit[obs$pmis == as.numeric(r) & obs$j15]), 3),
               list_size = stats::median(obs$n_k[obs$pmis == as.numeric(r) & obs$j15]))
  }))
  utils::write.table(fit_ok, file.path(out, "fits_and_list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

options(width = 220)
show <- function(x) print(x, row.names = FALSE)
g1 <- if (misseg) names(groups) else "all wavelengths"
cat("\n== Fig. 2d,", paste(g1, collapse = ", "), "\n")
show(est[est$group %in% g1 & est$measure != "step follows noise" & est$observed != "1000 cells", ])
cat("\n== Differences\n")
show(dif[dif$group %in% g1 & dif$observed != "1000 cells", ])
cat("\n== Each arm's predicted change against the noise (Spearman)\n")
show(est[est$group %in% g1 & est$measure == "step follows noise" & est$observed != "1000 cells" &
           grepl("no selection", est$set), ])
if (misseg) {
  cat("\n== Fits that worked, rows with every forecast, share of cells on the list\n")
  show(fit_ok)
}
cat("\n== Landscape accuracy of the same fits (no selection)\n")
show(acc)

## figures ------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
cols <- c(KARST = "#D55E00", alfak = "#0072B2", `true landscape` = "#222222",
          `true landscape, forecast at 5e-5` = "#009E73", `sampling noise alone` = "#999999")
e <- est[est$measure %in% c("raw", "noise held fixed") & grepl("no selection", est$set), ]
e$arm_l <- factor(arm_label[e$arm], levels = names(cols))
e$daughter <- sub(" [|].*$", "", e$set)
if (step == 1L) {
  e$cells <- as.integer(sub(" cells", "", e$observed))
  e$measure <- factor(e$measure, levels = c("raw", "noise held fixed"))
  p <- ggplot(e[e$daughter == "15 days", ], aes(cells, rho, colour = arm_l)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = arm_l), alpha = 0.12, colour = NA) +
    geom_line() + geom_point(size = 1.6) +
    scale_x_log10(breaks = c(100, 200, 400, 1000)) +
    scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
    facet_grid(measure ~ group) +
    labs(x = "cells per sample in the observed change (parent and daughter)",
         y = "Spearman, observed vs predicted change",
         caption = paste("Round 6ae2a24 fits (1000 cells per time point), two time points, every minobs;",
                         "daughter nearest 15 days. Bands: 95% bootstrap over simulations.")) +
    theme_bw(base_size = 9) + theme(legend.position = "top")
  ggsave(file.path(out, "figures", "fig2d_truth_cells.pdf"), p, width = 10, height = 5)
} else {
  e$measure <- factor(e$measure, levels = c("raw", "noise held fixed"))
  e$observed <- factor(e$observed, levels = c("Salehi-sized samples", "1000 cells"),
                       labels = c("observed on the Salehi-sized samples", "observed on the 1000 cells"))
  p <- ggplot(e[e$daughter == "15 days", ], aes(measure, rho, colour = arm_l)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.25) +
    scale_colour_manual(values = cols, name = NULL) + facet_grid(observed ~ group) +
    labs(x = NULL, y = "Spearman, observed vs predicted change",
         caption = paste("Fits of inputs subsampled to the Salehi inputs' cell counts, two time points, every minobs;",
                         "daughter nearest 15 days. Bars: 95% bootstrap over simulations.")) +
    theme_bw(base_size = 9) + theme(legend.position = "top")
  ggsave(file.path(out, "figures", "fig2d_truth_salehi_cells.pdf"), p, width = 10, height = 5.5)
}
cat("\nwritten:", out, "\n")

## missegregation rates: both arms against the rate --------------------------
suffix <- if (which_sims != "all") paste0("_", which_sims) else ""
arm_files <- file.path(root, "analysis", paste0(c("beta_fixed", "beta_true"), suffix))
if (misseg && all(file.exists(file.path(arm_files, "fig2d.tsv")))) {
  fig_dir <- file.path(root, "analysis", paste0("figures", suffix))
  dir.create(fig_dir, showWarnings = FALSE)
  beta_label <- c(beta_fixed = "fitted and forecast at 5e-5 (as on Salehi)",
                  beta_true = "fitted and forecast at the true rate")
  both <- do.call(rbind, lapply(c("beta_fixed", "beta_true"), function(a) {
    x <- utils::read.delim(file.path(root, "analysis", paste0(a, suffix), "fig2d.tsv"), stringsAsFactors = FALSE)
    x$beta <- beta_label[[a]]
    x
  }))
  both$rate <- as.numeric(sub("^rate ", "", both$group))
  both$arm_l <- factor(arm_label[both$arm], levels = names(cols))
  both$beta <- factor(both$beta, levels = beta_label)
  # Salehi (round fits, two time points, every trajectory, no selection;
  # scripts/fig2d_noise.R) as reference lines
  sal <- NULL
  sal_file <- Sys.getenv("F2T_SALEHI_TSV", "")
  if (nzchar(sal_file) && file.exists(sal_file)) {
    sal <- utils::read.delim(sal_file, stringsAsFactors = FALSE)
    sal <- sal[sal$group == "2 time points" & sal$subset == "every trajectory" & sal$arm %in% names(arm_label), ]
    sal$arm_l <- factor(arm_label[sal$arm], levels = names(cols))
    sal$measure <- factor(sal$measure, levels = c("raw", "noise held fixed"))
  }
  for (ob in c("Salehi-sized samples", "1000 cells")) {
    e2 <- both[both$set == "15 days | no selection" & both$observed == ob & both$measure %in% c("raw", "noise held fixed"), ]
    e2$measure <- factor(e2$measure, levels = c("raw", "noise held fixed"))
    p <- ggplot(e2, aes(rate, rho, colour = arm_l)) +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = arm_l), alpha = 0.12, colour = NA) +
      geom_line() + geom_point(size = 1.4)
    if (!is.null(sal) && ob == "Salehi-sized samples") {
      p <- p + geom_hline(data = sal, aes(yintercept = rho, colour = arm_l), linetype = "dashed", linewidth = 0.4)
    }
    p <- p + scale_x_log10(breaks = c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4), labels = c("5e-5", "1e-4", "2e-4", "4e-4", "8e-4")) +
      scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
      scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
      facet_grid(measure ~ beta) +
      labs(x = "missegregation rate of the simulation (ABM_varyMisseg, wavelength 0.8)",
           y = "Spearman, observed vs predicted change",
           caption = paste0("Inputs subsampled to the Salehi inputs' cell counts, two time points, every minobs; observed on the ",
                            if (ob == "1000 cells") "1000 cells stored" else "Salehi-sized samples",
                            "; daughter nearest 15 days. Bands: 95% bootstrap over simulations.",
                            if (!is.null(sal) && ob == "Salehi-sized samples") "\nDashed: Salehi, round 6ae2a24 (two time points, every trajectory, no selection)." else "")) +
      theme_bw(base_size = 9) + theme(legend.position = "top")
    ggsave(file.path(fig_dir, if (ob == "1000 cells") "fig2d_misseg_1000_cells.pdf" else "fig2d_misseg.pdf"), p,
           width = 9, height = 5.5)
  }
  # landscape accuracy and the fits that worked
  acc2 <- do.call(rbind, lapply(c("beta_fixed", "beta_true"), function(a) {
    x <- utils::read.delim(file.path(root, "analysis", paste0(a, suffix), "accuracy.tsv"), stringsAsFactors = FALSE)
    x <- x[x$score == "r", ]
    rbind(data.frame(beta = beta_label[[a]], rate = as.numeric(sub("^rate ", "", x$group)), method = "KARST",
                     value = x$karst_median),
          data.frame(beta = beta_label[[a]], rate = as.numeric(sub("^rate ", "", x$group)), method = "alfak",
                     value = x$alfak_median))
  }))
  acc2$beta <- factor(acc2$beta, levels = beta_label)
  p <- ggplot(acc2, aes(rate, value, colour = method)) + geom_line() + geom_point(size = 1.4) +
    scale_x_log10(breaks = c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4), labels = c("5e-5", "1e-4", "2e-4", "4e-4", "8e-4")) +
    scale_colour_manual(values = cols[c("KARST", "alfak")], name = NULL) + facet_grid(. ~ beta) +
    labs(x = "missegregation rate of the simulation", y = "median r, fitted vs true fitness",
         caption = "The fits behind the Fig. 2d rows (no selection); every karyotype of the landscape.") +
    theme_bw(base_size = 9) + theme(legend.position = "top")
  ggsave(file.path(fig_dir, "accuracy_misseg.pdf"), p, width = 7, height = 3.2)
  # how spread out the samples are (fig2d_truth.R --stage=dispersion)
  dfile <- file.path(root, "dispersion.tsv")
  if (file.exists(dfile)) {
    d <- utils::read.delim(dfile, stringsAsFactors = FALSE)
    d$x <- ifelse(d$set == "Salehi", "Salehi", ifelse(d$set == "ABM main set", "main set\n5e-5",
                                                       format(d$pmis, scientific = TRUE)))
    d$x <- factor(d$x, levels = c("main set\n5e-5", format(c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4), scientific = TRUE), "Salehi"))
    long <- do.call(rbind, lapply(c("simpson", "top", "once", "new_in_daughter"), function(m) {
      data.frame(x = d$x, set = d$set, measure = m, value = d[[m]])
    }))
    mlab <- c(simpson = "Simpson diversity", top = "share of the most frequent karyotype",
              once = "share of cells whose karyotype occurs once", new_in_daughter = "share of daughter cells new to the parent")
    long$measure <- factor(mlab[long$measure], levels = mlab)
    p <- ggplot(long, aes(x, value, fill = set)) + geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
      geom_jitter(data = long[long$set == "Salehi", ], width = 0.2, height = 0, size = 0.5, alpha = 0.6) +
      scale_fill_manual(values = c(`ABM main set` = "#BBBBBB", ABM_varyMisseg = "#FFFFFF", Salehi = "#E6C200"), name = NULL) +
      facet_wrap(~measure, nrow = 1, scales = "free_y") +
      labs(x = "missegregation rate (ABM, wavelength 0.8) or Salehi", y = NULL,
           caption = paste("Parent (last sampling time before day 120; Salehi: the last two passages) and daughter",
                           "(two sampling times later), 300 cells each (10 draws). Salehi: distinct pairs with at least 300 cells",
                           "(points; the lines differ widely).")) +
      theme_bw(base_size = 9) + theme(legend.position = "top", axis.text.x = element_text(size = 7))
    ggsave(file.path(fig_dir, "dispersion_misseg.pdf"), p, width = 11, height = 3.4)
  }
  cat("written:", fig_dir, "\n")
}
