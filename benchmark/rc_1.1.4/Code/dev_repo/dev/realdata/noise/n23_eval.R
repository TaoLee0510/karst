#!/usr/bin/env Rscript
# N2 and N3 of dev/realdata/noise/PLAN.md (corrections of 2026-11-20 11:40).
#   N2: r_hat against the true rate on the simulations (e1-e4, C2; "ends"),
#       monotonicity, false noise flags, and N_e on e3 against the true bottleneck.
#   N3: Salehi, r_hat (1 division per day) as the forecast rate, nearest rate of
#       the decomposition study's forecast grid (fit x1); gain over x1 with
#       trajectory-cluster bootstrap; parallel-lineage agreement of r_hat.
# Writes realdata/noise/evaluation.rds.
#   Rscript dev/realdata/noise/n23_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "noise")
rd <- function(kind) readRDS(file.path(root, paste0("n1_", kind, ".rds")))
med <- function(x) stats::median(x, na.rm = TRUE)

## ---------------- N2 ----------------
sim <- rbind(rd("eset"), rd("c2"))
ends <- sim[sim$design == "ends", ]
ends$true_multiple <- ends$true_beta / 5e-5
ends$within2 <- is.finite(ends$r_hat) & ends$r_hat > 0 & abs(log2(ends$r_hat / ends$true_beta)) <= 1
ends$label <- ifelse(ends$kind == "c2", paste0("C2 x", round(ends$true_multiple)), paste0(ends$set, " x", ends$true_multiple))
ends$n_cells <- ifelse(ends$kind == "eset", sub("^.*_n([0-9]{4})_k.*$", "\\1", ends$dataset_id), NA)
by_label <- do.call(rbind, lapply(split(ends, ends$label), function(x)
  data.frame(label = x$label[1], n = nrow(x), within2 = mean(x$within2), median_multiple = med(x$multiple),
             q25 = stats::quantile(x$multiple, 0.25, na.rm = TRUE, names = FALSE),
             q75 = stats::quantile(x$multiple, 0.75, na.rm = TRUE, names = FALSE),
             zero = mean(x$r_hat == 0), saturated = mean(!is.finite(x$r_hat)), flagged = mean(x$flagged),
             u = med(x$u), n_new = med(x$n_new))))
cat("== N2: r_hat on the simulations (ends), by set and true rate ==\n"); print(by_label, digits = 3, row.names = FALSE)
m <- function(lbl) by_label$median_multiple[by_label$label == lbl]
crit_within <- mean(ends$within2)
crit_mono <- isTRUE(m("e1 x1") < m("e4 x4") && m("e4 x4") < m("e4 x16") && m("C2 x3") < m("C2 x10") && m("C2 x10") < m("C2 x30"))
crit_flag <- mean(ends$flagged)
e3 <- ends[ends$set == "e3" & ends$replay & is.finite(ends$true_bottleneck), ]
e3$ne_within3 <- is.finite(e3$ne) & abs(log(e3$ne / e3$true_bottleneck)) <= log(3)
crit_ne <- mean(e3$ne_within3)
N2 <- c(within2 = crit_within >= 0.70, monotone = crit_mono, flags = crit_flag <= 0.10, drift = crit_ne >= 0.70)
cat(sprintf("\nN2: within 2x %.3f (>= 0.70); monotone %s; flagged %.3f (<= 0.10); e3 N_e within 3x %.3f (>= 0.70, n = %d)\n",
            crit_within, crit_mono, crit_flag, crit_ne, nrow(e3)))
print(N2)
cat("\nby cells per time point (e sets):\n")
print(stats::aggregate(cbind(within2, r_hat_multiple = multiple) ~ n_cells + set, ends[ends$kind == "eset", ],
                       function(v) round(stats::median(v, na.rm = TRUE), 3)))
cat("\ne3 drift: N_e / true bottleneck by bottleneck class:\n")
e3$b_class <- sub("^.*_b([0-9]{5})_.*$", "\\1", e3$dataset_id)
print(stats::aggregate(cbind(ratio = ne / true_bottleneck, within3 = ne_within3) ~ b_class, e3,
                       function(v) round(stats::median(v[is.finite(v)]), 3)))
adj <- sim[sim$design == "adjacent", ]
adj$within2 <- is.finite(adj$r_hat) & adj$r_hat > 0 & abs(log2(adj$r_hat / adj$true_beta)) <= 1
cat("\nadjacent intervals (report only): within 2x by set:\n")
print(stats::aggregate(cbind(within2, multiple) ~ set, adj, function(v) round(stats::median(v[is.finite(v)]), 3)))

## ---------------- N3 ----------------
sal <- rd("salehi")
salehi_a <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)
xlab <- function(m) ifelse(m < 1, paste0("x1_", round(1 / m)), paste0("x", m))
nearest <- function(mult) ifelse(!is.finite(mult), max(grid), ifelse(mult <= 0, min(grid),
                                 grid[vapply(mult, function(v) which.min(abs(log(grid) - log(v))), integer(1))]))
fc_dir <- file.path(dev_results, "realdata", "beta", "decomp")
read_fc <- function(f) {
  files <- list.files(file.path(fc_dir, paste0("F_", xlab(f))), pattern = "__x1__m10[.]Rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(p) { x <- readRDS(p); if (identical(x$status, "ok")) cbind(x$rows[c("id", "uid", "days", "ll_fc_s_e2")], F = f) else NULL }))
}
fc <- do.call(rbind, lapply(grid, read_fc))
fc <- fc[fc$days <= 60, ]
base <- fc[fc$F == 1, ]
fc$gain <- fc$ll_fc_s_e2 - base$ll_fc_s_e2[match(paste(fc$id, fc$uid), paste(base$id, base$uid))]
fc$line <- salehi_a$line[match(fc$id, salehi_a$id)]
fc$group <- ifelse(fc$line %in% c("SA609", "SA039"), "dev", "confirm")
set.seed(20261120)
boot_ci <- function(g, ids, reps = 2000L) {
  u <- unique(ids); by_id <- split(g, ids)
  b <- replicate(reps, { s <- sample(u, length(u), replace = TRUE); mean(unlist(by_id[s], use.names = FALSE)) })
  stats::quantile(b, c(0.025, 0.975), names = FALSE)
}
n3 <- function(rho) {
  s <- sal[sal$rho == rho, ]
  pick <- stats::setNames(nearest(s$multiple), s$dataset_id)
  x <- fc[fc$id %in% names(pick), ]
  x <- x[x$F == pick[x$id], ]
  est <- do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
    y <- x[x$group == gr, ]; ci <- boot_ci(y$gain, y$id)
    data.frame(rho = rho, group = gr, passages = nrow(y), trajectories = length(unique(y$id)), gain = mean(y$gain), lo = ci[1], hi = ci[2])
  }))
  picks <- stats::aggregate(F ~ line, x[!duplicated(x$id), ], function(v) paste(names(table(v)), table(v), sep = ":", collapse = " "))
  list(est = est, picks = picks)
}
p1 <- n3(1); p2 <- n3(0.141)
cat("\n== N3: r_hat (nearest forecast rate) as the forecast rate; gain over x1 per later passage ==\n")
print(rbind(p1$est, p2$est), digits = 3, row.names = FALSE)
cat("\npicks by line (1 division per day):\n"); print(p1$picks, row.names = FALSE)
cat("\npicks by line (0.141 per day):\n"); print(p2$picks, row.names = FALSE)
# parallel lineages (as B and beta_eval.R's S2), r_hat at 1 division per day
s1 <- sal[sal$rho == 1, ]
sa <- salehi_a[salehi_a$id %in% s1$dataset_id, ]
sa$r <- s1$r_hat[match(sa$id, s1$dataset_id)]
keys <- split(sa, list(sa$line, sa$treat_first, sa$treat_last), drop = TRUE)
agree <- unlist(lapply(keys, function(g) {
  if (nrow(g) < 2) return(NULL)
  cmb <- utils::combn(nrow(g), 2)
  vapply(seq_len(ncol(cmb)), function(j) {
    a <- g[cmb[1, j], ]; b <- g[cmb[2, j], ]
    if (length(intersect(c(a$first_uid, a$last_uid), c(b$first_uid, b$last_uid)))) return(NA)
    if (a$r == 0 && b$r == 0) return(TRUE)
    if (a$r == 0 || b$r == 0 || !is.finite(a$r) || !is.finite(b$r)) return(FALSE)
    abs(log2(a$r / b$r)) <= 1
  }, logical(1))
}))
pl <- mean(agree, na.rm = TRUE)
e1 <- p1$est
N3 <- c(confirm = e1$gain[e1$group == "confirm"] >= 0.5 * 0.4531 && e1$lo[e1$group == "confirm"] > 0,
        dev = e1$gain[e1$group == "dev"] >= 0, lineages = pl >= 0.70)
cat(sprintf("\nparallel lineages within 2x: %.3f (%d pairs)\n", pl, sum(!is.na(agree))))
print(N3)
cat("\n== Salehi r_hat by line (1 / day), u, n_new, N_e (report) ==\n")
s1$line <- s1$set
print(stats::aggregate(cbind(multiple, u, n_new, new_d1, new_d2, new_d3plus, ne, flagged) ~ line, s1,
                       function(v) round(stats::median(v[is.finite(v)]), 3)))
nb <- rd("numbat")
cat("\n== NUMBAT (0.141 / day; report) ==\n")
print(stats::aggregate(cbind(multiple, u, n_new, flagged) ~ set, nb, function(v) round(stats::median(v[is.finite(v)]), 3)))
dev_save_rds(list(ends = ends, by_label = by_label, N2 = N2, N2_values = c(within2 = crit_within, flagged = crit_flag, ne = crit_ne),
                  e3 = e3, N3 = N3, n3_rho1 = p1, n3_rho0141 = p2, lineages = pl, salehi = sal, numbat = nb),
             file.path(root, "evaluation.rds"))
