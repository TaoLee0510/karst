#!/usr/bin/env Rscript
# decomp/PLAN.md: the three leave-one-trajectory-out picks J (L = F), L (F = x1)
# and F (L = x1) on the Salehi later passages, their gains over (x1, x1) with
# trajectory-cluster bootstrap intervals, the pre-registered verdict, and the
# report-only tables. Writes realdata/beta/decomp/evaluation.rds.
#   Rscript dev/realdata/beta/decomp/decomp_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "beta")
# PLAN correction 07:53: the main analysis uses x1/4..x128 (28 combinations); x256 is report only
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
salehi_a <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
group_of <- function(line) ifelse(line %in% c("SA609", "SA039"), "dev", "confirm")

## ---------- read ----------
files <- list.files(file.path(root, "decomp"), pattern = "[.]Rds$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("/test/", files)]
recs <- lapply(files, readRDS)
status <- vapply(recs, function(x) x$status, character(1))
cat("forecast tasks:", length(recs), "\n"); print(table(status))
rows <- do.call(rbind, lapply(recs[status == "ok"], function(x) {
  r <- x$rows
  r$L <- x$task$L; r$F <- x$task$F
  r[c("id", "uid", "days", "L", "F", "ll_fc_s_e2", "ll_p0_s_e2", "support_size")]
}))
rows <- rows[rows$days <= 60 & is.finite(rows$ll_fc_s_e2), ]
rows_all <- rows   # with x256, for the report-only table
rows <- rows[rows$L <= max(grid) & rows$F <= max(grid), ]
rows$line <- salehi_a$line[match(rows$id, salehi_a$id)]
rows$group <- group_of(rows$line)
rows$key <- paste(rows$id, rows$uid)
rows$combo <- paste(rows$L, rows$F)
n_combo <- tapply(rows$combo, rows$key, function(x) length(unique(x)))
full <- names(n_combo)[n_combo == 28L]
cat("later passages with all 28 combinations:", length(full), "of", length(n_combo), "\n")
support_check <- tapply(rows$support_size, rows$key, function(x) length(unique(x)))
cat("passages whose support differs between combinations:", sum(support_check > 1), "\n")
rows <- rows[rows$key %in% full, ]
base <- rows[rows$L == 1 & rows$F == 1, c("key", "ll_fc_s_e2")]
rows$diff <- rows$ll_fc_s_e2 - base$ll_fc_s_e2[match(rows$key, base$key)]

sweeps <- list(J = function(k) rows$L == k & rows$F == k,
               L = function(k) rows$L == k & rows$F == 1,
               F = function(k) rows$L == 1 & rows$F == k)
# diff of every passage on one sweep, as a passage x rate matrix
sweep_matrix <- function(s) {
  m <- sapply(grid, function(k) { x <- rows[sweeps[[s]](k), ]; x$diff[match(full_keys, x$key)] })
  dimnames(m) <- list(full_keys, paste0("x", grid)); m
}
full_keys <- sort(unique(rows$key))
info <- unique(rows[c("key", "id", "line", "group")])
info <- info[match(full_keys, info$key), ]
M <- lapply(names(sweeps), sweep_matrix); names(M) <- names(sweeps)

## ---------- report only: curves ----------
curves <- do.call(rbind, lapply(names(M), function(s) do.call(rbind, lapply(split(seq_along(full_keys), info$line), function(i) {
  rbind(data.frame(sweep = s, line = info$line[i[1]], stat = "mean", n = length(i), t(colMeans(M[[s]][i, , drop = FALSE]))),
        data.frame(sweep = s, line = info$line[i[1]], stat = "median", n = length(i),
                   t(apply(M[[s]][i, , drop = FALSE], 2, stats::median))))
}))))
cat("\n== gain over (x1, x1) per later passage, by sweep, line and rate ==\n")
print(curves, digits = 2, row.names = FALSE)

## ---------- leave-one-trajectory-out picks ----------
# the pick maximises the mean over the other trajectories' passages (as E4 in beta_eval.R: mean over rows)
loo_pick_rows <- function(s) {
  ids <- unique(info$id)
  stats::setNames(vapply(ids, function(i) {
    line <- info$line[info$id == i][1]
    others <- info$line == line & info$id != i
    if (!any(others)) return(NA_real_)
    grid[which.max(colMeans(M[[s]][others, , drop = FALSE]))]
  }, numeric(1)), ids)
}
picks <- lapply(names(M), loo_pick_rows); names(picks) <- names(M)
gain_rows <- function(s, k_of_id) {
  k <- k_of_id[info$id]
  M[[s]][cbind(seq_along(full_keys), match(k, grid))]
}
set.seed(20261120)
boot_ci <- function(g, ids, reps = 2000L) {
  u <- unique(ids)
  by_id <- split(g, ids)
  b <- replicate(reps, { s <- sample(u, length(u), replace = TRUE); mean(unlist(by_id[s], use.names = FALSE)) })
  stats::quantile(b, c(0.025, 0.975), names = FALSE)
}
est <- do.call(rbind, lapply(names(M), function(s) do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
  g <- gain_rows(s, picks[[s]])
  i <- info$group == gr & is.finite(g)
  ci <- boot_ci(g[i], info$id[i])
  data.frame(estimator = s, group = gr, passages = sum(i), trajectories = length(unique(info$id[i])),
             gain = mean(g[i]), lo = ci[1], hi = ci[2])
}))))
cat("\n== leave-one-trajectory-out picks: gain over (x1, x1), per later passage, 95% trajectory bootstrap ==\n")
print(est, digits = 3, row.names = FALSE)
pick_table <- do.call(rbind, lapply(names(picks), function(s) {
  p <- data.frame(id = names(picks[[s]]), pick = picks[[s]])
  p$line <- salehi_a$line[match(p$id, salehi_a$id)]
  stats::aggregate(pick ~ line, p, function(x) paste(names(table(x)), table(x), sep = ":", collapse = " ")) |>
    transform(estimator = s)
}))
cat("\n== picks by line ==\n"); print(pick_table, row.names = FALSE)

## ---------- verdict (PLAN) ----------
e <- function(s, gr, col) est[est$estimator == s & est$group == gr, col]
A <- e("L", "confirm", "gain") > 0 && e("L", "confirm", "lo") > 0 && e("L", "dev", "gain") >= 0
B <- !A && e("F", "confirm", "gain") >= 0.75 * e("J", "confirm", "gain")
C <- A && e("F", "confirm", "gain") > 0 && e("F", "confirm", "lo") > 0
verdict <- if (C) "丙 两者都有" else if (A) "甲 拟合的 β 本身有贡献" else if (B) "乙 主要是预测更分散" else "丁 不确定"
cat("\nA (fitting rate itself):", A, " B (forecast spread):", B, " C (both):", C, "\nVERDICT:", verdict, "\n")

## ---------- report only: interaction, edge picks ----------
inter <- do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
  k <- picks$J[info$id]
  j <- match(k, grid)
  idx <- cbind(seq_along(full_keys), j)
  x <- M$J[idx] - M$L[idx] - M$F[idx]
  i <- info$group == gr & is.finite(x)
  data.frame(group = gr, J = mean(M$J[idx][i]), L_same_k = mean(M$L[idx][i]), F_same_k = mean(M$F[idx][i]),
             interaction = mean(x[i]))
}))
cat("\n== at J's pick: J gain, L and F gains at the same rate, interaction ==\n"); print(inter, digits = 3, row.names = FALSE)
edge <- vapply(picks, function(p) mean(p == max(grid), na.rm = TRUE), numeric(1))
cat("\nshare of trajectories picking the grid's top (x128):\n"); print(round(edge, 3))

## ---------- report only: x256 ----------
ra <- rows_all
ra$line <- salehi_a$line[match(ra$id, salehi_a$id)]
ra$key <- paste(ra$id, ra$uid)
b0 <- ra[ra$L == 1 & ra$F == 1, c("key", "ll_fc_s_e2")]
ra$diff <- ra$ll_fc_s_e2 - b0$ll_fc_s_e2[match(ra$key, b0$key)]
combos256 <- list("L=F=x256" = c(256, 256), "L=x256,F=x1" = c(256, 1), "L=x1,F=x256" = c(1, 256), "L=F=x128" = c(128, 128))
x256 <- do.call(rbind, lapply(names(combos256), function(nm) {
  lf <- combos256[[nm]]
  x <- ra[ra$L == lf[1] & ra$F == lf[2] & is.finite(ra$diff), ]
  do.call(rbind, lapply(split(x, x$line), function(y)
    data.frame(combo = nm, line = y$line[1], trajectories = length(unique(y$id)), passages = nrow(y), mean_gain = mean(y$diff))))
}))
cat("\n== x256 (report only): gain over (x1, x1) where available ==\n"); print(x256, digits = 3, row.names = FALSE)

## ---------- report only: E1-E3 on the 10-rate grid (horizon trajectories) ----------
fit_files <- list.files(file.path(root, "fits", "salehi"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
fr <- do.call(rbind, lapply(fit_files, function(f) {
  x <- readRDS(f); ok <- identical(x$status, "ok")
  data.frame(id = x$id, line = x$line, multiple = x$multiple, ok = ok,
             E1 = if (ok) x$criteria$objective else NA, E2 = if (ok) x$criteria$data_nll else NA,
             E3 = if (ok) x$criteria$cv_nll else NA)
}))
fr <- fr[fr$id %in% unique(info$id) & fr$multiple %in% grid, ]
complete <- tapply(fr$ok, fr$id, function(x) sum(x) == length(grid))
fr <- fr[fr$id %in% names(complete)[complete], ]
crit11 <- do.call(rbind, lapply(c("E1", "E2", "E3"), function(cn) do.call(rbind, lapply(split(fr, fr$line), function(y) {
  tot <- tapply(y[[cn]], y$multiple, sum)
  single <- vapply(split(y, y$id), function(z) z$multiple[which.min(z[[cn]])], numeric(1))
  data.frame(criterion = cn, line = y$line[1], n = length(single), pooled = as.numeric(names(tot))[which.min(tot)],
             single_median = stats::median(single), at_top = sum(single == max(grid)))
}))))
cat("\n== E1-E3 on the 10-rate grid (report only) ==\n"); print(crit11, row.names = FALSE)

dev_save_rds(list(rows = rows, curves = curves, picks = picks, estimates = est, verdict = verdict,
                  flags = c(A = A, B = B, C = C), interaction = inter, edge = edge, criteria = crit11, x256 = x256),
             file.path(root, "decomp", "evaluation.rds"))
