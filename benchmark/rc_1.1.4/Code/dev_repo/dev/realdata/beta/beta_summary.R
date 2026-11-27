#!/usr/bin/env Rscript
# Report tables of beta2 (dev/realdata/beta/PLAN.md, "只报告"): where each
# estimator's picks fall against the truth, the accuracy curve over the grid
# (the oracle rate), the criteria curves, the Salehi and NUMBAT picks.
# Reads realdata/beta/evaluation.rds and the fit records; prints, and writes
# realdata/beta/summary.rds.
#   Rscript dev/realdata/beta/beta_summary.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "beta")
ev <- readRDS(file.path(root, "evaluation.rds"))
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
crit <- c(E1 = "objective", E2 = "data_nll", E3 = "cv_nll")
med <- function(x) stats::median(x, na.rm = TRUE)

sim <- ev$sim_rows
grid_rows <- sim[!sim$sim, ]
complete <- tapply(grid_rows$ok, grid_rows$dataset_id, function(x) sum(x) == length(grid))
g <- grid_rows[grid_rows$dataset_id %in% names(complete)[complete], ]
g$truth_label <- ifelse(g$kind == "c2", paste0("c2 x", round(g$true_multiple)), paste0(g$set, " x", g$true_multiple))
g$family <- ifelse(g$kind == "c2", sub("^r1_b[^_]+_", "", g$set), g$set)

# 1. picks (single version) against the truth, all landscapes, by set and truth
picks <- do.call(rbind, lapply(names(crit), function(e) {
  do.call(rbind, lapply(split(g, g$dataset_id), function(x) {
    v <- x[[crit[[e]]]]
    if (!any(is.finite(v))) return(NULL)
    data.frame(estimator = e, dataset_id = x$dataset_id[1], set = x$set[1], kind = x$kind[1], truth_label = x$truth_label[1],
               split = x$split[1], true_multiple = x$true_multiple[1], pick = x$multiple[which.min(v)],
               oracle_d1 = x$multiple[which.min(x$d1)])
  }))
}))
pick_table <- function(p) {
  t <- table(p$truth_label, factor(p$pick, levels = grid))
  cbind(as.data.frame.matrix(t), n = rowSums(t))
}
cat("== picks of each estimator (single version, all landscapes), rows = set and true multiple ==\n")
for (e in names(crit)) { cat("\n", e, crit[[e]], "\n"); print(pick_table(picks[picks$estimator == e, ])) }
cat("\n== oracle rate (lowest d1 RMSE over the grid) ==\n")
orc <- picks[picks$estimator == "E1", ]
orc$pick <- orc$oracle_d1
print(pick_table(orc))

# 2. accuracy curve: median over datasets of RMSE(m) / RMSE(truth or nearest grid rate)
ratio_curve <- do.call(rbind, lapply(split(g, g$truth_label), function(x) {
  do.call(rbind, lapply(c("d0", "d1", "d2"), function(d) {
    ref <- tapply(seq_len(nrow(x)), x$dataset_id, function(i) {
      y <- x[i, ]
      at <- y$multiple[which.min(abs(log(y$multiple / y$true_multiple[1])))]
      y[[d]] / y[[d]][y$multiple == at]
    }, simplify = FALSE)
    m <- do.call(rbind, ref)
    data.frame(truth_label = x$truth_label[1], stratum = d, t(apply(m, 2, med)), check.names = FALSE)
  }))
}))
names(ratio_curve)[-(1:2)] <- paste0("x", grid)
cat("\n== median RMSE(m) / RMSE(true rate; nearest grid rate for C2), by set ==\n")
print(ratio_curve, digits = 3, row.names = FALSE)

# 3. criteria curves: median over datasets of criterion(m) - min over m, by set and truth
crit_curve <- do.call(rbind, lapply(names(crit), function(e) {
  do.call(rbind, lapply(split(g, g$truth_label), function(x) {
    m <- do.call(rbind, tapply(seq_len(nrow(x)), x$dataset_id, function(i) {
      v <- x[[crit[[e]]]][i]; v - min(v, na.rm = TRUE)
    }, simplify = FALSE))
    data.frame(estimator = e, truth_label = x$truth_label[1], t(apply(m, 2, med)), check.names = FALSE)
  }))
}))
names(crit_curve)[-(1:2)] <- paste0("x", grid)
cat("\n== median (criterion(m) - its minimum over the grid), by set ==\n")
print(crit_curve, digits = 3, row.names = FALSE)

# 4. Salehi and NUMBAT picks, single version, by line / cancer
real_picks <- function(kind, by) {
  files <- list.files(file.path(root, "fits", kind), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  r <- do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    ok <- identical(x$status, "ok")
    data.frame(unit = x[[by]], id = if (kind == "salehi") x$id else paste(x$cancer, x$sample), multiple = x$multiple, ok = ok,
               objective = if (ok) x$criteria$objective else NA, data_nll = if (ok) x$criteria$data_nll else NA,
               cv_nll = if (ok) x$criteria$cv_nll else NA)
  }))
  full <- tapply(r$ok, r$id, function(x) sum(x) == length(grid))
  r <- r[r$id %in% names(full)[full], ]
  out <- do.call(rbind, lapply(names(crit), function(e) {
    p <- do.call(rbind, lapply(split(r, r$id), function(y) {
      v <- y[[crit[[e]]]]
      if (!any(is.finite(v))) return(NULL)
      data.frame(unit = y$unit[1], id = y$id[1], pick = y$multiple[which.min(v)])
    }))
    pooled <- do.call(rbind, lapply(split(r, r$unit), function(y) {
      tot <- tapply(y[[crit[[e]]]], y$multiple, sum)
      data.frame(unit = y$unit[1], pooled = as.numeric(names(tot))[which.min(tot)])
    }))
    s <- do.call(rbind, lapply(split(p, p$unit), function(y)
      data.frame(unit = y$unit[1], n = nrow(y), median_pick = med(y$pick),
                 picks = paste(names(table(y$pick)), table(y$pick), sep = ":", collapse = " "))))
    cbind(estimator = e, merge(s, pooled, by = "unit"))
  }))
  list(rows = r, picks = out)
}
salehi <- real_picks("salehi", "line")
numbat <- real_picks("numbat", "cancer")
cat("\n== Salehi picks by line (single: distribution; pooled: one rate per line) ==\n")
print(salehi$picks, row.names = FALSE)
cat("\n== NUMBAT picks by cancer ==\n")
print(numbat$picks, row.names = FALSE)

# 5. C2: RMSE at beta_sim against the grid's best and against x1
c2 <- sim[sim$kind == "c2", ]
c2_tab <- do.call(rbind, lapply(split(c2, c2$set), function(x) {
  do.call(rbind, lapply(c("d0", "d1", "d2"), function(d) {
    by_id <- split(x, x$dataset_id)
    r <- do.call(rbind, lapply(by_id, function(y) {
      if (!all(y$ok)) return(NULL)
      c(bsim = y[[d]][y$sim], x1 = y[[d]][!y$sim & y$multiple == 1], best = min(y[[d]][!y$sim]))
    }))
    data.frame(set = x$set[1], stratum = d, n = nrow(r), bsim_over_x1 = med(r[, "bsim"] / r[, "x1"]),
               best_over_x1 = med(r[, "best"] / r[, "x1"]))
  }))
}))
cat("\n== C2: median RMSE ratio at beta_sim and at the grid's best rate, against x1 ==\n")
print(c2_tab, digits = 3, row.names = FALSE)

# 6. horizon (E4 / S1): per-cell log score difference from x1, by line and rate
hz <- ev$horizon
hz_curve <- NULL
if (!is.null(hz) && nrow(hz)) {
  base <- hz[hz$multiple == 1, c("id", "uid", "ll_fc_s_e2")]
  m <- merge(hz, base, by = c("id", "uid"), suffixes = c("", "_x1"))
  m$diff <- m$ll_fc_s_e2 - m$ll_fc_s_e2_x1
  hz_curve <- stats::aggregate(diff ~ line + multiple, m, mean)
  hz_curve <- stats::reshape(hz_curve, idvar = "line", timevar = "multiple", direction = "wide")
  cat("\n== Salehi horizon: mean per-cell log score (ll_fc_s_e2) minus that at x1, by line and rate (<= 60 days) ==\n")
  print(hz_curve, digits = 3, row.names = FALSE)
}

dev_save_rds(list(picks = picks, ratio_curve = ratio_curve, crit_curve = crit_curve, salehi = salehi$picks,
                  numbat = numbat$picks, c2 = c2_tab, horizon_curve = hz_curve),
             file.path(root, "summary.rds"))
