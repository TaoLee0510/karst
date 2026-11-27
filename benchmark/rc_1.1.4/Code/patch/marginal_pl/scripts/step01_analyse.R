# Steps 0 and 1 of PLAN.md (development half only).
#   Rscript scripts/step01_analyse.R OUT_DIR SIM_TASKS RES_DIR
args <- commandArgs(TRUE); out_dir <- args[[1]]; tasks <- utils::read.delim(args[[2]], stringsAsFactors = FALSE); res_dir <- args[[3]]
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
ws <- sprintf("w%.4f", c(1/22, 1/8, 1/3, 1))
rows <- list()
add <- function(x, family, unit, group, id, n_set) {
  if (!identical(x$status, "ok")) return(NULL)
  r <- data.frame(family = family, unit = unit, group = group, id = id, n_set = n_set,
                  pl_truth = x$pl[["truth"]], pl_karst = x$pl[["karst"]], pl_neutral = x$pl[["neutral"]], stringsAsFactors = FALSE)
  for (a in c("current", ws)) {
    for (m in c("sp_ou", "sp_un", "sp_in", "rmse_in")) r[[paste0(a, ".", m)]] <- x$arms[[a]]$acc[[m]]
    if (a != "current") { r[[paste0(a, ".b_sp_ou")]] <- x$arms[[a]]$acc_b[["sp_ou"]]; r[[paste0(a, ".conv")]] <- x$fits[[a]]$convergence }
  }
  r
}
for (f in list.files(file.path(out_dir, "own_dev"), "[.]rds$", full.names = TRUE)) {
  x <- readRDS(f); rows[[length(rows) + 1L]] <- add(x, "own", x$landscape_id, x$wavelength_class, x$task$dataset_id, if (is.null(x$n_sets)) NA else x$n_sets[["ou"]])
}
for (f in list.files(out_dir, "[.]rds$", recursive = TRUE, full.names = TRUE)) {
  if (!grepl("/paper_", f)) next
  x <- readRDS(f); t <- x$task
  if (t$half != "dev") next
  grp <- if (t$source == "abm") paste0("main w", x$info$w) else paste0("vM ", format(x$info$pmis, scientific = TRUE), " ", t$arm)
  ns <- if (is.null(x$n_sets)) NA else sum(x$n_sets[c("O", "unseen")], na.rm = TRUE)
  rows[[length(rows) + 1L]] <- add(x, "paper", paste(t$source, t$sim_id), grp, paste(t$source, t$arm, t$sim_id, t$design), ns)
}
d <- do.call(rbind, rows)
utils::write.table(d, file.path(res_dir, "dev_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("development datasets: own", sum(d$family == "own"), " paper", sum(d$family == "paper"), "\n")

## step 0 ---------------------------------------------------------------------------------------
d$tn <- d$pl_truth - d$pl_neutral; d$tk <- d$pl_truth - d$pl_karst; d$kn <- d$pl_karst - d$pl_neutral
s0 <- do.call(rbind, lapply(c("own", "paper"), function(fm) { x <- d[d$family == fm, ]
  data.frame(family = fm, datasets = nrow(x), truth_minus_neutral = stats::median(x$tn), truth_minus_karst = stats::median(x$tk),
             karst_minus_neutral = stats::median(x$kn), share_truth_above_karst = mean(x$tk > 0),
             S0a = stats::median(x$tn) > 0, S0b = stats::median(x$tk) > 0 && stats::median(x$tk) >= 0.2 * stats::median(x$tn)) }))
cat("\n== step 0: PL per t1 cell (medians)\n"); print(s0, row.names = FALSE, digits = 4)
cat("\n== step 0 by group (medians)\n"); print(stats::aggregate(cbind(tn, tk, kn) ~ family + group, d, stats::median), row.names = FALSE, digits = 3)
utils::write.table(s0, file.path(res_dir, "step0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(s0$S0a & s0$S0b)) {
  v <- sprintf("STOP at step 0: S0a %s, S0b %s", paste(s0$family, s0$S0a, collapse = " / "), paste(s0$family, s0$S0b, collapse = " / "))
  writeLines(v, file.path(res_dir, "verdict.txt")); cat("\n", v, "\n"); quit(save = "no")
}

## step 1 ---------------------------------------------------------------------------------------
cw <- function(x, v) mean(tapply(v, x$unit, mean, na.rm = TRUE), na.rm = TRUE)
e <- d[is.finite(d$n_set) & d$n_set >= 5, ]
g <- do.call(rbind, lapply(ws, function(w) do.call(rbind, lapply(c("own", "paper"), function(fm) { x <- e[e$family == fm, ]
  ok <- is.finite(x[[paste0(w, ".rmse_in")]]) & is.finite(x$current.rmse_in) & x$current.rmse_in > 0
  data.frame(w = w, family = fm, datasets = nrow(x), gain_sp_ou = cw(x, x[[paste0(w, ".sp_ou")]] - x$current.sp_ou),
             gain_sp_un = cw(x, x[[paste0(w, ".sp_un")]] - x$current.sp_un), gain_sp_in = cw(x, x[[paste0(w, ".sp_in")]] - x$current.sp_in),
             rmse_in_ratio = cw(x[ok, ], x[[paste0(w, ".rmse_in")]][ok] / x$current.rmse_in[ok]),
             gain_b_sp_ou = cw(x, x[[paste0(w, ".b_sp_ou")]] - x$current.sp_ou), not_converged = sum(x[[paste0(w, ".conv")]] != 0, na.rm = TRUE)) }))))
cat("\n== step 1: cluster-weighted gains over current (datasets with >= 5 in O + unseen)\n"); print(g, row.names = FALSE, digits = 4)
utils::write.table(g, file.path(res_dir, "step1_grid.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
mins <- tapply(g$gain_sp_ou, g$w, min)
best <- names(mins)[mins >= max(mins) - 0.005]; best <- best[which.min(as.numeric(sub("w", "", best)))]
gb <- g[g$w == best, ]
pass <- all(gb$gain_sp_ou >= 0.05) && all(gb$rmse_in_ratio <= 1.02)
cat("\nchosen w:", best, "\n")
cat("\n== by group at the chosen w (mean per dataset)\n")
e$gain <- e[[paste0(best, ".sp_ou")]] - e$current.sp_ou
print(stats::aggregate(cbind(current.sp_ou, gain) ~ family + group, e, mean), row.names = FALSE, digits = 3)
v <- if (pass) sprintf("PASS step 1: w = %s (own gain %.3f, paper gain %.3f; rmse_in ratio %.3f / %.3f)", best, gb$gain_sp_ou[1], gb$gain_sp_ou[2], gb$rmse_in_ratio[1], gb$rmse_in_ratio[2]) else
  sprintf("STOP at step 1: w = %s (own gain %.3f, paper gain %.3f; rmse_in ratio %.3f / %.3f)", best, gb$gain_sp_ou[1], gb$gain_sp_ou[2], gb$rmse_in_ratio[1], gb$rmse_in_ratio[2])
writeLines(c(v, paste("w", best)), file.path(res_dir, "verdict.txt")); cat("\n", v, "\n")
