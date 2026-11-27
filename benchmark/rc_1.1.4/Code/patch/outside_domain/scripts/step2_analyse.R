# Step 2 of PLAN.md: the candidate rules (B1, B2, B1B2) against KARST as it
# is, judged on the confirmation half (even replicates); the development half
# is reported beside it. Run only after step 1's branch is "FITNESS".
#   Rscript scripts/step2_analyse.R OUT_DIR TASK_FILE RESULT_DIR
args <- commandArgs(TRUE)
out_dir <- args[[1]]; tasks <- utils::read.delim(args[[2]], stringsAsFactors = FALSE); res_dir <- args[[3]]
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
variants <- c("B1", "B2", "B1B2")
files <- list.files(out_dir, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
fc <- list(); a1 <- list()
for (f in files) {
  r <- readRDS(f)
  if (!identical(r$status, "ok")) next
  t <- r$task
  half <- tasks$half[tasks$source == t$source & tasks$arm == t$arm & tasks$sim_id == t$sim_id & tasks$design == t$design][1]
  key <- data.frame(source = t$source, arm_set = t$arm, sim_id = t$sim_id, design = t$design, half = half, w = r$info$w, pmis = r$info$pmis)
  ob <- r$observed_shares
  x <- r$forecasts
  x$obs_halo_change <- ob["daughter", "halo_r4"] - ob["parent", "halo_r4"]
  x$pred_halo_change <- x$halo_r4 - ob["parent", "halo_r4"]
  fc[[length(fc) + 1L]] <- cbind(key[rep(1, nrow(x)), ], x[, c("arm", "obs_halo_change", "pred_halo_change", "angle", "overlap", "w1")], row.names = NULL)
  a <- r$a1[r$a1$kind == "observed", ]
  if (nrow(a)) a1[[length(a1) + 1L]] <- cbind(key[rep(1, nrow(a)), ], a[, c("dmajor", "truth_rel", "rel_current", paste0("rel_", variants))], row.names = NULL)
}
fc <- do.call(rbind, fc); a1 <- do.call(rbind, a1)
options(width = 220)
judge <- function(h) {
  F <- fc[fc$half == h, ]; A <- a1[a1$half == h, ]
  do.call(rbind, lapply(variants, function(v) {
    # C1: bias of the observed out-of-domain karyotypes, and their rank correlation with the truth
    b_cur <- stats::median(A$rel_current - A$truth_rel); b_v <- stats::median(A[[paste0("rel_", v)]] - A$truth_rel)
    rho_cur <- suppressWarnings(stats::cor(A$rel_current, A$truth_rel, method = "spearman"))
    rho_v <- suppressWarnings(stats::cor(A[[paste0("rel_", v)]], A$truth_rel, method = "spearman"))
    # C2: level-aligned forecasts, paired per input
    cur <- F[F$arm == "current_aligned", ]; var <- F[F$arm == paste0(v, "_aligned"), ]
    m <- merge(cur, var, by = c("source", "arm_set", "sim_id", "design"), suffixes = c("_cur", "_var"))
    err_cur <- abs(m$pred_halo_change_cur - m$obs_halo_change_cur); err_var <- abs(m$pred_halo_change_var - m$obs_halo_change_cur)
    ok <- err_cur > 1e-9
    halo_ratio <- stats::median(err_var[ok] / err_cur[ok])
    w1_ratio <- stats::median(m$w1_var / m$w1_cur, na.rm = TRUE)
    data.frame(half = h, variant = v, a1_n = nrow(A), bias_current = b_cur, bias_variant = b_v, spearman_current = rho_cur, spearman_variant = rho_v,
               C1 = abs(b_v) <= abs(b_cur) / 2 && rho_v >= rho_cur - 0.02, inputs = nrow(m), halo_error_ratio = halo_ratio, w1_ratio = w1_ratio,
               C2 = halo_ratio <= 0.8 && w1_ratio <= 1.00)
  }))
}
res <- rbind(judge("dev"), judge("confirm"))
utils::write.table(res, file.path(res_dir, "step2_judgement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("== C1 and C2 (C3, the own benchmark, is judged separately)\n"); print(res, row.names = FALSE, digits = 3)
grp <- function(x) ifelse(x$source == "abm", paste0("main w", x$w), paste0("vM ", format(x$pmis, scientific = TRUE), " ", x$arm_set))
fc$group <- grp(fc)
cat("\n== confirmation half: median angle / W1 / halo-change error by group and arm (aligned arms)\n")
cf <- fc[fc$half == "confirm" & grepl("aligned$|^truth", fc$arm), ]
cf$halo_err <- abs(cf$pred_halo_change - cf$obs_halo_change)
print(reshape(stats::aggregate(cbind(angle, w1, halo_err) ~ group + arm, cf, stats::median), idvar = "group", timevar = "arm", direction = "wide"),
      row.names = FALSE, digits = 3)
