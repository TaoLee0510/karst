# Step 1 of PLAN.md: diagnosis on the simulations, from the step-1 arms only
# (truth, truth_true_beta, current_native, current_aligned; A1's rel_current).
# The branch is decided on the development half (odd replicates).
#   Rscript scripts/step1_analyse.R OUT_DIR TASK_FILE RESULT_DIR
args <- commandArgs(TRUE)
out_dir <- args[[1]]; tasks <- utils::read.delim(args[[2]], stringsAsFactors = FALSE); res_dir <- args[[3]]
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(out_dir, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
step1_arms <- c("truth", "truth_true_beta", "current_native", "current_aligned")
rows <- list(); a1 <- list(); status <- list()
for (f in files) {
  r <- readRDS(f)
  t <- r$task
  key <- data.frame(source = t$source, arm_set = t$arm, sim_id = t$sim_id, design = t$design, stringsAsFactors = FALSE)
  status[[length(status) + 1L]] <- cbind(key, status = r$status, same_input = if (is.null(r$same_input)) NA else r$same_input)
  if (!identical(r$status, "ok")) next
  info <- r$info
  half <- tasks$half[tasks$source == t$source & tasks$arm == t$arm & tasks$sim_id == t$sim_id & tasks$design == t$design][1]
  ob <- r$observed_shares
  fc <- r$forecasts[r$forecasts$arm %in% step1_arms, ]
  for (i in seq_len(nrow(fc))) {
    row <- cbind(key, half = half, w = info$w, pmis = info$pmis, beta_fit = r$beta_fit, horizon = r$horizon, arm = fc$arm[i])
    for (g in c("inside", "halo", "island")) for (rr in c(3, 4, 6)) {
      col <- paste0(g, "_r", rr)
      row[[paste0("obs_", col)]] <- ob["daughter", col] - ob["parent", col]
      row[[paste0("pred_", col)]] <- fc[[col]][i] - ob["parent", col]
      row[[paste0("parent_", col)]] <- ob["parent", col]
    }
    row$angle <- fc$angle[i]; row$overlap <- fc$overlap[i]; row$w1 <- fc$w1[i]; row$parent_rate_mean <- fc$parent_rate_mean[i]
    rows[[length(rows) + 1L]] <- row
  }
  a <- r$a1
  if (nrow(a)) a1[[length(a1) + 1L]] <- cbind(key[rep(1, nrow(a)), ], half = half, w = info$w, pmis = info$pmis,
                                              a[, c("kind", "dmajor", "truth_rel", "rel_current")], row.names = NULL)
}
st <- do.call(rbind, status); d <- do.call(rbind, rows); a1 <- do.call(rbind, a1)
utils::write.table(st, file.path(res_dir, "status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(d, file.path(res_dir, "forecast_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 220)
cat("inputs:", nrow(st), " ok:", sum(st$status == "ok"), " input check failed:", sum(st$same_input %in% FALSE), "\n")
print(table(st$status, paste(st$source, st$arm_set)))

## branch, on the development half --------------------------------------------------
wide <- function(x, col) {
  w <- reshape(x[, c("source", "arm_set", "sim_id", "design", "half", "w", "pmis", "arm", paste0("obs_", col), paste0("pred_", col))],
               idvar = c("source", "arm_set", "sim_id", "design", "half", "w", "pmis", paste0("obs_", col)), timevar = "arm", direction = "wide")
  names(w) <- sub(paste0("^pred_", col, "[.]"), "", names(w))
  names(w)[names(w) == paste0("obs_", col)] <- "obs"
  w
}
branch <- function(x, col = "halo_r4") {
  w <- wide(x, col)
  up <- w$obs > 0
  data.frame(inputs = nrow(w), observed_up = sum(up),
             karst_native_down_given_up = mean(w$current_native[up] < 0), karst_aligned_down_given_up = mean(w$current_aligned[up] < 0),
             truth_up_given_up = mean(w$truth[up] > 0),
             truth_true_beta_up_given_up = if ("truth_true_beta" %in% names(w)) mean(w$truth_true_beta[up & !is.na(w$truth_true_beta)] > 0) else NA,
             median_obs = stats::median(w$obs), median_karst_native = stats::median(w$current_native), median_truth = stats::median(w$truth))
}
dev <- d[d$half == "dev", ]
b <- branch(dev)
cat("\n== branch (development half, halo = out of the domain, <= 4 steps from a major)\n"); print(b, row.names = FALSE, digits = 3)
verdict <- if (b$karst_native_down_given_up < 0.5) "NOT REPRODUCED: stop" else
  if (b$truth_up_given_up >= 0.7) "FITNESS: go to step 2" else "FLUX / LIST: stop, the forecast rate is the user's decision"
cat("verdict:", verdict, "\n")

## descriptive -------------------------------------------------------------------------
grp <- function(x) ifelse(x$source == "abm", paste0("main w", x$w), paste0("vM ", format(x$pmis, scientific = TRUE), " ", x$arm_set))
for (col in c("halo_r4", "inside_r4", "island_r4", "halo_r3", "halo_r6")) {
  cat("\n==", col, ": all inputs, by group (both halves)\n")
  dd <- d; dd$group <- grp(dd)
  out <- do.call(rbind, lapply(split(dd, dd$group), function(x) cbind(group = x$group[1], branch(x, col))))
  print(out, row.names = FALSE, digits = 3)
}
cat("\n== Fig. 2c-like angle and W1 to the observed daughter (median), by group and arm\n")
dd <- d; dd$group <- grp(dd)
print(reshape(stats::aggregate(cbind(angle, w1) ~ group + arm, dd, stats::median), idvar = "group", timevar = "arm", direction = "wide"), row.names = FALSE, digits = 3)
cat("\n== parent mean division rate (truth vs KARST native), median by group\n")
print(reshape(stats::aggregate(parent_rate_mean ~ group + arm, dd[dd$arm %in% c("truth", "current_native"), ], stats::median),
              idvar = "group", timevar = "arm", direction = "wide"), row.names = FALSE, digits = 3)

## A1 ------------------------------------------------------------------------------------
if (!is.null(a1)) {
  a1$band <- cut(a1$dmajor, c(0, 2, 4, Inf), labels = c("2", "3-4", ">4"))
  a1$bias <- a1$rel_current - a1$truth_rel
  a1$group <- grp(a1)
  utils::write.table(a1, file.path(res_dir, "a1_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n== A1: KARST relative fitness minus true relative fitness, out of the domain (both halves)\n")
  s <- do.call(rbind, lapply(split(a1, list(a1$kind, a1$band), drop = TRUE), function(x) data.frame(kind = x$kind[1], band = x$band[1], n = nrow(x),
    bias_median = stats::median(x$bias), truth_median = stats::median(x$truth_rel), karst_median = stats::median(x$rel_current),
    spearman = suppressWarnings(stats::cor(x$rel_current, x$truth_rel, method = "spearman")))))
  print(s, row.names = FALSE, digits = 3)
  cat("by group, observed out-of-domain karyotypes:\n")
  o <- a1[a1$kind == "observed", ]
  print(do.call(rbind, lapply(split(o, o$group), function(x) data.frame(group = x$group[1], n = nrow(x), bias_median = stats::median(x$bias),
    truth_median = stats::median(x$truth_rel), karst_median = stats::median(x$rel_current),
    spearman = suppressWarnings(stats::cor(x$rel_current, x$truth_rel, method = "spearman"))))), row.names = FALSE, digits = 3)
}
writeLines(c(verdict), file.path(res_dir, "verdict.txt"))
