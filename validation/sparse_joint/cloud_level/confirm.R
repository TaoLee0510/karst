#!/usr/bin/env Rscript
# Follow-up step 1: confirmation of the switching rule on new landscapes
# (sets dev2, short2, diverse2; README, "Follow-up"). The rule was fixed
# before these fits were made:
#   fit with balance_level = "free" and the history at the interval
#   (level_dt); if its fitted level exceeds 1, the fit is
#   initial_state = "free", otherwise the default (base).
# Criteria, as for the first round: on dev2 core and on short2, mean RMSE at
# distances 0, 1 and 2 at most 3% above base and direction accuracy at most
# 0.01 below; on diverse2 not worse by the same rules, per mechanism (`role`:
# rate, smooth, call), every one passing. Up to 2026-11-18 diverse2 was
# judged pooled, which is now only reported. Bootstrap over landscapes, 2000
# draws, seed 1.
#
#   Rscript validation/sparse_joint/cloud_level/confirm.R
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
out_dir <- file.path(script_dir, "results", "report_confirm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s <- utils::read.csv(file.path(script_dir, "results", "scores_confirm.csv"), stringsAsFactors = FALSE)
cat("Fits:\n")
print(table(s$set, s$method, s$status))
s <- s[s$status == "ok", ]
s$group <- ifelse(s$set == "dev2", paste("dev2", s$role), s$set)
threshold <- 1
level <- s[s$method == "level_dt", c("dataset_id", "balance_level")]
names(level)[2] <- "level"
base <- merge(s[s$method == "base", ], level, by = "dataset_id")
free <- s[s$method == "free", ]
free <- free[match(base$dataset_id, free$dataset_id), ]
keep <- !is.na(free$dataset_id) & is.finite(base$level)
base <- base[keep, ]
free <- free[keep, ]
switch_to_free <- base$level > threshold
rule <- base
columns <- c(paste0("rmse_d", 0:2), paste0("direction_d", 0:2), "var_pred", "change_pred", "spearman_d1")
columns <- intersect(columns, names(base))
rule[switch_to_free, columns] <- free[switch_to_free, columns]

boot <- function(index, landscapes, statistic, draws = 2000L) {
  groups <- split(index, landscapes)
  set.seed(1)
  values <- vapply(seq_len(draws), function(i) statistic(unlist(groups[sample.int(length(groups), replace = TRUE)],
                                                                use.names = FALSE)), numeric(1))
  stats::quantile(values, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
arms <- list(rule = rule, free = free)
rows <- list()
# The groups, and diverse2 by mechanism.
diverse <- which(base$set == "diverse2")
group_rows <- split(seq_len(nrow(base)), base$group)
if (length(diverse)) group_rows <- c(group_rows, split(diverse, paste("diverse2", base$role[diverse])))
for (g in names(group_rows)) {
  u <- group_rows[[g]]
  for (arm in names(arms)) {
    a <- arms[[arm]]
    for (d in 0:2) {
      r <- paste0("rmse_d", d)
      k <- paste0("direction_d", d)
      ratio <- function(i) mean(a[[r]][i]) / mean(base[[r]][i])
      diff <- function(i) mean(a[[k]][i], na.rm = TRUE) - mean(base[[k]][i], na.rm = TRUE)
      ci_r <- boot(u, base$landscape_id[u], ratio)
      ci_d <- boot(u, base$landscape_id[u], diff)
      rows[[length(rows) + 1L]] <- data.frame(group = g, arm = arm, distance = d, n = length(u),
        switched = mean(switch_to_free[u]), base_rmse = mean(base[[r]][u]), rmse_ratio = ratio(u),
        low = ci_r[1L], high = ci_r[2L], base_direction = mean(base[[k]][u], na.rm = TRUE), direction_diff = diff(u),
        direction_low = ci_d[1L], direction_high = ci_d[2L], stringsAsFactors = FALSE)
    }
  }
}
table <- do.call(rbind, rows)
utils::write.csv(table, file.path(out_dir, "rule.csv"), row.names = FALSE)
show <- table
num <- vapply(show, is.numeric, logical(1))
show[num] <- lapply(show[num], function(v) signif(v, 3))
cat("\n== Rule (switch to free when the fitted level > 1) and free, against base ==\n")
print(show, row.names = FALSE)
judge <- function(g, label = g) {
  t <- table[table$group == g & table$arm == "rule", ]
  if (!nrow(t)) return(NULL)
  data.frame(group = label, worst_rmse_ratio = max(t$rmse_ratio), worst_direction_diff = min(t$direction_diff),
             pass = max(t$rmse_ratio) <= 1.03 && min(t$direction_diff) >= -0.01)
}
mechanisms <- if (length(diverse)) sort(unique(paste("diverse2", base$role[diverse]))) else character()
per_mechanism <- do.call(rbind, lapply(mechanisms, judge))
criterion2 <- if (!is.null(per_mechanism)) {
  data.frame(group = "diverse2, every mechanism", worst_rmse_ratio = max(per_mechanism$worst_rmse_ratio),
             worst_direction_diff = min(per_mechanism$worst_direction_diff), pass = all(per_mechanism$pass))
}
verdict <- rbind(judge("dev2 core"), judge("short2"), criterion2, per_mechanism,
                 judge("diverse2", "diverse2 pooled (reported)"), judge("dev2 stress"))
utils::write.csv(verdict, file.path(out_dir, "verdict.csv"), row.names = FALSE)
cat("\n== Verdict (diverse2 pooled and dev2 stress reported, not criteria) ==\n")
print(verdict, row.names = FALSE)
cat("\nSwitched by group and level class:\n")
print(round(tapply(switch_to_free, paste(base$group, base$wavelength_class), mean), 2))
