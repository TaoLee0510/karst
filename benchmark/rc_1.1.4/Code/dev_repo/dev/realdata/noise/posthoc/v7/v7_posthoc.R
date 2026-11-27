# Post hoc descriptions for report_v7.md (not criteria). Development half only;
# the confirmation half stays unread. Paper simulations: report only.
source("dev/R/common.R"); source("dev/realdata/noise/v7/common_v7.R")
ev <- readRDS(v7_path("calibration.rds")); D <- ev$data
auc <- function(s, y) { s[is.na(s)] <- -Inf; r <- rank(s); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n0 * n1) }
dv <- D[D$half == "dev", ]; dv$l2 <- log2(dv$S_B); dv$grp <- ifelse(dv$family == "c2", "C2", dv$design)
A <- function(x, s = "S_B") { x <- x[x$m == 1 | x$m >= 4, ]; if (length(unique(x$m >= 4)) < 2) return(NA); auc(x[[s]], x$m >= 4) }
cat("== AUC (x1 vs >= x4), development half ==\n")
cat("by design:\n"); print(round(sapply(split(dv, dv$grp), A), 3))
cat("S_A by design:\n"); print(round(sapply(split(dv, dv$grp), A, s = "S_A"), 3))
cat("by depth pair:\n"); print(round(sapply(split(dv, dv$pair), A), 3))
cat("by wavelength:\n"); print(round(sapply(split(dv, dv$wavelength), A), 3))
cat("by bottleneck:\n"); print(round(sapply(split(dv, dv$bottleneck), A), 3))
cat("by interval:\n"); print(round(sapply(split(dv, dv$dt), A), 3))
cat("\n== median log2 S_B, development half, by design and multiple ==\n")
print(round(tapply(dv$l2, list(dv$grp, dv$m), stats::median), 2))
cat("by wavelength and multiple:\n"); print(round(tapply(dv$l2, list(dv$wavelength, dv$m), stats::median), 2))
cat("by depth pair and multiple:\n"); print(round(tapply(dv$l2, list(dv$pair, dv$m), stats::median), 2))
cat("\n== spread: IQR of log2 S_B within design x multiple (development) ==\n")
print(round(tapply(dv$l2, list(dv$grp, dv$m), stats::IQR), 2))
# paired: the same landscape, design, window and depth pair at x1 and x16 / x4
key <- function(x) sub("_x[0-9]{2}_", "_", x$dataset_id)
dv$key <- key(dv)
pair_share <- function(hi) {
  a <- dv[dv$m == 1, ]; b <- dv[dv$m == hi, ]
  j <- match(a$key, b$key); ok <- !is.na(j)
  data.frame(pairs = sum(ok), higher = mean(b$l2[j[ok]] > a$l2[ok]), median_diff = stats::median(b$l2[j[ok]] - a$l2[ok]))
}
cat("\n== paired (same landscape, design, window, depth pair), development half ==\n")
print(rbind(x4 = pair_share(4), x16 = pair_share(16)), digits = 3)
cat("paired x16 vs x1 by design:\n")
a <- dv[dv$m == 1, ]; b <- dv[dv$m == 16, ]; j <- match(a$key, b$key)
print(round(tapply(b$l2[j] > a$l2, a$grp, mean, na.rm = TRUE), 3))
cat("paired x16 vs x1 by depth pair:\n"); print(round(tapply(b$l2[j] > a$l2, a$pair, mean, na.rm = TRUE), 3))
cat("\n== between-landscape spread at x1 vs the x1 -> x16 shift (development, median over designs) ==\n")
lm_ <- tapply(dv$l2[dv$m == 1], dv$landscape_id[dv$m == 1], stats::median); print(round(range(lm_), 2))
pp <- D[D$family == "paper", ]; cat("\npaper simulations, AUC S_B:", round(A(pp), 3), " S_A:", round(A(pp, "S_A"), 3), "\n")
cat("S_A not computable, development, by design and multiple:\n"); print(tapply(is.na(dv$S_A), list(dv$grp, dv$m), mean) |> round(2))
