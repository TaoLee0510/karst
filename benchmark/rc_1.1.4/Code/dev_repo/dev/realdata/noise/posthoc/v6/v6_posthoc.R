# Post hoc descriptions for report_v6.md (not criteria).
source("dev/R/common.R")
ev <- readRDS(file.path(dev_results, "realdata", "noise", "v6", "calibration.rds"))
D <- ev$data; D$l2 <- log2(D$S_B)
cat("== median log2 S_B by family and multiple (both halves) ==\n")
print(round(tapply(D$l2, list(D$family, D$m), stats::median), 2))
cat("\n== slope of log2 S_B on log2 multiple, per family (all datasets) ==\n")
for (fm in c("e", "C2", "paper")) { x <- D[D$family == fm, ]; f <- stats::lm(l2 ~ log2(m), x)
  cat(sprintf("%-6s n %4d  c %.2f  d %.3f  sigma %.2f  R2 %.2f\n", fm, nrow(x), coef(f)[1], coef(f)[2], summary(f)$sigma, summary(f)$r.squared)) }
x <- D[D$family == "e", ]; f <- stats::lm(l2 ~ log2(m) + factor(set == "e2"), x)
cat("\n== x1 e sets: share raised >= 2x by bottleneck / depth / interval (both halves) ==\n")
e1 <- D[D$family == "e" & D$m == 1, ]
sc <- readRDS; meta <- do.call(rbind, lapply(list.files(file.path(dev_results, "realdata", "noise", "v6", "stats"), pattern = "^eset", full.names = TRUE), function(f)
  do.call(rbind, lapply(readRDS(f), function(r) data.frame(set = r$set, dataset_id = r$dataset_id, bottleneck = r$bottleneck)))))
e1$bottleneck <- meta$bottleneck[match(paste(e1$set, e1$dataset_id), paste(meta$set, meta$dataset_id))]
e1$r2 <- e1$beta_used >= 1e-4 * (1 - 1e-9)
print(round(tapply(e1$r2, e1$bottleneck, mean), 2)); print(round(tapply(e1$r2, e1$depth, mean), 2)); print(round(tapply(e1$r2, e1$dt, mean), 2))
cat("\n== paper: per window, median log2 S_B at x1 and x16 ==\n")
p <- D[D$family == "paper", ]; print(round(tapply(p$l2, list(p$window, p$m), stats::median), 2))
cat("\n== real data: S_B (log2) quantiles ==\n")
R <- ev$real; R$l2 <- log2(R$S_B)
print(do.call(rbind, lapply(split(R, R$line), function(x) data.frame(line = x$line[1], n = nrow(x), q25 = quantile(x$l2, .25), med = median(x$l2), q75 = quantile(x$l2, .75), identical_pairs_zero = sum(x$I0 + x$I1 == 0)))), digits = 3, row.names = FALSE)
