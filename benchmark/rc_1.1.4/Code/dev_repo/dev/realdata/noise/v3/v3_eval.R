#!/usr/bin/env Rscript
# Criteria of dev/realdata/noise/PLAN_v3.md on the new simulations, at the
# enabling thresholds T = 300 and 500 cells per time point, the choice of T,
# and the report-only tables. Writes realdata/noise/v3/evaluation.rds.
#   Rscript dev/realdata/noise/v3/v3_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "noise", "v3", "analysis")
files <- list.files(root, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
recs <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  base <- data.frame(kind = x$kind, set = x$set, dataset_id = x$dataset_id, variant = x$variant, truth = x$truth,
                     depth = x$depth, bottleneck = x$bottleneck, dt = x$dt, status = x$status, stringsAsFactors = FALSE)
  if (!identical(x$status, "ok")) return(cbind(base, n_abnormal = NA, excess = NA, n_edges = NA, p = NA, verdict = "failed"))
  cbind(base, x$summary[c("n_abnormal", "excess", "n_edges", "p", "verdict")])
}))
cat("analyses:", nrow(recs), " failed:", sum(recs$status != "ok"), "\n")
print(table(recs$kind, recs$variant))
samp <- function(x) if (nrow(x)) mean(x$verdict == "sampling") else NA_real_
e3 <- recs[recs$kind == "eset" & recs$set == "e3", ]
sel_sets <- recs$variant == "none" & ((recs$kind == "eset" & recs$set %in% c("e1", "e2", "e4")) | recs$kind == "c2")

criteria_at <- function(T) {
  e3n <- e3[e3$variant == "none" & e3$depth >= T, ]
  s200 <- samp(e3n[e3n$bottleneck == 200, ]); s2k <- samp(e3n[e3n$bottleneck == 2000, ]); s20k <- samp(e3n[e3n$bottleneck == 20000, ])
  c1 <- s20k
  c2 <- isTRUE(s200 > s2k && s2k > s20k && (s200 - s20k) >= 0.30)
  c3 <- samp(recs[sel_sets & recs$depth >= T, ])
  c4 <- samp(e3[e3$variant == "dm20" & e3$bottleneck == 20000 & e3$depth >= T, ])
  data.frame(T = T, c1_nodrift_sampling = c1, c1 = c1 <= 0.10, s200 = s200, s2000 = s2k, s20000 = s20k, c2 = c2,
             c3_selection_sampling = c3, c3 = c3 <= 0.10, c4_dm20_sampling = c4, c4 = c4 >= 0.70,
             pass = c1 <= 0.10 && c2 && c3 <= 0.10 && c4 >= 0.70)
}
crit <- rbind(criteria_at(300), criteria_at(500))
cat("\n== criteria (PLAN_v3) ==\n"); print(crit, digits = 3, row.names = FALSE)
chosen <- if (crit$pass[crit$T == 300]) 300 else if (crit$pass[crit$T == 500]) 500 else NA
cat("\nENABLING THRESHOLD:", if (is.na(chosen)) "none (not passed)" else chosen, "\n")

lev <- c("no_anomaly", "undetermined", "selection", "sampling", "failed")
shares <- function(x) { t <- prop.table(table(factor(x$verdict, levels = lev))); round(as.numeric(t), 3) }
tab <- do.call(rbind, lapply(split(recs, list(recs$kind, recs$set, recs$variant, recs$depth, recs$bottleneck), drop = TRUE), function(x)
  data.frame(kind = x$kind[1], set = x$set[1], variant = x$variant[1], depth = x$depth[1], bottleneck = x$bottleneck[1],
             n = nrow(x), t(setNames(shares(x), lev)))))
tab <- tab[order(tab$kind, tab$set, tab$variant, tab$bottleneck, tab$depth), ]
cat("\n== verdict shares by set, injection, bottleneck and cells per time point ==\n"); print(tab, row.names = FALSE)
dev_save_rds(list(records = recs, criteria = crit, threshold = chosen, table = tab),
             file.path(dev_results, "realdata", "noise", "v3", "evaluation.rds"))
