#!/usr/bin/env Rscript
# Report-only breakdowns of the v3 evaluation (not criteria of PLAN_v3.md):
# the criteria at T = 1000, where the "sampling" verdicts of criterion 3 come
# from (set and bottleneck), criterion 1 by interval length, and region
# sampling by bottleneck. Reads realdata/noise/v3/evaluation.rds.
#   Rscript dev/realdata/noise/v3/v3_extra.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 200)
ev <- readRDS(file.path(dev_results, "realdata", "noise", "v3", "evaluation.rds"))
recs <- ev$records
samp <- function(x) if (nrow(x)) mean(x$verdict == "sampling") else NA_real_
share <- function(x, v) if (nrow(x)) mean(x$verdict == v) else NA_real_
recs$drift_class <- ifelse(recs$bottleneck <= 300, "strong (200/300)", "weak (2000-20000)")

cat("== criteria at T = 1000 (report only) ==\n")
e3 <- recs[recs$kind == "eset" & recs$set == "e3", ]
sel <- recs$variant == "none" & ((recs$kind == "eset" & recs$set %in% c("e1", "e2", "e4")) | recs$kind == "c2")
at <- function(T) {
  e3n <- e3[e3$variant == "none" & e3$depth >= T, ]
  data.frame(T = T,
             c1 = samp(e3n[e3n$bottleneck == 20000, ]),
             s200 = samp(e3n[e3n$bottleneck == 200, ]), s2000 = samp(e3n[e3n$bottleneck == 2000, ]),
             s20000 = samp(e3n[e3n$bottleneck == 20000, ]),
             c3 = samp(recs[sel & recs$depth >= T, ]),
             c3_weak_drift = samp(recs[sel & recs$depth >= T & recs$bottleneck > 300, ]),
             c4 = samp(e3[e3$variant == "dm20" & e3$bottleneck == 20000 & e3$depth >= T, ]))
}
print(rbind(at(300), at(500), at(1000)), digits = 3, row.names = FALSE)

cat("\n== criterion 3: sampling share of the selection simulations by group (no injection) ==\n")
s <- recs[sel, ]
s$group <- ifelse(s$kind == "c2", "C2", s$set)
tab3 <- do.call(rbind, lapply(split(s, list(s$group, s$drift_class, s$depth), drop = TRUE), function(x)
  data.frame(group = x$group[1], drift = x$drift_class[1], depth = x$depth[1], n = nrow(x),
             sampling = samp(x), selection = share(x, "selection"))))
tab3 <- tab3[order(tab3$group, tab3$drift, tab3$depth), ]
print(tab3, digits = 3, row.names = FALSE)
cat("\nC2 by true rate (no injection, depth >= 300):\n")
c2 <- recs[recs$kind == "c2" & recs$variant == "none", ]
print(do.call(rbind, lapply(split(c2, list(c2$truth, c2$bottleneck, c2$depth), drop = TRUE), function(x)
  data.frame(beta = x$truth[1], bottleneck = x$bottleneck[1], depth = x$depth[1], n = nrow(x),
             no_anomaly = share(x, "no_anomaly"), sampling = samp(x)))), digits = 3, row.names = FALSE)

cat("\n== criterion 1 (e3, bottleneck 20000, no injection) by interval between the two time points ==\n")
c1 <- e3[e3$variant == "none" & e3$bottleneck == 20000, ]
print(do.call(rbind, lapply(split(c1, list(c1$dt, c1$depth), drop = TRUE), function(x)
  data.frame(dt = x$dt[1], depth = x$depth[1], n = nrow(x), sampling = samp(x), undetermined = share(x, "undetermined")))),
  digits = 3, row.names = FALSE)

cat("\n== region sampling: selection and sampling shares by bottleneck (depth >= 300) ==\n")
rg <- recs[recs$variant %in% c("reg05", "reg10") & recs$depth >= 300, ]
rg$group <- ifelse(rg$kind == "c2", "C2", rg$set)
print(do.call(rbind, lapply(split(rg, list(rg$group, rg$variant, rg$bottleneck), drop = TRUE), function(x)
  data.frame(group = x$group[1], variant = x$variant[1], bottleneck = x$bottleneck[1], n = nrow(x),
             selection = share(x, "selection"), sampling = samp(x)))), digits = 3, row.names = FALSE)

cat("\n== 200 cells per time point (e sets) ==\n")
d2 <- recs[recs$kind == "eset" & recs$depth == 200, ]
print(do.call(rbind, lapply(split(d2, list(d2$set, d2$variant), drop = TRUE), function(x)
  data.frame(set = x$set[1], variant = x$variant[1], n = nrow(x), no_anomaly = share(x, "no_anomaly"),
             undetermined = share(x, "undetermined"), selection = share(x, "selection"), sampling = samp(x)))),
  digits = 3, row.names = FALSE)
