#!/usr/bin/env Rscript
# P4 of dev/realdata/noise/PLAN_v2.md: verdict shares on the simulations, the
# three pre-registered criteria, and the report-only tables.
#   Rscript dev/realdata/noise/v2_eval.R
# Writes realdata/noise/v2_evaluation.rds.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "noise", "v2")
dirs <- list.dirs(root, recursive = FALSE)
recs <- do.call(rbind, lapply(dirs, function(d) {
  files <- list.files(d, pattern = "[.]rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    base <- data.frame(kind = x$kind, set = x$set, dataset_id = x$dataset_id, variant = x$variant, rho = x$rho,
                       truth = x$truth, status = x$status, stringsAsFactors = FALSE)
    if (!identical(x$status, "ok")) return(cbind(base, n0 = NA, n1 = NA, dt = NA, n_testable = NA, n_edges = NA,
                                                  n_abnormal = NA, n_rise = NA, n_fall = NA, excess = NA, moran = NA,
                                                  p = NA, near_share = NA, verdict = "failed"))
    cbind(base, x$summary[c("n0", "n1", "dt", "n_testable", "n_edges", "n_abnormal", "n_rise", "n_fall", "excess",
                            "moran", "p", "near_share", "verdict")])
  }))
}))
cat("analyses:", nrow(recs), " failed:", sum(recs$status != "ok"), "\n")
sims <- recs[recs$kind %in% c("eset", "c2"), ]
sims$family <- ifelse(sims$kind == "c2", paste0("C2 ", sims$set), sims$set)
sims$verdict <- factor(sims$verdict, levels = c("no_anomaly", "undetermined", "selection", "sampling", "failed"))
share <- function(x) round(prop.table(table(x$verdict)), 3)
tab <- do.call(rbind, lapply(split(sims, list(sims$family, sims$variant), drop = TRUE), function(x)
  data.frame(family = x$family[1], variant = x$variant[1], n = nrow(x), t(as.matrix(share(x))),
             median_abnormal = stats::median(x$n_abnormal, na.rm = TRUE), median_excess = stats::median(x$excess, na.rm = TRUE),
             median_edges = stats::median(x$n_edges, na.rm = TRUE), check.names = FALSE)))
tab <- tab[order(tab$variant, tab$family), ]
cat("\n== verdict shares by set and injected noise ==\n"); print(tab, row.names = FALSE)

v <- function(x) mean(x$verdict == "sampling")
c1 <- v(sims[sims$set == "e3" & sims$variant == "none", ])
sel <- sims[sims$variant == "none" & (sims$set %in% c("e1", "e2", "e4") | sims$kind == "c2"), ]
c2v <- v(sel)
c3 <- v(sims[sims$set == "e3" & sims$variant == "dm20", ])
crit <- c(neutral_le_0.10 = c1 <= 0.10, selection_le_0.10 = c2v <= 0.10, dm20_ge_0.70 = c3 >= 0.70)
cat(sprintf("\nP4: e3 sampling %.3f (<= 0.10); selection sims sampling %.3f (<= 0.10); e3 + phi 20 sampling %.3f (>= 0.70)\n", c1, c2v, c3))
print(crit)
cat("\nselection sims by family (sampling share):\n")
print(stats::aggregate(cbind(sampling = verdict == "sampling", selection = verdict == "selection") ~ family, sel, mean))
cat("\nby cells per time point (e sets, none):\n")
es <- sims[sims$kind == "eset", ]
es$cells <- sub("^.*_n([0-9]{4})_k.*$", "\\1", es$dataset_id)
print(stats::aggregate(cbind(sampling = verdict == "sampling", selection = verdict == "selection",
                             no_anomaly = verdict == "no_anomaly") ~ set + variant + cells, es, mean))
dev_save_rds(list(records = recs, table = tab, criteria = crit, values = c(neutral = c1, selection = c2v, dm20 = c3)),
             file.path(dev_results, "realdata", "noise", "v2_evaluation.rds"))
