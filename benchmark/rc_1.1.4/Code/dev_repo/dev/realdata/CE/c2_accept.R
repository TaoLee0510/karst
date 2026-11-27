#!/usr/bin/env Rscript
# Acceptance of a C2 round (dev/realdata/PLAN.md, phase C): dataset medians of
# the last-passage singleton share (counts), the domain coverage at the worst
# time point and the second-time-point Pearson dispersion (ce_base fits with
# cd342db, successful fits), per setting, against the ranges of regimes H and M
# (or W: ploidy-4 share). Among the settings inside all ranges of a regime, the
# one closest to its target (sum of |log(setting / target)| over the three
# medians) is chosen.
#   Rscript dev/realdata/CE/c2_accept.R ROUND
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
round <- commandArgs(TRUE)[[1L]]
s <- c2_settings()
s <- s[as.character(s$round) == round, ]

ranges <- list(H = list(singleton = c(0.34, 1.0), coverage = c(0.10, 0.40), dispersion = c(0.86, 2.5),
                        target = c(singleton = 0.64, coverage = 0.19, dispersion = 1.65)),
               M = list(singleton = c(0.10, 0.50), coverage = c(0.55, 0.92), dispersion = c(2.0, 6.0),
                        target = c(singleton = 0.25, coverage = 0.85, dispersion = 4.18)))

per_dataset <- do.call(rbind, lapply(s$setting, function(setting) {
  files <- list.files(file.path(c2_root, "grid", setting, "datasets"), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    d <- readRDS(f)
    st <- c2_count_stats(d)
    fit_path <- file.path(c2_root, "fits", "ce_base", setting, paste0(d$dataset_id, ".rds"))
    fit <- if (file.exists(fit_path)) readRDS(fit_path) else NULL
    st$setting <- setting
    st$fit <- if (is.null(fit)) "missing" else fit$status
    st$fail_class <- if (!is.null(fit) && fit$status != "ok") fit$condition_class else NA_character_
    st$coverage_min <- if (!is.null(fit) && fit$status == "ok") min(fit$diagnostics$domain_coverage$fraction) else NA_real_
    st$dispersion <- if (!is.null(fit) && fit$status == "ok") fit$dispersion else NA_real_
    st$n_major <- if (!is.null(fit) && fit$status == "ok") fit$diagnostics$n_major else NA_real_
    st
  }))
}))

med <- function(x) stats::median(x, na.rm = TRUE)
summary <- do.call(rbind, lapply(split(per_dataset, per_dataset$setting), function(x) {
  data.frame(setting = x$setting[[1L]], datasets = nrow(x), fitted = sum(x$fit == "ok"),
             failed = sum(x$fit == "failed"), missing = sum(x$fit == "missing"),
             singleton = med(x$singleton_share_last), coverage = med(x$coverage_min), dispersion = med(x$dispersion),
             karyotypes = med(x$karyotypes), n_major = med(x$n_major), shared = med(x$shared),
             ploidy4_first = med(x$ploidy4_first), ploidy4_last = med(x$ploidy4_last),
             rising = mean(x$ploidy4_last - x$ploidy4_first >= 0.1), stringsAsFactors = FALSE)
}))
summary <- summary[match(s$setting, summary$setting), ]
inside <- function(v, r) is.finite(v) && v >= r[[1L]] && v <= r[[2L]]
choice <- list()
for (regime in names(ranges)) {
  r <- ranges[[regime]]
  ok <- vapply(seq_len(nrow(summary)), function(i) {
    inside(summary$singleton[i], r$singleton) && inside(summary$coverage[i], r$coverage) &&
      inside(summary$dispersion[i], r$dispersion)
  }, logical(1))
  distance <- vapply(seq_len(nrow(summary)), function(i) {
    sum(abs(log(c(summary$singleton[i], summary$coverage[i], summary$dispersion[i]) / r$target)))
  }, numeric(1))
  summary[[paste0("in_", regime)]] <- ok
  summary[[paste0("distance_", regime)]] <- distance
  choice[[regime]] <- if (any(ok)) summary$setting[ok][which.min(distance[ok])] else NA_character_
}
summary$in_W <- with(summary, is.finite(ploidy4_first) & ploidy4_first >= 0.02 & ploidy4_first <= 0.5 & rising >= 0.25)
choice$W <- if (any(summary$in_W & s$ploidy_aware[match(summary$setting, s$setting)] == 1)) {
  summary$setting[summary$in_W & s$ploidy_aware[match(summary$setting, s$setting)] == 1][1L]
} else NA_character_
options(width = 220)
print(summary, digits = 3, row.names = FALSE)
cat("\nchosen:", paste(names(choice), unlist(choice), sep = " = ", collapse = "; "), "\n")
dev_save_rds(list(round = round, per_dataset = per_dataset, summary = summary, choice = choice, ranges = ranges),
             file.path(c2_root, paste0("accept_", round, ".rds")))
