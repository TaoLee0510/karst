#!/usr/bin/env Rscript
# Identity check of dev/realdata/noise/PLAN_v4.md step 2: on the datasets of
# tasks/check.tsv, karst_v4 without noise_phi (T0), karst_v4 with
# noise_phi = 1 (T2_phi1) and karst_6ae2a24 (T0__karst_6ae2a24) must give
# identical objectives and predictions.
#   Rscript dev/realdata/noise/v4/check_identity.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
root <- file.path(dev_results, "realdata", "noise", "v4", "fits")
tasks <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "noise", "v4", "tasks", "check.tsv"), stringsAsFactors = FALSE)
tasks <- unique(tasks[c("kind", "set", "variant", "dataset_id")])
get <- function(method, t) {
  f <- file.path(root, method, sprintf("%s_%s_%s", t$kind, t$set, t$variant), paste0(t$dataset_id, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}
out <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
  t <- tasks[i, ]
  a <- get("T0__karst_6ae2a24", t); b <- get("T0", t); c <- get("T2_phi1", t)
  cmp <- function(x, y) {
    if (is.null(x) || is.null(y) || !identical(x$status, "ok") || !identical(y$status, "ok")) return(NA_real_)
    max(abs(x$objective - y$objective), abs(x$predictions$relative_fitness - y$predictions$relative_fitness),
        abs(x$predictions$relative_sd - y$predictions$relative_sd), na.rm = TRUE)
  }
  data.frame(dataset_id = t$dataset_id, set = t$set, base_vs_v4 = cmp(a, b), v4_vs_phi1 = cmp(b, c),
             identical_base_v4 = !is.null(a) && !is.null(b) && identical(a$predictions, b$predictions),
             identical_v4_phi1 = !is.null(b) && !is.null(c) && identical(b$predictions, c$predictions))
}))
print(out, row.names = FALSE)
cat("ALL IDENTICAL:", all(out$identical_base_v4 & out$identical_v4_phi1), "\n")
