#!/usr/bin/env Rscript
# D3 (WGD) of dev/realdata/PLAN.md, section C–E: identity check, W2 and W3
# (W1 when the C2 W datasets exist). Reads the local copies:
#   dev/results/phase2/fits/{ce_base,d3_none,d3_w4,d3_w3}/<kind>/
#   dev/results/realdata/CE/d3/horizon/ and dev/results/realdata/B/horizon/
#   Rscript dev/realdata/CE/d3_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
fits <- file.path(dev_results, "phase2", "fits")
read_fit <- function(method, kind, id) {
  f <- file.path(fits, method, kind, paste0(id, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

## identity: D3 without wgd against the cd342db base
ident <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "CE", "tasks", "d3_ident.tsv"), stringsAsFactors = FALSE)
cat("== identity (d3 without wgd vs cd342db base) ==\n")
for (i in seq_len(nrow(ident))) {
  a <- read_fit("d3_none", ident$kind[i], ident$dataset_id[i]); b <- read_fit("ce_base", ident$kind[i], ident$dataset_id[i])
  if (is.null(a) || is.null(b)) { cat(ident$dataset_id[i], ": missing\n"); next }
  if (!identical(a$status, "ok") || !identical(b$status, "ok")) { cat(ident$dataset_id[i], ":", a$status, b$status, "\n"); next }
  same <- identical(a$near_nodes, b$near_nodes) && identical(a$far_nodes, b$far_nodes) && identical(a$objective, b$objective)
  diff <- max(abs(a$near_nodes$pred_relative - b$near_nodes$pred_relative), abs(a$near_nodes$model_sd - b$near_nodes$model_sd))
  cat(sprintf("%s %s: identical %s, max abs difference %.3g\n", ident$kind[i], ident$dataset_id[i], same, diff))
}

## W2: no-WGD dev GRF (dev split), d0-d2 RMSE ratio medians <= 1.01
w2 <- function(method, kind = "dev") {
  ids <- sub("\\.rds$", "", list.files(file.path(fits, method, kind)))
  rows <- do.call(rbind, lapply(ids, function(id) {
    a <- read_fit(method, kind, id); b <- read_fit("ce_base", kind, id)
    if (is.null(a) || is.null(b) || !identical(a$status, "ok") || !identical(b$status, "ok")) return(NULL)
    m <- merge(a$near_metrics[a$near_metrics$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse")],
               b$near_metrics[b$near_metrics$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse")], by = "stratum")
    data.frame(dataset_id = id, stratum = m$stratum, ratio = m$rmse.x / m$rmse.y)
  }))
  if (is.null(rows)) return(NULL)
  failed <- sum(vapply(ids, function(id) !identical(read_fit(method, kind, id)$status, "ok"), logical(1)))
  med <- tapply(rows$ratio, rows$stratum, stats::median, na.rm = TRUE)
  list(n = length(unique(rows$dataset_id)), failed = failed, median = med, pass = all(med <= 1.01))
}
cat("\n== W2: dev GRF (dev split), RMSE ratio D3/base, median over datasets ==\n")
W2 <- list()
for (m in c("d3_w4", "d3_w3")) {
  if (!dir.exists(file.path(fits, m, "dev"))) { cat(m, ": no fits yet\n"); next }
  W2[[m]] <- w2(m)
  if (is.null(W2[[m]])) { cat(m, ": no paired fits yet\n"); next }
  cat(sprintf("%s: %d datasets (%d failed); d0 %.4f, d1 %.4f, d2 %.4f -> %s\n", m, W2[[m]]$n, W2[[m]]$failed,
              W2[[m]]$median[["d0"]], W2[[m]]$median[["d1"]], W2[[m]]$median[["d2"]], if (W2[[m]]$pass) "pass" else "FAIL"))
}

## W3: SA906a eligible ends trajectories, minobs 10, paired ll_fc_e2 (D3 vs B's cd342db base)
read_rows <- function(dir, pattern) {
  fs <- list.files(dir, pattern = pattern, full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) { r <- readRDS(f); if (identical(r$status, "ok")) r$rows else NULL }))
}
tasks <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "CE", "tasks", "d3_salehi.tsv"), stringsAsFactors = FALSE)
ids <- unique(tasks$id)
base <- read_rows(file.path(dev_results, "realdata", "B", "horizon"), "__ends__karst__m10\\.Rds$")
base <- base[base$id %in% ids, ]
cat("\n== W3: SA906a, minobs 10, ll_fc_e2 (D3 - base), paired on (trajectory, descendant) ==\n")
W3 <- list()
d3dir <- file.path(dev_results, "realdata", "CE", "d3", "horizon")
for (m in c("d3_w4", "d3_w3")) {
  rows <- if (dir.exists(d3dir)) read_rows(d3dir, paste0("__ends__", m, "__m10\\.Rds$")) else NULL
  if (is.null(rows)) { cat(m, ": no rows yet\n"); next }
  p <- merge(rows[c("id", "uid", "days", "ll_fc_e2", "ll_p0_e2", "coverage")], base[c("id", "uid", "ll_fc_e2", "coverage")],
             by = c("id", "uid"), suffixes = c("", ".base"))
  d <- p$ll_fc_e2 - p$ll_fc_e2.base
  W3[[m]] <- list(rows = nrow(p), trajectories = length(unique(p$id)), mean = mean(d), sd = stats::sd(d),
                  better = mean(d > 0), pass = mean(d) >= 0, coverage = stats::median(p$coverage),
                  coverage_base = stats::median(p$coverage.base))
  cat(sprintf("%s: %d rows from %d trajectories (base %d rows); mean difference %.4f (sd %.4f), D3 better in %.2f -> %s; landscape coverage median %.3f (base %.3f)\n",
              m, nrow(p), W3[[m]]$trajectories, nrow(base), W3[[m]]$mean, W3[[m]]$sd, W3[[m]]$better,
              if (W3[[m]]$pass) "pass" else "FAIL", W3[[m]]$coverage, W3[[m]]$coverage_base))
}

## W1 (C2 W): only when the C2 datasets and their fits exist
c2 <- file.path(dev_results, "realdata", "CE", "c2")
cat("\n== W1 (C2 W) ==\n")
if (!dir.exists(c2)) cat("C2 datasets not available: W1 pending\n")
dev_save_rds(list(W2 = W2, W3 = W3), file.path(dev_results, "realdata", "CE", "d3", "evaluation.rds"))
