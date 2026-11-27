#!/usr/bin/env Rscript
# beta2 of dev/realdata/beta/PLAN.md: picks of E1-E4 (single and pooled), the
# criteria N, A1, A2, S1, S2 on the development split, then on the judging
# split for the estimators that pass development. Writes realdata/beta/evaluation.rds.
#   Rscript dev/realdata/beta/beta_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "beta")
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
crit_names <- c(E1 = "objective", E2 = "data_nll", E3 = "cv_nll")

read_records <- function(dir) {
  files <- list.files(dir, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  lapply(files, readRDS)
}

## ---------- simulations ----------
sim_rows <- do.call(rbind, lapply(read_records(file.path(root, "fits")), function(r) {
  if (!r$kind %in% c("eset", "c2")) return(NULL)
  base <- data.frame(kind = r$kind, set = r$set, dataset_id = r$dataset_id, multiple = r$multiple,
                     sim = identical(r$label, "bsim"), true_multiple = r$true_beta / 5e-5, ok = identical(r$status, "ok"),
                     stringsAsFactors = FALSE)
  if (!base$ok) return(cbind(base, objective = NA, data_nll = NA, cv_nll = NA, d0 = NA, d1 = NA, d2 = NA))
  m <- r$metrics
  cbind(base, objective = r$criteria$objective, data_nll = r$criteria$data_nll, cv_nll = r$criteria$cv_nll,
        d0 = m$rmse[m$stratum == "d0"], d1 = m$rmse[m$stratum == "d1"], d2 = m$rmse[m$stratum == "d2"])
}))
sim_rows$landscape <- as.integer(sub("^(mt|c2)_L([0-9]+)_.*$", "\\2", sim_rows$dataset_id))
sim_rows$split <- ifelse(sim_rows$landscape <= 6, "dev", "test")
sim_rows$cell <- ifelse(sim_rows$kind == "eset", sub("^mt_L[0-9]+_", "", sim_rows$dataset_id), sim_rows$set)
sim_rows$group <- ifelse(sim_rows$kind == "eset", paste(sim_rows$set, sim_rows$cell, sim_rows$true_multiple), sim_rows$set)
# the e-sets reuse dataset ids across e1-e4: identify a dataset by kind/set/id
sim_rows$dataset_id <- paste(sim_rows$kind, sim_rows$set, sim_rows$dataset_id, sep = "/")
grid_rows <- sim_rows[!sim_rows$sim, ]
cat("simulation fits:", nrow(grid_rows), "ok", sum(grid_rows$ok), "\n")
print(stats::aggregate(ok ~ set, grid_rows, function(x) round(mean(x), 3)))

# datasets with every grid rate fitted (a failure at every rate is a domain-coverage error of the data)
complete <- stats::aggregate(ok ~ dataset_id, grid_rows, function(x) sum(x) == length(grid))
complete_ids <- complete$dataset_id[complete$ok]
cat("datasets with all 9 rates:", length(complete_ids), "of", length(unique(grid_rows$dataset_id)), "\n")

pick_single <- function(rows, column) {
  do.call(rbind, lapply(split(rows, rows$dataset_id), function(x) {
    x <- x[is.finite(x[[column]]), ]
    if (!nrow(x)) return(NULL)
    data.frame(dataset_id = x$dataset_id[1], pick = x$multiple[which.min(x[[column]])])
  }))
}
pick_pooled <- function(rows, column, group = "group") {
  do.call(rbind, lapply(split(rows, rows[[group]]), function(x) {
    tot <- tapply(x[[column]], x$multiple, sum)
    p <- as.numeric(names(tot))[which.min(tot)]
    data.frame(dataset_id = unique(x$dataset_id), pick = p)
  }))
}
rmse_at <- function(rows, ids, multiples, d) {
  key <- paste(rows$dataset_id, rows$multiple)
  rows[[d]][match(paste(ids, multiples), key)]
}
sim_eval <- function(picks, rows) {
  info <- unique(rows[c("dataset_id", "set", "kind", "cell", "group", "true_multiple")])
  p <- merge(picks, info, by = "dataset_id")
  for (d in c("d0", "d1", "d2")) {
    p[[paste0(d, "_pick")]] <- rmse_at(rows, p$dataset_id, p$pick, d)
    p[[paste0(d, "_fixed")]] <- rmse_at(rows, p$dataset_id, rep(1, nrow(p)), d)
    p[[paste0("r_", d)]] <- p[[paste0(d, "_pick")]] / p[[paste0(d, "_fixed")]]   # = against the truth for e1-e3
  }
  p$within2 <- abs(log2(p$pick / p$true_multiple)) <= 1
  p
}
criteria_sims <- function(p, pooled) {
  med <- function(x) stats::median(x, na.rm = TRUE)
  e1 <- p[p$set == "e1", ]
  cells <- stats::aggregate(cbind(r_d0, r_d1, r_d2) ~ cell, e1, med)
  n_ok <- sum(cells$r_d0 <= 1.03 & cells$r_d1 <= 1.03 & cells$r_d2 <= 1.03)
  cell_med <- c(med(cells$r_d0), med(cells$r_d1), med(cells$r_d2))
  e23 <- sapply(c("e2", "e3"), function(s) { x <- p[p$set == s, ]; c(med(x$r_d0), med(x$r_d1), med(x$r_d2)) })
  N <- n_ok >= 18 && all(cell_med <= 1.02) && all(e23 <= 1.03)
  alt <- p[p$set == "e4" | p$kind == "c2", ]
  A1_value <- if (pooled) mean(tapply(alt$within2, alt$group, function(x) x[1])) else mean(alt$within2)
  A1 <- A1_value >= (if (pooled) 0.80 else 0.70)
  truths <- stats::aggregate(cbind(r_d0, r_d1, r_d2) ~ kind + true_multiple, alt, med)
  truths$ok <- ifelse(truths$true_multiple >= 4, truths$r_d0 < 1 & truths$r_d1 < 1 & truths$r_d2 < 1,
                      truths$r_d0 <= 1.02 & truths$r_d1 <= 1.02 & truths$r_d2 <= 1.02)
  A2 <- all(truths$ok)
  list(N = N, A1 = A1, A2 = A2, e1_cells_ok = n_ok, e1_cell_median = cell_med, e23 = e23, A1_value = A1_value,
       truths = truths, within_by_truth = stats::aggregate(within2 ~ kind + true_multiple, alt, mean))
}

## ---------- Salehi ----------
salehi_a <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
sal_rows <- do.call(rbind, lapply(read_records(file.path(root, "fits", "salehi")), function(r) {
  data.frame(id = r$id, line = r$line, group = r$group, multiple = r$multiple, ok = identical(r$status, "ok"),
             objective = if (identical(r$status, "ok")) r$criteria$objective else NA,
             data_nll = if (identical(r$status, "ok")) r$criteria$data_nll else NA,
             cv_nll = if (identical(r$status, "ok")) r$criteria$cv_nll else NA, stringsAsFactors = FALSE)
}))
# forecast_score.R writes list(task, status, rows) to <id>__ends__<method>__m10.Rds;
# the methods are the benchmark's beta labels (x1_4, x1_2, x1, ..., x64).
hz_files <- list.files(file.path(root, "horizon"), pattern = "[.]Rds$", full.names = TRUE)
hz_status <- if (length(hz_files)) vapply(hz_files, function(f) readRDS(f)$status, character(1)) else character(0)
if (length(hz_status)) { cat("horizon tasks by status:\n"); print(table(hz_status)) }
hz <- if (length(hz_files)) do.call(rbind, lapply(hz_files, function(f) { x <- readRDS(f); if (identical(x$status, "ok")) x$rows else NULL })) else NULL
x_multiple <- function(label) ifelse(grepl("^x1_", label), 1 / as.numeric(sub("^x1_", "", label)), as.numeric(sub("^x", "", label)))
if (!is.null(hz)) {
  hz$multiple <- x_multiple(hz$method)
  hz <- hz[hz$days <= 60 & is.finite(hz$ll_fc_s_e2), ]
  hz$line <- salehi_a$line[match(hz$id, salehi_a$id)]
}
salehi_picks <- function(column, pooled) {
  x <- sal_rows[sal_rows$ok & is.finite(sal_rows[[column]]), ]
  full <- tapply(x$multiple, x$id, length) == length(grid)
  x <- x[x$id %in% names(full)[full], ]
  if (!pooled) return(do.call(rbind, lapply(split(x, x$id), function(y) data.frame(id = y$id[1], pick = y$multiple[which.min(y[[column]])]))))
  do.call(rbind, lapply(split(x, x$line), function(y) {
    tot <- tapply(y[[column]], y$multiple, sum)
    data.frame(id = unique(y$id), pick = as.numeric(names(tot))[which.min(tot)])
  }))
}
e4_picks <- function() {
  if (is.null(hz)) return(NULL)
  do.call(rbind, lapply(unique(hz$id), function(i) {
    others <- hz[hz$line == hz$line[hz$id == i][1] & hz$id != i, ]
    if (!nrow(others)) return(NULL)
    s <- tapply(others$ll_fc_s_e2, others$multiple, mean)
    data.frame(id = i, pick = as.numeric(names(s))[which.max(s)])
  }))
}
salehi_criteria <- function(picks, group, single) {
  ids <- salehi_a$id[ifelse(salehi_a$line %in% c("SA609", "SA039"), "dev", "confirm") == group]
  p <- picks[picks$id %in% ids, ]
  S1_value <- NA
  if (!is.null(hz)) {
    h <- hz[hz$id %in% p$id, ]
    at <- h[h$multiple == p$pick[match(h$id, p$id)], c("id", "uid", "ll_fc_s_e2")]
    base <- h[h$multiple == 1, c("id", "uid", "ll_fc_s_e2")]
    m <- merge(at, base, by = c("id", "uid"), suffixes = c("_pick", "_fixed"))
    S1_value <- mean(m$ll_fc_s_e2_pick - m$ll_fc_s_e2_fixed)
  }
  S2_value <- NA
  if (single) {
    sa <- salehi_a[salehi_a$id %in% p$id, ]
    keys <- split(sa, list(sa$line, sa$treat_first, sa$treat_last), drop = TRUE)
    agree <- unlist(lapply(keys, function(g) {
      if (nrow(g) < 2) return(NULL)
      cmb <- utils::combn(nrow(g), 2)
      vapply(seq_len(ncol(cmb)), function(j) {
        a <- g[cmb[1, j], ]; b <- g[cmb[2, j], ]
        if (length(intersect(c(a$first_uid, a$last_uid), c(b$first_uid, b$last_uid)))) return(NA)
        abs(log2(p$pick[p$id == a$id] / p$pick[p$id == b$id])) <= 1
      }, logical(1))
    }))
    S2_value <- mean(agree, na.rm = TRUE)
  }
  list(S1 = is.finite(S1_value) && S1_value >= 0, S1_value = S1_value,
       S2 = if (single) is.finite(S2_value) && S2_value >= 0.70 else NA, S2_value = S2_value,
       picks_by_line = stats::aggregate(pick ~ line, merge(p, salehi_a[c("id", "line")], by = "id"),
                                        function(x) paste(names(table(x)), table(x), sep = ":", collapse = " ")))
}

## ---------- run ----------
results <- list()
for (e in names(crit_names)) for (pooled in c(FALSE, TRUE)) {
  name <- paste0(e, if (pooled) "_pooled" else "_single")
  res <- list()
  for (sp in c("dev", "test")) {
    # the judging split only for an estimator that passed development (PLAN, beta2)
    if (sp == "test" && (!isTRUE(res$dev$pass) || Sys.getenv("BETA_DEV_ONLY") == "1")) break
    rows <- grid_rows[grid_rows$split == sp & grid_rows$dataset_id %in% complete_ids, ]
    picks <- if (pooled) pick_pooled(rows, crit_names[[e]]) else pick_single(rows, crit_names[[e]])
    p <- sim_eval(picks, rows)
    cs <- criteria_sims(p, pooled)
    ss <- salehi_criteria(salehi_picks(crit_names[[e]], pooled), if (sp == "dev") "dev" else "confirm", !pooled)
    res[[sp]] <- list(sims = cs, salehi = ss, picks = p,
                      pass = cs$N && cs$A1 && cs$A2 && isTRUE(ss$S1) && (pooled || isTRUE(ss$S2)))
  }
  res$verdict <- isTRUE(res$dev$pass) && isTRUE(res$test$pass)
  results[[name]] <- res
}
e4 <- e4_picks()
if (!is.null(e4)) {
  results$E4_pooled <- list(dev = salehi_criteria(e4, "dev", FALSE))
  results$E4_pooled$dev$pass <- isTRUE(results$E4_pooled$dev$S1)
  if (results$E4_pooled$dev$pass && Sys.getenv("BETA_DEV_ONLY") != "1") {
    results$E4_pooled$test <- salehi_criteria(e4, "confirm", FALSE)
    results$E4_pooled$test$pass <- isTRUE(results$E4_pooled$test$S1)
  }
  results$E4_pooled$verdict <- isTRUE(results$E4_pooled$dev$pass) && isTRUE(results$E4_pooled$test$pass)
}

summ <- do.call(rbind, lapply(names(results), function(n) {
  r <- results[[n]]
  row <- function(sp) {
    x <- r[[sp]]
    if (is.null(x)) return(data.frame(N = NA, e1_cells = NA, A1 = NA, A1_value = NA, A2 = NA, S1 = NA, S1_value = NA,
                                      S2 = NA, S2_value = NA, pass = NA))
    if (is.null(x$sims)) return(data.frame(N = NA, e1_cells = NA, A1 = NA, A1_value = NA, A2 = NA,
                                           S1 = x$S1, S1_value = x$S1_value, S2 = NA, S2_value = NA, pass = x$pass))
    data.frame(N = x$sims$N, e1_cells = x$sims$e1_cells_ok, A1 = x$sims$A1, A1_value = x$sims$A1_value, A2 = x$sims$A2,
               S1 = x$salehi$S1, S1_value = x$salehi$S1_value, S2 = x$salehi$S2, S2_value = x$salehi$S2_value, pass = x$pass)
  }
  rbind(cbind(estimator = n, split = "dev", row("dev")), cbind(estimator = n, split = "test", row("test")))
}))
print(summ, digits = 3, row.names = FALSE)
cat("\nverdicts:\n"); print(vapply(results, function(r) isTRUE(r$verdict), logical(1)))
dev_save_rds(list(results = results, summary = summ, sim_rows = sim_rows, sal_rows = sal_rows, horizon = hz),
             file.path(root, "evaluation.rds"))
