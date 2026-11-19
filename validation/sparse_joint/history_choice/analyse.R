#!/usr/bin/env Rscript
# History choice study: the rules applied to the candidates of every fit, and
# their accuracy against the truth.
#
#   Rscript validation/sparse_joint/history_choice/analyse.R RESULTS_DIR
#
# Reads RESULTS_DIR/dev/*.rds and RESULTS_DIR/abm/*.rds (dev_task.R,
# abm_task.R) and writes RESULTS_DIR/report.md and tables. Rules, each the
# lowest value among the candidates that have one:
#   current  nlminb's objective (karst 1.1.1);
#   A        the objective at the posterior mode for a candidate whose mode
#            was found by the second search, nlminb's for the others;
#   B        the objective at the posterior mode for every candidate;
#   oracle   the candidate closest to the truth (mean RMSE over distances 0-2
#            on the development benchmark, highest Spearman correlation over
#            all domain states on the ABM), for scale only.
# Criteria agreed on 2026-11-19, on the fits whose choice B changes:
#   development benchmark: median RMSE ratio B / current at most 1.00 over
#   distances 0, 1 and 2 together and at most 1.02 at each, and median 95%
#   coverage at most 0.01 lower; ABM: median difference of Spearman correlations B - current
#   at least 0 (all domain states), and B not significantly worse by a sign
#   test (one-sided p >= 0.05).
args <- commandArgs(TRUE)
results <- args[[1]]

read_part <- function(part) {
  files <- list.files(file.path(results, part), pattern = "\\.rds$", full.names = TRUE)
  lapply(files, readRDS)
}
rules <- c("current", "A", "B")
dev_tasks <- utils::read.delim(file.path(results, "dev_tasks.tsv"), stringsAsFactors = FALSE)
role_of <- function(dataset_id) dev_tasks$role[match(dataset_id, dev_tasks$dataset_id)]

choose <- function(candidates) {
  values <- list(current = candidates$objective,
                 A = ifelse(candidates$restart %in% TRUE, candidates$posterior_value, candidates$objective),
                 B = candidates$posterior_value)
  vapply(values, function(v) {
    v[!is.finite(candidates$objective)] <- NA
    if (!any(is.finite(v))) NA_real_ else candidates$history[[which.min(v)]]
  }, numeric(1))
}

value_gap <- function(candidates) {
  ok <- is.finite(candidates$objective) & is.finite(candidates$posterior_value)
  data.frame(gap = candidates$posterior_value[ok] - candidates$objective[ok],
             restart = candidates$restart[ok] %in% TRUE, convergence = candidates$convergence[ok])
}

sign_test <- function(better, worse) {
  n <- better + worse
  if (n == 0L) return(NA_real_)
  stats::binom.test(worse, n, alternative = "greater")$p.value
}

lines <- c("# History choice study", "",
           sprintf("Results: `%s`, analysed %s.", normalizePath(results), format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "")
emit <- function(...) lines <<- c(lines, ...)
table_md <- function(d, digits = 3) {
  d <- as.data.frame(d)
  fmt <- vapply(d, function(x) if (is.numeric(x)) TRUE else FALSE, logical(1))
  for (j in which(fmt)) d[[j]] <- ifelse(is.na(d[[j]]), "NA", formatC(d[[j]], digits = digits, format = "fg"))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- |", ncol(d)), collapse = "")),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
}

## Development benchmark -----------------------------------------------------
dev <- read_part("dev")
dev_ok <- Filter(function(r) identical(r$status, "ok"), dev)
emit("## Development benchmark (two time points)", "",
     sprintf("%d datasets read, %d fitted, %d failed; kept candidate identical to the fit returned: %d of %d.",
             length(dev), length(dev_ok), length(dev) - length(dev_ok),
             sum(vapply(dev_ok, function(r) isTRUE(r$identical_kept), logical(1))), length(dev_ok)), "")
if (length(dev_ok)) {
  dev_rows <- do.call(rbind, lapply(dev_ok, function(r) {
    chosen <- choose(r$candidates)
    s <- r$scores
    mean_rmse <- tapply(s$rmse, s$history, mean)
    oracle <- as.numeric(names(mean_rmse))[which.min(mean_rmse)]
    kept <- r$candidates$history[r$candidates$kept][[1L]]
    do.call(rbind, lapply(c(rules, "oracle"), function(rule) {
      h <- if (rule == "oracle") oracle else chosen[[rule]]
      sc <- s[abs(s$history - h) < 1e-9, ]
      data.frame(dataset_id = r$dataset_id, scenario = r$scenario, rule = rule, history = h,
                 kept_matches_current = abs(kept - chosen[["current"]]) < 1e-9,
                 stratum = sc$stratum, rmse = sc$rmse, rmse_zero = sc$rmse_zero, spearman = sc$spearman,
                 coverage95 = sc$coverage95, stringsAsFactors = FALSE)
    }))
  }))
  utils::write.csv(dev_rows, file.path(results, "dev_rows.csv"), row.names = FALSE)
  gaps <- do.call(rbind, lapply(dev_ok, function(r) value_gap(r$candidates)))
  emit(sprintf("Candidates: %d; posterior-mode objective differs from nlminb's by more than 1e-3 at %d (%.1f%%), found by the second search at %d.",
               nrow(gaps), sum(abs(gaps$gap) > 1e-3), 100 * mean(abs(gaps$gap) > 1e-3), sum(gaps$restart)),
       sprintf("Current rule reproduces the fit's own choice at %d of %d datasets.",
               sum(dev_rows$kept_matches_current[dev_rows$rule == "current" & dev_rows$stratum == "d0"]),
               sum(dev_rows$rule == "current" & dev_rows$stratum == "d0")), "")
  wide <- function(rule) dev_rows[dev_rows$rule == rule, c("dataset_id", "scenario", "stratum", "history", "rmse",
                                                           "spearman", "coverage95")]
  cur <- wide("current")
  for (rule in c("A", "B", "oracle")) {
    other <- wide(rule)
    m <- merge(cur, other, by = c("dataset_id", "scenario", "stratum"), suffixes = c("_cur", "_new"))
    changed <- m[abs(m$history_cur - m$history_new) > 1e-9, ]
    n_changed <- length(unique(changed$dataset_id))
    emit(sprintf("### Rule %s against the current rule", rule), "",
         sprintf("Choice changed at %d of %d datasets (core %d, stress %d).", n_changed, length(unique(m$dataset_id)),
                 length(unique(changed$dataset_id[role_of(changed$dataset_id) == "core"])),
                 length(unique(changed$dataset_id[role_of(changed$dataset_id) == "stress"]))), "")
    if (n_changed) {
      per <- do.call(rbind, lapply(split(changed, changed$stratum), function(d) {
        ratio <- d$rmse_new / d$rmse_cur
        data.frame(distance = d$stratum[1], datasets = nrow(d),
                   median_rmse_ratio = stats::median(ratio), q25 = stats::quantile(ratio, 0.25),
                   q75 = stats::quantile(ratio, 0.75), better = sum(ratio < 1), worse = sum(ratio > 1),
                   sign_p_worse = sign_test(sum(ratio < 1), sum(ratio > 1)),
                   median_spearman_diff = stats::median(d$spearman_new - d$spearman_cur, na.rm = TRUE),
                   median_coverage_diff = stats::median(d$coverage95_new - d$coverage95_cur, na.rm = TRUE))
      }))
      emit(table_md(per))
      if (rule == "B") {
        dev_b <- per
        dev_b_overall <- stats::median(changed$rmse_new / changed$rmse_cur)
        emit(sprintf("Median RMSE ratio over distances 0-2 together: %.3f.", dev_b_overall), "")
      }
    } else if (rule == "B") dev_b <- NULL
  }
  dev_rows$role <- role_of(dev_rows$dataset_id)
  pooled <- do.call(rbind, lapply(split(dev_rows, list(dev_rows$rule, dev_rows$stratum, dev_rows$role),
                                        drop = TRUE), function(d) {
    data.frame(rule = d$rule[1], distance = d$stratum[1], set = d$role[1],
               datasets = nrow(d), mean_rmse = mean(d$rmse), mean_rmse_zero = mean(d$rmse_zero),
               median_spearman = stats::median(d$spearman, na.rm = TRUE),
               mean_coverage = mean(d$coverage95, na.rm = TRUE))
  }))
  pooled <- pooled[order(pooled$set, pooled$distance, match(pooled$rule, c(rules, "oracle"))), ]
  emit("### All datasets, by rule", "", table_md(pooled))
  # Every dataset whose choice A or B changes, with the RMSE ratio to the
  # current rule's choice averaged over distances 0-2.
  per_dataset <- do.call(rbind, lapply(split(dev_rows, dev_rows$dataset_id), function(d) {
    h <- vapply(rules, function(rule) d$history[d$rule == rule][[1L]], numeric(1))
    rmse <- vapply(rules, function(rule) mean(d$rmse[d$rule == rule]), numeric(1))
    data.frame(dataset = d$dataset_id[1], set = d$role[1], current = h[["current"]], A = h[["A"]], B = h[["B"]],
               rmse_ratio_A = rmse[["A"]] / rmse[["current"]], rmse_ratio_B = rmse[["B"]] / rmse[["current"]])
  }))
  per_dataset <- per_dataset[per_dataset$A != per_dataset$current | per_dataset$B != per_dataset$current, ]
  if (nrow(per_dataset)) {
    emit("### Datasets whose choice A or B changes", "",
         "History chosen by each rule; mean RMSE over distances 0-2 relative to the current rule's choice.", "",
         table_md(per_dataset, digits = 4))
  }
} else dev_b <- NULL

## ABM -----------------------------------------------------------------------
abm <- read_part("abm")
abm_ok <- Filter(function(r) identical(r$status, "ok"), abm)
emit("## ALFA-K paper benchmark: ABM and missegregation simulations", "",
     sprintf("%d tasks read, %d fitted, %d failed; kept candidate identical to the fit returned: %d of %d.",
             length(abm), length(abm_ok), length(abm) - length(abm_ok),
             sum(vapply(abm_ok, function(r) isTRUE(r$identical_kept), logical(1))), length(abm_ok)), "")
failed_abm <- Filter(function(r) !identical(r$status, "ok"), abm)
if (length(failed_abm)) {
  emit("Failed:", "", paste0("- ", vapply(failed_abm, function(r) paste(r$task$source, r$task$dataset, r$task$design,
                                                                        r$task$method, r$condition_class), ""), collapse = NULL), "")
}
abm_b <- NULL
if (length(abm_ok)) {
  abm_rows <- do.call(rbind, lapply(abm_ok, function(r) {
    chosen <- choose(r$candidates)
    s <- r$scores
    all_rho <- s[s$group == "all", ]
    oracle <- if (any(is.finite(all_rho$spearman))) all_rho$history[which.max(all_rho$spearman)] else NA_real_
    kept <- r$candidates$history[r$candidates$kept][[1L]]
    restart_kept <- vapply(c(rules, "oracle"), function(rule) {
      h <- if (rule == "oracle") oracle else chosen[[rule]]
      any(r$candidates$restart[abs(r$candidates$history - h) < 1e-9] %in% TRUE)
    }, logical(1))
    do.call(rbind, lapply(c(rules, "oracle"), function(rule) {
      h <- if (rule == "oracle") oracle else chosen[[rule]]
      sc <- s[abs(s$history - h) < 1e-9, ]
      data.frame(key = paste(r$task$source, r$task$dataset, r$task$design, r$task$method),
                 set = r$task$set, source = r$task$source, design = r$task$design, method = r$task$method,
                 rule = rule, history = h, rescued_chosen = restart_kept[[rule]],
                 kept_matches_current = abs(kept - chosen[["current"]]) < 1e-9,
                 n_candidates = nrow(r$candidates), n_failed = sum(!is.finite(r$candidates$objective)),
                 n_rescued = sum(r$candidates$restart %in% TRUE),
                 group = sc$group, n = sc$n, pearson = sc$pearson, spearman = sc$spearman,
                 stringsAsFactors = FALSE)
    }))
  }))
  utils::write.csv(abm_rows, file.path(results, "abm_rows.csv"), row.names = FALSE)
  empty <- setdiff(unique(abm_rows$group), unique(abm_rows$group[abm_rows$n > 0]))
  if (length(empty)) {
    abm_rows <- abm_rows[!abm_rows$group %in% empty, ]
    emit(sprintf("Group left out, no state in any fit: %s (a fit's domain is the majors and their one-step neighbours).",
                 paste(empty, collapse = ", ")), "")
  }
  gaps <- do.call(rbind, lapply(abm_ok, function(r) value_gap(r$candidates)))
  first <- abm_rows[abm_rows$rule == "current" & abm_rows$group == "all", ]
  emit(sprintf("Candidates: %d; posterior-mode objective differs from nlminb's by more than 1e-3 at %d (%.1f%%), found by the second search at %d; candidates whose mode was not found at all: %d.",
               nrow(gaps), sum(abs(gaps$gap) > 1e-3), 100 * mean(abs(gaps$gap) > 1e-3), sum(gaps$restart),
               sum(first$n_failed)),
       sprintf("Current rule reproduces the fit's own choice at %d of %d tasks.", sum(first$kept_matches_current),
               nrow(first)), "")
  gap_summary <- data.frame(candidates = c("all", "second search", "first search"),
                            n = c(nrow(gaps), sum(gaps$restart), sum(!gaps$restart)),
                            median_gap = c(stats::median(gaps$gap), stats::median(gaps$gap[gaps$restart]),
                                           stats::median(gaps$gap[!gaps$restart])),
                            share_above_1e3 = c(mean(gaps$gap > 1e-3), mean(gaps$gap[gaps$restart] > 1e-3),
                                                mean(gaps$gap[!gaps$restart] > 1e-3)),
                            share_below_1e3 = c(mean(gaps$gap < -1e-3), mean(gaps$gap[gaps$restart] < -1e-3),
                                                mean(gaps$gap[!gaps$restart] < -1e-3)))
  emit("Objective at the posterior mode minus nlminb's (ABM candidates):", "", table_md(gap_summary))
  for (set in c("random", "failed")) {
    rows <- abm_rows[abm_rows$set == set, ]
    if (!nrow(rows)) next
    emit(sprintf("### Set %s (%d tasks)", set, length(unique(rows$key))), "")
    rescued <- rows[rows$group == "all", ]
    emit(sprintf("Tasks with a candidate rescued by the second search: %d; with a candidate whose mode was not found: %d.",
                 sum(rescued$rule == "current" & rescued$n_rescued > 0), sum(rescued$rule == "current" & rescued$n_failed > 0)),
         sprintf("Rescued candidate chosen: current %d, A %d, B %d.", sum(rescued$rule == "current" & rescued$rescued_chosen),
                 sum(rescued$rule == "A" & rescued$rescued_chosen), sum(rescued$rule == "B" & rescued$rescued_chosen)), "")
    cur <- rows[rows$rule == "current", ]
    for (rule in c("A", "B", "oracle")) {
      other <- rows[rows$rule == rule, ]
      m <- merge(cur, other, by = c("key", "group"), suffixes = c("_cur", "_new"))
      changed <- m[abs(m$history_cur - m$history_new) > 1e-9 & !is.na(m$history_new), ]
      n_changed <- length(unique(changed$key))
      emit(sprintf("Rule %s against the current rule: choice changed at %d of %d tasks.", rule, n_changed,
                   length(unique(m$key))), "")
      if (n_changed) {
        per <- do.call(rbind, lapply(split(changed, changed$group), function(d) {
          diff <- d$spearman_new - d$spearman_cur
          data.frame(group = d$group[1], tasks = nrow(d), median_spearman_diff = stats::median(diff, na.rm = TRUE),
                     mean_spearman_diff = mean(diff, na.rm = TRUE), better = sum(diff > 0, na.rm = TRUE),
                     worse = sum(diff < 0, na.rm = TRUE),
                     sign_p_worse = sign_test(sum(diff > 0, na.rm = TRUE), sum(diff < 0, na.rm = TRUE)),
                     median_pearson_diff = stats::median(d$pearson_new - d$pearson_cur, na.rm = TRUE))
        }))
        per <- per[match(c("all", "major", "shell1", "shell2"), per$group), ]
        per <- per[!is.na(per$group), ]
        emit(table_md(per))
        if (rule == "B" && set == "random") abm_b <- per
        if (rule == "B") {
          by_stratum <- as.data.frame(table(paste(changed$design_cur[changed$group == "all"],
                                                  changed$method_cur[changed$group == "all"])))
          names(by_stratum) <- c("design / method", "tasks changed")
          emit("Changed by design and method:", "", table_md(by_stratum))
        }
      }
    }
    pooled <- do.call(rbind, lapply(split(rows, list(rows$rule, rows$group), drop = TRUE), function(d) {
      data.frame(rule = d$rule[1], group = d$group[1], tasks = nrow(d),
                 median_spearman = stats::median(d$spearman, na.rm = TRUE),
                 mean_spearman = mean(d$spearman, na.rm = TRUE), mean_pearson = mean(d$pearson, na.rm = TRUE))
    }))
    pooled <- pooled[order(match(pooled$group, c("all", "major", "shell1", "shell2")),
                           match(pooled$rule, c(rules, "oracle"))), ]
    emit("All tasks of the set, by rule:", "", table_md(pooled))
    # Every task whose choice A or B changes.
    all_states <- rows[rows$group == "all", ]
    per_task <- do.call(rbind, lapply(split(all_states, all_states$key), function(d) {
      pick <- function(rule, column) d[[column]][d$rule == rule][[1L]]
      data.frame(task = d$key[1], current = pick("current", "history"), A = pick("A", "history"),
                 B = pick("B", "history"), oracle = pick("oracle", "history"),
                 current_rescued = pick("current", "rescued_chosen"),
                 spearman_current = pick("current", "spearman"),
                 diff_A = pick("A", "spearman") - pick("current", "spearman"),
                 diff_B = pick("B", "spearman") - pick("current", "spearman"),
                 diff_oracle = pick("oracle", "spearman") - pick("current", "spearman"),
                 stringsAsFactors = FALSE)
    }))
    per_task <- per_task[per_task$A != per_task$current | per_task$B != per_task$current, ]
    if (nrow(per_task)) {
      emit("Tasks whose choice A or B changes (history chosen by each rule; Spearman correlation over all domain states, and each rule's difference from the current rule):",
           "", table_md(per_task, digits = 3))
    }
  }
}

## Criteria --------------------------------------------------------------------
emit("## Criteria for rule B (agreed 2026-11-19)", "")
if (is.null(dev_b)) {
  emit("- Development benchmark: B changed no choice; criterion not applicable.")
} else {
  ok_rmse <- dev_b_overall <= 1.00 && all(dev_b$median_rmse_ratio <= 1.02)
  ok_cov <- all(is.na(dev_b$median_coverage_diff) | dev_b$median_coverage_diff >= -0.01)
  emit(sprintf("- Development benchmark, median RMSE ratio at most 1.00 over distances 0-2 and at most 1.02 at each: %s (overall %.3f; %s).",
               if (ok_rmse) "met" else "not met", dev_b_overall,
               paste(sprintf("%s %.3f", dev_b$distance, dev_b$median_rmse_ratio), collapse = ", ")),
       sprintf("- Development benchmark, median coverage at most 0.01 lower: %s.", if (ok_cov) "met" else "not met"))
}
if (is.null(abm_b)) {
  emit("- ABM (random set): B changed no choice; criterion not applicable.")
} else {
  all_row <- abm_b[abm_b$group == "all", ]
  emit(sprintf("- ABM (random set), median Spearman difference at least 0 (all domain states): %s (%.4f).",
               if (isTRUE(all_row$median_spearman_diff >= 0)) "met" else "not met", all_row$median_spearman_diff),
       sprintf("- ABM (random set), not significantly worse by a sign test: %s (better %d, worse %d, one-sided p %.3g).",
               if (is.na(all_row$sign_p_worse) || all_row$sign_p_worse >= 0.05) "met" else "not met",
               all_row$better, all_row$worse, all_row$sign_p_worse))
}
writeLines(lines, file.path(results, "report.md"))
cat(lines, sep = "\n")
