# Scoreboard aggregation.
#
# Datasets are averaged within landscape and scenario first, then landscapes
# are weighted equally. Paired comparisons resample landscapes, the independent
# unit; patients and scenarios of one landscape are not independent replicates.

stt_metric_columns <- c("rmse", "rmse_zero", "median_abs_error", "median_abs_error_zero",
                        "extreme_fraction", "bias", "error_sd", "truth_mean", "truth_sd",
                        "skill", "slope", "spearman", "auc_beneficial",
                        "edge_rmse", "edge_sign_accuracy", "coverage95", "interval_score95",
                        "prediction_fraction")

stt_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# A failed fit is scored as the null predictor, whose rank correlation is taken
# as zero (score.R, stt_rank_correlation()). Metrics tables scored before that
# rule have NA there instead, which would leave failed fits out of a method's
# mean Spearman correlation; this sets them to what the scoring now gives.
stt_null_rank <- function(metrics, null_methods = "zero") {
  no_ordering <- metrics$status != "ok" | metrics$method %in% null_methods | metrics$prediction_fraction == 0
  fix <- is.na(metrics$spearman) & no_ordering & metrics$n_nodes >= 3 & metrics$truth_sd > 0
  metrics$spearman[which(fix)] <- 0
  metrics
}

# One row per landscape x scenario x method x distance.
stt_landscape_scenario_means <- function(metrics) {
  key <- c("landscape_id", "wavelength_class", "scenario", "method", "stratum")
  stats::aggregate(metrics[stt_metric_columns], metrics[key], stt_mean)
}

stt_pooled_scoreboard <- function(cells, metrics) {
  per_landscape <- stats::aggregate(cells[stt_metric_columns],
                                    cells[c("landscape_id", "method", "stratum")], stt_mean)
  board <- stats::aggregate(per_landscape[stt_metric_columns],
                            per_landscape[c("method", "stratum")], stt_mean)
  per_fit <- metrics[!duplicated(metrics[c("dataset_id", "method")]), ]
  runtime <- data.frame(method = sort(unique(per_fit$method)), stringsAsFactors = FALSE)
  runtime$failure_fraction <- as.numeric(tapply(per_fit$status != "ok", per_fit$method, mean)[runtime$method])
  runtime$median_seconds <- as.numeric(tapply(per_fit$seconds, per_fit$method, stats::median)[runtime$method])
  merge(board, runtime, by = "method")
}

# Paired difference of a metric between a candidate and a baseline, by
# distance, with a landscape bootstrap interval. Only landscape x scenario
# cells that both methods have are used, so a method fitted on part of the
# benchmark is compared on that part.
stt_paired_comparison <- function(cells, candidate, baseline, metric = "rmse",
                                  n_boot = 2000L, seed = 20261110L) {
  key <- c("landscape_id", "scenario")
  rows <- lapply(c("d0", "d1", "d2"), function(d) {
    pick <- function(method) {
      part <- cells[cells$method == method & cells$stratum == d, c(key, metric)]
      part[is.finite(part[[metric]]), ]
    }
    both <- merge(pick(candidate), pick(baseline), by = key, suffixes = c(".candidate", ".baseline"))
    if (!nrow(both)) return(NULL)
    per_landscape <- stats::aggregate(both[paste0(metric, c(".candidate", ".baseline"))],
                                      both["landscape_id"], mean)
    if (nrow(per_landscape) < 3L) return(NULL)
    a <- per_landscape[[2L]]
    b <- per_landscape[[3L]]
    difference <- a - b
    set.seed(seed)
    boot <- replicate(n_boot, mean(sample(difference, replace = TRUE)))
    data.frame(
      candidate = candidate, baseline = baseline, metric = metric, stratum = d,
      n_landscapes = nrow(per_landscape),
      baseline_mean = mean(b), candidate_mean = mean(a),
      difference = mean(difference),
      ci_low = unname(stats::quantile(boot, 0.025)), ci_high = unname(stats::quantile(boot, 0.975)),
      relative_change = mean(difference) / mean(b),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Scenario x distance cells: candidate RMSE relative to the baseline.
stt_safety_cells <- function(cells, candidate, baseline, metric = "rmse") {
  cells <- cells[cells$stratum %in% c("d0", "d1", "d2"), ]
  by_scenario <- stats::aggregate(cells[metric], cells[c("scenario", "method", "stratum")], stt_mean)
  a <- by_scenario[by_scenario$method == candidate, c("scenario", "stratum", metric)]
  b <- by_scenario[by_scenario$method == baseline, c("scenario", "stratum", metric)]
  names(a)[3L] <- "candidate_value"
  names(b)[3L] <- "baseline_value"
  out <- merge(a, b, by = c("scenario", "stratum"))
  out$ratio <- out$candidate_value / out$baseline_value
  cbind(candidate = candidate, baseline = baseline, out, stringsAsFactors = FALSE)
}

stt_markdown_table <- function(x, digits = 3L) {
  x[] <- lapply(x, function(col) if (is.numeric(col)) formatC(col, digits = digits, format = "f") else as.character(col))
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  rule <- paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|")
  body <- apply(x, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, rule, body)
}

# The acceptance rule for a candidate, evaluated on the three distances: RMSE
# below the no-information answer and Spearman correlation above the original
# ALFA-K (`acceptance$alfak`), each with the landscape-bootstrap 95% interval
# of the paired difference excluding zero.
#
# Predicting no fitness differences is the bar for error, because the original
# ALFA-K is far above it (a few unconstrained fits dominate its RMSE). Rank
# agreement is compared with the original ALFA-K, since it is what it does
# carry and it does not depend on the reference.
#
# A criterion is evaluated only when the candidate and the baseline are
# compared on every landscape of the split (and on at least three). Otherwise,
# for instance when the baseline was not fitted on the split or on some of its
# landscapes, its status is "not evaluated" and `passed` is NA: it never counts
# as passed.
stt_acceptance <- function(cells, candidate, acceptance) {
  methods <- unique(cells$method)
  n_split <- length(unique(cells$landscape_id))
  criterion <- function(label, baseline, metric, lower_is_better) {
    paired <- if (all(c(candidate, baseline) %in% methods)) {
      stt_paired_comparison(cells, candidate, baseline, metric)
    }
    do.call(rbind, lapply(c("d0", "d1", "d2"), function(d) {
      row <- if (is.null(paired)) NULL else paired[paired$stratum == d, , drop = FALSE]
      found <- !is.null(row) && nrow(row) == 1L
      n <- if (found) row$n_landscapes else {
        common <- function(method) {
          unique(cells$landscape_id[cells$method == method & cells$stratum == d & is.finite(cells[[metric]])])
        }
        length(intersect(common(candidate), common(baseline)))
      }
      note <- if (!baseline %in% methods) {
        sprintf("%s not scored on this split", baseline)
      } else if (!candidate %in% methods) {
        sprintf("%s not scored on this split", candidate)
      } else if (!found || n < n_split) {
        sprintf("compared on %d of %d landscapes", n, n_split)
      } else NA_character_
      evaluated <- is.na(note)
      passed <- if (!evaluated) NA else if (lower_is_better) row$ci_high < 0 else row$ci_low > 0
      value <- function(column) if (found) row[[column]] else NA_real_
      data.frame(
        candidate = candidate, criterion = label, stratum = d, n_landscapes = n,
        baseline_mean = value("baseline_mean"), candidate_mean = value("candidate_mean"),
        difference = value("difference"), ci_low = value("ci_low"), ci_high = value("ci_high"),
        status = if (!evaluated) "not evaluated" else if (passed) "passed" else "failed",
        passed = passed, note = note,
        stringsAsFactors = FALSE
      )
    }))
  }
  rbind(
    criterion("RMSE below the no-information answer", acceptance$zero, "rmse", TRUE),
    criterion(paste("Spearman above", acceptance$alfak), acceptance$alfak, "spearman", FALSE)
  )
}

# One line per candidate: passes only when every criterion was evaluated and
# passed.
stt_acceptance_summary <- function(verdict) {
  vapply(unique(verdict$candidate), function(candidate) {
    part <- verdict[verdict$candidate == candidate, ]
    failed <- part[part$status == "failed", ]
    open <- part[part$status == "not evaluated", ]
    describe <- function(x) paste(sprintf("%s at %s", x$criterion, x$stratum), collapse = "; ")
    verdict_text <- if (nrow(failed)) {
      sprintf("fails %d of %d criteria (%s)", nrow(failed), nrow(part), describe(failed))
    } else if (nrow(open)) {
      sprintf("passes %d of %d criteria; %d not evaluated (%s: %s), so it does not pass the rule",
              sum(part$status == "passed"), nrow(part), nrow(open), describe(open),
              paste(unique(open$note), collapse = "; "))
    } else {
      sprintf("passes every criterion (%d of %d)", nrow(part), nrow(part))
    }
    sprintf("- `%s` %s.", candidate, verdict_text)
  }, character(1), USE.NAMES = FALSE)
}

stt_stage_report <- function(config, split) {
  metrics <- utils::read.csv(file.path(config$output_dir, sprintf("metrics_%s.csv", split)),
                             stringsAsFactors = FALSE)
  # Methods scored in earlier rounds but no longer configured are left out.
  metrics <- metrics[metrics$method %in% names(config$methods), ]
  metrics <- stt_null_rank(metrics, config$acceptance$zero)
  # The verdict is made on the core scenarios. Stress scenarios break an
  # assumption on purpose and are tabulated after it.
  stress <- metrics[metrics$scenario %in% config$scenarios$scenario[config$scenarios$role == "stress"], ]
  metrics <- metrics[!metrics$scenario %in% stress$scenario, ]
  cells <- stt_landscape_scenario_means(metrics)
  board <- stt_pooled_scoreboard(cells, metrics)
  # A method scored on fewer datasets is averaged over fewer of them, so its
  # pooled row is not comparable; the paired tables use common cells only.
  board$n_datasets <- as.integer(tapply(metrics$dataset_id, metrics$method,
                                        function(x) length(unique(x)))[board$method])
  board <- board[order(board$stratum, board$method), ]
  methods <- unique(metrics$method)
  acceptance <- config$acceptance
  # Methods that are neither baselines nor named candidates (sensitivity
  # variants, for instance) are tabulated but not judged.
  candidates <- intersect(acceptance$candidates, methods)
  # Every method against the no-information answer and against the original
  # ALFA-K settings.
  pairs <- list()
  for (method in setdiff(methods, acceptance$zero)) {
    if (acceptance$zero %in% methods) pairs[[length(pairs) + 1L]] <- c(method, acceptance$zero)
    if (!method %in% acceptance$original) {
      for (original in intersect(acceptance$original, methods)) pairs[[length(pairs) + 1L]] <- c(method, original)
    }
  }
  # RMSE and median absolute error: lower is better. Spearman and AUC: higher
  # is better, and unlike the errors they do not depend on the reference.
  paired <- do.call(rbind, lapply(pairs, function(pair) {
    do.call(rbind, lapply(c("rmse", "median_abs_error", "spearman", "auc_beneficial"), function(metric) {
      stt_paired_comparison(cells, pair[[1L]], pair[[2L]], metric)
    }))
  }))
  if (is.null(paired)) stop("No paired comparison is available yet.", call. = FALSE)
  verdict <- do.call(rbind, lapply(candidates, function(candidate) stt_acceptance(cells, candidate, acceptance)))
  # Cohort methods against the single-dataset model of the same initial state
  # (and the cohort with balance also against the cohort without it), on the
  # cohort scenario.
  cohort_cells <- cells[cells$scenario == config$cohort_scenario, ]
  cohort_paired <- do.call(rbind, lapply(intersect(names(acceptance$cohort), methods), function(method) {
    do.call(rbind, lapply(intersect(acceptance$cohort[[method]], methods), function(baseline) {
      do.call(rbind, lapply(c("rmse", "bias", "spearman", "coverage95"), function(metric) {
        stt_paired_comparison(cohort_cells, method, baseline, metric)
      }))
    }))
  }))
  report_dir <- file.path(config$output_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cells, file.path(report_dir, sprintf("cells_%s.csv", split)), row.names = FALSE)
  utils::write.csv(board, file.path(report_dir, sprintf("scoreboard_%s.csv", split)), row.names = FALSE)
  utils::write.csv(paired, file.path(report_dir, sprintf("paired_%s.csv", split)), row.names = FALSE)
  if (!is.null(verdict)) {
    utils::write.csv(verdict, file.path(report_dir, sprintf("acceptance_%s.csv", split)), row.names = FALSE)
  }
  stress_table <- if (nrow(stress)) {
    part <- stress[stress$stratum %in% c("d0", "d1", "d2"), ]
    table <- stats::aggregate(part[c("rmse", "bias", "spearman", "coverage95")],
                              part[c("scenario", "stratum", "method")], stt_mean)
    table$n_datasets <- stats::aggregate(part$dataset_id, part[c("scenario", "stratum", "method")],
                                         function(x) length(unique(x)))$x
    utils::write.csv(table, file.path(report_dir, sprintf("stress_%s.csv", split)), row.names = FALSE)
    table[order(table$scenario, table$stratum, table$method), ]
  } else NULL
  show <- c("stratum", "method", "n_datasets", "rmse", "rmse_zero", "truth_sd", "median_abs_error",
            "slope", "spearman", "auc_beneficial", "edge_sign_accuracy", "coverage95",
            "extreme_fraction", "failure_fraction", "median_seconds")
  lines <- c(
    sprintf("# Sparse two-timepoint scoreboard (%s landscapes)", split), "",
    sprintf("%d datasets, %d landscapes, %d core scenarios. `rmse_zero` is the error of predicting no fitness differences; `truth_sd` is the error of knowing only the mean fitness of the stratum.",
            length(unique(metrics$dataset_id)), length(unique(metrics$landscape_id)), length(unique(metrics$scenario))),
    "", "## Pooled over scenarios", "", stt_markdown_table(board[show]),
    "", "## Acceptance", "",
    sprintf(paste("Rule: at distances 0, 1 and 2, RMSE below `%s` and Spearman correlation above `%s`, each with",
                  "the landscape-bootstrap 95%% interval of the paired difference excluding zero. A criterion whose",
                  "baseline is missing, or that compares on fewer than all %d landscapes, is not evaluated and does",
                  "not count as passed."),
            acceptance$zero, acceptance$alfak, length(unique(metrics$landscape_id))),
    "",
    if (is.null(verdict)) "No candidate method has been scored yet." else {
      c(stt_acceptance_summary(verdict), "",
        stt_markdown_table(verdict, 4L))
    },
    "", "## Paired differences (candidate minus baseline, landscape bootstrap 95% interval)", "",
    stt_markdown_table(paired[c("candidate", "baseline", "metric", "stratum", "n_landscapes", "baseline_mean",
                                "candidate_mean", "difference", "ci_low", "ci_high", "relative_change")], 4L),
    if (!is.null(cohort_paired)) {
      utils::write.csv(cohort_paired, file.path(report_dir, sprintf("cohort_%s.csv", split)), row.names = FALSE)
      c("", sprintf("## Cohort fits against the single-dataset fit of the same initial state (scenario %s; failed fits scored as the no-information answer)", config$cohort_scenario), "",
        stt_markdown_table(cohort_paired[c("candidate", "baseline", "metric", "stratum", "n_landscapes", "baseline_mean",
                                           "candidate_mean", "difference", "ci_low", "ci_high")], 4L))
    },
    if (!is.null(stress_table)) c("", "## Stress scenarios (not part of the verdict)", "", stt_markdown_table(stress_table))
  )
  path <- file.path(report_dir, sprintf("scoreboard_%s.md", split))
  writeLines(lines, path)
  cat(lines, sep = "\n")
  invisible(path)
}
