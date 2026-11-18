# karst's missegregation rate per input (round ac36fd1): the candidate of
# lowest objective, as fit_sparse_joint(beta = "auto") picks it.
#
#   Rscript scripts/beta_picks.R TASK_DIR
#
# TASK_DIR holds beta_tasks.tsv (scripts/make_beta_tasks.R). Each candidate is
# read from BENCH_RESULTS/candidates/<source>/<dataset>/<design>/<label>/, the
# 5e-5 candidate of an input without its own task from fits/<...>/karst (the
# reference round's fit, linked in). The objective is the fit's negative log
# marginal likelihood (run.Rds `objective`; for a reused fit, which predates
# that field, the lowest finite objective of its history profile, which is
# the chosen history's: scripts/beta_check.R checks this on refits). As
# karst's sparse_joint_choose_candidate(): the lowest finite objective, the
# first in the grid's order (1/16 ... 64) on a tie, whether or not its
# optimisation was accepted. A fit that failed has no objective; an input
# whose candidates all failed has no pick (karst stops: a failed fit).
#
# Writes to BENCH_RESULTS/tables/:
#   beta_candidates.tsv  one row per input and candidate: status (fit ok,
#                        failed, missing = no run.Rds), objective, history,
#                        convergence, fit and output seconds, peak memory
#   beta_picks.tsv       one row per input: the data's rate (pm_data; the
#                        true rate of a missegregation simulation), the
#                        Salehi context (PDX or in vitro: lines SA039 and
#                        SA906), candidates fitted / failed / missing, the
#                        pick (label, multiple, rate), the time beta = "auto"
#                        takes (every candidate's fit plus the picked one's
#                        outputs) and whether the 5e-5 fit was reused
#   beta_candidate_scores.tsv  ABM and missegregation: landscape scores of
#                        every candidate and minobs (score_landscape.Rds or
#                        score.Rds; Pearson r, Spearman rho, rescaled R^2)

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
tasks <- utils::read.delim(file.path(args[[1]], "beta_tasks.tsv"), stringsAsFactors = FALSE,
                           colClasses = c(minobs = "character"))
key <- c("source", "dataset", "design")
inputs <- unique(tasks[c(key, "minobs", "mem")])
has_x1 <- unique(tasks[tasks$beta_label == "x1", key])
inputs$reuse_x1 <- !paste(inputs$source, inputs$dataset, inputs$design) %in%
  paste(has_x1$source, has_x1$dataset, has_x1$design)
cand_dir <- function(i, label) {
  if (label == "x1" && inputs$reuse_x1[i]) {
    return(file.path(paths$results, "fits", inputs$source[i], inputs$dataset[i], inputs$design[i], "karst"))
  }
  file.path(paths$results, "candidates", inputs$source[i], inputs$dataset[i], inputs$design[i], label)
}
labels <- beta_label(karst_beta_multiples)

## Candidates -------------------------------------------------------------------------------
read_candidate <- function(i, k) {
  label <- labels[k]
  dir <- cand_dir(i, label)
  minobs <- strsplit(inputs$minobs[i], ",", fixed = TRUE)[[1]]
  runs <- lapply(minobs, function(m) {
    p <- file.path(dir, paste0("minobs_", m), "run.Rds")
    if (file.exists(p)) readRDS(p) else NULL
  })
  have <- !vapply(runs, is.null, logical(1))
  base <- data.frame(source = inputs$source[i], dataset = inputs$dataset[i], design = inputs$design[i], label = label,
                     multiple = karst_beta_multiples[k], beta = karst_beta_default * karst_beta_multiples[k],
                     reused = label == "x1" && inputs$reuse_x1[i], stringsAsFactors = FALSE)
  if (!any(have)) {
    return(cbind(base, status = "missing", outputs_ok = 0L, objective = NA_real_, objective_from = NA_character_,
                 history = NA_real_, convergence = NA_real_, fit_seconds = NA_real_, output_seconds = NA_real_,
                 fit_peak_mb = NA_real_, task_peak_mb = NA_real_, condition_class = NA_character_))
  }
  r <- runs[have][[1]]
  fit_ok <- !is.null(r$fit_seconds)
  objective <- r$objective
  from <- "run"
  if (fit_ok && is.null(objective)) {
    prof <- r$history_profile$objective
    objective <- if (length(prof) && any(is.finite(prof))) min(prof[is.finite(prof)]) else NA_real_
    from <- "history_profile"
  }
  peaks <- vapply(minobs, function(m) {
    p <- file.path(dir, paste0("minobs_", m), "task_peak.Rds")
    if (file.exists(p)) readRDS(p)$peak_mb else NA_real_
  }, numeric(1))
  cbind(base, status = if (fit_ok) "ok" else "failed",
        outputs_ok = sum(vapply(runs[have], function(x) identical(x$status, "ok"), logical(1))),
        objective = if (fit_ok) first_or_na(objective) else NA_real_, objective_from = if (fit_ok) from else NA_character_,
        history = first_or_na(r$history), convergence = first_or_na(r$convergence),
        fit_seconds = if (fit_ok) r$fit_seconds else r$seconds,
        output_seconds = if (fit_ok) sum(vapply(runs[have], function(x) x$seconds, numeric(1))) else 0,
        fit_peak_mb = first_or_na(r$fit_peak_mb), task_peak_mb = if (any(is.finite(peaks))) max(peaks, na.rm = TRUE) else NA_real_,
        condition_class = first_or_na(r$condition_class, NA_character_))
}
cand <- do.call(rbind, lapply(seq_len(nrow(inputs)), function(i) {
  do.call(rbind, lapply(seq_along(labels), function(k) read_candidate(i, k)))
}))

## Picks -----------------------------------------------------------------------------------
meta <- utils::read.csv(paths$salehi_meta)
pdx_of <- stats::setNames(meta$PDX_id, paste(meta$datasetname, meta$timepoint, sep = "_"))
pm_data <- function(source, dataset) {
  pm <- ifelse(source == "abm", bench_settings$abm$pm, bench_settings$salehi$pm)
  mis <- source == "misseg"
  pm[mis] <- vapply(dataset[mis], function(d) abm_sim_info(d)$pmis, numeric(1))
  pm
}
picks <- do.call(rbind, lapply(split(cand, cand[key], drop = TRUE), function(d) {
  d <- d[order(d$multiple), ]
  usable <- which(is.finite(d$objective))
  k <- if (length(usable)) usable[[which.min(d$objective[usable])]] else NA_integer_
  chosen <- if (is.na(k)) d[d$label == "x1", ] else d[k, ]
  data.frame(d[1, key], candidates_ok = sum(d$status == "ok"), candidates_failed = sum(d$status == "failed"),
             candidates_missing = sum(d$status == "missing"), reused_x1 = any(d$reused),
             pick_label = if (is.na(k)) NA_character_ else d$label[k],
             pick_multiple = if (is.na(k)) NA_real_ else d$multiple[k], pick_beta = if (is.na(k)) NA_real_ else d$beta[k],
             pick_objective = if (is.na(k)) NA_real_ else d$objective[k],
             x1_objective = d$objective[d$label == "x1"],
             pick_outputs_ok = chosen$outputs_ok,
             auto_seconds = sum(d$fit_seconds, na.rm = TRUE) + chosen$output_seconds,
             candidate_peak_mb_max = suppressWarnings(max(c(d$fit_peak_mb, d$task_peak_mb), na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
picks$candidate_peak_mb_max[!is.finite(picks$candidate_peak_mb_max)] <- NA
picks$pm_data <- pm_data(picks$source, picks$dataset)
picks$context <- ifelse(picks$source != "salehi", NA_character_,
                        ifelse(pdx_of[sub("_l_.*$", "", picks$dataset)] %in% c("SA039", "SA906"), "in vitro", "PDX"))
picks$group <- ifelse(grepl("all$", picks$design), "all time points", "2 time points")
picks <- picks[order(picks$source, picks$dataset, picks$design), ]
out <- file.path(paths$results, "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.table(cand, file.path(out, "beta_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(picks, file.path(out, "beta_picks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Landscape scores of every candidate (ABM, missegregation) ----------------------------------
scores <- do.call(rbind, lapply(which(inputs$source %in% c("abm", "misseg")), function(i) {
  do.call(rbind, lapply(labels, function(label) {
    dir <- cand_dir(i, label)
    do.call(rbind, lapply(strsplit(inputs$minobs[i], ",", fixed = TRUE)[[1]], function(m) {
      d <- file.path(dir, paste0("minobs_", m))
      p <- c(file.path(d, "score_landscape.Rds"), file.path(d, "score.Rds"))
      p <- p[file.exists(p)]
      if (!length(p)) return(NULL)
      s <- readRDS(p[1])$summary
      data.frame(inputs[i, key], label = label, minobs = m, r = s$r, rho = s$rho, R = s$R, Rxv = s$Rxv,
                 stringsAsFactors = FALSE)
    }))
  }))
}))
if (!is.null(scores)) {
  utils::write.table(scores, file.path(out, "beta_candidate_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

## Console ---------------------------------------------------------------------------------
options(width = 200)
cat("inputs:", nrow(picks), " candidates:", nrow(cand), "\n")
print(table(status = cand$status, label = factor(cand$label, labels)))
cat("objective from the history profile (reused 5e-5 fits):", sum(cand$objective_from %in% "history_profile"), "\n")
cat("inputs with a missing candidate:", sum(picks$candidates_missing > 0), "  with no pick (every candidate failed):",
    sum(is.na(picks$pick_label)), "\n")
lv <- factor(picks$pick_label, labels)
cat("\n== Picks by source and design group (multiple of 5e-5)\n")
print(table(paste(picks$source, picks$group), lv))
m <- picks[picks$source == "misseg" & !is.na(picks$pick_beta), ]
if (nrow(m)) {
  cat("\n== Missegregation simulations: pick by true rate\n")
  print(table(true_rate = m$pm_data, pick = factor(m$pick_label, labels)))
  ratio <- abs(log2(m$pick_beta / m$pm_data))
  cat("within 4x of the true rate:", round(mean(ratio <= 2 + 1e-9), 3), "  median |log2(pick / true)|:",
      round(stats::median(ratio), 2), "\n")
}
sal <- picks[picks$source == "salehi", ]
if (nrow(sal)) {
  cat("\n== Salehi: pick by context and design\n")
  print(table(paste(sal$context, sal$design), factor(sal$pick_label, labels)))
}
cat("\n== Time of beta = \"auto\" per input (s): median by source and group\n")
print(stats::aggregate(auto_seconds ~ source + group, picks, function(v) round(stats::median(v))))
