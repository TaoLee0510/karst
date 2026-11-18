# Run time and peak memory of every fit, one row per input and method.
#
#   Rscript scripts/collect_resources.R [METHOD=REGISTRY ...]
#
# METHOD=REGISTRY names the job registry whose tasks ran each method's fits
# in this round's results (default alfak=8327282, the round alfak's fits
# are linked from, and karst=jobs, this round's BENCH_RESULTS/jobs.tsv);
# tasks of the other method in a registry are ignored. karst's other
# multi-time-point methods (karst_pairwise, karst_state_space) are counted as
# karst: same registry, one fit per input.
#
# Reads run.Rds under BENCH_RESULTS/fits (fit time; karst's peak memory up to
# the end of its fit, from the multi-timepoint round on) and
# BENCH_RESULTS/tables/sacct.tsv (hpc/sacct_jobs.sh: each array task's peak
# memory from Slurm's accounting) and writes BENCH_RESULTS/tables/
# resources.Rds:
#   source, dataset, design, method  the input and the method
#   minobs_fits, failed              minobs folders and how many did not end
#                                    with status "ok" (a killed task leaves
#                                    no run.Rds: status "missing")
#   fit_seconds                      time to produce alfak's outputs for every
#                                    minobs: alfak, the sum of its fits (one
#                                    per minobs); karst, its one fit plus the
#                                    outputs derived for each minobs. Forecasts
#                                    and scoring are not included.
#   task_peak_mb                     the largest peak memory of the input's
#                                    array tasks (alfak: one task per minobs;
#                                    karst: one task), over every attempt; a
#                                    task includes its forecasts and scoring
#   tasks, attempts                  array tasks found in the accounting and
#                                    their attempts (a task already done is
#                                    skipped and its attempt uses little
#                                    memory, hence the maximum)
#   karst_fit_peak_mb               karst only: peak memory of the process up
#                                    to the end of the fit (input and fit)
#   process_peak_mb                  the task process's own peak (fit,
#                                    forecasts, scores; task_peak.Rds, from
#                                    the multi-timepoint round on)
#   peak_mb                          the largest of the three: Slurm samples
#                                    memory every 30 s or so and misses tasks
#                                    shorter than that (they show the batch
#                                    shell's few MB)

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
tab <- file.path(paths$results, "tables")
bench_dir <- Sys.getenv("BENCH_DIR", getwd())
registry_of <- c(alfak = "8327282", karst = "jobs")
for (a in commandArgs(TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  stopifnot(length(kv) == 2L, kv[1] %in% names(registry_of))
  registry_of[[kv[1]]] <- kv[2]
}

## Fit time ----------------------------------------------------------------------
dirs <- Sys.glob(file.path(paths$results, "fits", "*", "*", "*", "*", "minobs_*"))
runs <- do.call(rbind, lapply(dirs, function(d) {
  p <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
  path <- file.path(d, "run.Rds")
  r <- if (file.exists(path)) readRDS(path) else list(status = "missing", seconds = NA_real_)
  tp <- file.path(d, "task_peak.Rds")
  data.frame(source = p[5], dataset = p[4], design = p[3], method = p[2], minobs = sub("minobs_", "", p[1]),
             status = r$status, seconds = first_or_na(r$seconds), fit_seconds = first_or_na(r$fit_seconds),
             fit_peak_mb = first_or_na(r$fit_peak_mb),
             process_peak_mb = if (file.exists(tp)) first_or_na(readRDS(tp)$peak_mb) else NA_real_,
             stringsAsFactors = FALSE)
}))
max_or_na <- function(v) if (any(is.finite(v))) max(v, na.rm = TRUE) else NA_real_
key <- c("source", "dataset", "design", "method")
per_input <- do.call(rbind, lapply(split(runs, runs[key], drop = TRUE), function(x) {
  # karst: one fit shared by the minobs folders; a failed fit records its
  # time as `seconds` of every folder and no `fit_seconds`.
  shared <- x$fit_seconds[!is.na(x$fit_seconds)]
  t <- if (is_karst_method(x$method[1])) {
    if (length(shared)) shared[1] + sum(x$seconds[!is.na(x$fit_seconds)], na.rm = TRUE) else x$seconds[1]
  } else sum(x$seconds, na.rm = TRUE)
  data.frame(x[1, key], minobs_fits = nrow(x), failed = sum(x$status != "ok"),
             missing = sum(x$status == "missing"), fit_seconds = t,
             karst_fit_peak_mb = if (is_karst_method(x$method[1])) max_or_na(x$fit_peak_mb) else NA_real_,
             process_peak_mb = max_or_na(x$process_peak_mb), stringsAsFactors = FALSE)
}))

## Peak memory ---------------------------------------------------------------------
acct_path <- file.path(tab, "sacct.tsv")
if (file.exists(acct_path)) {
  acct <- utils::read.delim(acct_path, stringsAsFactors = FALSE, colClasses = c(registry = "character", array_index = "character"))
  acct <- acct[nzchar(acct$array_index), ]
  task_tables <- lapply(stats::setNames(nm = unique(acct$task_file)), function(tf) {
    path <- if (startsWith(tf, "/")) tf else file.path(bench_dir, tf)
    utils::read.delim(path, stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
  })
  acct <- do.call(rbind, lapply(split(acct, acct$task_file), function(a) {
    tt <- task_tables[[a$task_file[1]]]
    row <- tt[match(as.integer(a$array_index), tt$task), c("source", "dataset", "design", "method", "minobs")]
    # task tables of the alfak2 era (hpc/jobs/8327282.tsv points to them, by
    # absolute path) name alfak's tasks after the package
    row$method[row$method %in% "alfakR"] <- "alfak"
    cbind(a, row)
  }))
  method_registry <- unname(registry_of[ifelse(is_karst_method(acct$method), "karst", acct$method)])
  acct <- acct[!is.na(acct$method) & acct$registry == method_registry, ]
  acct$task_id <- paste(acct$task_file, acct$array_index)
  mem <- do.call(rbind, lapply(split(acct, acct[key], drop = TRUE), function(x) data.frame(
    x[1, key], task_peak_mb = max(x$maxrss_kb, na.rm = TRUE) / 1024, tasks = length(unique(paste(x$minobs))),
    attempts = nrow(x), oom_attempts = sum(x$state == "OUT_OF_MEMORY"), stringsAsFactors = FALSE)))
  per_input <- merge(per_input, mem, by = key, all.x = TRUE)
} else {
  message("No ", acct_path, ": peak memory from the accounting is left out (run hpc/sacct_jobs.sh).")
  per_input$task_peak_mb <- NA_real_
  per_input$tasks <- NA_integer_
  per_input$attempts <- NA_integer_
  per_input$oom_attempts <- NA_integer_
}
# The peak used: the largest of the accounting's, the task process's own and
# the fit's (the accounting misses tasks shorter than its sampling interval).
per_input$peak_mb <- apply(per_input[c("task_peak_mb", "process_peak_mb", "karst_fit_peak_mb")], 1, max_or_na)
saveRDS(per_input, file.path(tab, "resources.Rds"))
cat("inputs:", nrow(per_input), " with accounting:", sum(!is.na(per_input$task_peak_mb)),
    " with the process's own peak:", sum(!is.na(per_input$process_peak_mb)),
    " below 100 MB (likely missed by the accounting):", sum(per_input$peak_mb < 100, na.rm = TRUE), "\n")
print(stats::aggregate(cbind(fit_seconds, peak_mb) ~ source + method + design, per_input,
                       function(v) round(stats::median(v), 1), na.action = stats::na.pass))
