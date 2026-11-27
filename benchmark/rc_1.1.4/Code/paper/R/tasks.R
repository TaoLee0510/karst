# The round's task tables and the fit folders they promise.
#
# A task killed for time or memory before its fit ends leaves no run.Rds.
# karst writes its minobs folders only after the fit (R/fit.R), so a killed
# karst task leaves no folder at all, while alfak's folder appears with its
# first output: without a check against the task tables, karst's killed fits
# would vanish from the tables and alfak's would count as failed. The
# collect steps (scripts/collect_abm.R, collect_salehi.R) therefore check
# every fit folder of the round's task tables first. A missing output stops
# them; with ALLOW_MISSING=1 they list it (tables/missing_<name>.tsv) and go
# on, a fit without run.Rds then counting as a failed fit (status
# "missing").
#   BENCH_TASKS  the round's task folder (scripts/make_tasks.R), default
#                tasks/<commit> for results_<commit> and for the folder of a
#                variant of the default, results_<commit>_free (R/fit.R),
#                relative to this folder; "none" skips the check (a results
#                folder whose fits are all linked)
# A variant's fits (karst_free) go to a results folder of their own: the
# check of that folder takes the variant's tasks only, the check of the
# round's folder every other task.

# The variant of the default whose results folder BENCH_RESULTS is (e.g.
# "karst_free" for results_<commit>_free), NA for the round's own folder.
bench_variant <- function(paths = bench_paths()) {
  v <- names(karst_variant_methods)
  hit <- v[endsWith(basename(paths$results), sub("^karst", "", v))]
  if (length(hit)) hit[[1L]] else NA_character_
}

bench_task_dir <- function(paths = bench_paths()) {
  dir <- Sys.getenv("BENCH_TASKS", "")
  if (nzchar(dir)) return(dir)
  commit <- sub("^results_", "", basename(paths$results))
  variant <- bench_variant(paths)
  if (!is.na(variant)) commit <- sub(paste0(sub("^karst", "", variant), "$"), "", commit)
  file.path("tasks", commit)
}

# One row per fit folder the tables tasks_<source>.tsv of the round expect
# (source, dataset, design, method, minobs, dir). NULL with BENCH_TASKS=none.
# A round that refits a check subset has only that subset in its tables; the
# folders it links from another round are not checked.
round_expected_fits <- function(sources, paths = bench_paths()) {
  dir <- bench_task_dir(paths)
  if (identical(dir, "none")) return(NULL)
  files <- file.path(dir, paste0("tasks_", sources, ".tsv"))
  files <- files[file.exists(files)]
  if (!length(files)) {
    stop("No task table tasks_", paste(sources, collapse = "/tasks_"), ".tsv in ", dir, ": set BENCH_TASKS to the ",
         "round's task folder, or BENCH_TASKS=none to collect without the check")
  }
  tasks <- do.call(rbind, lapply(files, function(f) {
    utils::read.delim(f, stringsAsFactors = FALSE, colClasses = c(minobs = "character"))[
      c("source", "dataset", "design", "method", "minobs")]
  }))
  variant <- bench_variant(paths)
  tasks <- if (is.na(variant)) tasks[!tasks$method %in% names(karst_variant_methods), , drop = FALSE] else
    tasks[tasks$method == variant, , drop = FALSE]
  m <- strsplit(tasks$minobs, ",", fixed = TRUE)
  rows <- tasks[rep(seq_len(nrow(tasks)), lengths(m)), c("source", "dataset", "design", "method")]
  rows$minobs <- unlist(m)
  rows <- unique(rows)
  rownames(rows) <- NULL
  rows$dir <- file.path(paths$results, "fits", rows$source, rows$dataset, rows$design, rows$method,
                        paste0("minobs_", rows$minobs))
  rows
}

# Checks the outputs of the expected fit folders: run.Rds, and the files
# `need(row, run)` names as missing from a folder whose run.Rds exists (`run`
# its content). Writes tables/missing_<name>.tsv when anything is missing
# and stops, unless ALLOW_MISSING=1. Returns the expected fit folders
# without run.Rds (for the caller's run table).
check_round_outputs <- function(expected, name, need = function(row, run) character(), paths = bench_paths()) {
  if (is.null(expected)) {
    cat("Task check skipped (BENCH_TASKS=none)\n")
    return(NULL)
  }
  missing <- lapply(seq_len(nrow(expected)), function(i) {
    run_path <- file.path(expected$dir[i], "run.Rds")
    if (!file.exists(run_path)) return("run.Rds")
    need(expected[i, ], readRDS(run_path))
  })
  n <- lengths(missing)
  out <- expected[n > 0, , drop = FALSE]
  out$missing <- vapply(missing[n > 0], paste, character(1), collapse = ",")
  path <- file.path(paths$results, "tables", paste0("missing_", name, ".tsv"))
  if (nrow(out)) {
    utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (file.exists(path)) {
    file.remove(path)
  }
  cat("Task check (", bench_task_dir(paths), "): ", nrow(expected), " fit folders expected, ", nrow(out),
      " with outputs missing\n", sep = "")
  if (nrow(out)) {
    print(utils::head(out[setdiff(names(out), "dir")], 20), row.names = FALSE)
    if (Sys.getenv("ALLOW_MISSING", "0") != "1") {
      stop(nrow(out), " fit folders of the round's tasks have outputs missing (", path, "): rerun their tasks, or ",
           "set ALLOW_MISSING=1 to collect a fit without run.Rds as failed")
    }
  }
  out[grepl("(^|,)run[.]Rds", out$missing), , drop = FALSE]
}
