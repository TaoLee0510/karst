# Shared helpers of the dev autorun (dev/PLAN.md). Sourced by every dev script.
#
# Paths: DEV_REPO is the karst working tree that holds dev/ (local clone or
# the HPC copy in /share/lab_crd/taoli/Project/karst_dev/repo); DEV_RESULTS
# holds results (dev/results locally, git-ignored; /share/.../karst_dev/results
# on HPC); DEV_INPUTS holds the read-only inputs of the e sets.

dev_repo <- local({
  env <- Sys.getenv("DEV_REPO", "")
  if (nzchar(env)) normalizePath(env) else {
    here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
    candidate <- normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(candidate, "dev", "PLAN.md"))) candidate else normalizePath(getwd())
  }
})
dev_results <- Sys.getenv("DEV_RESULTS", file.path(dev_repo, "dev", "results"))
dev_inputs <- Sys.getenv("DEV_INPUTS", file.path(dev_results, "inputs"))
# The alfak2-era inputs of e1-e4 (read only). Locally only what was copied.
dev_backup <- Sys.getenv("DEV_BACKUP",
  file.path(Sys.getenv("RC_HOME"), "data", "dev", "backup"))  # rc freeze: was the alfak2-era backup

# The validation scripts, sourced unchanged: series_methods/setup.R sources
# multi_timepoint/setup.R, which sources short_interval/setup.R and the
# development benchmark (stt_*).
dev_source_validation <- function() {
  if (exists("sm_config", mode = "function")) return(invisible())
  sm_script_dir <<- file.path(dev_repo, "validation", "sparse_joint", "series_methods")
  mt_script_dir <<- file.path(dev_repo, "validation", "sparse_joint", "multi_timepoint")
  si_script_dir <<- file.path(dev_repo, "validation", "sparse_joint", "short_interval")
  source(file.path(sm_script_dir, "setup.R"))
  invisible()
}

# Folder of one e set's inputs: datasets/, truth/, dataset_index.csv.
dev_set_inputs <- function(set) {
  local <- file.path(dev_inputs, set)
  if (dir.exists(local)) return(local)
  if (set == "e1") file.path(dev_backup, "multi_timepoint", "results", "series") else
    file.path(dev_backup, "series_methods", "results", "inputs", set)
}

# Truth folder of a set: e2 and e3 share e1's truth populations (their own
# truth/ is a link to a folder of the alfak2 era that no longer exists).
dev_truth_dir <- function(set) {
  if (set %in% c("e2", "e3")) file.path(dev_set_inputs("e1"), "truth") else file.path(dev_set_inputs(set), "truth")
}

# The generating configuration of a set (abm settings, design), with its
# output folder moved to dev results so that nothing is written to inputs.
dev_set_config <- function(set) {
  dev_source_validation()
  config <- if (set == "e1") mt_config(file.path(dev_results, "abm_cache", "e1"), "series") else {
    cfg <- sm_config(set)
    cfg$output_dir <- file.path(dev_results, "abm_cache", set)
    cfg
  }
  config$input_dir <- dev_set_inputs(set)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  config
}

dev_load_abm <- function() {
  if (exists("dev_abm_cpp", mode = "function")) return(invisible())
  # one cache per R process: tasks compiling into a shared folder on one node
  # collide ("no such symbol sourceCpp_..._dev_abm_cpp")
  cache <- file.path(tempdir(), "dev_cpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(dev_repo, "dev", "abm", "abm_dev.cpp"), cacheDir = cache, env = globalenv())
  invisible()
}

# Run the dev ABM from `initial` (named counts). Returns the populations at the
# recorded days and the passages (day, cells before, survival, cells after).
dev_run_abm <- function(initial, abm, record_days, seed, beta, culling_survival,
                        centroids = matrix(0, 1, length(strsplit(names(initial)[1], ".", fixed = TRUE)[[1]])),
                        wavelength = 1, constant_fitness = 0, landscape_type = 0L,
                        landscape_par = numeric(), landscape_ref = integer()) {
  dev_load_abm()
  record_steps <- as.integer(round(sort(unique(record_days)) / abm$delta_t))
  raw <- dev_abm_cpp(
    initial_karyotypes = matrix(as.integer(do.call(rbind, strsplit(names(initial), ".", fixed = TRUE))),
                                nrow = length(initial)),
    initial_counts = as.numeric(initial), centroids = centroids, wavelength = wavelength, beta = beta,
    delta_t = abm$delta_t, n_steps = max(record_steps), record_steps = record_steps,
    max_population = abm$max_pop, culling_survival = culling_survival,
    min_cn = abm$min_cn, max_cn = abm$max_cn, seed = seed, constant_fitness = constant_fitness,
    landscape_type = as.integer(landscape_type), landscape_par = as.numeric(landscape_par),
    landscape_ref = as.integer(landscape_ref))
  populations <- lapply(stats::setNames(record_steps, as.character(record_steps * abm$delta_t)), function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
  passages <- as.data.frame(raw$passages)
  passages$day <- passages$step * abm$delta_t
  list(populations = populations, passages = passages)
}

# The window of a dataset of e1-e4 replayed with the dev ABM (same seed, so the
# same trajectory), with its passages. `seed_offset` != 0 gives an independent
# replicate of the same window from the same start.
dev_window_item <- function(config, dataset) {
  plan <- if (config$set %in% c("series", "heldout")) mt_window_plan(config) else sm_window_plan(config)
  for (item in plan) {
    same <- item$landscape$landscape_id == dataset$landscape_id && item$bottleneck == dataset$bottleneck &&
      item$spacing == dataset$spacing && (is.null(item$start) || item$start == dataset$t0)
    if (same) return(item)
  }
  stop("No window of the plan matches ", dataset$dataset_id, call. = FALSE)
}

dev_replay_window <- function(config, dataset, truth, seed_offset = 0, record_index = NULL) {
  item <- dev_window_item(config, dataset)
  days <- dataset$t0 + item$spacing * (seq_len(config$max_samples) - 1L)
  start <- stt_snapshot(truth, dataset$t0)
  constant <- 0
  if (identical(config$set, "e3")) {
    fitness <- stt_grf_fitness(names(start), truth$centroids, truth$wavelength)
    constant <- sum(start * pmax(fitness, 0)) / sum(start)
  }
  record <- days[-1L] - dataset$t0
  if (!is.null(record_index)) record <- record[record_index - 1L]
  run <- dev_run_abm(start, config$abm, record_days = record, seed = item$seed + seed_offset,
                     beta = truth$beta, culling_survival = item$bottleneck / config$nb,
                     centroids = truth$centroids, wavelength = truth$wavelength, constant_fitness = constant)
  list(start = start, populations = run$populations, passages = run$passages, item = item)
}

# Passages between the first and the last time point of a k-sample dataset:
# the cells that survived each one (`after`).
dev_passages_between <- function(passages, last_day) {
  passages[passages$day > 0 & passages$day <= last_day, , drop = FALSE]
}

dev_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(object, temporary)
  file.rename(temporary, path)
}

dev_task_row <- function(task_file, index) {
  tasks <- utils::read.delim(task_file, stringsAsFactors = FALSE)
  tasks[as.integer(index), , drop = FALSE]
}

# rc freeze: a path stored by an as-run dev study (inside a saved record) maps
# into the freeze: the as-run dev results root to DEV_RESULTS, the alfak2-era
# backup to DEV_BACKUP (data/dev/backup).
rc_as_run_path <- function(p) {
  kdr <- "/share/lab_crd/taoli/Project/karst_dev/results/"
  bk <- "/share/lab_crd/taoli/Project/alfak2_results_backup_2026-11-17/validation/sparse_joint/"
  ifelse(startsWith(p, kdr), file.path(dev_results, substring(p, nchar(kdr) + 1L)),
         ifelse(startsWith(p, bk), file.path(dev_backup, substring(p, nchar(bk) + 1L)), p))
}
