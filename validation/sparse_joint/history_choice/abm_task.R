#!/usr/bin/env Rscript
# History choice study, ALFA-K paper benchmark ABM / missegregation
# simulations: one fit task (source, dataset, design, method) fitted as the
# benchmark fits it (benchmark/paper/R/fit.R: fit_karst()), every candidate
# history scored against the true fitness of the simulation's landscape.
#
#   Rscript .../history_choice/abm_task.R TASK_FILE ROW OUT_DIR HC_DIR
#
# Run from benchmark/paper (its R/ helpers and bench_paths()) with karst
# 1.1.1 or later on the library path. TASK_FILE has columns source, dataset,
# design, method, pm; HC_DIR is this directory.
args <- commandArgs(TRUE)
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
key <- paste(task$source, task$dataset, task$design, task$method, sep = "__")
out_file <- file.path(out_dir, paste0(key, ".rds"))
if (file.exists(out_file)) quit(save = "no")
for (f in list.files("R", full.names = TRUE)) source(f)
here <- args[[4]]
source(file.path(here, "capture.R"))
cat("karst", as.character(utils::packageVersion("karst")), find.package("karst"), "\n")
stopifnot(utils::packageVersion("karst") >= "1.1.1")

paths <- bench_paths()
raw <- if (task$source == "abm") paths$abm_raw else paths$misseg_raw
sim <- read_abm_sim(task$dataset, raw)
info <- abm_sim_info(task$dataset)
yi <- abm_input(sim, abm_designs[abm_designs$design == task$design, ])
x <- as.matrix(yi$x)
storage.mode(x) <- "double"
times <- input_times(yi)
extra <- list()
if (task$method %in% names(karst_series_methods)) extra$series <- karst_series_methods[[task$method]]

hc_install()
result <- hc_fit_with_candidates(function() {
  do.call(karst::fit_sparse_joint, c(list(x, dt = diff(times), beta = task$pm), extra))
})

score_fit <- function(fit) {
  s <- fit$summary
  truth <- grf_fitness(s$karyotype, info$w, sim$true_landscape)
  groups <- list(all = rep(TRUE, nrow(s)), major = s$major, shell1 = !s$major & s$shell == 1,
                 shell2 = !s$major & s$shell == 2)
  do.call(rbind, lapply(names(groups), function(g) {
    use <- groups[[g]] & is.finite(s$relative_fitness)
    ok <- sum(use) >= 3L && stats::sd(truth[use]) > 0 && stats::sd(s$relative_fitness[use]) > 0
    data.frame(group = g, n = sum(use),
               pearson = if (ok) stats::cor(s$relative_fitness[use], truth[use]) else NA_real_,
               spearman = if (ok) stats::cor(s$relative_fitness[use], truth[use], method = "spearman") else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

if (identical(result$status, "ok")) {
  scores <- lapply(seq_len(nrow(result$candidates)), function(i) {
    fit <- result$fits[[i]]
    if (is.null(fit)) return(NULL)
    cbind(history = result$candidates$history[[i]], score_fit(fit))
  })
  result$scores <- do.call(rbind, scores)
  result$fits <- NULL
  result$fit <- NULL
}
result$task <- task
result$karst_version <- as.character(utils::packageVersion("karst"))
temporary <- paste0(out_file, ".tmp", Sys.getpid())
saveRDS(result, temporary)
file.rename(temporary, out_file)
cat(key, result$status, if (!is.null(result$candidates)) nrow(result$candidates), "candidates",
    if (!is.null(result$identical_kept)) paste("kept identical:", result$identical_kept), round(result$seconds), "s\n")
