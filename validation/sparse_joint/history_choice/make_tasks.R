#!/usr/bin/env Rscript
# Task tables of the history choice study.
#
#   Rscript .../history_choice/make_tasks.R KARST_DIR ROUND FAILED_TSV OUT_DIR [N_PER_STRATUM]
#
# dev_tasks.tsv: the 168 scoreboard datasets of the development benchmark
#   (96 core, 72 stress), one task each.
# abm_tasks.tsv: fit tasks of the ALFA-K paper benchmark's ABM and
#   missegregation simulations, from the task tiers of ROUND
#   (benchmark/paper/tasks/ROUND/tasks_{abm,misseg}_mem*.tsv), with the memory
#   of their tier:
#   - set "failed": the ABM and missegregation fit tasks of FAILED_TSV (source,
#     dataset, design, method; the tasks whose fit failed in round 2c60af1
#     with karst_inner_failure, where candidates rescued by the second search
#     are);
#   - set "random": N_PER_STRATUM (default 34) others drawn at random (seed 1)
#     from each design x method stratum (w2/w4/w8_ends karst; w4/w8_all
#     karst, karst_pairwise, karst_state_space), every tier included.
args <- commandArgs(TRUE)
karst_dir <- args[[1]]
round <- args[[2]]
failed_tsv <- args[[3]]
out_dir <- args[[4]]
n_per <- if (length(args) >= 5L) as.integer(args[[5]]) else 34L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

index <- utils::read.csv(file.path(karst_dir, "benchmark", "results", "sparse_two_timepoint", "dataset_index_dev.csv"),
                         stringsAsFactors = FALSE)
dev <- index[index$scoreboard, c("dataset_id", "scenario", "role")]
utils::write.table(dev, file.path(out_dir, "dev_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tier_dir <- file.path(karst_dir, "benchmark", "paper", "tasks", round)
tier_files <- Sys.glob(file.path(tier_dir, "tasks_*_mem*.tsv"))
tier_files <- tier_files[grepl("tasks_(abm|misseg)_mem", basename(tier_files))]
tasks <- do.call(rbind, lapply(tier_files, function(f) {
  d <- utils::read.delim(f, stringsAsFactors = FALSE)
  d$mem <- sub("^.*_mem([0-9]+G)\\.tsv$", "\\1", basename(f))
  d[, c("source", "dataset", "design", "method", "pm", "mem")]
}))
key <- function(d) paste(d$source, d$dataset, d$design, d$method)
failed <- utils::read.delim(failed_tsv, header = FALSE, col.names = c("source", "dataset", "design", "method"),
                            stringsAsFactors = FALSE)
failed <- failed[failed$source %in% c("abm", "misseg"), ]
in_failed <- key(tasks) %in% key(failed)
stopifnot(sum(in_failed) == nrow(failed))
chosen_failed <- tasks[in_failed, ]
chosen_failed$set <- "failed"
rest <- tasks[!in_failed, ]
set.seed(1)
strata <- split(rest, paste(rest$design, rest$method))
random <- do.call(rbind, lapply(strata, function(d) d[sort(sample.int(nrow(d), min(n_per, nrow(d)))), ]))
random$set <- "random"
abm <- rbind(chosen_failed, random)
rownames(abm) <- NULL
utils::write.table(abm, file.path(out_dir, "abm_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("dev tasks:", nrow(dev), "\n")
cat("abm tasks:", nrow(abm), "\n")
print(table(abm$set, abm$mem))
print(table(paste(abm$design, abm$method), abm$set))
