# Split task tables into memory tiers, one Slurm array per tier.
#
#   Rscript scripts/split_tasks.R TASK_FILE ...
#
# Writes, next to each TASK_FILE, <name>_mem<tier>.tsv with the rows of that
# tier renumbered from 1 (the original number is kept in `orig_task`; the
# outputs do not depend on the task number), and tiers.tsv listing every
# tier file with its memory and time limit, for hpc/submit.sh.
#
# The tiers come from the test round (4 ABM simulations, misseg at the
# highest rate 8e-4, 3 Salehi trajectories; 32 GB per task). Peak memory grows
# with the number of time points alfak fits and falls with minobs:
#   karst, two time points, ABM and misseg         < 2 GB    ->   4 GB,  6 h
#   karst, two time points, Salehi; all time points    8 time points x 1000 cells: 0.8 GB
#                                               and 8.5 min for the fit
#                                               (karst multi-timepoint)
#                                                         ->   8 GB, 12 h
#   alfak, two time points (ends designs)          < 9.3 GB  ->  16 GB,  6 h
#   alfak, all time points, otherwise              < 9.7 GB  ->  32 GB, 12 h
#   alfak, ABM w4_all minobs 5, w8_all minobs 10
#                                               < 15 GB   ->  64 GB, 12 h
#   alfak, ABM w8_all minobs 5, misseg w8_all, Salehi all minobs 5
#                                 30 GB, or killed at 32 GB -> 128 GB, 12 h
# karst's other multi-time-point methods (c7f809c on; the multi-timepoint
# session's simulations: pairwise up to 6 GB and 2x the time of anchored,
# state space up to 10.5 GB and 11x, over 6 h for some 8-time-point fits):
#   karst_pairwise                                       ->  16 GB, 12 h
#   karst_state_space                                    ->  32 GB, 12 h
# 12 h is qos xxlarge's limit; the add-hours monitor extends running tasks.

memory_tier <- function(tasks) {
  all <- tasks$design %in% c("w4_all", "w8_all", "all")
  heaviest <- (tasks$design == "w8_all" & (tasks$minobs == "5" | tasks$source == "misseg")) |
    (tasks$design == "all" & tasks$minobs == "5")
  heavy <- (tasks$source == "abm" & tasks$design == "w4_all" & tasks$minobs == "5") |
    (tasks$source == "abm" & tasks$design == "w8_all" & tasks$minobs == "10")
  ifelse(tasks$method == "karst_pairwise", "16G", ifelse(tasks$method == "karst_state_space", "32G",
  ifelse(tasks$method == "karst", ifelse(all | tasks$source == "salehi", "8G", "4G"),
         ifelse(!all, "16G", ifelse(heaviest, "128G", ifelse(heavy, "64G", "32G"))))))
}
time_limit <- c("4G" = "06:00:00", "8G" = "12:00:00", "16G" = "12:00:00", "32G" = "12:00:00",
                "64G" = "12:00:00", "128G" = "12:00:00")

files <- commandArgs(TRUE)
if (!length(files)) stop("Usage: Rscript scripts/split_tasks.R TASK_FILE ...")
listing <- list()
for (path in files) {
  tasks <- utils::read.delim(path, stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
  tier <- memory_tier(tasks)
  for (m in intersect(names(time_limit), unique(tier))) {
    part <- tasks[tier == m, ]
    part <- data.frame(task = seq_len(nrow(part)), part[setdiff(names(part), "task")], orig_task = part$task)
    out <- file.path(dirname(path), paste0(sub("\\.tsv$", "", basename(path)), "_mem", m, ".tsv"))
    utils::write.table(part, out, sep = "\t", quote = FALSE, row.names = FALSE)
    listing[[length(listing) + 1L]] <- data.frame(file = out, tasks = nrow(part), mem = m, time = time_limit[[m]])
  }
}
listing <- do.call(rbind, listing)
utils::write.table(listing, file.path(dirname(files[1]), "tiers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(transform(listing, file = basename(file)), row.names = FALSE)
cat("Total:", sum(listing$tasks), "tasks\n")
