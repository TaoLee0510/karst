# Task tables of karst's candidate missegregation rates (the round of
# ac36fd1, fit_sparse_joint(beta = "auto")): one task per input and candidate
# rate, one fit each.
#
#   Rscript scripts/make_beta_tasks.R OUT_DIR REFERENCE_RESULTS
#
# Inputs: every karst input of the benchmark, as in scripts/make_tasks.R
# (ABM and missegregation w2/w4/w8_ends and w4/w8_all; Salehi ends and all),
# "all" fitted by the default method, anchored. Candidates: 5e-5 times
# karst_beta_multiples (R/fit.R; the grid of beta = "auto"). The candidate at
# 5e-5 is the reference round's karst fit wherever that took 5e-5 (ABM,
# Salehi, and the missegregation simulations at 5e-5); it is not refitted,
# and scripts/beta_picks.R reads it from fits/<...>/karst, linked from
# REFERENCE_RESULTS (results_05a8855). The other missegregation simulations
# get it as a task. A task writes to BENCH_RESULTS/candidates/ (run_task.R,
# column beta_label) and runs no forecasts.
#
# Memory: per input from the reference round's peak of its karst task (fit,
# outputs and forecasts; REFERENCE_RESULTS/tables/resources.Rds), below
# 2.5 GB -> 4 GB, below 5 GB -> 8 GB, otherwise 16 GB (a peak not measured:
# 8 GB for two time points, 16 GB otherwise). The domain, and with it the size
# of the model, does not depend on the rate. Time: 6 h for 4 GB, 12 h above.
#
# Writes to OUT_DIR: beta_tasks.tsv (every task), beta_mem<tier>.tsv (one
# array each, rows renumbered, `orig_task` kept), beta_check.tsv (the 5e-5
# candidate refitted for 12 inputs whose reference fit is reused, as
# beta_label x1_check: are they identical, and does the reused fit's
# objective equal the profile's?) and beta_tiers.tsv.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
out <- args[[1]]
reference <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
s <- bench_settings

## Inputs (karst tasks of make_tasks.R) -----------------------------------------------
inputs <- function(source, dataset, designs, minobs, n0, nb, pm) {
  merge(data.frame(source = source, dataset = dataset, pm_data = pm, n0 = n0, nb = nb,
                   minobs = paste(minobs, collapse = ","), stringsAsFactors = FALSE),
        data.frame(design = designs, stringsAsFactors = FALSE))
}
sims <- sub("\\.Rds$", "", list.files(paths$abm_raw, pattern = "\\.Rds$"))
misseg_sims <- sub("\\.Rds$", "", list.files(paths$misseg_raw, pattern = "\\.Rds$"))
traj <- salehi_trajectories(paths)
x <- rbind(
  inputs("abm", sims, abm_designs$design, s$minobs, s$abm$n0, s$abm$nb, s$abm$pm),
  do.call(rbind, lapply(misseg_sims, function(id) {
    inputs("misseg", id, abm_designs$design, s$misseg$minobs, s$misseg$n0, s$misseg$nb, abm_sim_info(id)$pmis)
  })),
  inputs("salehi", traj$traj_id, "ends", s$minobs, s$salehi$n0, s$salehi$nb, s$salehi$pm),
  inputs("salehi", traj$traj_id[traj$length > 2], "all", s$minobs, s$salehi$n0, s$salehi$nb, s$salehi$pm))
# the reference round's karst fit is the 5e-5 candidate where it took 5e-5
x$reuse_x1 <- abs(x$pm_data / karst_beta_default - 1) < 1e-9
cat("inputs:", nrow(x), " (5e-5 candidate reused:", sum(x$reuse_x1), ")\n")
print(table(x$source, x$design))

## Memory tier ---------------------------------------------------------------------------
res <- readRDS(file.path(reference, "tables", "resources.Rds"))
res <- res[res$method == "karst", c("source", "dataset", "design", "peak_mb")]
x <- merge(x, res, by = c("source", "dataset", "design"), all.x = TRUE)
all_design <- grepl("all$", x$design)
x$mem <- ifelse(is.na(x$peak_mb), ifelse(all_design, "16G", "8G"),
                ifelse(x$peak_mb < 2500, "4G", ifelse(x$peak_mb < 5000, "8G", "16G")))
cat("inputs without a measured peak:", sum(is.na(x$peak_mb)), "\n")
time_limit <- c("4G" = "06:00:00", "8G" = "12:00:00", "16G" = "12:00:00")

## Tasks ---------------------------------------------------------------------------------
cand <- merge(x, data.frame(beta_multiple = karst_beta_multiples))
cand <- cand[!(cand$beta_multiple == 1 & cand$reuse_x1), ]
cand$beta_label <- beta_label(cand$beta_multiple)
cand$pm <- karst_beta_default * cand$beta_multiple
cand$method <- "karst"
cols <- c("source", "dataset", "design", "method", "minobs", "n0", "nb", "pm", "beta_multiple", "beta_label", "mem")
cand <- cand[order(cand$source, cand$dataset, cand$design, cand$beta_multiple), cols]
write_tsv <- function(d, name) {
  d <- cbind(task = seq_len(nrow(d)), d)
  utils::write.table(d, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file.path(out, name))
}
write_tsv(cand, "beta_tasks.tsv")
cat("candidate tasks:", nrow(cand), "\n")
print(table(cand$source, cand$beta_label))
listing <- list()
for (m in names(time_limit)) {
  part <- cand[cand$mem == m, ]
  if (!nrow(part)) next
  part <- data.frame(task = seq_len(nrow(part)), part, orig_task = match(
    paste(part$source, part$dataset, part$design, part$beta_label),
    paste(cand$source, cand$dataset, cand$design, cand$beta_label)))
  file <- file.path(out, paste0("beta_mem", m, ".tsv"))
  utils::write.table(part, file, sep = "\t", quote = FALSE, row.names = FALSE)
  listing[[length(listing) + 1L]] <- data.frame(file = file, tasks = nrow(part), mem = m, time = time_limit[[m]])
}

## Check of the reused 5e-5 fits -----------------------------------------------------------
# Two inputs per source and design group; Salehi "all": one trajectory where
# anchored left out time points (a sweep) and one where it did not.
# `at`: positions in the sorted inputs (0 first, 1 last), to spread the check
# over wavelengths and trajectories
pick <- function(src, design, at) {
  d <- x[x$source == src & x$design == design & x$reuse_x1, ]
  d[order(d$dataset), ][unique(1L + round(at * (nrow(d) - 1L))), ]
}
sal_runs <- readRDS(file.path(reference, "tables", "salehi_runs.Rds"))
sweep <- sort(unique(sal_runs$fi[sal_runs$method == "karst" & sal_runs$design == "all" &
                                   sal_runs$dropped_samples %in% seq_len(100)]))
no_sweep <- sort(setdiff(x$dataset[x$source == "salehi" & x$design == "all"], sweep))
check <- rbind(pick("abm", "w2_ends", 0), pick("abm", "w4_ends", 0.3), pick("abm", "w8_ends", 0.55),
               pick("abm", "w4_all", 0.8), pick("abm", "w8_all", 1),
               pick("misseg", "w8_ends", 0.5), pick("misseg", "w8_all", 1),
               pick("salehi", "ends", c(0, 1)),
               x[x$source == "salehi" & x$design == "all" & x$dataset %in% c(sweep[1], no_sweep[length(no_sweep)]), ])
check$beta_multiple <- 1
check$beta_label <- "x1_check"
check$pm <- karst_beta_default
check$method <- "karst"
check$mem <- "16G"
check <- check[cols]
write_tsv(check, "beta_check.tsv")
listing[[length(listing) + 1L]] <- data.frame(file = file.path(out, "beta_check.tsv"), tasks = nrow(check), mem = "16G",
                                              time = "06:00:00")
listing <- do.call(rbind, listing)
utils::write.table(listing, file.path(out, "beta_tiers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(transform(listing, file = basename(file)), row.names = FALSE)
cat("check inputs:\n")
print(check[c("source", "dataset", "design")], row.names = FALSE)
cat("Total:", sum(listing$tasks), "tasks\n")
