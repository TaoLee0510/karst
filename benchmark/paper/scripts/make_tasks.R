# Task tables for the HPC arrays.
#
#   Rscript scripts/make_tasks.R OUT_DIR [METHOD ...]
#
# Writes tasks_abm.tsv, tasks_misseg.tsv, tasks_salehi.tsv and tasks_test.tsv
# (a small subset of the three) to OUT_DIR. One row is one array task: a
# method fitted to one input (data set and design) at the listed minobs
# values, then its forecasts and scores (scripts/run_task.R). alfak gets one
# task per minobs; karst fits once per input and derives every minobs. Both
# methods get every design: first and last time point ("ends") and every time point
# ("all"; karst from its multi-timepoint version on). METHOD (alfak,
# karst, karst_pairwise, karst_state_space) keeps the tasks of those
# methods only, e.g. karst for a new round whose alfak fits are linked from
# an earlier one; the default is alfak and karst. karst_pairwise and
# karst_state_space (fit_sparse_joint(series = ...), c7f809c on) get the
# designs with more than two time points only. karst_free (R/fit.R) gets
# every design; its tasks go to a results folder of their own. A new round writes its
# tables to a folder of its own (tasks/<commit>), so that the tables of
# earlier rounds, which their job records refer to, stay as they were.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
out <- args[[1]]
methods <- if (length(args) > 1L) args[-1L] else c("alfak", "karst")
stopifnot(all(methods %in% c("alfak", "karst", names(karst_series_methods), names(karst_variant_methods))))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
s <- bench_settings

expand <- function(source, dataset, designs, ends, minobs, n0, nb, pm) {
  base <- merge(data.frame(dataset = dataset, pm = pm, stringsAsFactors = FALSE),
                data.frame(design = designs, ends = ends, stringsAsFactors = FALSE))
  alfak <- merge(base, data.frame(minobs = as.character(minobs), stringsAsFactors = FALSE))
  alfak$method <- "alfak"
  karst <- base
  karst$minobs <- rep(paste(minobs, collapse = ","), nrow(karst))
  karst$method <- rep("karst", nrow(karst))
  series <- do.call(rbind, lapply(names(karst_series_methods), function(m) {
    x <- karst[!karst$ends, , drop = FALSE]
    x$method <- rep(m, nrow(x))
    x
  }))
  variants <- do.call(rbind, lapply(names(karst_variant_methods), function(m) {
    x <- karst
    x$method <- rep(m, nrow(x))
    x
  }))
  tasks <- rbind(alfak, karst, series, variants)
  tasks <- tasks[tasks$method %in% methods, ]
  tasks$source <- source
  tasks$n0 <- n0
  tasks$nb <- nb
  tasks[order(tasks$method, tasks$dataset, tasks$design, tasks$minobs),
        c("source", "dataset", "design", "method", "minobs", "n0", "nb", "pm")]
}

sims <- sub("\\.Rds$", "", list.files(paths$abm_raw, pattern = "\\.Rds$"))
abm <- expand("abm", sims, abm_designs$design, abm_designs$ends, s$minobs, s$abm$n0, s$abm$nb, s$abm$pm)

misseg_sims <- sub("\\.Rds$", "", list.files(paths$misseg_raw, pattern = "\\.Rds$"))
misseg <- do.call(rbind, lapply(misseg_sims, function(id) {
  expand("misseg", id, abm_designs$design, abm_designs$ends, s$misseg$minobs, s$misseg$n0, s$misseg$nb,
         abm_sim_info(id)$pmis)
}))

traj <- salehi_trajectories(paths)
salehi <- rbind(
  expand("salehi", traj$traj_id, "ends", TRUE, s$minobs, s$salehi$n0, s$salehi$nb, s$salehi$pm),
  expand("salehi", traj$traj_id[traj$length > 2], "all", FALSE, s$minobs, s$salehi$n0, s$salehi$nb,
         s$salehi$pm)
)

write_tasks <- function(x, name) {
  x <- cbind(task = seq_len(nrow(x)), x)
  utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(name, nrow(x), "tasks\n")
  print(table(x$method, x$design))
}
write_tasks(abm, "tasks_abm.tsv")
write_tasks(misseg, "tasks_misseg.tsv")
write_tasks(salehi, "tasks_salehi.tsv")

# Test subset: one simulation per wavelength (all designs, both methods), one
# missegregation simulation at the highest rate, and three Salehi
# trajectories of different lengths.
test_sims <- c("w_0p2_m_5e-05_rep_01", "w_0p4_m_5e-05_rep_01", "w_0p8_m_5e-05_rep_01", "w_1p6_m_5e-05_rep_01")
test_traj <- c(traj$traj_id[traj$length == 2][1], traj$traj_id[traj$length == 4][1],
               traj$traj_id[traj$length == max(traj$length)][1])
test <- rbind(abm[abm$dataset %in% test_sims, ],
              misseg[misseg$dataset == misseg_sims[grepl("_m_8e-04_rep_01$", misseg_sims)], ],
              salehi[salehi$dataset %in% test_traj, ])
write_tasks(test, "tasks_test.tsv")
