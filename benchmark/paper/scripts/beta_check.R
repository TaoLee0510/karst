# Can the reference round's karst fits stand in for the 5e-5 candidate?
# (round ac36fd1; scripts/make_beta_tasks.R, beta_check.tsv)
#
#   Rscript scripts/beta_check.R TASK_DIR
#
# For each input of TASK_DIR/beta_check.tsv, the 5e-5 fit refitted with this
# round's karst (candidates/<...>/x1_check) against the reused fit
# (fits/<...>/karst, linked from the reference round):
#   identical     landscape.Rds, bootstrap_res.Rds and xval.Rds identical for
#                 every minobs
#   objective     the refit's objective (run.Rds) equals the lowest finite
#                 objective of the reused fit's history profile, which
#                 scripts/beta_picks.R uses for the reused fits
# Writes BENCH_RESULTS/tables/beta_check.tsv; passes if every row passes.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
check <- utils::read.delim(file.path(args[[1]], "beta_check.tsv"), stringsAsFactors = FALSE,
                           colClasses = c(minobs = "character"))
rows <- lapply(seq_len(nrow(check)), function(i) {
  rel <- file.path(check$source[i], check$dataset[i], check$design[i])
  new <- file.path(paths$results, "candidates", rel, "x1_check")
  old <- file.path(paths$results, "fits", rel, "karst")
  ms <- strsplit(check$minobs[i], ",", fixed = TRUE)[[1]]
  same <- vapply(ms, function(m) {
    all(vapply(c("landscape.Rds", "bootstrap_res.Rds", "xval.Rds"), function(f) {
      a <- file.path(new, paste0("minobs_", m), f)
      b <- file.path(old, paste0("minobs_", m), f)
      file.exists(a) && file.exists(b) && identical(readRDS(a), readRDS(b))
    }, logical(1)))
  }, logical(1))
  rn <- file.path(new, paste0("minobs_", ms[1]), "run.Rds")
  ro <- file.path(old, paste0("minobs_", ms[1]), "run.Rds")
  obj_new <- if (file.exists(rn)) first_or_na(readRDS(rn)$objective) else NA_real_
  prof <- if (file.exists(ro)) readRDS(ro)$history_profile$objective else NULL
  obj_old <- if (length(prof) && any(is.finite(prof))) min(prof[is.finite(prof)]) else NA_real_
  data.frame(check[i, c("source", "dataset", "design")], identical = all(same), objective_refit = obj_new,
             objective_reused = obj_old, objective_equal = isTRUE(all.equal(obj_new, obj_old, tolerance = 1e-10)),
             stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
res$pass <- res$identical & res$objective_equal
dir.create(file.path(paths$results, "tables"), recursive = TRUE, showWarnings = FALSE)
utils::write.table(res, file.path(paths$results, "tables", "beta_check.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
options(width = 200, digits = 12)
print(res, row.names = FALSE)
cat("\nReused 5e-5 fits stand in for the candidate:", all(res$pass), "\n")
