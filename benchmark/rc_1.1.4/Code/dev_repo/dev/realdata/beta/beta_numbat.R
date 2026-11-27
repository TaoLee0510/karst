#!/usr/bin/env Rscript
# beta1 of dev/realdata/beta/PLAN.md on the NUMBAT cancers (report only): one
# input with the diploid row restored, first and last time point, birth rate
# 0.141 per day, fitted at every rate of the grid, with the three criteria.
#   Rscript dev/realdata/beta/beta_numbat.R TASK_FILE ID   (columns: cancer, sample)
# Output: DEV_RESULTS/realdata/beta/fits/numbat/<cancer>/<sample>/<b...>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
source(file.path(dev_repo, "dev", "realdata", "beta", "R", "beta.R"))
Sys.setenv(OMP_NUM_THREADS = "1")
out_dir <- file.path(beta_results, "fits", "numbat", task$cancer, task$sample)
r <- rd_restore(task$cancer, task$sample)
if (r$status != "ok") quit(save = "no")
counts <- rd_ends(r$restored)
dt <- diff(range(r$times))
bounds <- rd_cn_bounds(task$cancer)
seed <- beta_seed(paste(task$cancer, task$sample))
for (m in beta_grid) {
  label <- beta_label(m)
  path <- file.path(out_dir, paste0(label, ".rds"))
  if (file.exists(path)) next
  set.seed(seed)
  run <- beta_fit_timed(list(counts, dt = dt, beta = beta_default * m, birth_rate = rd_rates[rd_central],
                             min_cn = bounds[["min_cn"]], max_cn = bounds[["max_cn"]]))
  key <- list(kind = "numbat", cancer = task$cancer, sample = task$sample, beta = beta_default * m, label = label,
              multiple = m, seconds = run$seconds, karst = version, cells = colSums(counts), dt = dt)
  record <- if (inherits(run$fit, "error")) c(key, beta_failure(run$fit)) else
    c(key, list(status = "ok", criteria = beta_criteria(run$fit, seed)))
  dev_save_rds(record, path)
}
