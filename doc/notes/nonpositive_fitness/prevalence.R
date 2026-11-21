# Prevalence of non-positive GRF fitness in the own benchmark's held-out data.
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/prevalence.R
# Notes: doc/notes/nonpositive_fitness.md. KARST_OWN_ROUND names the own benchmark's
# round (default 6ae2a24); results are read from and written to round_<round>/ of this
# folder. OWN_RESULTS is that round's held-out results (datasets/, nodes/, fits/,
# dataset_index_heldout.csv); the default is the copy made by
# doc/report/build/fetch_own.sh. NF_HERE is this folder (default
# doc/notes/nonpositive_fitness), NF_CORES the workers (default SLURM_CPUS_PER_TASK, else 3).
here <- Sys.getenv("NF_HERE", "doc/notes/nonpositive_fitness")
round <- Sys.getenv("KARST_OWN_ROUND", "6ae2a24")
own <- Sys.getenv("OWN_RESULTS", file.path("doc/report/data", paste0("own_", round)))
out_dir <- file.path(here, paste0("round_", round))
dir.create(out_dir, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("NF_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "3")))
source("benchmark/sparse_two_timepoint/landscape.R")
ids <- sub("\\.rds$", "", list.files(file.path(own, "nodes"), pattern = "^heldout_.*\\.rds$"))
rows <- lapply(ids, function(id) {
  d <- readRDS(file.path(own, "datasets", paste0(id, ".rds"))); tr <- d$truth
  f <- function(k) stt_grf_fitness(k, tr$centroids, tr$wavelength)
  f0 <- f(names(tr$frequency_t0)); f1 <- f(names(tr$frequency_t1))
  cnt <- d$counts; fo <- f(rownames(cnt))
  nd <- readRDS(file.path(own, "nodes", paste0(id, ".rds"))); fn <- f(nd$karyotype)
  data.frame(id = id, lam = tr$wavelength, scenario = d$scenario,
             mass_t0 = sum(tr$frequency_t0[f0 <= 0]), mass_t1 = sum(tr$frequency_t1[f1 <= 0]),
             kary_t0 = mean(f0 <= 0), kary_t1 = mean(f1 <= 0),
             obs_kary = mean(fo <= 0), obs_cells = sum(cnt[fo <= 0, ]) / sum(cnt),
             d0 = mean(fn[nd$distance == 0] <= 0), d1 = mean(fn[nd$distance == 1] <= 0), d2 = mean(fn[nd$distance == 2] <= 0),
             fmed_obs = stats::weighted.mean(fo, rowSums(cnt)), fmin_d1 = min(fn[nd$distance == 1]))
})
r <- do.call(rbind, rows)
saveRDS(r, file.path(out_dir, "prevalence.rds"))
q <- function(x) sprintf("%.3g (%.3g–%.3g) max %.3g", median(x), quantile(x, .25), quantile(x, .75), max(x))
for (l in sort(unique(r$lam))) { s <- r[r$lam == l, ]; cat("\nλ =", l, " datasets", nrow(s), "\n")
  for (v in c("mass_t0","mass_t1","kary_t1","obs_kary","obs_cells","d0","d1","d2","fmed_obs")) cat(sprintf("  %-9s %s\n", v, q(s[[v]]))) }
cat("\nall\n"); for (v in c("mass_t1","obs_cells","d0","d1","d2")) cat(sprintf("  %-9s %s\n", v, q(r[[v]])))
