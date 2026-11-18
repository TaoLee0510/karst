# Salehi et al. passaging data, as the paper's 02_run_alfak.R uses it.
#
# Each input file is one trajectory: consecutive passages of one lineage. The
# paper fits all of them ("all"); here each is also reduced to its first and
# last passage ("ends"). Every first + middle passage pair of a trajectory is
# the first + last pair of a shorter trajectory, so "ends" over all
# trajectories covers both splits.

salehi_trajectories <- function(paths = bench_paths()) {
  files <- list.files(paths$salehi_inputs, pattern = "\\.Rds$")
  id <- sub("\\.Rds$", "", files)
  length <- as.integer(sub(".*_l_(\\d+)_d1_.*", "\\1", id))
  data.frame(traj_id = id, length = length, stringsAsFactors = FALSE)
}

salehi_input <- function(traj_id, ends, paths = bench_paths()) {
  yi <- readRDS(file.path(paths$salehi_inputs, paste0(traj_id, ".Rds")))
  if (ends && ncol(yi$x) > 2L) yi$x <- yi$x[, c(1L, ncol(yi$x)), drop = FALSE]
  yi
}

# Time between the trajectory's last passage and its furthest daughter
# passage in days (02_run_alfak.R), or NULL when nothing follows.
salehi_forecast_time <- function(traj_id, lineages, meta) {
  lin <- lineages[[traj_id]]
  if (is.null(lin) || length(lin$dec1) == 0) return(NULL)
  parent <- utils::tail(lin$ids, 1)
  parent_p <- as.numeric(sub("X", "", meta$timepoint[meta$uid %in% parent]))
  daughter_p <- as.numeric(sub("X", "", meta$timepoint[meta$uid %in% lin$dec1]))
  days_per_p <- ifelse(meta$PDX_id[meta$uid %in% parent] %in% c("SA906", "SA039"), 5, 15)
  (max(daughter_p) - parent_p) * days_per_p
}

# Forecast replicates for one fit directory, from the last passage of the
# input (02_run_alfak.R, with seed r for replicate r).
forecast_salehi_dir <- function(fit_dir, yi, sim_time, pm, nb, settings = bench_settings$forecast_salehi) {
  lscape <- readRDS(file.path(fit_dir, "landscape.Rds"))
  x0 <- rep(0, nrow(lscape))
  names(x0) <- lscape$k
  last <- yi$x[rownames(yi$x) %in% lscape$k, ncol(yi$x)]
  x0[names(last)] <- last
  # No cell of the last passage in the landscape: nothing to start from (the
  # paper's 02_run_alfak.R stops here and leaves no forecast).
  if (sum(x0) <= 0) stop("no cell of the last passage is in the landscape")
  x0 <- x0 / sum(x0)
  out <- file.path(fit_dir, "abm_preds")
  dir.create(out, showWarnings = FALSE)
  start <- proc.time()[["elapsed"]]
  for (rep in seq_len(settings$reps)) {
    path <- file.path(out, sprintf("rep_%d.Rds", rep))
    if (file.exists(path)) next
    y <- alfakR::predict_evo(lscape = lscape, p = pm, times = seq(0, sim_time, by = 1), x0 = x0,
                             prediction_type = "ABM", abm_pop_size = settings$abm_pop_size,
                             abm_delta_t = settings$abm_delta_t, abm_max_pop = nb,
                             abm_record_interval = settings$abm_record_interval,
                             abm_culling_survival = settings$abm_culling_survival, abm_seed = rep)
    saveRDS(y, path)
  }
  proc.time()[["elapsed"]] - start
}
