# ABM simulations of the paper: identifiers, true fitness and inputs.

# "w_0p8_m_5e-05_rep_01" -> wavelength 0.8, missegregation rate 5e-05, rep "01"
abm_sim_info <- function(sim_id) {
  parts <- strsplit(sim_id, "_", fixed = TRUE)[[1]]
  data.frame(sim_id = sim_id, w = as.numeric(sub("p", ".", parts[2], fixed = TRUE)),
             pmis = as.numeric(parts[4]), abmrep = parts[6], stringsAsFactors = FALSE)
}

read_abm_sim <- function(sim_id, dir) {
  sim <- readRDS(file.path(dir, paste0(sim_id, ".Rds")))
  # ABM_CELLS (scripts/fig2d_truth.R, step 2): a table of cell counts; every
  # time point of the simulation is then subsampled to a count drawn from it,
  # so that the fits, forecasts and scores see Salehi-sized samples.
  cells <- Sys.getenv("ABM_CELLS", "")
  if (nzchar(cells)) sim$abm_output$x <- abm_subsample_x(sim$abm_output$x, sim_id, utils::read.delim(cells)$cells)
  sim
}

# Seed of a simulation ("w_0p8_m_5e-05_rep_07"): fixed by its wavelength,
# missegregation rate and replicate, so that a subsample is the same in every
# job. The main ABM set (rate 5e-5) keeps the seeds it had before the rates of
# ABM_varyMisseg were added (docs/fig2d_misseg_6ae2a24.md); that set's 5e-5
# simulations, different simulations under the same names, share them.
abm_sim_seed <- function(sim_id) {
  info <- abm_sim_info(sim_id)
  rate <- match(info$pmis, c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4))
  stopifnot(!is.na(rate))
  9100000L + 100000L * (rate - 1L) + 1000L * match(info$w, c(0.2, 0.4, 0.8, 1.6)) + as.integer(info$abmrep)
}

# Each time point (column) of the count matrix `x` subsampled without
# replacement to a number of cells drawn from `cells` (at most the cells the
# time point has), with the simulation's seed (abm_sim_seed()). Rows and
# columns are kept; the class of `x` (data frame or matrix) too.
abm_subsample_x <- function(x, sim_id, cells, seed = abm_sim_seed(sim_id)) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  n <- sample(cells, ncol(x), replace = TRUE)
  for (j in seq_len(ncol(x))) {
    counts <- as.numeric(x[, j])
    n_j <- min(n[j], sum(counts))
    idx <- rep.int(seq_along(counts), counts)
    x[, j] <- tabulate(idx[sample.int(length(idx), n_j)], length(counts))
  }
  x
}

# True fitness of karyotype labels on a Gaussian-random-field landscape
# (utils_karyo.R getf of the paper, for many karyotypes at once).
grf_fitness <- function(k, wavelength, centroids) {
  kmat <- do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
  scalef <- 1 / (pi * sqrt(nrow(centroids)))
  vapply(seq_len(nrow(kmat)), function(i) {
    d <- sqrt(colSums((t(centroids) - kmat[i, ])^2))
    sum(sin(d / wavelength) * scalef)
  }, numeric(1))
}

# Sampling times of the window of a design: the last `window` sampling times
# before `last_time`, or only the first and last of them.
abm_window_times <- function(times, window, ends, last_time = 120) {
  w <- utils::tail(times[times < last_time], window)
  if (ends) w[c(1L, length(w))] else w
}

# The input of a design in the form alfakR::alfak() takes (S02 of the paper:
# the columns of the window, rows unchanged).
abm_input <- function(sim, design, last_time = 120) {
  x <- sim$abm_output$x
  times <- as.numeric(colnames(x))
  keep <- abm_window_times(times, design$window, design$ends, last_time)
  list(x = x[, times %in% keep, drop = FALSE], dt = sim$abm_output$dt)
}
