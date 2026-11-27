# ABM simulations of the paper: identifiers, true fitness and inputs.

# "w_0p8_m_5e-05_rep_01" -> wavelength 0.8, missegregation rate 5e-05, rep "01"
abm_sim_info <- function(sim_id) {
  parts <- strsplit(sim_id, "_", fixed = TRUE)[[1]]
  data.frame(sim_id = sim_id, w = as.numeric(sub("p", ".", parts[2], fixed = TRUE)),
             pmis = as.numeric(parts[4]), abmrep = parts[6], stringsAsFactors = FALSE)
}

read_abm_sim <- function(sim_id, dir) readRDS(file.path(dir, paste0(sim_id, ".Rds")))

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
