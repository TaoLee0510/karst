# Paths and settings of the benchmark.
#
# Paths come from environment variables so that the same code runs on the
# HPC (inside the SIF) and on a workstation:
#   PAPER_DATA     data folder of the ALFA-K paper (raw/ and processed/)
#   BENCH_RESULTS  where this benchmark writes its fits, forecasts and scores
#                  (default results_<commit of hpc/karst.lock>, as
#                  hpc/env.sh)
#
# Settings: the algorithm settings of alfak follow PANcanKFLs (nboot 500,
# correct_efflux TRUE, passage_times NULL, minobs 5/10/20); the settings tied
# to a data set (n0, nb, pm) follow the paper for that data set.

bench_paths <- function() {
  data <- Sys.getenv("PAPER_DATA", "/share/lab_crd/M010_ALFAK_2023/ALFA-K/data")
  results <- Sys.getenv("BENCH_RESULTS", "")
  if (!nzchar(results)) {
    lock <- readLines("hpc/karst.lock")
    commit <- sub("^KARST_COMMIT=", "", grep("^KARST_COMMIT=", lock, value = TRUE))
    results <- file.path("/share/lab_crd/taoli/Project/karst/benchmark/paper", paste0("results_", commit))
  }
  list(data = data,
       abm_raw = file.path(data, "raw", "ABM"),
       misseg_raw = file.path(data, "raw", "ABM_varyMisseg"),
       salehi_inputs = file.path(data, "processed", "salehi", "alfak_inputs"),
       salehi_lineages = file.path(data, "processed", "salehi", "lineages.Rds"),
       salehi_meta = file.path(data, "raw", "salehi", "metadata.csv"),
       paper_abm = file.path(data, "processed", "ABM"),
       paper_misseg = file.path(data, "processed", "ABM_vary_misseg"),
       paper_salehi = file.path(data, "processed", "salehi", "alfak_outputs"),
       results = results)
}

bench_settings <- list(
  minobs = c(5, 10, 20),
  # alfak (PANcanKFLs) and the matching karst outputs
  nboot = 500,
  correct_efflux = TRUE,
  # data sets (paper): ABM S02/S03, missegregation S02b, Salehi 02_run_alfak.R
  abm = list(n0 = 2e5, nb = 2e7, pm = 5e-5, last_time = 120),
  misseg = list(n0 = 2e5, nb = 2e7, minobs = 10, last_time = 120),
  # Most Salehi inputs contain karyotypes called with 0 copies of a
  # chromosome; alfak takes them as they are, so karst's copy-number lower
  # bound is lowered to admit them (the ABM never produces them).
  salehi = list(n0 = 1e5, nb = 1e7, pm = 5e-5, karst_args = list(min_cn = 0)),
  # forecasts: ABM as S03, Salehi as 02_run_alfak.R; replicate r uses seed r
  forecast_abm = list(reps = 5, horizon = 10, abm_pop_size = 2e5, abm_max_pop = 2e7, abm_delta_t = 0.1,
                      abm_culling_survival = 0.01, abm_record_interval = 10),
  # (the Salehi forecasts take abm_max_pop = nb, as 02_run_alfak.R)
  forecast_salehi = list(reps = 5, abm_pop_size = 1e5, abm_delta_t = 0.1,
                         abm_culling_survival = 0.01, abm_record_interval = 10)
)

# Time windows of the ABM data. Each paper window (the last 2, 4 or 8 passages
# before day 120) is fitted with all its time points ("all", alfak only, as in the
# paper) and with its first and last time point ("ends", both methods). The two
# fits of the 2-passage window coincide, so it has only "ends".
abm_designs <- data.frame(design = c("w2_ends", "w4_ends", "w8_ends", "w4_all", "w8_all"),
                          window = c(2L, 4L, 8L, 4L, 8L),
                          ends = c(TRUE, TRUE, TRUE, FALSE, FALSE),
                          stringsAsFactors = FALSE)
