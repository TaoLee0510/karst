# Paths and settings of the benchmark.
#
# Paths come from environment variables (rc freeze: Code/rc_env.sh and
# Code/paper/hpc/env.sh set them from RC_HOME, the rc_1.1.4 folder):
#   PAPER_DATA     the inputs: RC_HOME/data/paper (raw/ and processed/salehi)
#   PAPER_REF      the paper's published ALFA-K outputs, read as reference by
#                  some analyses: RC_HOME/data/reference/paper (processed/)
#   BENCH_RESULTS  where this benchmark writes its fits, forecasts and scores
#                  (a folder under RC_RESULTS; never results_rc itself)
#
# Settings: the algorithm settings of alfak follow PANcanKFLs (nboot 500,
# correct_efflux TRUE, passage_times NULL, minobs 5/10/20); the settings tied
# to a data set (n0, nb, pm) follow the paper for that data set.

bench_paths <- function() {
  rc_home <- Sys.getenv("RC_HOME", "")
  if (!nzchar(rc_home)) stop("RC_HOME is not set (source Code/rc_env.sh)", call. = FALSE)
  data <- Sys.getenv("PAPER_DATA", file.path(rc_home, "data", "paper"))
  ref <- Sys.getenv("PAPER_REF", file.path(rc_home, "data", "reference", "paper"))
  results <- Sys.getenv("BENCH_RESULTS", "")
  if (!nzchar(results)) stop("BENCH_RESULTS is not set (Code/paper/hpc/env.sh)", call. = FALSE)
  # ALIGN=karst (R/align.R): version B of the tests writes and reads a round
  # of its own inside the round, <round>/align_karst/; `results_base` is the
  # round, whose fits it is made from.
  results_base <- results
  if (align_mode() == "karst") results <- file.path(results, "align_karst")
  list(data = data, ref = ref, results_base = results_base,
       abm_raw = file.path(data, "raw", "ABM"),
       misseg_raw = file.path(data, "raw", "ABM_varyMisseg"),
       salehi_inputs = file.path(data, "processed", "salehi", "alfak_inputs"),
       salehi_lineages = file.path(data, "processed", "salehi", "lineages.Rds"),
       salehi_meta = file.path(data, "raw", "salehi", "metadata.csv"),
       paper_abm = file.path(ref, "processed", "ABM"),
       paper_misseg = file.path(ref, "processed", "ABM_vary_misseg"),
       paper_salehi = file.path(ref, "processed", "salehi", "alfak_outputs"),
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
  # (ABM_FORECAST_HORIZON: fewer sampling times, scripts/fig2d_truth.R)
  forecast_abm = list(reps = 5, horizon = as.integer(Sys.getenv("ABM_FORECAST_HORIZON", "10")),
                      abm_pop_size = 2e5, abm_max_pop = 2e7, abm_delta_t = 0.1,
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
