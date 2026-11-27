# Roots of the patch studies (C4 forecast_rate, C6 outside_domain, C7
# marginal_pl) in the rc freeze, sourced by their scripts in place of the
# as-run paths into the karst clone and the patch folders. RC_HOME and
# RC_UPSTREAM come from Code/rc_env.sh (Code/hpc/rc_job.sh passes them in).
rc_home <- Sys.getenv("RC_HOME"); if (!nzchar(rc_home)) stop("RC_HOME is not set", call. = FALSE)
rc_up <- Sys.getenv("RC_UPSTREAM", Sys.getenv("RC_RESULTS"))
rc_round <- file.path(rc_up, "part1_paper", "round_6ae2a24")                 # was results_6ae2a24
rc_rate_grid <- file.path(rc_up, "part1_paper", "salehi_rate_grid")          # was results_6ae2a24/beta_grid
rc_fig2d_truth <- file.path(rc_up, "part4_rate_choice", "fig2d_truth")       # was results_6ae2a24/fig2d_truth
rc_fig2d_misseg <- file.path(rc_up, "part4_rate_choice", "fig2d_misseg")     # was results_6ae2a24/fig2d_misseg
rc_gate_vm <- file.path(rc_up, "part4_rate_choice", "dispersion_gate", "results", "step1", "vm")
rc_shared_parent <- file.path(rc_round, "analysis", "fig2c_two_timepoints", "shared_full_parent")  # Salehi maps, daughter angles
rc_paper_data <- Sys.getenv("PAPER_DATA", file.path(rc_home, "data", "paper"))
rc_abm_R <- file.path(rc_home, "Code", "paper", "R", "abm.R")               # = the staged fig2d_misseg/code/R/abm.R
rc_own_code <- file.path(rc_home, "Code", "own", "benchmark", "sparse_two_timepoint")
# a task table's dataset folder: relative to RC_HOME in the freeze (e.g. data/own_benchmark_dev)
rc_resolve <- function(p) if (startsWith(p, "/")) p else file.path(rc_home, p)
