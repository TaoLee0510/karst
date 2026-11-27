# Step 0, Salehi (descriptive; PLAN.md): PL per t1 cell of KARST's current
# rates against neutral, on the 34 windows of ../outside_rank step 3.
# rc freeze: run as Rscript Code/patch/marginal_pl/scripts/salehi_pl.R (or through Code/hpc/rc_job.sh, RC_ARGS);
# the 35 window fits are the copies of ../outside_rank/data/salehi_fits kept with this study's frozen results;
# output to RC_RESULTS/part4_rate_choice/marginal_pl/results/step01
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "..", "rc_paths.R"))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "rank.R")); source(file.path(here, "R", "mpl.R"))
k <- or_local_karst()
fits <- list.files(file.path(rc_up, "part4_rate_choice", "marginal_pl", "upstream", "outside_rank_salehi_fits"), full.names = TRUE)
r <- do.call(rbind, lapply(fits, function(f) {
  w <- sub("[.]Rds$", "", basename(f)); if (w == "SA000_X6_l_6_d1_4_d2_0") return(NULL)
  fit <- readRDS(f); S <- mpl_setup(k, fit)
  data.frame(window = w, family = if (grepl("^SA906", w)) "SA906" else if (grepl("^SA609|^SA000", w)) "SA609" else "other",
             t0_karyotypes = length(S$n0), t0_in_domain = sum(S$inside), share_t0_cells_in_domain = sum(S$n0[S$inside]) / sum(S$n0),
             pl_karst = mpl_pl_rates(S, S$self * exp(S$lb)), pl_neutral = mpl_pl_rates(S, rep(0, length(S$n0))))
}))
r$karst_minus_neutral <- r$pl_karst - r$pl_neutral
out_dir <- file.path(Sys.getenv("RC_RESULTS"), "part4_rate_choice", "marginal_pl", "results", "step01"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(r, file.path(out_dir, "salehi_pl.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats::aggregate(cbind(share_t0_cells_in_domain, karst_minus_neutral) ~ family, r, stats::median), digits = 3, row.names = FALSE)
cat("KARST above neutral in", sum(r$karst_minus_neutral > 0), "of", nrow(r), "windows\n")
