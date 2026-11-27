# Free disk space after the missegregation-rate round (ac36fd1) is composed:
# delete landscape_posterior_samples.Rds (about 17 MB or more per minobs) in
# every candidate folder that is not used by an arm. Kept: the picked
# candidate of each input (karst_beta_auto), the 5e-5 candidate of the
# missegregation simulations (karst_beta_fixed), and every other file of
# every candidate (run.Rds, landscape.Rds, bootstrap_res.Rds, xval.Rds,
# karst_fit.Rds, scores). No script of the benchmark reads the posterior
# samples; the reused 5e-5 fits of the reference round are not touched (they
# are not under candidates/).
#
#   Rscript scripts/prune_beta_candidates.R [--delete]
#
# Without --delete, lists what would be deleted (count and size) only.
# Reads BENCH_RESULTS/tables/beta_picks.tsv (scripts/beta_picks.R) and
# refuses to run while an input has a missing candidate or no
# karst_beta_auto link (run scripts/compose_beta.R first).

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
delete <- identical(commandArgs(TRUE)[1], "--delete")
picks <- utils::read.delim(file.path(paths$results, "tables", "beta_picks.tsv"), stringsAsFactors = FALSE)
stopifnot(!any(picks$candidates_missing > 0))
links <- file.path(paths$results, "fits", picks$source, picks$dataset, picks$design, "karst_beta_auto")
stopifnot(all(nzchar(Sys.readlink(links))))
root <- file.path(paths$results, "candidates")
keep <- unique(c(
  file.path(root, picks$source, picks$dataset, picks$design, ifelse(is.na(picks$pick_label), "x1", picks$pick_label)),
  with(picks[picks$source == "misseg", ], file.path(root, source, dataset, design, "x1"))))
files <- Sys.glob(file.path(root, "*", "*", "*", "*", "minobs_*", "landscape_posterior_samples.Rds"))
cand <- dirname(dirname(files))
drop <- files[!cand %in% keep]
size <- sum(file.size(drop), na.rm = TRUE)
cat(sprintf("posterior-sample files under candidates/: %d; in candidates kept for an arm: %d; to delete: %d (%.1f GB)\n",
            length(files), sum(cand %in% keep), length(drop), size / 1024^3))
if (delete) {
  ok <- file.remove(drop)
  cat("deleted:", sum(ok), " failed:", sum(!ok), "\n")
} else {
  cat("dry run; pass --delete to delete\n")
}
