# Per-scenario view of a compare_<candidate>.rds table written by compare.R.
for (f in c("config.R")) source(file.path("benchmark/sparse_two_timepoint", f))
args <- commandArgs(TRUE); candidate <- args[1]; config <- stt_config()
rows <- readRDS(file.path(config$output_dir, "report", paste0("compare_", candidate, ".rds")))
main <- rows[rows$stratum %in% c("d0", "d1", "d2"), ]
cand <- main[main$method == candidate, ]; zero <- main[main$method == "zero", ]; alfak <- main[main$method == "alfak_minobs20", ]
key <- function(x) paste(x$dataset_id, x$stratum)
cand$ratio_zero <- cand$rmse / zero$rmse[match(key(cand), key(zero))]
cand$spearman_alfak <- alfak$spearman[match(key(cand), key(alfak))]
options(width = 200)
cat("datasets per scenario:\n"); print(table(cand$scenario[cand$stratum == "d0"]))
cat("\n== RMSE relative to the no-information answer (mean over datasets) ==\n"); print(round(tapply(cand$ratio_zero, list(cand$scenario, cand$stratum), mean), 2))
cat("\n== worst dataset per scenario: RMSE relative to the no-information answer ==\n"); print(round(tapply(cand$ratio_zero, list(cand$scenario, cand$stratum), max), 2))
cat("\n== Spearman, candidate (alfak, minobs = 20) ==\n")
a <- tapply(cand$spearman, list(cand$scenario, cand$stratum), mean, na.rm = TRUE); b <- tapply(cand$spearman_alfak, list(cand$scenario, cand$stratum), mean, na.rm = TRUE)
print(noquote(matrix(sprintf("%.2f (%.2f)", a, b), nrow = nrow(a), dimnames = dimnames(a))))
cat("\n== coverage of 95% intervals ==\n"); print(round(tapply(cand$coverage95, list(cand$scenario, cand$stratum), mean, na.rm = TRUE), 2))
d0 <- cand[cand$stratum == "d0", ]
cat("\n== fit status / seconds ==\n"); print(table(d0$scenario, d0$status)); print(round(tapply(d0$seconds, d0$scenario, median)))
cat("\nnot converged or Hessian not positive definite:", sum(grepl("conv=1|pdHess=FALSE", d0$detail)), "of", nrow(d0), "\n")
