
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/rescore_summary.R
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
r <- readRDS(file.path(out_dir, "rescore.rds"))
lab <- c(karst_sparse_joint_balance = "KARST", alfak_minobs5 = "alfak5", alfak_minobs10 = "alfak10", alfak_minobs20 = "alfak20")
r$m <- lab[r$method]
cell <- aggregate(r[c("rmse", "rmse_zero", "spearman")], r[c("m", "convention", "wavelength_class", "landscape_id", "stratum")], mean, na.rm = TRUE)
pool <- aggregate(cell[c("rmse", "rmse_zero", "spearman")], cell[c("m", "convention", "wavelength_class", "stratum")], mean, na.rm = TRUE)
pool$ratio <- pool$rmse / pool$rmse_zero
w <- reshape(pool[c("m", "convention", "wavelength_class", "stratum", "ratio", "spearman")], idvar = c("m", "wavelength_class", "stratum"), timevar = "convention", direction = "wide")
w <- w[order(w$wavelength_class, w$stratum, w$m), ]
w[-(1:3)] <- round(w[-(1:3)], 3)
cat("Pooled RMSE ratio vs null and mean Spearman, by convention\n")
print(w[w$wavelength_class == "w0.2", ], row.names = FALSE)
cat("\nOther wavelengths: largest change of ratio or Spearman between conventions\n")
o <- w[w$wavelength_class != "w0.2", ]
o$dr <- pmax(abs(o$ratio.death - o$ratio.clip), abs(o$ratio.exclude - o$ratio.clip)); o$ds <- pmax(abs(o$spearman.death - o$spearman.clip), abs(o$spearman.exclude - o$spearman.clip), na.rm = TRUE)
print(aggregate(cbind(dr, ds) ~ wavelength_class, o, max))
nn <- aggregate(cbind(nonpos, n_nodes) ~ wavelength_class + stratum, r[r$convention == "clip" & r$m == "KARST", ], sum); nn$share <- round(nn$nonpos / nn$n_nodes, 4); print(nn[nn$nonpos > 0, ])
# per landscape at 0.2
c2 <- cell[cell$wavelength_class == "w0.2", ]; c2$ratio <- c2$rmse / c2$rmse_zero
print(reshape(c2[c2$stratum != "d0", c("m","convention","landscape_id","stratum","ratio")], idvar=c("m","landscape_id","stratum"), timevar="convention", direction="wide"), row.names=FALSE, digits=3)
