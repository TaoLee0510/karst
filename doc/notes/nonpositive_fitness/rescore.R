# Rescore the own benchmark's held-out core fits with three conventions for the
# true rate of karyotypes with GRF fitness f <= 0:
#   clip     (current, KARST ABM: no division, cells survive)  rate = log1p(max(f,0) dt)/dt
#   death    (dies at rate |f|, the signed scale the paper scores on) rate = log1p(f dt)/dt
#   exclude  (lethal at once, as in the paper ABM: no finite rate)   f <= 0 nodes not scored
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/rescore.R
# Notes: doc/notes/nonpositive_fitness.md. OWN_RESULTS is the own benchmark's held-out
# results (datasets/, nodes/, fits/, dataset_index_heldout.csv); the default is the
# copy made by doc/report/build/fetch_own.sh.
here <- "doc/notes/nonpositive_fitness"
own <- Sys.getenv("OWN_RESULTS", "doc/report/data/own_e7e32f8")
for (p in c("landscape.R", "config.R", "score.R")) source(file.path("benchmark/sparse_two_timepoint", p))
idx <- read.csv(file.path(own, "dataset_index_heldout.csv"), stringsAsFactors = FALSE)
idx <- idx[idx$scoreboard & idx$role == "core", ]
methods <- c("karst_sparse_joint_balance", "alfak_minobs5", "alfak_minobs10", "alfak_minobs20")
one <- function(id) {
  d <- readRDS(file.path(own, "datasets", paste0(id, ".rds"))); tr <- d$truth
  nodes <- readRDS(file.path(own, "nodes", paste0(id, ".rds")))
  f <- stt_grf_fitness(nodes$karyotype, tr$centroids, tr$wavelength); dt <- tr$abm_delta_t
  conv <- list(clip = list(rate = log1p(pmax(f, 0) * dt) / dt, keep = rep(TRUE, length(f))),
               death = list(rate = log1p(f * dt) / dt, keep = rep(TRUE, length(f))),
               exclude = list(rate = log1p(pmax(f, 0) * dt) / dt, keep = f > 0))
  do.call(rbind, lapply(methods, function(m) {
    path <- file.path(own, "fits", paste0(id, "__", m, ".rds")); if (!file.exists(path)) return(NULL)
    fit <- readRDS(path)
    pred <- if (is.null(fit$prediction)) data.frame(karyotype = character(), fitness = numeric()) else fit$prediction
    do.call(rbind, lapply(names(conv), function(cn) {
      k <- conv[[cn]]$keep
      s <- stt_score_predictions(nodes[k, , drop = FALSE], conv[[cn]]$rate[k], pred, 5)$metrics
      s <- s[s$stratum %in% c("d0", "d1", "d2"), c("stratum", "n_nodes", "rmse", "rmse_zero", "spearman")]
      cbind(dataset_id = id, method = m, convention = cn, s, nonpos = sapply(s$stratum, function(st) sum(f[nodes$distance == as.integer(sub("d", "", st))] <= 0)))
    }))
  }))
}
r <- do.call(rbind, parallel::mclapply(idx$dataset_id, one, mc.cores = 3))
r <- merge(r, idx[c("dataset_id", "landscape_id", "wavelength_class", "scenario")], by = "dataset_id")
saveRDS(r, file.path(here, "rescore.rds"))
# sanity: clip reproduces the benchmark's metrics
m <- read.csv(file.path(own, "metrics_heldout.csv"), stringsAsFactors = FALSE)
chk <- merge(r[r$convention == "clip", ], m[c("dataset_id", "method", "stratum", "rmse", "spearman")], by = c("dataset_id", "method", "stratum"))
cat("check rows", nrow(chk), "max |rmse diff|", max(abs(chk$rmse.x - chk$rmse.y), na.rm = TRUE), "max |spearman diff|", max(abs(chk$spearman.x - chk$spearman.y), na.rm = TRUE), "\n")
