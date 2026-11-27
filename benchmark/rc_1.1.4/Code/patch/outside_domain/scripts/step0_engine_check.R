# Step 0 of PLAN.md: the deterministic engine against codex's five-seed ABM
# ensembles on the 35 Salehi windows (shared full parent, list E), for the
# KARST native and ALFA-K rate maps.
#   Rscript scripts/step0_engine_check.R [CORES]
# Reads data/salehi_maps/map_XX.rds (codex: keys, allocated, rates, day) and
# the panel's forecasts.rds / daughter_angles.tsv (read only). Writes
# results/step0/engine_check.tsv and prints the verdict.
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "engine.R"))
cores <- as.integer(commandArgs(TRUE)[1]); if (is.na(cores)) cores <- 3L
k <- od_karst(use_package = FALSE)
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
data_root <- rc_paper_data
panel <- rc_shared_parent
pairs <- read.delim(file.path(panel, "daughter_angles.tsv"), stringsAsFactors = FALSE)
abm <- readRDS(file.path(panel, "forecasts.rds"))
meta <- read.csv(file.path(data_root, "raw/salehi/metadata.csv"), stringsAsFactors = FALSE)
inputs_dir <- file.path(data_root, "processed/salehi/alfak_inputs"); inputs <- list.files(inputs_dir)
maps <- list.files(file.path(rc_shared_parent, "maps"), full.names = TRUE)  # = the as-run data/salehi_maps copy
out_dir <- file.path(Sys.getenv("RC_RESULTS"), "part4_rate_choice", "outside_domain", "results", "step0"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

one <- function(path) {
  z <- readRDS(path)
  started <- proc.time()[["elapsed"]]
  raw <- readRDS(file.path(inputs_dir, paste0(z$window, ".Rds")))$x
  x <- as.matrix(raw[, c(1L, ncol(raw)), drop = FALSE]); storage.mode(x) <- "double"; x <- x[rowSums(x) > 0, , drop = FALSE]
  dom <- od_domain(k, x)
  eng <- engine_build(z$keys)
  n0 <- numeric(length(z$keys)); n0[match(names(z$allocated), z$keys)] <- z$allocated
  stopifnot(!anyNA(match(names(z$allocated), z$keys)))
  par <- z$parent_counts / sum(z$parent_counts)
  c_par <- od_centre(par)
  g_keys <- od_group(od_distance(k, dom, z$keys))
  rows <- list()
  for (m in c("karst_native", "alfak")) {
    f <- engine_run(eng, n0, z[[m]], 5e-5, 0.1, steps = z$day * 10L)[, 1]
    f <- f[f > 0]
    a <- abm[[z$window]][[m]]$mean_freq
    sd_ <- od_shares(f, g_keys[match(names(f), z$keys)])
    sa <- od_shares(a, g_keys[match(names(a), z$keys)])
    for (j in which(pairs$window == z$window)) {
      p <- pairs[j, ]
      fj <- paste(meta[meta$uid == p$daughter_uid, c("datasetname", "timepoint")], collapse = "_")
      tru <- readRDS(file.path(inputs_dir, inputs[grepl(fj, inputs)][1]))$x
      d <- stats::setNames(as.numeric(tru[, ncol(tru)]), rownames(tru)); d <- d[d > 0]
      obs <- od_centre(d) - c_par
      rows[[length(rows) + 1L]] <- data.frame(window = z$window, daughter = p$daughter, method = m, states = length(z$keys),
        angle_det = od_angle(od_centre(f) - c_par, obs), angle_abm = od_angle(od_centre(a) - c_par, obs),
        angle_codex = if (m == "alfak") p$alfak_angle else p$karst_native_angle,
        det_vs_abm = od_angle(od_centre(f) - c_par, od_centre(a) - c_par),
        inside_det = sd_[["inside"]], inside_abm = sa[["inside"]], halo_det = sd_[["halo"]], halo_abm = sa[["halo"]],
        island_det = sd_[["island"]], island_abm = sa[["island"]], seconds = proc.time()[["elapsed"]] - started)
    }
  }
  do.call(rbind, rows)
}
res <- do.call(rbind, parallel::mclapply(maps, function(p) tryCatch(one(p), error = function(e) { message(basename(p), ": ", conditionMessage(e)); NULL }),
                                         mc.cores = cores, mc.preschedule = FALSE))
utils::write.table(res, file.path(out_dir, "engine_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
cat("rows:", nrow(res), " windows:", length(unique(res$window)), "\n")
cat("check: ABM angles reproduce codex's (max |diff|):", signif(max(abs(res$angle_abm - res$angle_codex)), 3), "\n")
v <- do.call(rbind, lapply(split(res, res$method), function(d) data.frame(method = d$method[1], n = nrow(d),
  angle_diff_median = median(abs(d$angle_det - d$angle_abm)), angle_diff_max = max(abs(d$angle_det - d$angle_abm)),
  det_vs_abm_median = median(d$det_vs_abm),
  share_diff_median_pp = 100 * median(abs(c(d$inside_det - d$inside_abm, d$halo_det - d$halo_abm, d$island_det - d$island_abm))),
  halo_diff_median_pp = 100 * median(abs(d$halo_det - d$halo_abm)))))
v$pass <- v$angle_diff_median <= 3 & v$share_diff_median_pp <= 2
print(v, row.names = FALSE, digits = 3)
cat("\nStep 0", if (all(v$pass)) "PASSES" else "FAILS", "\n")
