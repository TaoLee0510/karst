# Study B (PLAN.md), Salehi (descriptive): one of codex's windows (shared full
# parent, list E, deterministic engine as ../outside_rank step 3). Forecast
# rate g = x1..x64 of 5e-5 from the round's fit (f = x1), the beta grid's fits
# at f = g, the mixtures from f = x1, and ALFA-K at g = x1 for reference.
#   Rscript scripts/salehi_task.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: window, map (path of the window's map, a copy of
# ../outside_domain/data/salehi_maps). Writes OUT_DIR/salehi/<window>.rds.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "engine.R")); source(file.path(here, "R", "fr.R"))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_file <- file.path(args[[3]], "salehi", paste0(task$window, ".rds"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
started <- proc.time()[["elapsed"]]
k <- od_karst(use_package = TRUE)
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
paper_data <- rc_paper_data
pairs <- utils::read.delim(file.path(rc_shared_parent, "daughter_angles.tsv"), stringsAsFactors = FALSE)  # = the as-run data/ copy
meta <- utils::read.csv(file.path(paper_data, "raw", "salehi", "metadata.csv"), stringsAsFactors = FALSE)
inputs_dir <- file.path(paper_data, "processed", "salehi", "alfak_inputs"); inputs <- list.files(inputs_dir)
grid <- c(1, 2, 4, 8, 16, 32, 64); base_rate <- 5e-5
z <- readRDS(file.path(rc_shared_parent, "maps", basename(task$map))); w <- z$window  # = the as-run data/salehi_maps copy
fit_path <- function(m) if (m == 1) file.path(rc_round, "fits", "salehi", w, "ends", "karst", "minobs_10", "karst_fit.Rds") else
  file.path(rc_rate_grid, sprintf("x%02d", m), "fits", "salehi", w, "ends", "karst", "minobs_10", "karst_fit.Rds")
fits <- list()
for (m in grid) { p <- fit_path(m); if (file.exists(p)) { f <- readRDS(p); if (inherits(f, "karst_sparse_joint_fit")) fits[[as.character(m)]] <- f } }
info <- list(task = task, karst = as.character(utils::packageVersion("karst")), fits_found = names(fits))
if (is.null(fits[["1"]])) { save(c(info, list(status = "no x1 fit"))); quit(save = "no") }

raw <- readRDS(file.path(inputs_dir, paste0(w, ".Rds")))$x
x <- as.matrix(raw[, c(1L, ncol(raw)), drop = FALSE]); storage.mode(x) <- "double"; x <- x[rowSums(x) > 0, , drop = FALSE]
observed <- rownames(x)
dom <- od_domain(k, x)
dist_e <- od_distance(k, dom, z$keys)
dist_of <- function(labels) od_distance(k, dom, labels)
cnK <- od_parse(z$keys)
div_rate <- function(f) k$predict.karst_sparse_joint_fit(f, z$keys)$fitness_mean /
  (2 * k$sparse_joint_self_probability(cnK, f$data$beta) - 1)
rates <- lapply(fits, div_rate)
reproduce <- max(abs(rates[["1"]] / z$karst_native - 1))
eng <- engine_build(z$keys)
n0 <- numeric(length(z$keys)); n0[match(names(z$allocated), z$keys)] <- z$allocated
par <- z$parent_counts
run <- function(r, g) engine_run(eng, n0, r, base_rate * g, 0.1, steps = z$day * 10L)[, 1]
fc <- list()
for (g in grid) fc[[sprintf("f1_g%d", g)]] <- run(rates[["1"]], g)
for (m in setdiff(names(rates), "1")) fc[[sprintf("f%s_g%s", m, m)]] <- run(rates[[m]], as.numeric(m))
fc[["alfak_g1"]] <- run(z$alfak, 1)
for (top in c(4, 8, 16)) fc[[sprintf("mix_1_%d", top)]] <- fr_mix(fc[sprintf("f1_g%d", grid[grid <= top])])

rows <- list()
for (j in which(pairs$window == w)) {
  p <- pairs[j, ]
  fj <- paste(meta[meta$uid == p$daughter_uid, c("datasetname", "timepoint")], collapse = "_")
  tru <- readRDS(file.path(inputs_dir, inputs[grepl(fj, inputs)][1]))$x
  d <- stats::setNames(as.numeric(tru[, ncol(tru)]), rownames(tru)); d <- d[d > 0]
  sc <- fr_score(fc, par, d, observed, dist_e, dist_of)
  rows[[length(rows) + 1L]] <- cbind(window = w, daughter = p$daughter, lineage = p$lineage, sc, stringsAsFactors = FALSE)
}
save(c(info, list(status = "ok", reproduce_max_rel = reproduce, beta = vapply(fits, function(f) f$data$beta, numeric(1)),
                  rows = do.call(rbind, rows), seconds = proc.time()[["elapsed"]] - started)))
