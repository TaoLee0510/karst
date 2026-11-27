# One own-benchmark dataset (PLAN.md): refit KARST 1.1.3 as the round does
# (as ../outside_rank/scripts/own_task.R), then the PL of the truth, KARST and
# neutral rates (step 0), the update at every w (D1a; D1b reported), and the
# scores: accuracy on O + unseen out-of-domain nodes, in-domain RMSE, and the
# round's d0-d2 scores (K2).
#   Rscript scripts/own_task.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: split (dev | heldout), dataset_id, round_dir.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "rank.R")); source(file.path(here, "R", "mpl.R"))
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
for (part in c("landscape.R", "score.R")) source(file.path(rc_own_code, part))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
task$round_dir <- rc_resolve(task$round_dir)
out_file <- file.path(args[[3]], paste0("own_", task$split), paste0(task$dataset_id, ".rds"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
started <- proc.time()[["elapsed"]]
k <- od_karst(use_package = TRUE)
grid_w <- c(1/22, 1/8, 1/3, 1)
dataset <- readRDS(file.path(task$round_dir, "datasets", paste0(task$dataset_id, ".rds")))
node_path <- file.path(task$round_dir, "nodes", paste0(task$dataset_id, ".rds"))
nodes <- if (file.exists(node_path)) readRDS(node_path) else stt_scoring_nodes(dataset$counts, 1L, 8L, 2000L, 2000L, dataset$seed)
meta <- list(task = task, landscape_id = dataset$landscape_id, wavelength_class = dataset$wavelength_class, scenario = dataset$scenario)
counts <- dataset$counts; storage.mode(counts) <- "integer"
set.seed(dataset$seed)
fit <- tryCatch(karst::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                                        min_cn = 1L, max_cn = 8L, initial_state = "balance"), error = function(e) e)
if (inherits(fit, "error")) { save(c(meta, list(status = "fit failed", error = conditionMessage(fit)))); quit(save = "no") }
pred0 <- stats::predict(fit, nodes$karyotype)
identity <- NA_real_
stored_path <- file.path(task$round_dir, "fits", paste0(task$dataset_id, "__karst_sparse_joint_balance.rds"))
if (task$split == "heldout" && file.exists(stored_path)) {
  st <- readRDS(stored_path)
  if (identical(st$status, "ok")) { m <- match(nodes$karyotype, st$prediction$karyotype)
    identity <- if (anyNA(m)) Inf else max(abs(st$prediction$fitness[m] - pred0$relative_fitness), na.rm = TRUE) }
}
truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
in_dom <- nodes$karyotype %in% fit$summary$karyotype
in_O <- nodes$karyotype %in% fit$dropped
major <- nodes$karyotype %in% fit$summary$karyotype[fit$summary$major]
set_ou <- in_O | (!in_dom & nodes$distance >= 1)          # O + unseen out-of-domain nodes
set_un <- !in_dom & !in_O & nodes$distance >= 1             # unseen out-of-domain only
set_in <- in_dom & !major & nodes$distance == 0             # observed in-domain non-majors
rank_cor <- function(p, t) { if (length(p) < 3L || stats::sd(t) == 0) return(NA_real_); if (stats::sd(p) == 0) return(0); stats::cor(p, t, method = "spearman") }
crmse <- function(p, t) if (length(p) < 2L) NA_real_ else sqrt(mean(((p - mean(p)) - (t - mean(t)))^2))
acc <- function(p) c(sp_ou = rank_cor(p[set_ou], truth[set_ou]), sp_un = rank_cor(p[set_un], truth[set_un]),
                     sp_in = rank_cor(p[set_in], truth[set_in]), rmse_in = crmse(p[set_in], truth[set_in]))
k2 <- function(p) { sc <- stt_score_predictions(nodes[, c("karyotype", "distance", "count_t0", "count_t1")], truth,
                                                  data.frame(karyotype = nodes$karyotype, fitness = p, fitness_sd = pred0$relative_sd))
  sc$metrics[sc$metrics$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse", "spearman", "coverage95")] }

S <- mpl_setup(k, fit)
pl <- c(truth = mpl_pl_rates(S, stt_true_growth_rate(S$lab0, dataset$truth)), karst = mpl_pl_rates(S, S$self * exp(S$lb)),
        neutral = mpl_pl_rates(S, rep(0, length(S$n0))))
arms <- list(current = list(acc = acc(pred0$relative_fitness), k2 = k2(pred0$relative_fitness)))
fits <- list()
for (w in grid_w) {
  r <- mpl_fit(S, w); nm <- sprintf("w%.4f", w); fits[[nm]] <- r[c("delta", "convergence", "pl", "pl0")]
  pa <- mpl_predict(k, fit, nodes$karyotype, r$delta, S$free)$relative_fitness
  pb <- mpl_predict(k, fit, nodes$karyotype, r$delta, S$free, rule_o = TRUE)$relative_fitness
  arms[[nm]] <- list(acc = acc(pa), k2 = k2(pa), acc_b = acc(pb))
}
save(c(meta, list(status = "ok", identity = identity, n_O = length(fit$dropped), n_t0 = length(S$n0), n_t0_inside = sum(S$inside),
                  n_sets = c(ou = sum(set_ou), un = sum(set_un), inside = sum(set_in)), pl = pl, arms = arms, fits = fits,
                  seconds = proc.time()[["elapsed"]] - started)))
