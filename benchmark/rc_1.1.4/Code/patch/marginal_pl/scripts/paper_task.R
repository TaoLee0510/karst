# Paper simulations (PLAN.md), from the round's stored KARST fits (no refit):
# the PL of the truth, KARST and neutral rates (step 0), the update at every w
# (D1a; D1b reported), and accuracy against the true landscape on O, on a
# sample of unseen out-of-domain karyotypes within two steps of the observed
# ones, and on the observed in-domain non-majors.
#   Rscript scripts/paper_task.R TASK_FILE OUT_DIR      (all rows, OR_CORES workers)
# TASK_FILE: ../outside_rank/tasks/sim_tasks.tsv (source, arm, sim_id, design, rep, half).
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "rank.R")); source(file.path(here, "R", "mpl.R"))
tasks <- utils::read.delim(args[[1]], stringsAsFactors = FALSE); out_root <- args[[2]]
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
paper_data <- rc_paper_data
source(rc_abm_R)
k <- od_karst(use_package = TRUE)
grid_w <- c(1/22, 1/8, 1/3, 1)
grf <- function(cn, w, centroids) {
  s <- 1 / (pi * sqrt(nrow(centroids))); acc <- numeric(nrow(cn))
  for (c in seq_len(nrow(centroids))) acc <- acc + sin(sqrt(rowSums(sweep(cn, 2L, centroids[c, ], "-")^2)) / w)
  acc * s
}
nb <- function(labels) {
  cn <- od_parse(labels); out <- character(0)
  for (i in seq_len(ncol(cn))) for (s in c(-1L, 1L)) {
    m <- cn; m[, i] <- m[, i] + s; ok <- m[, i] >= 1L
    out <- c(out, do.call(paste, c(as.data.frame(m[ok, , drop = FALSE]), sep = ".")))
  }
  unique(out)
}
rank_cor <- function(p, t) { if (length(p) < 3L || stats::sd(t) == 0) return(NA_real_); if (stats::sd(p) == 0) return(0); stats::cor(p, t, method = "spearman") }
crmse <- function(p, t) if (length(p) < 2L) NA_real_ else sqrt(mean(((p - mean(p)) - (t - mean(t)))^2))

one <- function(i) {
  t <- tasks[i, ]
  out_file <- file.path(out_root, paste0("paper_", t$source, "_", t$arm), paste0(t$sim_id, "_", t$design, ".rds"))
  if (file.exists(out_file)) return(out_file)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
  started <- proc.time()[["elapsed"]]
  fit_path <- if (t$source == "abm") file.path(rc_fig2d_truth, "step2", "fits", "abm", t$sim_id, t$design, "karst", "minobs_10", "karst_fit.Rds") else
    file.path(rc_fig2d_misseg, t$arm, "fits", "misseg", t$sim_id, t$design, "karst", "minobs_10", "karst_fit.Rds")
  meta <- list(task = t, fit_path = fit_path)
  if (!file.exists(fit_path)) { save(c(meta, list(status = "no fit"))); return(out_file) }
  fit <- readRDS(fit_path)
  if (!inherits(fit, "karst_sparse_joint_fit")) { save(c(meta, list(status = "fit failed"))); return(out_file) }
  info <- abm_sim_info(t$sim_id)
  sim <- readRDS(file.path(paper_data, "raw", if (t$source == "abm") "ABM" else "ABM_varyMisseg", paste0(t$sim_id, ".Rds")))
  net <- function(labels) { cn <- od_parse(labels); grf(cn, info$w, sim$true_landscape) * (2 * k$sparse_joint_self_probability(cn, fit$data$beta) - 1) }
  # nodes: O, observed in-domain non-majors, and up to 3,000 unseen out-of-domain karyotypes within two steps
  x <- as.matrix(fit$data$counts); if (is.null(rownames(x))) rownames(x) <- fit$data$labels
  observed <- rownames(x)[rowSums(x) > 0]
  dom <- fit$summary$karyotype; maj <- dom[fit$summary$major]
  one_s <- setdiff(nb(observed), observed); two_s <- setdiff(nb(one_s), c(observed, one_s))
  unseen <- setdiff(c(one_s, two_s), dom)
  set.seed(abm_sim_seed(t$sim_id) + 11L)
  if (length(unseen) > 3000L) unseen <- sample(unseen, 3000L)
  O <- fit$dropped; inside <- setdiff(intersect(observed, dom), maj)
  labels <- c(O, unseen, inside); grp <- rep(c("O", "unseen", "inside"), c(length(O), length(unseen), length(inside)))
  tru <- net(labels)
  acc <- function(p) c(sp_ou = rank_cor(p[grp != "inside"], tru[grp != "inside"]), sp_un = rank_cor(p[grp == "unseen"], tru[grp == "unseen"]),
                       sp_in = rank_cor(p[grp == "inside"], tru[grp == "inside"]), rmse_in = crmse(p[grp == "inside"], tru[grp == "inside"]))
  S <- mpl_setup(k, fit)
  pl <- c(truth = mpl_pl_rates(S, net(S$lab0)), karst = mpl_pl_rates(S, S$self * exp(S$lb)), neutral = mpl_pl_rates(S, rep(0, length(S$n0))))
  p0 <- k$predict.karst_sparse_joint_fit(fit, labels)$relative_fitness
  arms <- list(current = list(acc = acc(p0))); fits <- list()
  for (w in grid_w) {
    r <- mpl_fit(S, w); nm <- sprintf("w%.4f", w); fits[[nm]] <- r[c("delta", "convergence", "pl", "pl0")]
    arms[[nm]] <- list(acc = acc(mpl_predict(k, fit, labels, r$delta, S$free)$relative_fitness),
                       acc_b = acc(mpl_predict(k, fit, labels, r$delta, S$free, rule_o = TRUE)$relative_fitness))
  }
  save(c(meta, list(status = "ok", info = info, n_sets = table(grp), n_t0 = length(S$n0), n_t0_inside = sum(S$inside), free = S$free,
                    pl = pl, arms = arms, fits = fits, seconds = proc.time()[["elapsed"]] - started)))
  out_file
}
cores <- as.integer(Sys.getenv("OR_CORES", "1"))
res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) tryCatch(one(i), error = function(e) paste("ERROR", tasks$sim_id[i], tasks$design[i], conditionMessage(e))),
                          mc.cores = cores, mc.preschedule = FALSE)
err <- grep("^ERROR", unlist(res), value = TRUE)
cat("rows", nrow(tasks), " errors", length(err), "\n"); if (length(err)) print(utils::head(err, 20))
