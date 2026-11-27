# Kester (benchmark/rc_1.1.4/PLAN.md, part 3): score every fit. Run after all
# fit tasks, as one analysis job (OR_CORES workers).
#   Rscript kester_score.R TASK_FILE ID OUT_DIR      (array of one)
# Salehi-like inputs:
#   paper: the paper's Fig. 2c-style scores of the mean ABM forecast on the
#          method's landscape (kester_paper_scores());
#   engine: both landscapes on list E with the deterministic engine from the
#          parent (the last training passage) to the held-out passage, plus
#          persistence and neutral (flat rates) at the same rate.
# GBM-like inputs (cell holdout): kflInfer's score of the held-out cells at
#   each passage, for KARST's native expected composition, and at the later
#   passage for both landscapes with the engine, persistence and neutral.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (f in c("common.R", "engine.R", "support.R", "kester.R")) source(file.path(here, "R", f))
out <- args[[3]]
k <- od_karst(use_package = TRUE)
tasks <- utils::read.delim(file.path(out, "tasks.tsv"), stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
manifest <- utils::read.delim(file.path(out, "fold_manifest.tsv"), stringsAsFactors = FALSE)
canon <- function(v) k$format_karyotypes(k$parse_karyotypes(v))
dir.create(file.path(out, "scores"), showWarnings = FALSE)
dir.create(file.path(out, "compositions"), showWarnings = FALSE)
keep_comp <- function(f, d) { s <- f / sum(f); s[s >= 1e-4 | names(s) %in% names(d)] }

score_input <- function(id) {
  inp <- readRDS(file.path(out, "inputs", paste0(id, ".rds")))
  xt <- inp$yi$x
  E <- fr_list_e(unique(canon(rownames(xt))))
  eng <- engine_build(E)
  first <- if (inp$kind == "salehi") inp$parent else { v <- inp$train_full[, 1]; v[v > 0] }
  n0 <- numeric(length(E)); n0[match(names(first), E)] <- first
  steps <- as.integer(round(inp$horizon_days * 10))
  run <- function(r, pm) stats::setNames(engine_run(eng, n0, r, pm, 0.1, steps = steps)[, 1], E)
  rows <- list(); comps <- list()
  persist <- stats::setNames(n0 / sum(n0), E)
  for (m in kester_grid) {
    pm <- 5e-5 * m
    arms <- list(persistence = persist, neutral = run(rep(1, length(E)), pm))
    paper <- list()
    for (method in c("karst", "alfak")) for (mo in c(5, 10, 20)) {
      dir <- file.path(out, "fits", id, sprintf("x%02d", m), method, paste0("minobs_", mo))
      run_ok <- file.exists(file.path(dir, "run.Rds")) && identical(readRDS(file.path(dir, "run.Rds"))$status, "ok")
      lab <- sprintf("%s_m%d", method, mo)
      if (!run_ok) { rows[[length(rows) + 1L]] <- data.frame(input = id, multiple = m, arm = lab, type = "fit", status = "fit failed"); next }
      if (method == "karst" && is.null(arms$karst) && file.exists(file.path(dir, "karst_fit.Rds"))) {
        fit <- readRDS(file.path(dir, "karst_fit.Rds"))
        arms$karst <- tryCatch(run(kester_karst_rates(fit, E, k), pm), error = function(e) e)
        if (inp$kind == "gbm") {
          for (j in seq_len(ncol(xt))) {
            col <- paste0("expected_t", j - 1L)
            if (is.null(inp$test) || !col %in% names(fit$summary)) next
            pred <- stats::setNames(fit$summary[[col]], fit$summary$karyotype)
            obs <- stats::setNames(inp$test[, j], rownames(inp$test))
            sc <- kester_score_counts(obs, pred, stats::setNames(inp$train_full[, 1], rownames(inp$train_full)), names(pred))
            rows[[length(rows) + 1L]] <- cbind(data.frame(input = id, multiple = m, arm = "karst_native_expected", type = paste0("passage_", j), status = "ok"), t(sc))
          }
        }
      }
      if (method == "alfak") arms[[lab]] <- tryCatch(run(kester_alfak_on(dir, E, k), pm), error = function(e) e)
      if (inp$kind == "salehi") {
        fc <- file.path(dir, "forecast.Rds")
        if (file.exists(fc) && identical(readRDS(fc)$status, "ok")) {
          preds <- lapply(list.files(file.path(dir, "abm_preds"), full.names = TRUE), readRDS)
          paper[[lab]] <- kester_paper_scores(preds, inp$parent, inp$target)
        } else paper[[lab]] <- NULL
      }
    }
    for (a in names(paper)) rows[[length(rows) + 1L]] <- cbind(data.frame(input = id, multiple = m, arm = a, type = "paper_abm", status = "ok"), t(paper[[a]]))
    for (a in names(arms)) {
      f <- arms[[a]]
      if (inherits(f, "error")) { rows[[length(rows) + 1L]] <- data.frame(input = id, multiple = m, arm = a, type = "engine", status = paste("failed:", conditionMessage(f))); next }
      if (inp$kind == "salehi") {
        sc <- kester_engine_scores(f, inp$parent, inp$target, inp$eps)
        comps[[paste(m, a)]] <- keep_comp(f, inp$target)
      } else {
        if (is.null(inp$test)) next
        obs <- stats::setNames(inp$test[, 2], rownames(inp$test))
        sc <- kester_score_counts(obs, f, stats::setNames(inp$train_full[, 1], rownames(inp$train_full)), E)
        comps[[paste(m, a)]] <- keep_comp(f, obs[obs > 0])
      }
      rows[[length(rows) + 1L]] <- cbind(data.frame(input = id, multiple = m, arm = a, type = "engine", status = "ok"), t(sc))
    }
  }
  r <- do.call(rbind, lapply(rows, function(x) { for (c in setdiff(c("angle", "w1", "overlap", "log_score", "target_cells_outside_E", "target_cells", "wasserstein",
    "coverage_target", "total_cells", "scored_cells", "scored_coverage", "nll_per_cell", "baseline_nll_per_cell"), names(x))) x[[c]] <- NA_real_; x }))
  saveRDS(r, file.path(out, "scores", paste0(id, ".rds")))
  saveRDS(list(observed = if (inp$kind == "salehi") inp$target else inp$test, predicted = comps), file.path(out, "compositions", paste0(id, ".rds")))
  id
}
cores <- as.integer(Sys.getenv("OR_CORES", "4"))
res <- parallel::mclapply(manifest$input, function(id) tryCatch(score_input(id), error = function(e) paste("ERROR", id, conditionMessage(e))),
                          mc.cores = cores, mc.preschedule = FALSE)
print(unlist(res))
all <- do.call(rbind, lapply(list.files(file.path(out, "scores"), full.names = TRUE), readRDS))
utils::write.table(all, file.path(out, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("score rows", nrow(all), "\n")
