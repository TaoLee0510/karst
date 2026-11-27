#!/usr/bin/env Rscript
# Step 1 of dev/realdata/noise/PLAN_v6.md: the cloud-level statistics S_A and
# S_B of one v4 screen folder (v3 simulations, unmodified, rho 1), of a fifth
# of the paper's missegregation simulations (windows 2, 4, 8 sampling times
# before day 120, first and last; 1000 cells and a 300-cell subsample; rho 1),
# or of the real data (Salehi rho 1; NUMBAT rho 0.141).
#   Rscript dev/realdata/noise/v6/stats_task.R TASK_FILE ID   (columns: kind, tag, part)
# Output: DEV_RESULTS/realdata/noise/v6/stats/<tag>.rds, a list of records.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
rd_load_karst()
Sys.setenv(OMP_NUM_THREADS = "1")
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "cloud_level.R"))
v4 <- file.path(dev_results, "realdata", "noise", "v4")
out_file <- file.path(dev_results, "realdata", "noise", "v6", "stats", paste0(task$tag, ".rds"))
col <- function(x, j) { v <- stats::setNames(as.numeric(x[, j]), rownames(x)); v[v > 0] }
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
one <- function(id, x0, x1, dt, rho, meta, keep_counts = FALSE) {
  rec <- tryCatch(list(status = "ok", stats = cl_stats(x0, x1, dt, rho, seed_of(id))),
                  error = function(e) list(status = "failed", message = conditionMessage(e)))
  if (keep_counts) rec <- c(rec, list(x0 = x0, x1 = x1))
  c(list(dataset_id = id, dt = dt), meta, rec)
}
bench <- function() {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f, local = globalenv())
  bench_paths()
}

started <- proc.time()[["elapsed"]]
if (task$kind %in% c("eset", "c2")) {
  files <- list.files(file.path(v4, "screen", task$tag), pattern = "[.]rds$", full.names = TRUE)
  recs <- lapply(files, function(p) {
    scr <- readRDS(p)
    meta <- scr[c("kind", "set", "variant", "landscape_id", "half", "depth", "bottleneck", "truth")]
    if (!identical(scr$status, "ok")) return(c(list(dataset_id = scr$dataset_id, status = "no screen"), meta))
    one(scr$dataset_id, scr$x0, scr$x1, scr$dt, 1, meta)
  })
} else if (task$kind == "paper") {
  paths <- bench()
  sims <- sort(sub("[.]Rds$", "", list.files(paths$misseg_raw, pattern = "[.]Rds$")))
  sims <- sims[(seq_along(sims) - 1L) %% 5L + 1L == task$part]
  sub_cells <- function(v, n, seed) {
    set.seed(seed)
    t <- table(sample(rep(names(v), v), n))
    stats::setNames(as.numeric(t), names(t))
  }
  recs <- unlist(lapply(sims, function(sim_id) {
    sim <- read_abm_sim(sim_id, paths$misseg_raw)
    info <- abm_sim_info(sim_id)
    x <- sim$abm_output$x; times <- as.numeric(colnames(x))
    unlist(lapply(c(2L, 4L, 8L), function(w) {
      tt <- abm_window_times(times, w, TRUE, 120)
      j0 <- which(times == tt[1]); j1 <- which(times == tt[2])
      x0 <- col(x, j0); x1 <- col(x, j1)
      lapply(c(1000L, 300L), function(n) {
        id <- sprintf("%s_w%d_n%04d", sim_id, w, n)
        y0 <- if (n < sum(x0)) sub_cells(x0, n, seed_of(paste0(id, "_t0"))) else x0
        y1 <- if (n < sum(x1)) sub_cells(x1, n, seed_of(paste0(id, "_t1"))) else x1
        one(id, y0, y1, tt[2] - tt[1], 1,
            list(kind = "paper", set = "paper", sim_id = sim_id, landscape_id = sim_id, truth = info$pmis, window = w,
                 depth = n, t0 = tt[1], t1 = tt[2], wavelength = info$w), keep_counts = TRUE)
      })
    }), recursive = FALSE)
  }), recursive = FALSE)
} else if (task$kind == "salehi") {
  paths <- bench()
  sal <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "decomp", "tasks", "fit_ext.tsv"),
                           stringsAsFactors = FALSE)
  recs <- lapply(seq_len(nrow(sal)), function(i) {
    yi <- salehi_input(sal$id[i], ends = TRUE, paths)
    x <- as.matrix(yi$x); tt <- input_times(yi)
    one(sal$id[i], col(x, 1), col(x, ncol(x)), diff(tt[c(1L, length(tt))]), 1,
        list(kind = "salehi", line = sal$line[i], group = sal$group[i]))
  })
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  recs <- lapply(seq_len(nrow(nb)), function(i) {
    r <- rd_restore(nb$cancer[i], nb$sample[i])
    if (r$status != "ok") return(NULL)
    x <- rd_ends(r$restored)
    one(nb$sample[i], col(x, 1), col(x, 2), diff(range(r$times)), 0.141, list(kind = "numbat", cancer = nb$cancer[i]))
  })
  recs <- Filter(Negate(is.null), recs)
} else stop("unknown kind ", task$kind)
dev_save_rds(recs, out_file)
cat(task$tag, length(recs), "records,", sum(vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "ok,",
    round(proc.time()[["elapsed"]] - started), "s\n")
