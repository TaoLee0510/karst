#!/usr/bin/env Rscript
# Step 1 of dev/realdata/noise/PLAN_v7.md: S_A and S_B (R/cloud_level.R, rho 1)
# of the v7 datasets of one e-set landscape and multiple, or of one C2
# landscape, multiple and bottleneck; or, for the report only, of a fifth of
# the paper's missegregation simulations (as v6), of Salehi (rho 1) or of
# NUMBAT (rho 0.141), recomputed here rather than read from v6.
#   Rscript dev/realdata/noise/v7/stats_task.R TASK_FILE ID   (columns: kind, tag, landscape, multiple, bottleneck, part)
# Output: DEV_RESULTS/realdata/noise/v7/stats/<tag>.rds, a list of records.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
rd_load_karst()
Sys.setenv(OMP_NUM_THREADS = "1")
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "cloud_level.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
out_file <- v7_path("stats", paste0(task$tag, ".rds"))
if (file.exists(out_file)) quit(save = "no")
col <- function(x, j) { v <- stats::setNames(as.numeric(x[, j]), rownames(x)); v[v > 0] }
one <- function(id, x0, x1, dt, rho, meta) {
  rec <- tryCatch(list(status = "ok", stats = cl_stats(x0, x1, dt, rho, v7_seed(id) %% 100000L)),
                  error = function(e) list(status = "failed", message = conditionMessage(e)))
  c(list(dataset_id = id, dt = dt), meta, rec)
}
bench <- function() {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(tempdir(), "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f, local = globalenv())
  bench_paths()
}
started <- proc.time()[["elapsed"]]
if (task$kind %in% c("e", "c2")) {
  land <- v7_landscapes[task$landscape, ]
  pattern <- if (task$kind == "e") sprintf("^%s_e[123]_x%02d_", land$landscape_id, as.integer(task$multiple)) else
    sprintf("^%s_c2_x%02d_B%04d_", land$landscape_id, as.integer(task$multiple), as.integer(task$bottleneck))
  files <- list.files(v7_path("sims", task$kind, "datasets"), pattern = pattern, full.names = TRUE)
  expected <- if (task$kind == "e") 160L else 20L
  if (length(files) != expected) stop("expected ", expected, " datasets, found ", length(files))
  recs <- lapply(files, function(p) {
    d <- readRDS(v7_check_read(p))
    meta <- d[c("family", "design", "landscape_id", "wavelength", "half", "multiple", "beta", "interval",
                "bottleneck", "depth_t0", "depth_t1")]
    one(d$dataset_id, col(d$counts, 1), col(d$counts, 2), d$interval, 1, meta)
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
  seed_v6 <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L   # the v6 subsample
  recs <- unlist(lapply(sims, function(sim_id) {
    sim <- read_abm_sim(sim_id, paths$misseg_raw)
    info <- abm_sim_info(sim_id)
    x <- sim$abm_output$x; times <- as.numeric(colnames(x))
    unlist(lapply(c(2L, 4L, 8L), function(w) {
      tt <- abm_window_times(times, w, TRUE, 120)
      x0 <- col(x, which(times == tt[1])); x1 <- col(x, which(times == tt[2]))
      lapply(c(1000L, 300L), function(n) {
        id <- sprintf("%s_w%d_n%04d", sim_id, w, n)
        y0 <- if (n < sum(x0)) sub_cells(x0, n, seed_v6(paste0(id, "_t0"))) else x0
        y1 <- if (n < sum(x1)) sub_cells(x1, n, seed_v6(paste0(id, "_t1"))) else x1
        one(id, y0, y1, tt[2] - tt[1], 1, list(family = "paper", design = "paper", landscape_id = sim_id,
                                               multiple = info$pmis / v7_beta_default, beta = info$pmis,
                                               window = w, depth_t0 = n, depth_t1 = n, half = "paper"))
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
        list(family = "salehi", line = sal$line[i], group = sal$group[i]))
  })
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  recs <- lapply(seq_len(nrow(nb)), function(i) {
    r <- rd_restore(nb$cancer[i], nb$sample[i])
    if (r$status != "ok") return(NULL)
    x <- rd_ends(r$restored)
    one(paste(nb$cancer[i], nb$sample[i], sep = "/"), col(x, 1), col(x, 2), diff(range(r$times)), 0.141,
        list(family = "numbat", cancer = nb$cancer[i]))
  })
  recs <- Filter(Negate(is.null), recs)
} else stop("unknown kind ", task$kind)
dev_save_rds(recs, out_file)
cat(task$tag, length(recs), "records,", sum(vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "ok,",
    round(proc.time()[["elapsed"]] - started), "s\n")
