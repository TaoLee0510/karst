# KARST's cross-validation against targets from each karyotype's own counts
# (karst 1.1.4, as_alfak(xval_target = "own"); docs/xv_target_6ae2a24.md).
#
#   Rscript scripts/xv_own_target.R --index N --chunks C     array task N of C
#   Rscript scripts/xv_own_target.R --combine                join the chunks
#
# For every KARST fit (methods karst, karst_pairwise, karst_state_space; each
# source, input and design) the held-out predictions of its stored
# cross-validation (xval.Rds, one per minobs) are compared with targets from
# each tested karyotype's own counts (karst:::sparse_joint_own_fitness(); one
# KARST fit serves every minobs, so the targets are computed once per fit):
#   xv_own_target   resampled once with seed 1, as as_alfak() draws them
#   xv_own_point    without resampling (the sensitivity analysis)
# with Pearson's r beside each. The held-out predictions do not depend on what
# they are compared with, so only the targets are new. On the landscape's
# scale, as as_alfak() puts them: (own net rate - reference + level) / divisor,
# with the efflux divisor of each karyotype and the level read back from the
# stored landscape.Rds (its mean times the divisor minus the model's relative
# fitness; one value per fit, `level_spread` its largest deviation). Only
# karyotypes with at least one cell are tested (n_cells; n_at_bound of them
# ended at the lower bound of the birth rate, mostly karyotypes the draw left
# without cells).
# Every 25th fit (and every fit with CHECK_ALL=1) is also run through
# karst::as_alfak(xval_target = "own") in full, its folds refitted:
# check_landscape (largest difference of the landscape's mean from the
# stored one) and check_xv (difference of the score from xv_own_target) must
# be at rounding level.
# Simulations (abm, misseg): Spearman's rho with the true fitness, on the
# karyotypes KARST tested with cells, of the own target (rho_truth_own,
# rho_truth_point), of KARST's default target (rho_truth_model: its estimate
# from all the data) and of its held-out prediction (rho_truth_prediction);
# and of alfak's target (rho_truth_alfak: its xval_data observation, one
# bootstrap draw of its first-stage fitness) on every karyotype alfak tested
# in the fit of the same input, design and minobs, and on those of them
# KARST tested with cells (rho_truth_alfak_common).
# Writes BENCH_RESULTS/tables/xv_own_target/chunk_<N>.Rds; --combine writes
# tables/xv_own_target.Rds and prints a summary. Needs the karst with
# sparse_joint_own_fitness() (1.1.4 on): KARST_LIB in hpc/env.sh.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
# OWN_CHUNK_DIR: another folder for the chunks (a test run)
# (hpc/env.sh passes it empty when it is not set)
chunk_dir <- Sys.getenv("OWN_CHUNK_DIR")
if (!nzchar(chunk_dir)) chunk_dir <- file.path(paths$results, "tables", "xv_own_target")
sources <- c("abm", "misseg", "salehi")

r2r <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum(obs^2)
}
spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
}
run_ok <- function(d) {
  p <- file.path(d, "run.Rds")
  file.exists(p) && identical(readRDS(p)$status, "ok") && file.exists(file.path(d, "xval.Rds"))
}
all_groups <- function() {
  sort(unlist(lapply(sources, function(src) {
    Sys.glob(file.path(paths$results, "fits", src, "*", "*", karst_xv_methods))
  })))
}

if (identical(args[1], "--combine")) {
  files <- list.files(chunk_dir, pattern = "^chunk_\\d+\\.Rds$", full.names = TRUE)
  rows <- do.call(rbind, lapply(files, readRDS))
  n_groups <- length(all_groups())
  done <- unique(paste(rows$source, rows$dataset, rows$design, rows$method))
  cat("chunks", length(files), " fits with rows", length(done), " of", n_groups, "fit folders\n")
  key <- paste(rows$source, rows$dataset, rows$design, rows$method, rows$minobs)
  if (anyDuplicated(key)) stop("duplicated fits in the chunks")
  checked <- rows[!is.na(rows$check_xv), ]
  cat("full as_alfak() checks:", nrow(checked), " largest |check_xv|", signif(max(abs(checked$check_xv)), 3),
      " largest |check_landscape|", signif(max(checked$check_landscape), 3), "\n")
  cat("largest level_spread", signif(max(rows$level_spread, na.rm = TRUE), 3), "\n")
  if (nrow(checked) && (max(abs(checked$check_xv)) > 1e-8 || max(checked$check_landscape) > 1e-8)) {
    stop("as_alfak(xval_target = \"own\") does not reproduce the scores computed here")
  }
  if (max(rows$level_spread, na.rm = TRUE) > 1e-8) stop("the level read back from landscape.Rds is not one value")
  saveRDS(rows, file.path(paths$results, "tables", "xv_own_target.Rds"))
  print(table(source = rows$source, method = rows$method, status = rows$status))
  quit(save = "no")
}

stopifnot(identical(args[1], "--index"), identical(args[3], "--chunks"))
index <- as.integer(args[2])
chunks <- as.integer(args[4])
groups <- all_groups()
mine <- groups[(seq_along(groups) - 1L) %% chunks + 1L == index]
check_every <- if (identical(Sys.getenv("CHECK_ALL"), "1")) 1L else 25L
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(chunk_dir, sprintf("chunk_%03d.Rds", index))
cat("chunk", index, "of", chunks, ":", length(mine), "fits; karst", as.character(utils::packageVersion("karst")),
    find.package("karst"), "\n")

input_first_time <- function(src, dataset, design) {
  yi <- if (src == "salehi") salehi_input(dataset, ends = design == "ends", paths) else {
    sim <- read_abm_sim(dataset, if (src == "abm") paths$abm_raw else paths$misseg_raw)
    abm_input(sim, abm_designs[abm_designs$design == design, ])
  }
  input_times(yi)[1]
}

one_group <- function(gdir, check) {
  p <- rev(strsplit(gdir, "/", fixed = TRUE)[[1]])
  method <- p[1]; design <- p[2]; dataset <- p[3]; src <- p[4]
  mdirs <- Sys.glob(file.path(gdir, "minobs_*"))
  mdirs <- mdirs[vapply(mdirs, function(d) run_ok(d) && file.exists(file.path(d, "karst_fit.Rds")), logical(1))]
  if (!length(mdirs)) return(NULL)
  fit <- readRDS(file.path(mdirs[1], "karst_fit.Rds"))
  if (inherits(fit, "karst_sparse_joint_rate_fit")) fit <- fit$fits[[fit$point]]
  xv <- lapply(mdirs, function(d) readRDS(file.path(d, "xval.Rds")))
  tested <- unique(unlist(lapply(xv, function(x) x$xval_data$k)))
  started <- proc.time()[["elapsed"]]
  own <- karst:::sparse_joint_own_fitness(fit, tested, resample = TRUE, seed = 1L)
  point <- karst:::sparse_joint_own_fitness(fit, tested, resample = FALSE)
  seconds <- proc.time()[["elapsed"]] - started
  bound <- stats::setNames(attr(own, "at_bound"), tested)
  reference <- fit$trend$reference$fitness
  sim_based <- src %in% c("abm", "misseg")
  if (sim_based) {
    sim <- read_abm_sim(dataset, if (src == "abm") paths$abm_raw else paths$misseg_raw)
    w <- abm_sim_info(dataset)$w
  }
  settings <- bench_settings[[src]]
  rows <- lapply(seq_along(mdirs), function(i) {
    d <- mdirs[i]
    minobs <- sub("minobs_", "", basename(d))
    x <- xv[[i]]$xval_data
    x <- x[stats::complete.cases(x[c("observation", "prediction")]), ]
    k <- x$k
    head <- data.frame(source = src, dataset = dataset, design = design, method = method, minobs = minobs,
                       stringsAsFactors = FALSE)
    empty <- cbind(head, xv_own = as.numeric(xv[[i]]$R2R), xv_own_target = NA_real_, pearson_own_target = NA_real_,
                   xv_own_point = NA_real_, pearson_own_point = NA_real_, n_karst = length(k), n_cells = 0L,
                   n_at_bound = NA_integer_, level_spread = NA_real_, rho_truth_own = NA_real_,
                   rho_truth_point = NA_real_, rho_truth_model = NA_real_, rho_truth_prediction = NA_real_,
                   rho_truth_alfak = NA_real_, rho_truth_alfak_common = NA_real_, check_landscape = NA_real_,
                   check_xv = NA_real_, own_seconds = seconds, status = "fewer than 2 karyotypes with cells",
                   stringsAsFactors = FALSE)
    if (length(k) < 2) return(empty)
    lsc <- readRDS(file.path(d, "landscape.Rds"))
    divisor <- if (isTRUE(bench_settings$correct_efflux)) {
      2 * karst:::sparse_joint_self_probability(karst:::parse_karyotypes(k), fit$data$beta) - 1
    } else rep(1, length(k))
    relative <- stats::predict(fit, k)$relative_fitness
    lv <- lsc$mean[match(k, lsc$k)] * divisor - relative
    level <- stats::median(lv)
    t_own <- (as.numeric(own[k]) - reference + level) / divisor
    t_point <- (as.numeric(point[k]) - reference + level) / divisor
    ok <- is.finite(t_own) & is.finite(t_point)
    out <- empty
    out$n_cells <- sum(ok)
    out$level_spread <- max(abs(lv - level))
    if (sum(ok) < 2) return(out)
    out$status <- "ok"
    out$n_at_bound <- sum(bound[k[ok]], na.rm = TRUE)
    out$xv_own_target <- r2r(t_own[ok], x$prediction[ok])
    out$pearson_own_target <- suppressWarnings(stats::cor(t_own[ok], x$prediction[ok]))
    out$xv_own_point <- r2r(t_point[ok], x$prediction[ok])
    out$pearson_own_point <- suppressWarnings(stats::cor(t_point[ok], x$prediction[ok]))
    if (sim_based) {
      truth <- grf_fitness(k[ok], w, sim$true_landscape)
      out$rho_truth_own <- spearman(t_own[ok], truth)
      out$rho_truth_point <- spearman(t_point[ok], truth)
      out$rho_truth_model <- spearman(x$observation[ok], truth)
      out$rho_truth_prediction <- spearman(x$prediction[ok], truth)
      dR <- file.path(dirname(dirname(d)), "alfak", basename(d))
      if (run_ok(dR)) {
        a <- readRDS(file.path(dR, "xval.Rds"))$xval_data
        a <- a[stats::complete.cases(a[c("observation", "prediction")]), ]
        if (nrow(a) >= 3) {
          out$rho_truth_alfak <- spearman(a$observation, grf_fitness(a$k, w, sim$true_landscape))
          common <- intersect(a$k, k[ok])
          out$rho_truth_alfak_common <- spearman(a$observation[match(common, a$k)], truth[match(common, k[ok])])
        }
      }
    }
    if (check && i == 1L) {
      run <- readRDS(file.path(d, "run.Rds"))
      as_args <- list(fit, minobs = as.numeric(minobs), nboot = 2L, n0 = settings$n0, nb = settings$nb,
                      correct_efflux = bench_settings$correct_efflux, xval_target = "own")
      if (identical(run$freq_origin, "time0")) as_args$first_time <- input_first_time(src, dataset, design)
      full <- do.call(karst::as_alfak, as_args)
      m <- match(lsc$k, full$landscape$k)
      out$check_landscape <- max(abs(full$landscape$mean[m] - lsc$mean))
      out$check_xv <- full$xval$R2R - out$xv_own_target
    }
    out
  })
  do.call(rbind, rows)
}

rows <- list()
for (j in seq_along(mine)) {
  check <- (match(mine[j], groups) - 1L) %% check_every == 0L
  r <- tryCatch(one_group(mine[j], check), error = function(e) {
    message("FAILED ", mine[j], ": ", conditionMessage(e))
    stop(e)
  })
  if (!is.null(r)) rows[[length(rows) + 1L]] <- r
  if (j %% 10L == 0L) cat(j, "of", length(mine), format(Sys.time()), "\n")
}
rows <- do.call(rbind, rows)
saveRDS(rows, out_file)
cat("chunk", index, "done:", if (is.null(rows)) 0 else nrow(rows), "rows\n")
