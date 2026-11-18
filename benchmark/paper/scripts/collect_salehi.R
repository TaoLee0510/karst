# Collect the Salehi scores, the method-independent references of Fig. 2 and
# the novel-karyotype regressions of Fig. 3.
#
#   Rscript scripts/collect_salehi.R
#
# Writes to BENCH_RESULTS/tables/:
#   salehi_scores.Rds     the rows of 03 for every fit (scores/salehi/*.Rds)
#                         with design and method, plus the paper's own rows
#                         (alfak_outputs_proc) as method "alfak_paper",
#                         design "all"
#   salehi_reference.Rds  list(sisters, unrelated) of 04: forked sister
#                         passages and pairs of unrelated lineages
#   salehi_novel.Rds      05's logistic regression of novel karyotype
#                         emergence for every fit and daughter passage,
#                         one row per coefficient (all fits; 05's selection
#                         by CV score is applied in the analysis)
#   salehi_inputs.Rds     per trajectory and design (ends: first and last
#                         passage; all: every passage), before any fit: the
#                         share of each passage inside karst's model domain
#                         and whether the input is analysable
#                         (R/coverage.R; needs karst on the library path)
#   salehi_runs.Rds       status, message and error class of every fit
# Uses SLURM_CPUS_PER_TASK cores for the regressions.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
out <- file.path(paths$results, "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
meta <- utils::read.csv(paths$salehi_meta)
lineages <- readRDS(paths$salehi_lineages)
input_files <- list.files(paths$salehi_inputs)
read_input <- local({
  cache <- new.env()
  function(file) {
    if (is.null(cache[[file]])) cache[[file]] <- readRDS(file.path(paths$salehi_inputs, file))$x
    cache[[file]]
  }
})

## Inputs and fit status --------------------------------------------------------
inputs <- salehi_input_coverage(paths)
saveRDS(inputs, file.path(out, "salehi_inputs.Rds"))
cat("trajectories:", sum(inputs$design == "ends"), " not analysable (ends):",
    sum(inputs$analysable[inputs$design == "ends"] %in% FALSE), " (all):",
    sum(inputs$analysable[inputs$design == "all"] %in% FALSE), " not judged (all):",
    sum(is.na(inputs$analysable[inputs$design == "all"])), "\n")
runs <- do.call(rbind, lapply(Sys.glob(file.path(paths$results, "fits", "salehi", "*", "*", "*", "minobs_*")),
                              function(d) {
  p <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
  run <- file.path(d, "run.Rds")
  r <- if (file.exists(run)) readRDS(run) else list(status = "missing", message = NA_character_)
  data.frame(fi = p[4], design = p[3], method = p[2], min_obs = sub("minobs_", "", p[1]), status = r$status,
             message = first_or_na(r$message, NA_character_),
             condition_class = first_or_na(r$condition_class, NA_character_),
             fit_seconds = first_or_na(r$seconds), karst_fit_seconds = first_or_na(r$fit_seconds),
             samples = first_or_na(r$samples, NA_integer_), karst_fit_peak_mb = first_or_na(r$fit_peak_mb),
             intermediate_weight = first_or_na(r$intermediate_weight),
             intermediate_dispersion = first_or_na(r$intermediate_dispersion),
        series = first_or_na(r$series, NA_character_), drift = first_or_na(r$drift), samples_fitted = first_or_na(r$samples_fitted, NA_integer_),
             dropped_samples = if (is.null(r$dropped_samples)) NA_integer_ else NROW(r$dropped_samples),
             # the first value where a multi-time-point fit reports one per interval
             history = first_or_na(r$history), convergence = first_or_na(r$convergence),
             dt_convergence = first_or_na(r$history_profile$convergence),
             candidates = if (is.null(r$history_profile)) NA_integer_ else nrow(r$history_profile),
             candidates_unconverged = if (is.null(r$history_profile)) NA_integer_ else
               sum(r$history_profile$convergence != 0),
             coverage = if (is.null(r$domain_coverage)) NA_real_ else min(r$domain_coverage$fraction),
             stringsAsFactors = FALSE)
}))
saveRDS(runs, file.path(out, "salehi_runs.Rds"))

## Scores ---------------------------------------------------------------------
scores <- do.call(rbind, lapply(list.files(file.path(paths$results, "scores", "salehi"), full.names = TRUE), readRDS))
paper <- do.call(rbind, lapply(list.files(file.path(paths$data, "processed", "salehi", "alfak_outputs_proc"),
                                          full.names = TRUE), readRDS))
cols <- setdiff(names(scores), c("design", "method"))
# Columns beyond 03 (denoised baseline, entropies) are not in the paper's rows.
for (n in setdiff(cols, names(paper))) paper[[n]] <- NA_real_
scores <- rbind(scores, cbind(design = "all", method = "alfak_paper", paper[cols], stringsAsFactors = FALSE))
saveRDS(scores, file.path(out, "salehi_scores.Rds"))
cat("score rows:", nrow(scores), "\n")
print(table(scores$method, scores$design))

## References of Fig. 2c, g-h (04) ----------------------------------------------
pairs_list <- lapply(meta$uid, function(id) meta[meta$parent == id & !is.na(meta$parent), ])
pairs_list <- pairs_list[vapply(pairs_list, nrow, integer(1)) > 1]
pairs_df <- do.call(rbind, lapply(pairs_list, function(xi) {
  combos <- utils::combn(nrow(xi), 2)
  do.call(rbind, lapply(seq_len(ncol(combos)), function(ci) {
    i <- combos[1, ci]
    j <- combos[2, ci]
    fi <- input_files[grepl(paste(c(xi$datasetname[i], xi$timepoint[i], "l_2"), collapse = "_"), input_files)]
    fj <- input_files[grepl(paste(c(xi$datasetname[j], xi$timepoint[j], "l_2"), collapse = "_"), input_files)]
    if (length(fi) == 0 || length(fj) == 0) return(NULL)
    data.frame(fi = fi, fj = fj, on_treat = paste0(xi$on_treatment[c(i, j)], collapse = ""), pdx = xi$PDX_id[1],
               stringsAsFactors = FALSE)
  }))
}))
sisters <- do.call(rbind, lapply(seq_len(nrow(pairs_df)), function(i) {
  cbind(pairs_df[i, , drop = FALSE], salehi_sister_metrics(read_input(pairs_df$fi[i]), read_input(pairs_df$fj[i])))
}))
sisters$lineage <- c(SA1035 = "SA1035", SA532 = "SA532", SA609 = "SA609", SA906 = "p53 k.o", SA039 = "p53 w.t",
                     SA535 = "SA535")[sisters$pdx]

mm <- meta
mm$fi <- lapply(paste0(mm$datasetname, "_", mm$timepoint, "_l_2"), function(p) input_files[grepl(p, input_files)])
mm <- mm[vapply(mm$fi, length, integer(1)) > 0, ]
mm$fi <- vapply(mm$fi, `[`, character(1), 1)
combos <- utils::combn(seq_len(nrow(mm)), 2)
combos <- combos[, apply(combos, 2, function(ci) length(unique(mm$PDX_id[ci])) > 1), drop = FALSE]
unrelated <- do.call(rbind, lapply(seq_len(ncol(combos)), function(i) {
  j1 <- combos[1, i]
  j2 <- combos[2, i]
  data.frame(PDX1 = mm$PDX_id[j1], PDX2 = mm$PDX_id[j2],
             angle = salehi_pair_angle(read_input(mm$fi[j1]), read_input(mm$fi[j2])), stringsAsFactors = FALSE)
}))
saveRDS(list(sisters = sisters, unrelated = unrelated), file.path(out, "salehi_reference.Rds"))
cat("sister pairs:", nrow(pairs_df), " unrelated pairs:", nrow(unrelated), "\n")

## Novel karyotypes, Fig. 3 (05) ------------------------------------------------
ours <- Sys.glob(file.path(paths$results, "fits", "salehi", "*", "*", "*", "minobs_*"))
ours <- ours[file.exists(file.path(ours, "landscape.Rds")) & file.exists(file.path(ours, "run.Rds"))]
ours <- ours[vapply(ours, function(d) identical(readRDS(file.path(d, "run.Rds"))$status, "ok"), logical(1))]
info_ours <- do.call(rbind, lapply(ours, function(d) {
  p <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
  data.frame(dir = d, fi = p[4], design = p[3], method = p[2], min_obs = sub("minobs_", "", p[1]),
             stringsAsFactors = FALSE)
}))
theirs <- Sys.glob(file.path(paths$paper_salehi, "*", "minobs_*"))
theirs <- theirs[file.exists(file.path(theirs, "landscape.Rds")) & file.exists(file.path(theirs, "xval.Rds"))]
info_paper <- data.frame(dir = theirs, fi = basename(dirname(theirs)), design = "all", method = "alfak_paper",
                         min_obs = sub("minobs_", "", basename(theirs)), stringsAsFactors = FALSE)
fits <- rbind(info_ours, info_paper)
fits <- fits[fits$fi %in% names(lineages), ]

novel_rows <- function(i) {
  fit <- fits[i, ]
  decs <- salehi_daughter_inputs(fit$fi, lineages, meta)
  # A fit without daughters, or whose regressions all fail, gives a 0-row
  # table, so that NULL marks a worker that did not deliver.
  if (!length(decs)) return(data.frame())
  lscape <- readRDS(file.path(fit$dir, "landscape.Rds"))
  rows <- do.call(rbind, lapply(decs, function(dec) {
    tryCatch({
      test <- salehi_novel_data(lscape, readRDS(file.path(paths$salehi_inputs, paste0(dec, ".Rds")))$x)
      res <- summary(suppressWarnings(stats::glm(y ~ ., data = test, family = stats::binomial)))$coefficients
      data.frame(fit[c("fi", "design", "method", "min_obs")], dec_id = dec, ids = rownames(res),
                 Estimate = res[, 1], p = res[, 4], n = nrow(test), n_new = sum(test$y), row.names = NULL,
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  }))
  if (is.null(rows)) data.frame() else rows
}
novel <- parallel::mclapply(seq_len(nrow(fits)), novel_rows, mc.cores = cores)
# A worker killed (e.g. for memory) leaves NULL or try-error values: stop
# rather than drop fits.
lost <- vapply(novel, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(lost)) stop(sum(lost), " fits lost their regressions (a worker failed, e.g. out of memory); use fewer cores.")
novel <- do.call(rbind, novel)
saveRDS(novel, file.path(out, "salehi_novel.Rds"))
cat("fits for the regressions:", nrow(fits), " regressions:", length(unique(paste(novel$fi, novel$design,
    novel$method, novel$min_obs, novel$dec_id))), "\n")
