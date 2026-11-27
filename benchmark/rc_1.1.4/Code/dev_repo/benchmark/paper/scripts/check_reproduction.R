# Check that the scoring code reproduces the paper's ABM scores.
#
#   Rscript scripts/check_reproduction.R [--out DIR] [--index N | SIM_ID ...]
#
# For each fit of the paper (PAPER_DATA/processed/ABM/<sim>/minobs_*_ntp_*),
# the scores are computed again from landscape.Rds, xval.Rds and the forecast
# replicates, and compared with the paper's abm_preds_summary.Rds and with
# the matching row of ABM_summaries.Rds. `--index N` checks the N-th
# simulation in sorted order (one array task per simulation); without
# simulations every one is checked. With `--out`, the differences are saved
# to DIR/<sim>.Rds. Stops if any difference exceeds 1e-8 or the pattern of
# missing values differs.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
take <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) return(NULL)
  value <- args[[i + 1L]]
  args <<- args[-c(i, i + 1L)]
  value
}
out <- take("--out")
index <- take("--index")
all_sims <- sort(list.files(paths$paper_abm))
sims <- if (!is.null(index)) all_sims[as.integer(index)] else if (length(args)) args else all_sims
summaries <- readRDS(file.path(paths$data, "processed", "ABM_summaries.Rds"))$dfl
key <- function(d) paste(d$w, d$abmrep, d$minobs, d$ntp)
tables <- c("df_so", "df_sc", "df_de", "df_dw", "df_a")

rows <- list()
for (sim_id in sims) {
  sim <- read_abm_sim(sim_id, paths$abm_raw)
  for (dj in list.files(file.path(paths$paper_abm, sim_id))) {
    fit_dir <- file.path(paths$paper_abm, sim_id, dj)
    saved_path <- file.path(fit_dir, "abm_preds_summary.Rds")
    if (!file.exists(saved_path)) next
    parts <- strsplit(dj, "_", fixed = TRUE)[[1]]
    info <- data.frame(w = abm_sim_info(sim_id)$w, abmrep = abm_sim_info(sim_id)$abmrep,
                       minobs = parts[2], ntp = parts[4], stringsAsFactors = FALSE)
    mine <- abm_fit_summary(fit_dir, sim, info)
    saved <- readRDS(saved_path)
    d <- c(summary = max(abs(unlist(mine$summary) - unlist(saved$summary)), na.rm = TRUE),
           vapply(tables, function(n) max(abs(as.matrix(mine[[n]]) - as.matrix(saved[[n]])), na.rm = TRUE), 0))
    same_na <- identical(is.na(unlist(mine$summary)), is.na(unlist(saved$summary))) &&
      all(vapply(tables, function(n) identical(is.na(as.matrix(mine[[n]])), is.na(as.matrix(saved[[n]]))), TRUE))
    paper_row <- summaries[key(summaries) == key(info), names(mine$summary)]
    d <- c(d, table_row = if (nrow(paper_row) == 1L) max(abs(unlist(mine$summary) - unlist(paper_row)), na.rm = TRUE)
           else NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(sim_id = sim_id, fit = dj, t(d), same_na = same_na)
  }
}
rows <- do.call(rbind, rows)
if (is.null(rows)) {
  cat("No scored paper fits for", paste(sims, collapse = ", "), "\n")
  quit(save = "no", status = 0)
}
if (!is.null(out)) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(rows, file.path(out, paste0(if (length(sims) == 1L) sims else "all", ".Rds")))
}
numbers <- setdiff(names(rows), c("sim_id", "fit", "same_na"))
cat("Fits checked:", nrow(rows), "\n")
print(sapply(rows[numbers], max, na.rm = TRUE))
cat("Missing values in the same places in all fits:", all(rows$same_na), "\n")
bad <- rows[apply(rows[numbers], 1, function(v) any(v > 1e-8, na.rm = TRUE)) | !rows$same_na, ]
if (nrow(bad)) {
  print(bad)
  stop("Scores differ from the paper's.")
}
cat("All scores match the paper's.\n")
