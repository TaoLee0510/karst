# Compose karst's "auto" arm from BENCH_RESULTS/tables/auto_picks.tsv
# (scripts/auto_picks.R): for each input, fits/<source>/<dataset>/<design>/
# karst_auto becomes a relative symbolic link to the picked method's folder
# next to it, so that every collect, score and analysis step reads it as a
# method of its own. Run again after new picks: existing karst_auto links
# are replaced (anything else of that name is refused).
#
#   Rscript scripts/compose_auto.R
#
# Afterwards the Salehi scores of the auto arm must be computed again
# (scripts/score_salehi.R reads every fit folder of a trajectory), and the
# collect and analysis steps rerun.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
picks <- utils::read.delim(file.path(paths$results, "tables", "auto_picks.tsv"), stringsAsFactors = FALSE)
stopifnot(all(picks$pick %in% c("karst", names(karst_series_methods))))
made <- 0L
missing <- character()
for (i in seq_len(nrow(picks))) {
  base <- file.path(paths$results, "fits", picks$source[i], picks$dataset[i], picks$design[i])
  if (!dir.exists(file.path(base, picks$pick[i]))) {
    missing <- c(missing, file.path(base, picks$pick[i]))
    next
  }
  link <- file.path(base, "karst_auto")
  if (file.exists(link) || nzchar(Sys.readlink(link))) {
    if (!nzchar(Sys.readlink(link))) stop("Not a link, refusing to replace: ", link)
    file.remove(link)
  }
  stopifnot(file.symlink(picks$pick[i], link))
  made <- made + 1L
}
cat("karst_auto links:", made, " picked fit folder missing:", length(missing), "\n")
if (length(missing)) print(utils::head(missing, 20))
print(table(source = picks$source, pick = picks$pick))
