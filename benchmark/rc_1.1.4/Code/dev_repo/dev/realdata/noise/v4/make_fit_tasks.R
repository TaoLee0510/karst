#!/usr/bin/env Rscript
# Fit task tables of dev/realdata/noise/PLAN_v4.md from the screen records:
# noise class T0, T1, T2, T2r; higher_beta class T0, T1 (report only).
#   Rscript dev/realdata/noise/v4/make_fit_tasks.R
# Writes dev/realdata/noise/v4/tasks/fits_{pilot,dev,confirm}.tsv and check.tsv;
# prints the class counts.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
root <- file.path(dev_results, "realdata", "noise", "v4", "screen")
dirs <- list.files(root, pattern = "^(eset|c2)_")
recs <- do.call(rbind, lapply(dirs, function(dd) {
  do.call(rbind, lapply(list.files(file.path(root, dd), full.names = TRUE), function(f) {
    x <- readRDS(f)
    data.frame(kind = x$kind, set = x$set, variant = x$variant, dataset_id = x$dataset_id, landscape_id = x$landscape_id,
               half = x$half, depth = x$depth, status = x$status, class = if (identical(x$status, "ok")) x$class else NA,
               stringsAsFactors = FALSE)
  }))
}))
cat("screened:", nrow(recs), " failed:", sum(recs$status != "ok"), "\n")
print(table(paste(recs$kind, recs$set), recs$class, useNA = "ifany"))
methods <- list(noise = c("T0", "T1", "T2", "T2r"), higher_beta = c("T0", "T1"))   # T2r: correction 16:24
rows <- do.call(rbind, lapply(names(methods), function(cl) {
  r <- recs[recs$status == "ok" & recs$class == cl, ]
  if (!nrow(r)) return(NULL)
  g <- r[rep(seq_len(nrow(r)), each = length(methods[[cl]])), ]
  g$method <- rep(methods[[cl]], times = nrow(r))
  g
}))
rows$keep <- 0L
cols <- c("kind", "set", "variant", "dataset_id", "method", "keep")
out <- file.path(dev_repo, "dev", "realdata", "noise", "v4", "tasks")
pilot <- rows[rows$landscape_id %in% c("nz_L01", "nz_c2_L01") & rows$half == "dev", ]
write.table(pilot[cols], file.path(out, "fits_pilot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (h in c("dev", "confirm")) {
  x <- rows[rows$half == h, ]
  x <- x[order(x$method != "T0", x$kind, x$set, x$variant, x$dataset_id), ]
  write.table(x[cols], file.path(out, sprintf("fits_%s.tsv", h)), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(h, "fits:", nrow(x), "\n")
}
cat("pilot fits:", nrow(pilot), "\n")
# identity check: two e-set and two C2 datasets without injection, from the development half
ok <- recs[recs$status == "ok" & recs$variant == "none" & recs$half == "dev", ]
pick <- rbind(utils::head(ok[ok$kind == "eset" & ok$set == "e1", ], 1), utils::head(ok[ok$kind == "eset" & ok$set == "e4", ], 1),
              utils::head(ok[ok$kind == "c2" & ok$depth == 300, ], 1), utils::head(ok[ok$kind == "c2" & ok$depth == 1000, ], 1))
chk <- pick[rep(seq_len(nrow(pick)), each = 2), ]
chk$method <- rep(c("T0", "T2_phi1"), times = nrow(pick)); chk$keep <- 1L
write.table(chk[cols], file.path(out, "check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
base <- pick; base$method <- "T0"; base$keep <- 1L
write.table(base[cols], file.path(out, "check_base.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
