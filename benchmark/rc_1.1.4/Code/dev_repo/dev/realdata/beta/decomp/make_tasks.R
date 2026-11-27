#!/usr/bin/env Rscript
# decomp/PLAN.md: the 31 (L, F) combinations of each Salehi trajectory with
# later passages (L = fitting rate, F = forecast rate, multiples of 5e-5).
#   Rscript dev/realdata/beta/decomp/make_tasks.R   (writes decomp/tasks/decomp.tsv)
here <- file.path("dev", "realdata", "beta")
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)
xlabel <- function(m) ifelse(m < 1, paste0("x1_", round(1 / m)), paste0("x", m))   # benchmark/paper/R/fit.R
ids <- utils::read.delim(file.path(here, "decomp", "tasks", "fit_ext.tsv"), stringsAsFactors = FALSE)$id
old <- utils::read.delim(file.path(here, "tasks", "beta_horizon.tsv"), stringsAsFactors = FALSE)
combos <- unique(rbind(data.frame(L = grid, F = grid), data.frame(L = grid, F = 1), data.frame(L = 1, F = grid)))
stopifnot(nrow(combos) == 31L)
tasks <- merge(data.frame(id = ids), combos)
tasks <- tasks[order(tasks$id, tasks$L, tasks$F), ]
out <- data.frame(id = tasks$id, method = xlabel(tasks$L), minobs = 10L, design = "ends", root = old$root[1],
                  support_roots = old$support_roots[1], pm = 5e-5 * tasks$F, nb = old$nb[1],
                  out = paste0("F_", xlabel(tasks$F)), L = tasks$L, F = tasks$F, stringsAsFactors = FALSE)
utils::write.table(out, file.path(here, "decomp", "tasks", "decomp.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "tasks\n")
