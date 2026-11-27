# rc freeze: GATE_RESULTS is the gate's results folder (inputs read from it,
# outputs written to it); the as-run copy read and wrote
# validation/sparse_joint/dispersion_gate/results from the repository root.
G <- Sys.getenv("GATE_RESULTS"); stopifnot(nzchar(G))
s <- read.delim(file.path(G, "step0/stats.tsv"), stringsAsFactors = FALSE)
f <- read.delim(file.path(G, "step1/tasks/fits_v7.tsv"), stringsAsFactors = FALSE)
ids <- unique(f$dataset_id)
tau <- c("100" = 0.602, "300" = 0.4015)
s0 <- read.delim(file.path(G, "step0/s0.tsv"))
tau <- stats::setNames(s0$tau[s0$statistic == "once"], s0$n[s0$statistic == "once"])
print(tau)
s <- s[s$dataset_id %in% ids, ]
g <- do.call(rbind, lapply(split(s, s$dataset_id), function(x) {
  use <- if (all(x$depth_t0 >= 300 & x$depth_t1 >= 300)) 300 else 100
  r <- x[x$n == use, ]
  data.frame(dataset_id = r$dataset_id, source = r$source, multiple = r$multiple, wavelength = r$wavelength, n = use,
             once = r$once, raised = r$once >= tau[as.character(use)])
}))
print(with(g, table(source, multiple, raised)))
write.table(g, file.path(G, "step1/gate_v7_dev.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
r <- g[g$raised, ]
t <- do.call(rbind, lapply(c(1, 2, 4, 8), function(m) data.frame(family = ifelse(r$source == "v7_c2", "c2", "e"), dataset_id = r$dataset_id, multiple = m)))
write.table(t, file.path(G, "step1/tasks/pick_v7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("raised:", nrow(r), " pick tasks:", nrow(t), "\n")
