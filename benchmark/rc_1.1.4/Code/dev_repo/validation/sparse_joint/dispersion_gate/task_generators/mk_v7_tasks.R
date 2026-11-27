# rc freeze: GATE_RESULTS is the gate's results folder (inputs read from it,
# outputs written to it); the as-run copy read and wrote
# validation/sparse_joint/dispersion_gate/results from the repository root.
G <- Sys.getenv("GATE_RESULTS"); stopifnot(nzchar(G))
s <- read.delim(file.path(G, "step0/stats.tsv"), stringsAsFactors = FALSE)
s <- unique(s[grepl("^v7_", s$source), c("source", "dataset_id", "design", "wavelength", "multiple")])
L <- as.integer(sub("^v7_L([0-9]+)_.*$", "\\1", s$dataset_id))
stopifnot(!anyNA(L))
s$half <- ifelse(((L - 1) %% 4) < 2, "dev", "confirm")
s$depth <- sub("^.*_d([0-9]+-[0-9]+)$", "\\1", s$dataset_id)
d <- s[s$half == "dev" & s$depth %in% c("0200-0200", "0500-0700"), ]
cat("dev datasets:", nrow(d), "\n"); print(table(d$source, d$multiple)); print(table(d$wavelength, d$depth))
t <- do.call(rbind, lapply(c(1, 2, 4, 8, 16), function(m) data.frame(family = ifelse(d$source == "v7_c2", "c2", "e"), dataset_id = d$dataset_id, multiple = m)))
write.table(t, file.path(G, "step1/tasks/fits_v7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("tasks:", nrow(t), "\n")
