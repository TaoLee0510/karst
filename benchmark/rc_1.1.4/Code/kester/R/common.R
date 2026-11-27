# Copy of ../outside_domain/R/common.R (unchanged), so this study is self-contained.
# Shared helpers of the out-of-domain study (patch/outside_domain/PLAN.md).
# The karst package is used read-only: with `use_package = TRUE` its namespace
# (the round's library on the HPC); otherwise its R sources are sourced from
# the working tree into an environment, with R stand-ins for the two compiled
# label helpers (local runs without an installed karst).

od_karst <- function(use_package = requireNamespace("karst", quietly = TRUE),
                     repo = Sys.getenv("KARST_REPO", "")) {
  if (use_package) return(asNamespace("karst"))
  if (!nzchar(repo)) stop("karst is not installed and KARST_REPO is not set", call. = FALSE)
  env <- new.env()
  sys.source(file.path(repo, "R", "karyotype-utils.R"), envir = env)
  sys.source(file.path(repo, "R", "sparse-joint.R"), envir = env)
  env$parse_karyotypes <- function(labels) {
    m <- do.call(rbind, lapply(strsplit(as.character(labels), ".", fixed = TRUE), as.integer))
    storage.mode(m) <- "integer"
    m
  }
  env$format_karyotypes <- function(karyotypes) apply(as.matrix(karyotypes), 1L, paste, collapse = ".")
  env
}

od_parse <- function(labels) {
  m <- do.call(rbind, lapply(strsplit(as.character(labels), ".", fixed = TRUE), as.integer))
  storage.mode(m) <- "integer"
  m
}

# Domain of a two-time-point fit as fit_sparse_joint() builds it (defaults:
# major share 0.01, 3 cells, at most 5000 states), its majors, and each label's
# distance to the nearest major (0 inside the domain).
od_domain <- function(k, counts, min_cn = 0L, max_cn = NULL) {
  cn <- od_parse(rownames(counts))
  if (is.null(max_cn)) max_cn <- max(8L, max(cn))
  k$sparse_joint_domain(counts, 0.01, 3, 5000, min_cn, max_cn)
}
od_distance <- function(k, dom, labels) {
  d <- rep(0, length(labels))
  out <- !labels %in% dom$labels
  if (any(out)) d[out] <- k$sparse_joint_major_distance(od_parse(labels[out]), dom$major)
  d
}

# Group of each label: inside the domain, halo (out of the domain, at most
# `halo` steps from a major) or island (farther).
od_group <- function(distance, halo = 4) {
  ifelse(distance == 0, "inside", ifelse(distance <= halo, "halo", "island"))
}
od_shares <- function(freq, group) {
  s <- tapply(freq, factor(group, levels = c("inside", "halo", "island")), sum)
  s[is.na(s)] <- 0
  s / sum(freq)
}

# Count-weighted mean karyotype, change vectors and their angle (degrees).
od_centre <- function(p) drop(crossprod(p / sum(p), od_parse(names(p))))
od_angle <- function(x, y) acos(pmax(-1, pmin(1, sum(x * y) / sqrt(sum(x^2)) / sqrt(sum(y^2))))) * 180 / pi

# Exact-karyotype overlap between two frequency vectors named by karyotype.
od_overlap <- function(p, q) {
  u <- union(names(p), names(q))
  a <- p[u]; a[is.na(a)] <- 0; b <- q[u]; b[is.na(b)] <- 0
  sum(pmin(a / sum(a), b / sum(b)))
}

# Wasserstein distance (Euclidean ground distance, p = 1) between two frequency
# vectors named by karyotype. A deterministic forecast puts some mass on every
# reachable state: states below `floor` are dropped and the rest renormalised.
# revsimplex (networkflow can stall on large supports).
od_w1 <- function(p, q, floor = 1e-6) {
  p <- p[p >= floor]; q <- q[q > 0]
  transport::wasserstein(transport::wpp(od_parse(names(p)), mass = p / sum(p)),
                         transport::wpp(od_parse(names(q)), mass = q / sum(q)), method = "revsimplex")
}
