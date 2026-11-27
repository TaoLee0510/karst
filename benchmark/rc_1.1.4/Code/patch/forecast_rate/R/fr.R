# Helpers of study B (PLAN.md). Needs common.R (od_*).

# One-step neighbours with at least one copy of every chromosome.
fr_nb <- function(labels) {
  cn <- od_parse(labels); out <- character(0)
  for (i in seq_len(ncol(cn))) for (s in c(-1L, 1L)) {
    m <- cn; m[, i] <- m[, i] + s; ok <- m[, i] >= 1L
    out <- c(out, do.call(paste, c(as.data.frame(m[ok, , drop = FALSE]), sep = ".")))
  }
  unique(out)
}

# List E: the training karyotypes and their one- and two-step neighbours.
fr_list_e <- function(observed) {
  one <- setdiff(fr_nb(observed), observed)
  two <- setdiff(fr_nb(one), c(observed, one))
  c(observed, one, two)
}

# Mixture of forecasts (shares on the same list), equal weights.
fr_mix <- function(forecasts) {
  m <- Reduce(`+`, lapply(forecasts, function(f) f / sum(f))) / length(forecasts)
  names(m) <- names(forecasts[[1]])
  m
}

# Scores of each forecast against the observed sample.
#   fc: named list of forecasts on list E (shares, named by karyotype)
#   parent, daughter: counts named by karyotype
#   observed: the training karyotypes ("new" = not among them)
#   dist_e: distance of each karyotype of E from the nearest major (0 inside the domain)
#   dist_of: function(labels) giving the same distance for other labels
fr_score <- function(fc, parent, daughter, observed, dist_e, dist_of) {
  c_par <- od_centre(parent); obs_change <- od_centre(daughter) - c_par
  halo <- function(freq, d) od_shares(freq, od_group(d, 4))[["halo"]]
  h_par <- halo(parent, dist_of(names(parent)))
  obs_halo_change <- halo(daughter, dist_of(names(daughter))) - h_par
  obs_new <- sum(daughter[!names(daughter) %in% observed]) / sum(daughter)
  in_e <- sum(daughter[names(daughter) %in% names(fc[[1]])]) / sum(daughter)
  do.call(rbind, lapply(names(fc), function(a) {
    f <- fc[[a]]; keep <- f > 0; fr <- f[keep]
    pred_halo_change <- halo(fr, dist_e[keep]) - h_par
    pred_new <- sum(fr[!names(fr) %in% observed]) / sum(fr)
    data.frame(arm = a, angle = od_angle(od_centre(fr) - c_par, obs_change),
               w1 = tryCatch(od_w1(fr, daughter, floor = 1e-5), error = function(e) NA_real_),
               overlap = od_overlap(fr, daughter),
               obs_halo_change = obs_halo_change, pred_halo_change = pred_halo_change,
               halo_err = abs(pred_halo_change - obs_halo_change),
               obs_new = obs_new, pred_new = pred_new, new_err = abs(pred_new - obs_new),
               daughter_in_E = in_e, stringsAsFactors = FALSE)
  }))
}
