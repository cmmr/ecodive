# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Associations between the features of a count table: methods whose result
# is a `dist` over features (vlr(), rho_p(), phi_s(), sparcc()). Facts about
# each method live here, keyed by the result's `method` attribute, never as
# flags on the result (plans/README.md, "Result attributes and plots"):
#
#   name      Display name: the result's title (see ?ecodive_display).
#   id        The `method` attribute of its results.
#   alt_ids   Other names, space-separated (none yet).
#   kind      'similarity' (larger is closer: rho_p, sparcc) or
#             'dissimilarity' (vlr, phi_s).
#   lower,    The range of its values. Numerical approximations may stray
#   upper     outside it (sparcc() warns when they do).
#   fisher_z  TRUE for Pearson-type correlations, for which Fisher's z
#             transform is meaningful.
#   sparse    TRUE when zeros mean "no edge was selected" (none yet).
#   strict    TRUE when its own cutoff keeps |r| > threshold rather than >=
#             (none yet).
#
# A row here means "over features": dist_clust_order() clusters a
# similarity as 1 - x, and plot.ecodive_dist() labels the axes "features"
# and draws a similarity on a diverging scale.
FEATURE_ASSOCS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                     id,     alt_ids, kind,          lower, upper, fisher_z, sparse, strict
      Variance of Log Ratios,   vlr,           , dissimilarity,     0,   Inf,    FALSE,  FALSE,  FALSE
      Proportionality (rho_p),  rho_p,         , similarity,       -1,     1,    FALSE,  FALSE,  FALSE
      Proportionality (phi_s),  phi_s,         , dissimilarity,     0,   Inf,    FALSE,  FALSE,  FALSE
      SparCC Correlation,       sparcc,        , similarity,       -1,     1,     TRUE,  FALSE,  FALSE
  " )

  df$alt_ids <- as.character(df$alt_ids)
  return (df)
})


# The FEATURE_ASSOCS row (a list) of association method `id`, or NULL.
feature_assoc_info <- function (id) {
  i <- match(id, FEATURE_ASSOCS$id)
  if (is.na(i)) return (NULL)
  as.list(FEATURE_ASSOCS[i, , drop = FALSE])
}


# TRUE if `d` is a similarity over features, by its `method` attribute.
is_similarity <- function (d) {
  info <- feature_assoc_info(dist_metric_of(d))
  !is.null(info) && identical(info$kind, 'similarity')
}
