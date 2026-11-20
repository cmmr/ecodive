# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Clustering the samples of a distance matrix: the `clust` argument of
# beta_div(), its metric functions, vlr(), rho_p() and phi_s(), and the
# default order of plot.ecodive_dist().
#
#   dist_result()        labels a distance matrix with its metric, then
#                        applies `clust`: the end of every metric function.
#   validate_clust()     checks `clust` (validate_args() calls it).
#   dist_clust()         reorders a distance matrix by `clust`.
#   dist_clust_order()   the order of the samples that `clust` gives.
#   dist_reorder()       a distance matrix with its samples in a new order.


# The distance matrix `d` from metric `metric` (its ID, such as 'bray'),
# with the ID as its `method` attribute, as stats::dist() and
# vegan::vegdist() set it, and the metric's display name as its title (see
# ?ecodive_display). Then clustered by `clust`, from the metric function's
# frame `env`.
dist_result <- function (d, metric, env = parent.frame()) {
  attr(d, 'method') <- metric
  d <- set_display(d, title = METRICS$name[match(metric, METRICS$id)])
  dist_clust(d, get('clust', envir = env), get('cpus', envir = env))
}



# The values `clust` can take, besides NULL and FALSE: stats::hclust()'s
# linkage methods, 'pcoa', and 'none'.
DIST_CLUST <- c(
  'average', 'complete', 'single', 'ward.D', 'ward.D2', 'mcquitty',
  'median', 'centroid', 'pcoa', 'none' )


# Checks `clust` in `env`: NULL, FALSE (made 'none'), or one of DIST_CLUST.
validate_clust <- function (env = parent.frame()) {
  clust <- get('clust', envir = env, inherits = FALSE)
  if (is.null(clust)) return (invisible())
  if (isFALSE(clust)) {
    assign('clust', 'none', envir = env)
    return (invisible())
  }
  if (!is.character(clust) || length(clust) != 1 || !clust %in% DIST_CLUST)
    stop(
      '`clust` must be NULL, FALSE, or one of ',
      paste0('"', DIST_CLUST, '"', collapse = ', '), '.', call. = FALSE )
  invisible()
}



# `d`, with its samples reordered by the clustering `clust` (checked by
# validate_clust()): attribute `order` gives the input position of each
# sample, now in order, and `clust` in its `display` list names the method,
# for plot(). NULL leaves `d` unchanged; 'none' keeps the order and sets
# `clust`, so that plot() keeps it too.
dist_clust <- function (d, clust, cpus = 1L) {

  if (is.null(clust)) return (d)

  if (!identical(clust, 'none')) {
    if (anyNA(d))
      stop(
        '`clust` needs every distance, so it cannot be used with `pairs` ',
        '(or distances that are NA).', call. = FALSE )
    ord <- dist_clust_order(d, clust, cpus)
    d   <- dist_reorder(d, ord)
    attr(d, 'order') <- ord
  }

  set_display(d, clust = clust)
}


# The `clust` that `d` was made with (from its `display` list), or NULL.
dist_clust_of <- function (d) {
  cl <- get_display(d, 'clust')
  if (is.character(cl) && length(cl) == 1 && !is.na(cl)) cl
}


# The ID of the metric that made `d`: its `method` attribute, as ecodive's
# distance functions and stats::dist() set it, or ''.
dist_metric_of <- function (d) {
  m <- attr(d, 'method', exact = TRUE)
  if (is.character(m) && length(m) == 1 && !is.na(m)) m else ''
}



# The order of the samples of `d` (a vector of their positions) that `clust`
# gives: the leaf order of stats::hclust() with that linkage, or for 'pcoa',
# their order along the first axis of pcoa(). Similarities from rho_p() are
# clustered as the dissimilarities 1 - rho.
dist_clust_order <- function (d, clust, cpus = 1L) {

  n <- attr(d, 'Size', exact = TRUE)
  if (n < 3 || identical(clust, 'none')) return (seq_len(n))

  if (identical(dist_metric_of(d), 'rho_p')) {
    a <- attributes(d)
    d <- 1 - d
    attributes(d) <- a
  }

  if (identical(clust, 'pcoa')) {
    p <- pcoa(d, n_axes = 1L, cpus = cpus)
    if (ncol(p) == 0) return (seq_len(n))
    return (order(p[, 1]))
  }

  dd <- d
  class(dd) <- 'dist'  # hclust() keeps its labels; the class needs no more
  stats::hclust(dd, method = clust)$order
}



# `d` with its samples in the order `ord` (a permutation of their positions):
# element (i, j) of the result is element (ord[i], ord[j]) of `d`. Keeps
# every attribute, with `Labels` reordered; samples without labels are
# labelled with their positions in `d`, so that each stays identifiable.
# Filled one column at a time, so that it needs little memory besides the
# result.
dist_reorder <- function (d, ord) {

  n   <- attr(d, 'Size', exact = TRUE)
  out <- numeric(length(d))
  pos <- 0L

  for (j in seq_len(n - 1L)) {
    r  <- (j + 1L):n
    a  <- ord[r]
    b  <- ord[[j]]
    lo <- pmin(a, b)
    hi <- pmax(a, b)
    k  <- length(r)
    out[pos + seq_len(k)] <- d[(lo - 1) * (n - lo / 2) + (hi - lo)]
    pos <- pos + k
  }

  att <- attributes(d)
  att$Labels <- if (is.null(att$Labels)) as.character(ord) else att$Labels[ord]
  attributes(out) <- att
  out
}
