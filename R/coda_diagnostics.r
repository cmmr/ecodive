# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# variation_array() and zero_pattern(): two first looks at a composition.
# Clean room: written from Aitchison's (1986) definitions and the identities
# in plans/coda_diagnostic_plots.md; compositions and zCompositions (GPL)
# were run for their outputs only. The variances are vlr()'s C worker.



#' Variation Array
#'
#' Aitchison's variation array, the standard first summary of a
#' composition: for every pair of features, the variance of their log ratio
#' (above the diagonal) and its mean (below), with each feature's share of
#' the total variance. Small variances mark nearly proportional pairs.
#'
#' @inheritParams vlr
#'
#' @param clust  Whether and how to reorder the features: `NULL` (the
#'        default) or `FALSE` to keep their input order, or one of
#'        [stats::hclust()]'s linkage methods (`"average"`, `"complete"`,
#'        `"single"`, `"ward.D"`, `"ward.D2"`, `"mcquitty"`, `"median"`,
#'        `"centroid"`), or `"pcoa"` to order them along the first axis of
#'        [pcoa()], each of the distances \eqn{\sqrt{\tau_{ij}}}. The result's
#'        rows and columns are then reordered, and its `order` attribute
#'        gives each feature's position in the input. With `NULL`, [plot()]
#'        clusters the features itself (average linkage), without changing
#'        the result.
#'
#' @return A \eqn{D \times D}{D x D} numeric matrix of class
#'   `c("ecodive_variation_array", "matrix", "array")`, with the features as
#'   its row and column names (their positions when `counts` has no feature
#'   names). Element `[i, j]` is
#'
#'   * above the diagonal (\eqn{i < j}): \eqn{\tau_{ij}}, the variance of
#'     \eqn{\log(x_i / x_j)}, which is [vlr()];
#'   * below the diagonal (\eqn{i > j}): \eqn{\xi_{ij}}, the mean of
#'     \eqn{\log(x_i / x_j)}, **row over column**: positive when the row's
#'     feature is the more abundant, on average in log terms;
#'   * on the diagonal: `NA`.
#'
#'   Attributes:
#'
#'   * `clr_var`: each feature's clr variance, named;
#'   * `total_var`: the total variance, their sum;
#'   * `pseudocount`: the pseudocount added before taking logs (`0` for
#'     none);
#'   * `n`: the number of samples;
#'   * `order` (with `clust`): each feature's position in the input;
#'   * `method`: `"variation_array"`;
#'   * `display`: its title, `"Variation Array"`, and `clust` (see
#'     [ecodive_display]).
#'
#'   `print()` shows the first 10 features and the total variance;
#'   `plot()` draws the array (see [plot.ecodive_variation_array()]).
#'   Arithmetic, and subsetting with `[`, give a plain matrix.
#'
#' @details
#' With \eqn{\ell_{si} = \log(x_{si} + c)} for sample \eqn{s}, feature
#' \eqn{i} and pseudocount \eqn{c} (natural logs), over \eqn{n} samples and
#' \eqn{D} features:
#' \deqn{\tau_{ij} = \frac{1}{n - 1} \sum_{s=1}^{n} \left(\ell_{si} -
#' \ell_{sj} - \overline{\ell_i} + \overline{\ell_j}\right)^2, \qquad
#' \xi_{ij} = \overline{\ell_i} - \overline{\ell_j}}
#' with \eqn{\overline{\ell_i}} feature \eqn{i}'s mean log over the samples.
#' \eqn{\tau} is the same as [vlr()]'s, computed by the same code, and so
#' independent of each sample's total, and of the other features. The means
#' \eqn{\xi_{ij}} depend on the samples' totals only through their common
#' scale, which cancels in the ratio.
#'
#' **Total and clr variances.** The total variance of the composition is
#' \deqn{\mathrm{totvar} = \frac{1}{2D} \sum_{i,j} \tau_{ij}}
#' (Aitchison 1986), the sum of the variances of its clr coordinates. Each
#' feature's clr variance follows from the \eqn{\tau}s alone,
#' \deqn{\mathrm{var}(\mathrm{clr}_i) = \frac{1}{D}\Big(\sum_j \tau_{ij} -
#' \mathrm{totvar}\Big)}
#' since \eqn{\sum_j \tau_{ij} = D\,\mathrm{var}(\mathrm{clr}_i) +
#' \mathrm{totvar}} (the clr coordinates sum to 0). They are the
#' `clr_var` attribute, and they equal the column variances of
#' `normalize(counts, "clr")` with the same pseudocount. A feature with a
#' large share of the total variance drives the clr PCA (`pca(counts, norm
#' = "clr")`) and the Aitchison distances ([aitchison()]) between samples.
#' Two displays split the same total variance: the covariance biplot
#' (`plot(pca(...), scaling = "covariance")`), whose rays' lengths
#' approximate the clr standard deviations, and the CoDa-dendrogram
#' (Pawlowsky-Glahn & Egozcue 2011), by balance.
#'
#' **Sign convention.** Aitchison (1986) and Pawlowsky-Glahn & Egozcue
#' (2011) put the means below the diagonal; ecodive's are of the row's
#' feature over the column's, so that \eqn{\xi_{ij} = -\xi_{ji}}. compositions'
#' `variation()` gives \eqn{\tau} alone, as a full symmetric matrix.
#'
#' @inheritSection vlr Zeros
#'
#' @section Size:
#'   The result is a dense \eqn{D \times D}{D x D} matrix: 2,000 features
#'   take 32 MB, and more give a warning. Filter out rare features first;
#'   the table itself stays sparse.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman & Hall.
#'
#' Pawlowsky-Glahn, V., & Egozcue, J. J. (2011). Exploring compositional
#' data with the CoDa-dendrogram. *Austrian Journal of Statistics*, 40,
#' 103–113. \doi{10.17713/ajs.v40i1&2.202}
#'
#' Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., & Bähler,
#' J. (2015). Proportionality: a valid alternative to correlation for
#' relative data. *PLoS Computational Biology*, 11, e1004075.
#' \doi{10.1371/journal.pcbi.1004075}
#'
#' @seealso [vlr()], the variances alone as a `dist`; [zero_pattern()];
#'   [normalize()] (`"clr"`); [pca()]; [replace_zeros()].
#'
#' @export
#' @examples
#'     # Zeros replaced first
#'     x  <- replace_zeros(ex_counts, 'mult')
#'     va <- variation_array(x)
#'     va
#'
#'     # The upper triangle is vlr(); the clr variances sum to the total
#'     attr(va, 'total_var')
#'     sum(attr(va, 'clr_var'))
#'
#'     plot(va)
#'
#'     # Clustered when computed: the result is reordered
#'     va <- variation_array(x, clust = 'ward.D2')
#'     attr(va, 'order')
#'
variation_array <- function (
    counts, margin = 1L, pseudocount = NULL, clust = NULL, cpus = n_cpus() ) {

  what <- 'the variation array'
  validate_clust()
  info <- log_ratio_input(what, zero_ok = TRUE, dense_out = FALSE)

  if (info$n < 2L)
    stop('`counts` must have at least 2 samples for ', what, '.', call. = FALSE)

  D <- info$D
  if (D > 2000L)
    warning(
      sprintf(
        'The result is a dense %d x %d matrix (%.2f GB). Filter out rare features first.',
        D, D, as.double(D) * D * 8 / 1e9 ),
      call. = FALSE )

  labels <- if (is.null(info$features)) as.character(seq_len(D)) else info$features

  # vlr()'s worker, which also returns each feature's mean log.
  tau  <- .Call(
    C_proportionality, counts, margin, 4L, NULL, seq_len(D) - 1L,
    labels, pseudocount, cpus )
  mlog <- attr(tau, 'mean_log', exact = TRUE)
  attr(tau, 'mean_log') <- NULL

  ord <- seq_len(D)
  if (!is.null(clust) && !identical(clust, 'none'))
    ord <- dist_clust_order(sqrt(tau), clust, cpus)

  tm <- as.matrix(tau)
  dimnames(tm) <- NULL

  # Each feature's clr variance, from the row sums of tau.
  totvar  <- sum(tm) / (2 * D)
  clr_var <- (rowSums(tm) - totvar) / D

  res <- tm[ord, ord, drop = FALSE]
  lo  <- lower.tri(res)
  res[lo]  <- outer(mlog[ord], mlog[ord], '-')[lo]
  diag(res) <- NA_real_
  dimnames(res) <- list(labels[ord], labels[ord])

  attr(res, 'clr_var')     <- stats::setNames(clr_var[ord], labels[ord])
  attr(res, 'total_var')   <- totvar
  attr(res, 'pseudocount') <- pseudocount
  attr(res, 'n')           <- info$n
  if (!is.null(clust) && !identical(clust, 'none'))
    attr(res, 'order') <- ord
  attr(res, 'method')      <- 'variation_array'
  class(res) <- c('ecodive_variation_array', 'matrix', 'array')

  res <- set_display(res, title = 'Variation Array')
  if (!is.null(clust)) res <- set_display(res, clust = clust)
  res
}



# A variation array as a plain matrix with its dimnames.
va_plain <- function (x) {
  if (!inherits(x, 'ecodive_variation_array')) return (x)
  dn <- dimnames(x)
  d  <- dim(x)
  attributes(x) <- NULL
  dim(x) <- d
  dimnames(x) <- dn
  x
}


#' @export
Ops.ecodive_variation_array <- function (e1, e2) {
  e1 <- va_plain(e1)
  if (missing(e2)) return (get(.Generic)(e1))
  get(.Generic)(e1, va_plain(e2))
}


#' @export
Math.ecodive_variation_array <- function (x, ...) {
  get(.Generic)(va_plain(x), ...)
}


#' @export
print.ecodive_variation_array <- function (x, ...) {

  D     <- nrow(x)
  n     <- attr(x, 'n', exact = TRUE)
  k     <- min(D, 10L)
  shown <- seq_len(k)

  cat(sprintf(
    '# Variation array of %d features (%s samples), total variance %s\n',
    D, if (is.null(n)) '?' else n, plot_num(attr(x, 'total_var', exact = TRUE)) ))
  cat('# upper: var(log ratio); lower: mean log(row / column)\n')

  m <- signif(va_plain(x)[shown, shown, drop = FALSE], 3)
  print.default(m, na.print = '', ...)
  if (k < D)
    cat(sprintf('# ... and %d more features\n', D - k))

  cv <- attr(x, 'clr_var', exact = TRUE)
  tv <- attr(x, 'total_var', exact = TRUE)
  if (length(cv) == D && isTRUE(tv > 0)) {
    cat('# % of the total variance in each clr coordinate:\n')
    print.default(round(100 * cv[shown] / tv, 1), ...)
    if (k < D) cat(sprintf('# ... and %d more features\n', D - k))
  }

  invisible(x)
}



#' Plot a variation array
#'
#' Draws a [variation_array()] as two heatmaps on one grid: the variances
#' of the log ratios above the diagonal, their means below, and beside each
#' row the feature's share of the total variance.
#'
#' @param x  A variation array, from [variation_array()].
#'
#' @param digits  The number of decimal places of the values written in the
#'        cells, when there are at most 15 features. Default `2`.
#'
#' @param key  Whether to draw the two colour keys on the right. Default
#'        `TRUE`.
#'
#' @param ...  Passed on to [graphics::image()], such as `main`, `sub`, or
#'        `useRaster`.
#'
#' @details
#' The features run from the first at the top to the last at the bottom,
#' and from left to right in the same order. Their order is the one `x` was
#' made with (its `clust`), or else average-linkage clustering of the
#' distances \eqn{\sqrt{\tau_{ij}}}, so that nearly proportional features
#' sit together. `x` itself isn't changed.
#'
#' Cells above the diagonal are coloured by \eqn{\tau_{ij}}, from light
#' (proportional) to dark red; cells below by the mean log ratio of the
#' row's feature over the column's, blue for negative and red for positive,
#' symmetric about 0. The bars on the right are each feature's clr variance
#' as a percentage of the total variance ("% clr var", as in
#' Pawlowsky-Glahn & Egozcue 2011).
#'
#' @return `x`, invisibly.
#'
#' @seealso [variation_array()], [plot.ecodive_dist()] for [vlr()] alone.
#'
#' @export
#' @examples
#'     x <- replace_zeros(ex_counts, 'mult')
#'     plot(variation_array(x))
#'
plot.ecodive_variation_array <- function (x, digits = 2L, key = TRUE, ...) {

  args <- list(...)
  if (!isTRUE(key) && !isFALSE(key))
    stop('`key` must be TRUE or FALSE.', call. = FALSE)
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0)
    stop('`digits` must be a number of at least 0.', call. = FALSE)

  D <- nrow(x)
  if (length(D) != 1 || is.na(D) || D < 2 || ncol(x) != D)
    stop('`x` must be a variation array of at least 2 features.', call. = FALSE)

  m      <- va_plain(x)
  labels <- rownames(m)
  if (is.null(labels)) labels <- as.character(seq_len(D))

  # tau as a full symmetric matrix (from the upper triangle).
  up  <- upper.tri(m)
  lo  <- lower.tri(m)
  tau <- matrix(0, D, D)
  tau[up] <- m[up]
  tau     <- tau + t(tau)


  #--------------------------------------------------
  # The order: x's own when it was made with `clust`, else clustered here.
  #--------------------------------------------------
  cl  <- dist_clust_of(x)
  ord <- seq_len(D)
  how <- dist_plot_order_name(cl)
  if (is.null(cl) && all(is.finite(tau))) {
    d   <- stats::as.dist(sqrt(pmax(tau, 0)))
    ord <- dist_clust_order(d, 'average')
    how <- dist_plot_order_name('average')
  }

  # The reordered array: upper triangle from tau, lower from the means,
  # which reorder as a skew-symmetric matrix (xi_ji = -xi_ij).
  xi <- matrix(0, D, D)
  xi[lo] <- m[lo]
  xi     <- xi - t(xi)
  tau    <- tau[ord, ord, drop = FALSE]
  xi     <- xi[ord, ord, drop = FALSE]
  labels <- labels[ord]

  zu <- tau; zu[!up] <- NA
  zl <- xi;  zl[!lo] <- NA


  #--------------------------------------------------
  # Colour scales.
  #--------------------------------------------------
  tu  <- zu[is.finite(zu)]
  tl  <- zl[is.finite(zl)]
  top_u <- if (length(tu)) max(tu) else 0
  top_l <- if (length(tl)) max(abs(tl)) else 0
  if (!(top_u > 0)) top_u <- 1
  if (!(top_l > 0)) top_l <- 1
  pal_u <- grDevices::hcl.colors(64, 'YlOrRd', rev = TRUE)
  pal_l <- grDevices::hcl.colors(65, 'Blue-Red 3')
  brk_u <- seq(0, top_u, length.out = length(pal_u) + 1L)
  brk_l <- seq(-top_l, top_l, length.out = length(pal_l) + 1L)


  #--------------------------------------------------
  # Margins: the row labels on the left; the bars and keys on the right.
  #--------------------------------------------------
  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  mai <- graphics::par('mai')
  mai[2L] <- max(mai[2L], max(graphics::strwidth(labels, 'inch', cex = 0.7), 0, na.rm = TRUE) + 0.25)
  mai[4L] <- max(mai[4L], if (key) 2.1 else 1.1)
  graphics::par(mai = mai)

  # image() draws its first row at the bottom: rows reversed.
  rows <- rev(seq_len(D))
  rast <- D > 100 && identical(grDevices::dev.capabilities('rasterImage')$rasterImage, 'yes')

  n   <- attr(x, 'n', exact = TRUE)
  tv  <- attr(x, 'total_var', exact = TRUE)
  sub <- paste0(
    D, ' features',
    if (length(n) == 1) paste0(', ', n, ' samples'),
    if (length(tv) == 1) paste0(', total variance ', plot_num(tv)),
    if (!is.null(how)) paste0(', ', how),
    '\nupper: var(log ratio); lower: mean log(row / column)' )

  plot_do(
    graphics::image, c(list(seq_len(D), seq_len(D), t(zu[rows, , drop = FALSE])), args),
    list(
      col = pal_u, breaks = brk_u, axes = FALSE, xlab = '', ylab = '',
      useRaster = rast, main = plot_title(x, 'Variation Array'), sub = sub ))
  graphics::image(
    seq_len(D), seq_len(D), t(zl[rows, , drop = FALSE]),
    col = pal_l, breaks = brk_l, add = TRUE, useRaster = rast )

  graphics::axis(2, at = D - seq_len(D) + 1, labels = labels, las = 1, tick = FALSE, cex.axis = 0.7)
  graphics::box()

  # Each cell's value, for small arrays; white on the darkest cells.
  if (D <= 15) {
    xy <- expand.grid(r = seq_len(D), c = seq_len(D))
    xy <- xy[xy$r != xy$c, , drop = FALSE]
    v  <- ifelse(xy$r < xy$c, tau[cbind(xy$r, xy$c)], xi[cbind(xy$r, xy$c)])
    dk <- ifelse(xy$r < xy$c, v > 0.55 * top_u, abs(v) > 0.55 * top_l)
    ok <- is.finite(v)
    if (any(ok)) {
      r <- round(v[ok], digits)
      r[r == 0] <- 0  # not "-0.00"
      graphics::text(
        xy$c[ok], D - xy$r[ok] + 1, formatC(r, format = 'f', digits = digits),
        cex = if (D <= 8) 0.7 else 0.55, col = ifelse(dk[ok], 'white', 'black') )
    }
  }


  #--------------------------------------------------
  # The right margin: % clr var bars, then the keys.
  #--------------------------------------------------
  usr <- graphics::par('usr')
  pin <- graphics::par('pin')
  ui  <- diff(usr[1:2]) / pin[[1]]   # user units per inch, across

  cv <- attr(x, 'clr_var', exact = TRUE)
  if (length(cv) == D && isTRUE(tv > 0)) {
    pct <- 100 * as.numeric(cv)[ord] / tv
    pct[!is.finite(pct)] <- 0
    x0  <- usr[[2]] + 0.1 * ui
    w   <- 0.65 * ui
    len <- w * pmax(pct, 0) / max(pct, 1e-300)
    yc  <- D - seq_len(D) + 1
    hh  <- min(0.4, 0.4 * 15 / max(D, 15))
    graphics::rect(x0, yc - hh, x0 + len, yc + hh, col = 'grey60', border = NA, xpd = NA)
    graphics::segments(x0, usr[[3]], x0, usr[[4]], col = 'grey30', xpd = NA)
    if (D <= 40)
      graphics::text(
        x0 + len + 0.03 * ui, yc, formatC(pct, format = 'f', digits = 1),
        adj = c(0, 0.5), cex = 0.55, xpd = NA )
    graphics::text(x0, usr[[4]], '% clr var', adj = c(0, -0.6), cex = 0.75, xpd = NA)
  }

  if (key) {
    H  <- diff(usr[3:4])
    kx <- usr[[2]] + c(1.2, 1.35) * ui
    dist_plot_key(pal_u, brk_u, 'var(log ratio)', x = kx, y = usr[[3]] + c(0.55, 1) * H)
    dist_plot_key(pal_l, brk_l, 'mean log ratio', x = kx, y = usr[[3]] + c(0, 0.4) * H)
  }

  invisible(x)
}



#' Zero Patterns
#'
#' Which combinations of zero features occur among the samples, how often
#' each occurs, and how often each feature is zero: a look at the zeros
#' before choosing a pseudocount or a [replace_zeros()] method.
#'
#' @inheritParams vlr
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**.
#'
#' @return An integer vector with one value per sample, its pattern's id,
#'   named by the samples, of class `ecodive_zero_pattern`. Pattern 1 is
#'   the most frequent; patterns as frequent as each other are numbered by
#'   their first sample. Attributes:
#'
#'   * `patterns`: a logical matrix, patterns \eqn{\times} features, `TRUE`
#'     where the feature is zero, with the ids as row names;
#'   * `n_samples`: the number of samples with each pattern;
#'   * `zeros`: each feature's proportion of samples where it is zero,
#'     named;
#'   * `n`, `D`: the numbers of samples and features.
#'
#'   `print()` summarizes it, and `plot()` draws the patterns (see
#'   [plot.ecodive_zero_pattern()]).
#'
#' @details
#' A sample's **pattern** is the set of features that are zero in it. Two
#' samples share a pattern only when exactly the same features are zero in
#' both. Sparse tables are read without being made dense.
#'
#' **Reading it.** Zeros come from two sources:
#'
#' * **Count zeros**, features present but not sampled at the depth
#'   sequenced, are scattered: many patterns, each in few samples, with the
#'   rarest features zero most often. These are what [replace_zeros()] and
#'   pseudocounts are for.
#' * **Structural zeros**, features truly absent from some samples, come in
#'   blocks: a few frequent patterns where whole groups of features are
#'   zero together, often in the samples of one group. Replacing them
#'   invents abundances; test them as absent with [structural_zeros()]
#'   instead, or analyse the groups separately.
#'
#' Microbiome tables are sparse enough that nearly every sample has a
#' pattern of its own. The patterns are then little use, but the features'
#' proportions of zeros (the top bars of the plot) still show which
#' features are rare, and filtering those out first leaves fewer zeros to
#' replace.
#'
#' The `patterns` matrix has one row per pattern, so when every sample has
#' its own it is as large as the table, as a dense logical matrix.
#'
#' zCompositions' `zPatterns()` (Palarea-Albaladejo & Martín-Fernández
#' 2015) finds the same patterns and counts; its numbering of them isn't
#' documented, and may differ.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015). zCompositions
#' — R package for multivariate imputation of left-censored data under a
#' compositional approach. *Chemometrics and Intelligent Laboratory
#' Systems*, 143, 85–96. \doi{10.1016/j.chemolab.2015.02.019}
#'
#' @seealso [replace_zeros()], [structural_zeros()], [variation_array()].
#'
#' @export
#' @examples
#'     zp <- zero_pattern(ex_counts)
#'     zp
#'
#'     # The patterns, TRUE where a feature is zero
#'     attr(zp, 'patterns')[1:3, 1:5]
#'
#'     plot(zp)
#'
zero_pattern <- function (counts, margin = 1L) {

  validate_counts()
  validate_margin()

  vals <- stored_values(counts)
  if (!is.numeric(vals) && !is.logical(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  remove('vals')

  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  n <- as.integer(dims[[margin]])
  D <- as.integer(dims[[3L - margin]])
  if (n < 1L || D < 1L)
    stop('`counts` must have at least 1 sample and 1 feature.', call. = FALSE)
  snames <- if (length(dn[[margin]])) as.character(dn[[margin]])
  fnames <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])

  zp <- .Call(C_zero_pattern, counts, margin)

  # Renumber: the most frequent first, ties by first sample (which is the
  # C code's numbering).
  P     <- nrow(zp$patterns)
  size  <- tabulate(zp$group, P)
  perm  <- order(-size, seq_len(P))
  newid <- integer(P)
  newid[perm] <- seq_len(P)

  pat <- zp$patterns[perm, , drop = FALSE]
  dimnames(pat) <- list(as.character(seq_len(P)), fnames)

  res <- newid[zp$group]
  names(res) <- snames
  attr(res, 'patterns')  <- pat
  attr(res, 'n_samples') <- size[perm]
  attr(res, 'zeros')     <- stats::setNames((n - zp$n_nonzero) / n, fnames)
  attr(res, 'n')         <- n
  attr(res, 'D')         <- D
  class(res) <- 'ecodive_zero_pattern'
  res
}



#' @export
print.ecodive_zero_pattern <- function (x, ...) {

  pat   <- attr(x, 'patterns', exact = TRUE)
  ns    <- attr(x, 'n_samples', exact = TRUE)
  zeros <- attr(x, 'zeros', exact = TRUE)
  n     <- attr(x, 'n', exact = TRUE)
  D     <- attr(x, 'D', exact = TRUE)
  fn    <- colnames(pat)
  if (is.null(fn)) fn <- as.character(seq_len(D))

  plural <- function (k, what) sprintf('%d %s%s', k, what, if (k == 1) '' else 's')

  cat(sprintf(
    '# %s among %s and %s\n',
    plural(nrow(pat), 'zero pattern'), plural(n, 'sample'), plural(D, 'feature') ))
  cat(sprintf(
    '# %s never zero, %s always zero\n',
    plural(sum(zeros == 0), 'feature'), sum(zeros == 1) ))

  if (nrow(pat)) {
    z <- fn[pat[1L, ]]
    zs <- if (length(z) == 0) 'none'
          else if (length(z) <= 5) paste0("'", z, "'", collapse = ', ')
          else paste0(paste0("'", z[1:5], "'", collapse = ', '), sprintf(', and %d more', length(z) - 5))
    cat(sprintf(
      '# most frequent pattern (1): %s, zeros in %s\n', plural(ns[[1]], 'sample'), zs ))
  }

  invisible(x)
}



#' Plot zero patterns
#'
#' Draws a [zero_pattern()] result: the most frequent patterns as rows, the
#' features that are sometimes zero as columns, zeros dark; with each
#' feature's percentage of zeros above, and each pattern's percentage of
#' the samples on the right.
#'
#' @param x  A result of [zero_pattern()].
#'
#' @param n_patterns  The number of patterns to draw, the most frequent
#'        first. Default `30`. `Inf` or `TRUE` for all.
#'
#' @param n_features  The largest number of features to draw. Default `50`.
#'        `Inf` or `TRUE` for all.
#'
#' @param ...  Passed on to [graphics::image()]; `main` and `sub` replace
#'        the title and subtitle.
#'
#' @details
#' **Features drawn.** Only features that are zero in some samples but not
#' all: those never zero and those always zero say nothing about patterns,
#' and are counted in the subtitle instead. When more than `n_features`
#' remain, those whose proportion of zeros is closest to one half are
#' drawn, as the most informative. Columns are sorted by their proportion
#' of zeros, the fewest on the left.
#'
#' **Patterns drawn.** The `n_patterns` most frequent, top to bottom.
#' Patterns are defined on all the features, so two rows can look the same
#' when they differ only in features that aren't drawn; the subtitle says
#' how many are hidden.
#'
#' When nearly every sample has its own pattern, as in most microbiome
#' tables, the top bars (each feature's zeros) are the useful part.
#'
#' @return `x`, invisibly.
#'
#' @seealso [zero_pattern()], [replace_zeros()], [structural_zeros()].
#'
#' @export
#' @examples
#'     plot(zero_pattern(ex_counts), n_features = 20)
#'
plot.ecodive_zero_pattern <- function (x, n_patterns = 30, n_features = 50, ...) {

  args <- list(...)
  n_patterns <- plot_count(n_patterns, 'n_patterns')
  n_features <- plot_count(n_features, 'n_features')
  if (n_patterns < 1) stop('`n_patterns` must be at least 1.', call. = FALSE)
  if (n_features < 1) stop('`n_features` must be at least 1.', call. = FALSE)

  pat   <- attr(x, 'patterns', exact = TRUE)
  ns    <- attr(x, 'n_samples', exact = TRUE)
  zeros <- as.numeric(attr(x, 'zeros', exact = TRUE))
  n     <- attr(x, 'n', exact = TRUE)
  D     <- attr(x, 'D', exact = TRUE)
  if (!is.matrix(pat) || length(zeros) != D || length(ns) != nrow(pat))
    stop('`x` must be a result of zero_pattern().', call. = FALSE)
  fn <- colnames(pat)
  if (is.null(fn)) fn <- as.character(seq_len(D))
  P  <- nrow(pat)


  #--------------------------------------------------
  # Which features and patterns are drawn.
  #--------------------------------------------------
  never  <- sum(zeros == 0)
  always <- sum(zeros == 1)
  vary   <- which(zeros > 0 & zeros < 1)
  if (length(vary) > n_features)
    vary <- vary[order(abs(zeros[vary] - 0.5), vary)][seq_len(n_features)]
  feats  <- vary[order(zeros[vary], vary)]
  beyond <- sum(zeros > 0 & zeros < 1) - length(feats)
  k      <- length(feats)
  pp     <- seq_len(min(P, n_patterns))
  np     <- length(pp)

  hidden <- D - k
  parts  <- c(
    if (never)  paste(never, 'never zero'),
    if (always) paste(always, 'always zero'),
    if (beyond) paste(beyond, 'beyond `n_features`') )
  sub <- sprintf(
    '%d of %d patterns shown (%s%% of samples); %d feature%s hidden%s',
    np, P, plot_num(100 * sum(ns[pp]) / n), hidden, if (hidden == 1) '' else 's',
    if (hidden && length(parts)) paste0(' (', paste(parts, collapse = ', '), ')') else '' )
  main <- if (!is.null(args$main)) args$main else 'Zero patterns'
  if (!is.null(args$sub)) sub <- args$sub
  args$main <- NULL
  args$sub  <- NULL


  #--------------------------------------------------
  # Layout: the main panel, bars above and to the right; title and
  # subtitle in the outer margins.
  #--------------------------------------------------
  opar <- graphics::par(c('mfrow', 'mai', 'oma', 'cex', 'mgp'))
  on.exit(graphics::par(opar))

  graphics::layout(matrix(c(2, 0, 1, 3), 2, byrow = TRUE), widths = c(4, 1), heights = c(1, 4))
  graphics::par(cex = 0.9 * opar$cex, mgp = c(2, 0.5, 0))
  main <- strsplit(plot_fit_outer(main, 1.1, 2), '\n', fixed = TRUE)[[1]]
  sub  <- strsplit(plot_fit_outer(sub, 0.75, 1), '\n', fixed = TRUE)[[1]]
  graphics::par(oma = c(0.5 + 0.9 * length(sub), 0, 0.8 + 1.4 * length(main), 0))

  flab <- fn[feats]
  plab <- as.character(pp)
  fcex <- if (k) min(0.7, 0.85 * 30 / max(k, 30) * 1.4) else 0.7
  fcex <- max(fcex, 0.35)
  bot  <- if (k) max(graphics::strwidth(flab, 'inch', cex = fcex), 0) + 0.2 else 0.3
  left <- max(graphics::strwidth(plab, 'inch', cex = 0.7), 0) + 0.45
  ml   <- c(bot, left, 0.05, 0.05)
  graphics::par(mai = ml)

  rows <- rev(pp)   # most frequent at the top
  if (k) {
    z <- t(pat[rows, feats, drop = FALSE] * 1)
    plot_do(
      graphics::image, c(list(seq_len(k), seq_len(np), z), args),
      list(
        col = c('grey95', 'grey20'), breaks = c(-0.5, 0.5, 1.5),
        axes = FALSE, xlab = '', ylab = '' ))
    if (k <= 80 && np <= 80)
      graphics::abline(
        v = seq_len(k - 1) + 0.5, h = seq_len(np - 1) + 0.5, col = 'white', lwd = 0.5 )
    graphics::axis(1, at = seq_len(k), labels = flab, las = 2, tick = FALSE, cex.axis = fcex)
  }
  else {
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, 1.5), ylim = c(0.5, np + 0.5), xaxs = 'i', yaxs = 'i')
    graphics::text(1, (np + 1) / 2, 'No feature is zero in some samples but not in others.', cex = 0.8)
  }
  graphics::axis(2, at = seq_len(np), labels = rev(plab), las = 1, tick = FALSE, cex.axis = if (np <= 40) 0.7 else 0.5)
  graphics::mtext('pattern', side = 2, line = left / graphics::par('csi') - 1, cex = 0.7)
  graphics::box()

  # Above: % zeros of each feature drawn.
  graphics::par(mai = c(0.05, left, 0.1, 0.05))
  graphics::plot.new()
  ztop <- if (k) max(pretty(c(0, 100 * max(zeros[feats])), 3)) else 100
  graphics::plot.window(xlim = c(0.5, max(k, 1) + 0.5), ylim = c(0, ztop), xaxs = 'i', yaxs = 'i')
  if (k)
    graphics::rect(seq_len(k) - 0.4, 0, seq_len(k) + 0.4, 100 * zeros[feats], col = 'grey45', border = NA)
  graphics::axis(2, at = pretty(c(0, ztop), 3), las = 1, cex.axis = 0.6)
  graphics::mtext('% zero', side = 2, line = 2, cex = 0.7)

  # Right: % of samples with each pattern drawn.
  pct <- 100 * ns[rows] / n
  graphics::par(mai = c(bot, 0.05, 0.05, 0.25))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, max(pct, 1e-300) * 1.05), ylim = c(0.5, np + 0.5), xaxs = 'i', yaxs = 'i')
  graphics::rect(0, seq_len(np) - 0.4, pct, seq_len(np) + 0.4, col = 'grey45', border = NA)
  graphics::axis(1, cex.axis = 0.6)
  graphics::mtext('% samples', side = 1, line = 1.6, cex = 0.7)

  for (i in seq_along(main))
    graphics::mtext(main[[i]], side = 3, outer = TRUE, line = 0.4 + 1.4 * (length(main) - i), font = 2, cex = 1.1)
  for (i in seq_along(sub))
    graphics::mtext(sub[[i]], side = 1, outer = TRUE, line = 0.1 + 0.9 * (i - 1), cex = 0.75)

  invisible(x)
}


# A title for the outer margin, wrapped to the device's width at text size
# `cex` and `font`.
plot_fit_outer <- function (main, cex, font) {
  if (!is.character(main) || length(main) != 1 || is.na(main)) return (main)
  w <- graphics::strwidth(main, 'inches', cex = cex / graphics::par('cex'), font = font)  # mtext()'s cex is absolute
  avail <- 0.95 * graphics::par('din')[[1]]
  if (!isTRUE(w > avail)) return (main)
  n <- ceiling(w / avail)
  paste(strwrap(main, width = ceiling(nchar(main) / n) + 1L), collapse = '\n')
}
