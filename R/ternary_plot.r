# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Ternary Plot
#'
#' Draws three-part compositions in a ternary diagram: an equilateral
#' triangle whose vertices are the pure parts, where each point's distance
#' from an edge is its share of the part at the opposite vertex. It has two
#' modes:
#'
#' * **Samples mode** (`parts` given): one point per sample, at its
#'   three-part subcomposition, such as three phyla, or Bacteroidota,
#'   Bacillota and everything else. The points can be **centred**
#'   (`center = TRUE`) and given **confidence or prediction regions**
#'   (`ellipses`), per group.
#' * **Features mode** (`parts = NULL`): one point per feature, at its mean
#'   relative abundance in each of three sample groups (`groups`), sized by
#'   its overall abundance, to show where each feature is enriched.
#'
#' The plot is drawn at once, and the drawn coordinates are returned
#' invisibly, as [graphics::barplot()] returns its bars' positions.
#'
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features), or
#'        proportions; see **Input Types**.
#'
#' @param parts  Samples mode's three parts, in vertex order (bottom left,
#'        bottom right, top):
#'   * a vector selecting exactly 3 features, by name, 1-based index, or a
#'     logical vector with one value per feature; or
#'   * a list of 3 such selectors, whose features are summed into each
#'     part (amalgamation). One element may be `NA`, for every feature in
#'     neither of the other two.
#'
#'   The parts are labelled with the list's names, else the features'
#'   names (`"Other"` for an `NA` part, `"Part 2"` for an unnamed part of
#'   several features). No feature may be in two parts. Default `NULL`:
#'   features mode.
#'
#' @param groups  The samples' groups: the name of a column of `data`, or a
#'        vector with one value per sample (matched to the samples by name
#'        when both are named, as for [permanova()]). In samples mode they
#'        colour the points, with a legend, and give one ellipse each; a
#'        numeric vector gives a colour gradient, and no ellipses. In
#'        features mode, required, with exactly 3 values or levels (`NA`
#'        is in none): the vertices, in the order of a factor's levels, else
#'        sorted.
#'
#' @param data  A data frame of sample metadata, with one row per sample
#'        (matched by row name, else in order), when `groups` names one of
#'        its columns.
#'
#' @param pseudocount  Samples mode, with `center = TRUE` or `ellipses`:
#'        added to each of the three part sums before taking logs. Default
#'        `NULL`, by the rule every log-ratio method in ecodive shares: none
#'        when `counts` has no zeros, else half its smallest non-zero value
#'        (with a warning). It is added to the points too, so that they
#'        match the ellipses. Not used otherwise: zeros are then drawn on
#'        the edges.
#'
#' @param center  Samples mode: whether to centre the points, perturbing
#'        each by the inverse of their geometric mean so that it moves to
#'        the middle of the triangle (see **Centring**). The gridlines are
#'        perturbed alike, and still read in the original proportions.
#'        Default `FALSE`.
#'
#' @param ellipses  Samples mode: `"none"` (the default), `"confidence"`
#'        for a confidence region of each group's mean, or `"prediction"`
#'        for a region that holds a new sample of the group, at
#'        `conf_level` (see **Ellipses**). Without `groups`, one for all
#'        the samples.
#'
#' @param conf_level  The confidence level of `ellipses`, strictly between
#'        0 and 1. Default `0.95`.
#'
#' @param grid  Whether to draw gridlines at proportions 0.2, 0.4, 0.6 and
#'        0.8 of each part, labelled on the edge they meet. Default `TRUE`.
#'
#' @param labels  Features mode: how many features to name, the most
#'        abundant first: a number, `Inf` (or `TRUE`) for all, or `0` (or
#'        `FALSE`) for none. Default `10`. Not used in samples mode.
#'
#' @param legend  Where to put the legend: a keyword for
#'        [graphics::legend()], such as `"topright"` (the default) or
#'        `"topleft"`, or `NULL` for none.
#'
#' @param col  The points' colours. Samples mode: as for
#'        [plot.ecodive_ord()]; without `groups`, one colour or one per
#'        sample; with them, one per group (in order, or named by group),
#'        or the colours of a numeric gradient. Features mode: one colour
#'        for all, one per feature of `counts`, or colours named by feature
#'        to replace those features' defaults (to mark the results of
#'        [da_test()], say). With `col`, features mode draws no legend.
#'
#' @param pch  The points' symbols, as for `plot()`. Default `1` (open
#'        circles), or `16` (filled circles) with `groups` and in features
#'        mode.
#'
#' @param ...  Passed on to [graphics::plot.default()], which sets up the
#'        plot, such as `main` or `sub`; `cex`, `lwd` and `bg` are passed to
#'        [graphics::points()] instead. They replace the defaults.
#'
#' @details
#' **Coordinates.** Each point is a composition \eqn{(p_1, p_2, p_3)},
#' closed to sum 1, drawn at
#' \deqn{x = p_2 + p_3 / 2, \quad y = \frac{\sqrt 3}{2} p_3,}
#' so that part 1's vertex is at \eqn{(0, 0)}, part 2's at \eqn{(1, 0)} and
#' part 3's at \eqn{(1/2, \sqrt 3 / 2)}. A gridline is a line of constant
#' share of one part: the labels along the bottom edge read part 2 (rising
#' to the right), those along the right edge part 3 (rising upwards), and
#' those along the left edge part 1 (rising downwards).
#'
#' **Samples mode.** Each sample's counts are summed over each part's
#' features, and the three sums closed: \eqn{p_{ik} = y_{ik} / \sum_k
#' y_{ik}}. Summing features (amalgamation) is not a log-ratio operation:
#' the log ratios of a sum are not those of its members. But it is what a
#' plot of phyla means, and the other parts' proportions are unchanged by
#' it. Samples whose three parts are all zero have no point, and are left
#' out with a warning.
#'
#' **Centring.** Microbiome subcompositions crowd against one vertex or
#' edge, where the triangle's geometry hides their differences. Centring
#' (von Eynatten, Pawlowsky-Glahn & Egozcue 2002) perturbs every sample by
#' the inverse of the samples' centre \eqn{g}, the closed geometric mean of
#' each part,
#' \deqn{p'_i = C(p_{i1} / g_1, \; p_{i2} / g_2, \; p_{i3} / g_3),
#'   \quad g_k \propto \exp\left(\frac{1}{n} \sum_i \log p_{ik}\right),}
#' with \eqn{C} the closure. This moves the centre to the middle of the
#' triangle, and spreads the points, without changing any difference
#' between two samples' log ratios: in Aitchison's geometry, perturbation is
#' a translation. The gridlines are perturbed by the same \eqn{g^{-1}}, so
#' they become curves, which keep their ends on the edges and still mark
#' the **original** proportions. The subtitle gives \eqn{g}.
#'
#' **Ellipses.** In the isometric log-ratio (ilr) coordinates
#' \deqn{z_i = V \log p_i, \quad
#'   V = \begin{pmatrix} 1/\sqrt 2 & -1/\sqrt 2 & 0 \\
#'   1/\sqrt 6 & 1/\sqrt 6 & -2/\sqrt 6 \end{pmatrix}}
#' (Egozcue et al. 2003), each group of \eqn{m \ge 3} samples, with mean
#' \eqn{\bar z} and covariance \eqn{S} (divisor \eqn{m - 1}), gets the
#' ellipse \eqn{(z - \bar z)^\top S^{-1} (z - \bar z) = r^2}, with
#' \eqn{F} the `conf_level` quantile of the F distribution on 2 and \eqn{m
#' - 2} degrees of freedom and
#' * `"confidence"`: \eqn{r^2 = \frac{2 (m - 1)}{m (m - 2)} F}, Hotelling's
#'   \eqn{T^2} region for the mean;
#' * `"prediction"`: \eqn{r^2 = \frac{2 (m + 1)(m - 1)}{m (m - 2)} F}, a
#'   region for one new sample.
#'
#' These assume that the log ratios are multivariate normal. The ellipse
#' is drawn back in the triangle, \eqn{p = C(\exp(V^\top z))} at 200
#' points, centred as the points are. Straight lines in ilr coordinates are
#' curves in the triangle, so the regions are closed curves, not ellipses,
#' there. The choice of ilr basis does not change them. A group with fewer
#' than 3 samples, or whose log ratios are collinear, gets none, with one
#' warning naming them.
#'
#' ggtern's `geom_confidence_tern()` draws a different curve: the
#' \eqn{\chi^2_2} contour of the samples, \eqn{r^2 =} `qchisq(level, 2)`
#' whatever \eqn{m} (checked from its output). That is the large-sample
#' limit of the prediction region, smaller for small groups, and not a region
#' for the mean. `compositions::ellipses()` draws the same quadratic form
#' with a radius the caller gives.
#'
#' **Features mode.** Every sample is closed over all the features
#' (relative abundances). For feature \eqn{j}, \eqn{m_{jg}} is its mean
#' relative abundance over the samples of group \eqn{g}, and its point is
#' \eqn{C(m_{j1}, m_{j2}, m_{j3})}: its share of the three groups' means.
#' Features absent from all three groups are left out (the subtitle counts
#' them), as are samples with no counts, and those whose group is `NA`.
#' Each point's size is proportional to the square root of the feature's
#' mean relative abundance over the three groups' samples (`cex` from 0.4
#' to 3). A point takes its group's colour when at least half its share is
#' in that group ("mostly in" it), and is grey otherwise. This is a visual
#' convention, not a test: give `col` to colour by differential abundance
#' instead. The `labels` most abundant features are named.
#'
#' @inheritSection documentation Input Types
#'
#' @section Same names in other packages:
#' `vcd::ternaryplot()`, `Ternary::TernaryPlot()`,
#' `robCompositions::ternaryDiag()`, `compositions::plot.acomp()` and the
#' ggtern package draw ternary diagrams too. None of them is called
#' `ternary_plot()`, so none is masked.
#'
#' @return Invisibly, a data frame with one row per point (sample or
#'   feature): the three closed proportions as drawn (centred, with `center
#'   = TRUE), named after the parts or groups; `x` and `y`, the point's
#'   coordinates; and `group` (samples mode, with `groups`) or `abundance`
#'   (features mode: the mean relative abundance that sizes the point). Its
#'   attributes:
#'
#'   * `mode`: `"samples"` or `"features"`.
#'   * `center`: the centre \eqn{g} that was removed, named by part, or
#'     `NULL`.
#'   * `pseudocount`: the pseudocount added to the part sums (`0` for none).
#'   * `ellipses`: a list with one matrix per drawn region, named by group
#'     (`"all"` without `groups`): 200 points with the three proportions,
#'     `x` and `y` as columns. Empty without ellipses.
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman & Hall.
#'
#' von Eynatten, H., Pawlowsky-Glahn, V., & Egozcue, J. J. (2002).
#' Understanding perturbation on the simplex: a simple method to better
#' visualize and interpret compositional data in ternary diagrams.
#' *Mathematical Geology*, 34, 249-257. \doi{10.1023/A:1014826205533}
#'
#' Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal,
#' C. (2003). Isometric logratio transformations for compositional data
#' analysis. *Mathematical Geology*, 35, 279-300.
#' \doi{10.1023/A:1023818214614}
#'
#' Hamilton, N. E., & Ferry, M. (2018). ggtern: ternary diagrams using
#' ggplot2. *Journal of Statistical Software*, 87(3).
#' \doi{10.18637/jss.v087.c03}
#'
#' @seealso [composition_barplot()], [ilr()], [balance()], [log_ratio()]
#'
#' @export
#' @examples
#'     # Three genera of ex_counts (samples in rows).
#'     ternary_plot(ex_counts, parts = c('Streptococcus', 'Bacteroides', 'Corynebacterium'))
#'
#'     # Two genera against everything else, in 30 samples of two groups,
#'     # centred, with 95% confidence regions of each group's mean.
#'     n <- 30; p <- 12
#'     grp    <- rep(c('A', 'B'), each = n / 2)
#'     counts <- sapply(seq_len(p), function (j)
#'       round(1000 * exp(-j / 3) * ifelse(grp == 'A', 1, j / 4) *
#'             exp(sin(seq_len(n) * j) / 2)))
#'     dimnames(counts) <- list(paste0('S', 1:n), paste0('Genus', 1:p))
#'     res <- ternary_plot(
#'       counts, parts = list(G1 = 'Genus1', G2 = 'Genus2', Rest = NA),
#'       groups = grp, center = TRUE, ellipses = 'confidence' )
#'     head(res)
#'
#'     # Features mode: where each genus is most abundant, among 3 groups.
#'     site <- rep(c('Gut', 'Skin', 'Oral'), each = n / 3)
#'     ternary_plot(counts, groups = site)
ternary_plot <- function (
    counts, parts = NULL, groups = NULL, data = NULL, margin = 1L,
    pseudocount = NULL, center = FALSE, ellipses = 'none',
    conf_level = 0.95, grid = TRUE, labels = 10, legend = 'topright',
    col = NULL, pch = NULL, ... ) {

  args <- list(...)

  validate_counts()
  validate_margin()

  for (arg in c('center', 'grid'))
    if (!isTRUE(get(arg)) && !isFALSE(get(arg)))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)
  if (!is.character(ellipses) || length(ellipses) != 1L ||
      !ellipses %in% c('none', 'confidence', 'prediction'))
    stop('`ellipses` must be "none", "confidence" or "prediction".', call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)
  labels <- plot_count(labels, 'labels')
  pos    <- plot_legend_pos(legend)

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
  n      <- as.integer(dims[[margin]])
  D      <- as.integer(dims[[3L - margin]])
  snames <- if (length(dn[[margin]])) as.character(dn[[margin]])
  fnames <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])

  if (is.null(parts)) {
    if (center)
      stop('`center` applies to samples mode only (with `parts`): ',
           'in features mode, group means are often 0.', call. = FALSE)
    if (ellipses != 'none')
      stop('`ellipses` apply to samples mode only (with `parts`): ',
           'in features mode, group means are often 0.', call. = FALSE)
    res <- tp_features(counts, margin, n, D, snames, fnames, groups, data, col, pch)
  }
  else {
    res <- tp_samples(
      counts, margin, n, D, snames, fnames, parts, groups, data,
      pseudocount, center, ellipses, conf_level, col, pch )
  }

  tp_draw(res, grid, labels, pos, args)

  invisible(res$out)
}



# The vertices' barycentric to Cartesian map: rows of `P` (closed, 3
# columns) to an n x 2 matrix of x and y.
tp_xy <- function (P) {
  P <- matrix(P, ncol = 3L)
  cbind(x = P[, 2] + P[, 3] / 2, y = sqrt(3) / 2 * P[, 3])
}


# Rows of `m` closed to sum 1.
tp_close <- function (m) {
  m <- matrix(m, ncol = 3L)
  m / rowSums(m)
}


# The fixed orthonormal ilr basis: 2 x 3, rows summing to 0.
tp_basis <- function () {
  rbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6))
}



#--------------------------------------------------
# `parts`: a list of 3 index vectors (in the order given) and their
# labels.
#--------------------------------------------------
tp_parts <- function (parts, fnames, D) {

  # One selector, in the order given (feature_set() checks it and sorts).
  pick <- function (x, arg) {
    i <- feature_set(x, fnames, D, arg)
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) i
    else if (is.character(x)) match(x, fnames)
    else as.integer(x)
  }

  if (is.list(parts)) {

    if (length(parts) != 3L)
      stop(sprintf(
        'A list `parts` must have 3 elements, one per part (it has %d).', length(parts)),
        call. = FALSE )

    is_na <- vapply(parts, function (x) length(x) == 1L && is.atomic(x) && is.na(x), TRUE)
    if (sum(is_na) > 1L)
      stop('Only one element of `parts` may be NA (every feature in neither other part).',
           call. = FALSE)

    idx <- vector('list', 3L)
    for (k in which(!is_na)) idx[[k]] <- pick(parts[[k]], sprintf('parts[[%d]]', k))

    all_i <- unlist(idx)
    if (anyDuplicated(all_i)) {
      dup <- unique(all_i[duplicated(all_i)])
      stop(sprintf(
        '`parts` selects %s more than once: %s.',
        if (length(dup) == 1L) 'a feature' else 'features',
        if (is.null(fnames)) paste(utils::head(dup, 5), collapse = ', ') else name_list(fnames[dup]) ),
        call. = FALSE )
    }

    if (any(is_na)) {
      k    <- which(is_na)
      rest <- setdiff(seq_len(D), all_i)
      if (length(rest) == 0L)
        stop(sprintf(
          '`parts[[%d]]` is NA, for every feature in neither other part, but there are none.', k),
          call. = FALSE )
      idx[[k]] <- rest
    }

    lab <- vapply(seq_len(3L), function (k) {
      if (is_na[[k]]) 'Other'
      else if (length(idx[[k]]) == 1L) tp_feature_name(idx[[k]], fnames)
      else paste('Part', k)
    }, '')
    nm <- names(parts)
    if (!is.null(nm)) lab[!is.na(nm) & nzchar(nm)] <- nm[!is.na(nm) & nzchar(nm)]
  }

  else {
    i <- pick(parts, 'parts')
    if (length(i) != 3L)
      stop(sprintf(
        '`parts` must select exactly 3 features, or be a list of 3 sets of them (it selects %d).',
        length(i) ), call. = FALSE )
    idx <- as.list(i)
    lab <- tp_feature_name(i, fnames)
    nm  <- names(parts)
    if (!is.null(nm)) lab[!is.na(nm) & nzchar(nm)] <- nm[!is.na(nm) & nzchar(nm)]
  }

  list(idx = idx, labels = make.unique(as.character(lab)))
}


tp_feature_name <- function (i, fnames) {
  if (is.null(fnames)) paste('Feature', i) else fnames[i]
}



#--------------------------------------------------
# Samples mode: the points, the centre and the ellipses.
#--------------------------------------------------
tp_samples <- function (
    counts, margin, n, D, snames, fnames, parts, groups, data,
    pseudocount, center, ellipses, conf_level, col, pch ) {

  pp   <- tp_parts(parts, fnames, D)
  part <- rep(-1L, D)
  for (k in 1:3) part[pp$idx[[k]]] <- k - 1L

  g <- plot_groups(n, snames, groups, data, ref = '`counts`')

  Y    <- .Call(C_ternary_sums, counts, margin, part)
  keep <- which(rowSums(Y) > 0)
  if (length(keep) == 0L)
    stop('Every sample has 0 in all three parts, so there is nothing to draw.', call. = FALSE)
  if (length(keep) < n) {
    drop <- setdiff(seq_len(n), keep)
    warning(
      sprintf('Left out %d sample%s whose three parts are all zero: %s',
        length(drop), if (length(drop) == 1L) '' else 's',
        name_list(if (is.null(snames)) drop else snames[drop]) ),
      call. = FALSE )
  }
  Y <- Y[keep, , drop = FALSE]
  if (!is.null(g$f))      g$f      <- g$f[keep]
  if (!is.null(g$values)) g$values <- g$values[keep]
  if (!is.null(col) && is.null(g$f) && is.null(g$values) && length(col) == n)
    col <- col[keep]
  if (!is.null(pch) && is.null(g$f) && length(pch) == n)
    pch <- pch[keep]

  need_log <- center || ellipses != 'none'
  pc <- 0
  if (need_log || !is.null(pseudocount)) {
    pc1 <- log_pseudocount(
      counts, pseudocount, 'the log ratios of a ternary plot', constant_ok = TRUE )
    if (need_log) pc <- pc1
  }

  P <- tp_close(Y + pc)
  m <- nrow(P)

  # The centre, and the points' ilr coordinates.
  g_cen <- NULL
  if (need_log) {
    L <- log(P)
    if (center) {
      g_cen <- exp(colMeans(L))
      g_cen <- g_cen / sum(g_cen)
      names(g_cen) <- pp$labels
    }
  }
  Pd <- if (center) tp_close(sweep(P, 2L, g_cen, '/')) else P

  # The ellipses, one per group (or one for all).
  ell  <- list()
  ell_col <- character(0)
  look <- ord_plot_look(g, col, pch)

  if (ellipses != 'none') {
    Z <- L %*% t(tp_basis())
    if (!is.null(g$values)) {
      warning(
        'Numeric `groups` are drawn as a gradient, so no ellipses are drawn. ',
        'Give factor(groups) for one per value.', call. = FALSE )
    }
    else {
      if (is.null(g$f)) {
        sets <- list(all = seq_len(m))
        ecol <- if (is.null(col) || length(col) != 1L) graphics::par('fg') else col
        ecol <- list(all = ecol)
      }
      else {
        lev  <- levels(g$f)
        sets <- lapply(seq_along(lev), function (k) which(as.integer(g$f) == k))
        names(sets) <- ifelse(is.na(lev), 'NA', lev)
        used <- lengths(sets) > 0
        sets <- sets[used]
        ecol <- as.list(look$pal[used])
        names(ecol) <- names(sets)
      }

      none <- character(0)
      for (nm in names(sets)) {
        e <- tp_ellipse(Z[sets[[nm]], , drop = FALSE], ellipses, conf_level)
        if (is.null(e)) { none <- c(none, nm); next }
        Pe <- tp_close(exp(e %*% tp_basis()))
        if (center) Pe <- tp_close(sweep(Pe, 2L, g_cen, '/'))
        colnames(Pe) <- pp$labels
        ell[[nm]]     <- cbind(Pe, tp_xy(Pe))
        ell_col[[nm]] <- ecol[[nm]]
      }
      if (length(none))
        warning(
          sprintf('No ellipse for %s with fewer than 3 samples or collinear log ratios: %s.',
            if (length(none) == 1L) '1 group' else paste(length(none), 'groups'),
            name_list(none) ),
          call. = FALSE )
    }
  }

  xy  <- tp_xy(Pd)
  out <- data.frame(Pd, xy, check.names = FALSE)
  names(out)[1:3] <- pp$labels
  rownames(out) <- if (is.null(snames)) as.character(keep) else snames[keep]
  if (!is.null(groups)) out$group <- if (!is.null(g$values)) g$values else g$f
  attr(out, 'mode')        <- 'samples'
  attr(out, 'center')      <- g_cen
  attr(out, 'pseudocount') <- pc
  attr(out, 'ellipses')    <- ell

  # Subtitle.
  sub <- character(0)
  if (center)
    sub <- c(sub, sprintf(
      'Centred: perturbed by the inverse of the geometric mean (%s)',
      paste(vapply(g_cen, plot_num, ''), collapse = ', ') ))
  if (need_log && pc > 0)
    sub <- c(sub, sprintf('pseudocount %s', plot_num(pc)))
  if (!need_log && any(P == 0))
    sub <- c(sub, 'zeros drawn on the edges')
  if (length(keep) < n)
    sub <- c(sub, sprintf('%d sample%s with no counts in the parts left out',
                          n - length(keep), if (n - length(keep) == 1L) '' else 's'))
  if (ellipses != 'none')
    sub <- c(sub, sprintf('%s %s regions', plot_pct(conf_level), ellipses))

  list(
    out = out, xy = xy, vertex = pp$labels, center = g_cen,
    main = paste('Ternary plot of', paste(pp$labels, collapse = ', ')),
    sub = sub,
    col = look$col, pch = look$pch, cex = 1,
    ell = ell, ell_col = ell_col, label_idx = integer(0),
    legend = if (!is.null(g$f) || !is.null(g$values))
      list(labels = look$legend, col = look$legend_col, pch = look$legend_pch, title = g$title) )
}


# The ellipse of ilr points `Z` (m x 2), as 200 points (200 x 2), or NULL
# for fewer than 3 points or a singular covariance.
tp_ellipse <- function (Z, kind, conf_level) {

  m <- nrow(Z)
  if (m < 3L) return (NULL)
  zbar <- colMeans(Z)
  S    <- stats::cov(Z)
  ev   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(ev)) || !(ev[[1]] > 0) || ev[[2]] <= ev[[1]] * 1e-10) return (NULL)

  r2 <- tp_radius2(m, kind, conf_level)
  R  <- chol(S)
  th <- seq(0, 2 * pi, length.out = 201L)[-201L]
  U  <- cbind(cos(th), sin(th))
  sweep(sqrt(r2) * U %*% R, 2L, zbar, '+')
}


# The squared Mahalanobis radius of a region for m samples.
tp_radius2 <- function (m, kind, conf_level) {
  f <- stats::qf(conf_level, 2, m - 2)
  if (kind == 'confidence') 2 * (m - 1) / (m * (m - 2)) * f
  else                      2 * (m + 1) * (m - 1) / (m * (m - 2)) * f
}



#--------------------------------------------------
# Features mode: each feature's share of the three groups' means.
#--------------------------------------------------
tp_features <- function (counts, margin, n, D, snames, fnames, groups, data, col, pch) {

  if (is.null(groups))
    stop('Give `parts` (3 features, or 3 sets of them) for samples mode, ',
         'or `groups` with exactly 3 levels for features mode.', call. = FALSE)

  g <- plot_groups(n, snames, groups, data, ref = '`counts`')
  v <- if (!is.null(g$values)) g$values else g$f
  f <- if (is.factor(v)) v else factor(v)
  f <- factor(f, levels = levels(f)[!is.na(levels(f))])
  k <- nlevels(f)
  if (k != 3L)
    stop(sprintf(
      'Without `parts`, `groups` must have exactly 3 levels (it has %d): ', k),
      'each feature is placed by its mean abundance in the 3 groups.', call. = FALSE)
  lev <- levels(f)

  gi  <- ifelse(is.na(f), -1L, as.integer(f) - 1L)
  res <- .Call(C_ternary_means, counts, margin, as.integer(gi))

  empty <- which(res$sizes == 0L)
  if (length(empty))
    stop(sprintf(
      '%s %s no samples with counts above 0: %s.',
      if (length(empty) == 1L) 'Group' else 'Groups',
      if (length(empty) == 1L) 'has' else 'have', name_list(lev[empty]) ),
      call. = FALSE )
  zero <- which(gi >= 0L & !(res$totals > 0))
  if (length(zero))
    warning(
      sprintf('Left out %d sample%s whose counts are all zero: %s',
        length(zero), if (length(zero) == 1L) '' else 's',
        name_list(if (is.null(snames)) zero else snames[zero]) ),
      call. = FALSE )

  M    <- res$means
  tot  <- rowSums(M)
  keep <- which(tot > 0)
  if (length(keep) == 0L)
    stop('No feature has counts in the 3 groups, so there is nothing to draw.', call. = FALSE)
  P   <- M[keep, , drop = FALSE] / tot[keep]
  ab  <- res$overall[keep]
  nms <- tp_feature_name(keep, fnames)

  # Colours: the group's when it holds at least half the share, else muted.
  pal  <- plot_palette(3L)
  top  <- max.col(P, ties.method = 'first')
  most <- P[cbind(seq_along(top), top)] >= 0.5
  pcol <- ifelse(most, pal[top], PLOT_MUTED)
  legend <- list(
    labels = c(lev, 'No majority'), col = c(pal, PLOT_MUTED), pch = 16, title = 'Mostly in' )

  if (!is.null(col)) {
    if (!is.character(col) && !is.numeric(col))
      stop('`col` must be a vector of colours.', call. = FALSE)
    cn <- names(col)
    if (!is.null(cn) && all(nzchar(cn))) {
      all_names <- tp_feature_name(seq_len(D), fnames)
      bad <- setdiff(cn, all_names)
      if (length(bad))
        stop('`col` names features not in `counts`: ', name_list(bad), '.', call. = FALSE)
      hit <- match(cn, nms)
      pcol[hit[!is.na(hit)]] <- col[!is.na(hit)]
    }
    else if (length(col) == 1L) pcol[] <- col
    else if (length(col) == D)  pcol <- col[keep]
    else
      stop(sprintf(
        'In features mode, `col` must be one colour, one per feature (%d), or named by feature.', D),
        call. = FALSE )
    legend <- NULL
  }

  if (is.null(pch)) pch <- 16
  else if (length(pch) == D) pch <- pch[keep]

  # Sizes: cex 0.4 to 3, in proportion to the square root of the abundance.
  cex <- 0.4 + 2.6 * sqrt(ab / max(ab))

  colnames(P) <- lev
  xy  <- tp_xy(P)
  out <- data.frame(P, xy, abundance = ab, check.names = FALSE)
  rownames(out) <- make.unique(nms)
  attr(out, 'mode')        <- 'features'
  attr(out, 'center')      <- NULL
  attr(out, 'pseudocount') <- 0
  attr(out, 'ellipses')    <- list()

  sub <- c(
    sprintf('%d feature%s', length(keep), if (length(keep) == 1L) '' else 's'),
    'size: mean relative abundance' )
  if (length(keep) < D)
    sub <- c(sub, sprintf('%d absent from all 3 groups left out', D - length(keep)))

  list(
    out = out, xy = xy, vertex = lev, center = NULL,
    main = paste('Mean relative abundance in', paste(lev, collapse = ', ')),
    sub = sub, col = pcol, pch = pch, cex = cex,
    ell = list(), ell_col = character(0),
    label_idx = order(-ab, seq_along(ab)),
    legend = legend )
}



#--------------------------------------------------
# Drawing.
#--------------------------------------------------
tp_draw <- function (res, grid, labels, pos, args) {

  h   <- sqrt(3) / 2
  pt_args <- args[intersect(names(args), c('cex', 'lwd', 'bg'))]
  fr_args <- plot_drop(args, c('cex', 'lwd', 'bg', 'sub'))
  sub     <- if ('sub' %in% names(args)) args[['sub']] else res$sub

  plot_do(
    graphics::plot.default, c(list(res$xy[, 1], res$xy[, 2]), fr_args),
    list(
      type = 'n', axes = FALSE, xlab = '', ylab = '', asp = 1,
      xlim = c(-0.08, 1.08), ylim = c(-0.08, h + 0.06), main = res$main ))
  tp_subtitle(sub)

  if (grid) tp_grid(res$center, res$vertex)

  graphics::polygon(c(0, 1, 0.5), c(0, 0, h), border = graphics::par('fg'))

  # Vertex labels, outside the corners.
  vx <- c(0, 1, 0.5); vy <- c(0, 0, h)
  plot_label(vx[1:2], vy[1:2], res$vertex[1:2], col = graphics::par('fg'), pos = 1, font = 2, cex = 0.9)
  plot_label(vx[[3]], vy[[3]], res$vertex[[3]], col = graphics::par('fg'), pos = 3, font = 2, cex = 0.9)

  # Points: features largest first, so that small ones stay visible.
  ord <- if (length(res$label_idx)) res$label_idx else seq_len(nrow(res$xy))
  rep_n <- function (v) if (length(v) == nrow(res$xy)) v[ord] else v
  plot_do(
    graphics::points, c(list(res$xy[ord, 1], res$xy[ord, 2]), pt_args),
    list(col = rep_n(res$col), pch = rep_n(res$pch), cex = rep_n(res$cex)) )

  for (nm in names(res$ell))
    graphics::polygon(res$ell[[nm]][, 'x'], res$ell[[nm]][, 'y'],
                      border = res$ell_col[[nm]], col = NA, lwd = 1.5)

  # Feature names.
  if (labels > 0 && length(res$label_idx)) {
    i   <- utils::head(res$label_idx, labels)
    col <- if (length(res$col) == nrow(res$xy)) res$col[i] else res$col
    plot_label(res$xy[i, 1], res$xy[i, 2], rownames(res$out)[i], col = col, pos = 3, cex = 0.7)
  }

  if (!is.null(res$legend))
    plot_legend(pos, res$legend$labels, res$legend$col, res$legend$pch, res$legend$title)

  invisible()
}


# The subtitle's pieces, joined by "; " onto as few lines as fit the
# figure's width, under the plot. A user's `sub` is drawn as given.
tp_subtitle <- function (sub) {
  if (is.null(sub) || length(sub) == 0L || identical(sub, '')) return (invisible())
  if (is.character(sub) && length(sub) > 1L) {
    sub[[1]] <- plot_cap(sub[[1]])
    cex   <- 0.85
    avail <- 0.95 * graphics::par('fin')[[1]]
    lines <- sub[[1]]
    for (piece in sub[-1L]) {
      try_line <- paste(lines[[length(lines)]], piece, sep = '; ')
      if (graphics::strwidth(try_line, 'inches', cex = cex) <= avail) lines[[length(lines)]] <- try_line
      else lines <- c(lines, piece)
    }
    lines[-length(lines)] <- paste0(lines[-length(lines)], ';')
    graphics::mtext(lines, side = 1, line = 1 + 1.1 * (seq_along(lines) - 1L), cex = cex)
  }
  else {
    if (is.character(sub)) sub <- plot_cap(sub)
    graphics::mtext(sub, side = 1, line = 1, cex = 0.85)
  }
  invisible()
}


# Gridlines at 0.2, 0.4, 0.6 and 0.8 of each part, perturbed by the
# inverse of `center` (or not, for NULL), with their labels on the edge
# where the next part (cyclically) is 0. Centred gridlines can crowd
# together, so a label that would overlap one already drawn, or a vertex's
# label (`vertex`), is left out; its gridline and tick are still drawn.
tp_grid <- function (center, vertex) {

  h    <- sqrt(3) / 2
  s    <- seq(0, 1, length.out = 101L)
  out  <- list(c(-h, 0.5), c(0, -1), c(h, 0.5))   # outward normals: left, bottom, right edges
  pert <- function (P) if (is.null(center)) P else tp_close(sweep(P, 2L, center, '/'))
  cex  <- 0.65

  # Boxes (x0, x1, y0, y1) of the drawn labels: the vertices' first.
  box <- function (x, y, w, hh) c(x - w / 2, x + w / 2, y - hh / 2, y + hh / 2)
  vh  <- graphics::strheight('M', cex = 0.9, font = 2)
  vw  <- graphics::strwidth(vertex, cex = 0.9, font = 2)
  boxes <- list(
    box(0,   -vh * 1.5,      vw[[1]], vh * 1.6),
    box(1,   -vh * 1.5,      vw[[2]], vh * 1.6),
    box(0.5, h + vh * 1.5,   vw[[3]], vh * 1.6) )
  overlaps <- function (b)
    any(vapply(boxes, function (a) a[[1]] < b[[2]] && b[[1]] < a[[2]] &&
                                     a[[3]] < b[[4]] && b[[3]] < a[[4]], TRUE))

  for (k in 1:3) {
    k1 <- k %% 3L + 1L
    k2 <- k1 %% 3L + 1L
    nv <- out[[k]]
    for (t in c(0.2, 0.4, 0.6, 0.8)) {
      P <- matrix(0, length(s), 3L)
      P[, k]  <- t
      P[, k1] <- (1 - t) * (1 - s)
      P[, k2] <- (1 - t) * s
      xy <- tp_xy(pert(P))
      graphics::lines(xy[, 1], xy[, 2], col = 'grey85')
      # The end where part k1 is 0.
      e <- xy[length(s), ]
      graphics::segments(e[[1]], e[[2]], e[[1]] + 0.012 * nv[[1]], e[[2]] + 0.012 * nv[[2]], col = 'grey50')
      lab <- format(t)
      lx  <- e[[1]] + 0.04 * nv[[1]]
      ly  <- e[[2]] + 0.04 * nv[[2]]
      b   <- box(lx, ly, 1.2 * graphics::strwidth(lab, cex = cex), 1.4 * graphics::strheight(lab, cex = cex))
      if (overlaps(b)) next
      boxes[[length(boxes) + 1L]] <- b
      graphics::text(lx, ly, lab, col = 'grey40', cex = cex, xpd = NA)
    }
  }
  invisible()
}
