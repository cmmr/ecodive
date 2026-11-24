# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Composition Barplot
#'
#' Stacked bars of relative abundance, one per sample, drawn so that they
#' can be read: the samples ordered by their compositions (neighbouring bars
#' are similar), the most abundant features on the common baseline, and,
#' with `feature_groups`, one hue per higher-level group (such as a phylum)
#' with shades for its most abundant members (microshades' idea; Dahl et al.
#' 2022). An optional heatmap of centred log ratios below the bars gives
#' every drawn feature its own aligned row.
#'
#' The bars show **proportions only**: each sample closed to 1. They are not
#' a compositional analysis, and they say nothing about how much of anything
#' a sample held. For absolute loads, scale the counts with [qmp_scale()]
#' and draw them with [graphics::barplot()]; for compositional readings of
#' the same samples, see the heatmap (`heatmap = TRUE`), [aitchison()]
#' distances with [pcoa()], the covariance biplot of [pca()] with `norm =
#' "clr"` (`plot(..., scaling = "covariance")`), and the planned
#' `variation_array()` and `ternary_plot()`.
#'
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features), or
#'        proportions; see **Input Types**. Samples whose counts are all
#'        zero are left out, with a warning.
#'
#' @param groups  Draws the samples in panels, one per group, separated by
#'        a gap and labelled beneath: the name of a column of `data`, or a
#'        vector with one value per sample (matched by name when both are
#'        named, as for [permanova()]). Every distinct value is a group, in
#'        the order of a factor's levels or else sorted; `NA` is a group of
#'        its own. Samples are never averaged. Default: `NULL`, no panels.
#'
#' @param data  A data frame of sample metadata, with one row per sample, when
#'        `groups` names one of its columns.
#'
#' @param feature_groups  Each feature's higher-level group, such as its
#'        phylum: a vector with one value per feature, matched to the
#'        features by name when its values are named by every feature of
#'        `counts`, else by position. `NA` is the group `"Unassigned"`.
#'        Default: `NULL`, no grouping.
#'
#' @param n_features  How many features to draw: without `feature_groups`,
#'        the `n_features` with the largest mean proportion (default `12`);
#'        with them, the `n_features` most abundant members of each drawn
#'        group (default `4`, microshades' maximum). A whole number of at
#'        least 1, or `Inf` for all.
#'
#' @param n_feature_groups  With `feature_groups`, how many groups to draw
#'        in colour: those with the largest mean proportion. The others are
#'        pooled as `"Other"`. A whole number of at least 1, or `Inf`.
#'        Default: `5`.
#'
#' @param clust  How to order the samples (within each group, with
#'        `groups`), from their [aitchison()] distances: `"pcoa"` (the
#'        default) along the first axis of [pcoa()], which is the first
#'        principal component of their centred log ratios; one of
#'        [stats::hclust()]'s linkage methods (`"average"`, `"complete"`,
#'        `"single"`, `"ward.D"`, `"ward.D2"`, `"mcquitty"`, `"median"` or
#'        `"centroid"`) for its leaf order; or `"none"`, `FALSE` or `NULL`
#'        for their input order. The same values as [beta_div()]'s `clust`.
#'
#' @param heatmap  Also draw, below the bars, a heatmap of the drawn
#'        features' centred log ratios (see Details). Default: `FALSE`.
#'
#' @param pseudocount  Added to every count before taking logs, for the
#'        ordering and the heatmap only: the bars show the proportions
#'        themselves. Default: `NULL`, by the rule every log-ratio method in
#'        ecodive shares: none when `counts` has no zeros, else half its
#'        smallest non-zero value. Unlike elsewhere, no warning is given,
#'        since the pseudocount only orders the bars and colours the heatmap.
#'
#' @param legend  Where to draw the legend, outside the plot: `"right"` (the
#'        default), `"bottom"`, or `NULL` (or `FALSE`) for none.
#'
#' @param col  Colours to replace the defaults: one per drawn segment, in
#'        stack order (bottom to top), or named by segment (see Value) to
#'        replace only those.
#'
#' @param ...  Passed on to [graphics::barplot()], such as `main`, `ylab`
#'        or `cex.names`. `sub` replaces the subtitle.
#'
#' @details
#' **Which features.** Each sample's counts are divided by their sum, and
#' the features are ranked by their mean proportion over the samples (ties
#' in their order in `counts`).
#'
#' * Without `feature_groups`: the `n_features` most abundant features are
#'   drawn, and the rest summed into `"Other"`.
#' * With `feature_groups`: the `n_feature_groups` groups with the largest
#'   mean proportion (the sum of their members') are drawn; within each,
#'   its `n_features` most abundant members, and the rest of the group as
#'   `"Other <group>"`. The remaining groups are pooled as `"Other"`.
#'   microshades ranks by the sum over samples instead, which ranks the
#'   same.
#'
#' Features whose proportion is zero in every sample are never drawn, and
#' segments that would hold nothing are left out.
#'
#' **Stack order**, from the bottom: the drawn groups from most to least
#' abundant, within each its members from most to least abundant and then
#' `"Other <group>"`; `"Other"` on top. Without groups, the features from
#' most to least abundant, then `"Other"`. The legend lists the segments
#' from the top down in this order (the most abundant first), each group's
#' under its name in bold.
#'
#' **Colours.** Group \eqn{g} of the \eqn{G} drawn has hue \eqn{h_g = 15 +
#' 360 (g - 1) / G}; its \eqn{k} drawn members are shaded
#' `hcl(h_g, c = 55, l = seq(35, 80, length.out = k + 1))`, darkest for
#' the most abundant, and its `"Other <group>"` takes the lightest. Without
#' groups, the features take [grDevices::hcl.colors()]'s `"Dark 3"`
#' palette. `"Other"` is `"grey80"`. These follow microshades' scheme (one
#' hue per group, shades from dark to light, the group's "Other" lightest),
#' re-derived in HCL rather than copied. Distinct hues with a luminance ramp
#' are reasonably robust to colour-vision deficiencies; for strict needs,
#' give `col`.
#'
#' **Sample order.** By `clust`, from the samples' [aitchison()] distances
#' (with `pseudocount`), as [beta_div()] orders them (`clust = "pcoa"` is
#' `order(pcoa(aitchison(counts))[, 1])`). With `groups`, the samples are
#' then sorted by group, keeping that order within each group. The
#' subtitle says how the samples are ordered, and how much of the total
#' `"Other"` holds (its mean proportion over the samples).
#'
#' **Heatmap.** With `heatmap = TRUE`, the drawn features' centred log
#' ratios: each sample's clr (\eqn{\log(x + c)} minus the sample's mean of
#' those logs over **every** feature, with \eqn{c} the pseudocount), minus
#' the feature's mean clr over the drawn samples. Red is above the
#' feature's average, blue below, on a scale symmetric about 0. Rows are in
#' the legend's order and columns line up with the bars. This is the
#' compositional (CoDa-consistent) reading of the same samples: a
#' feature's row compares it with the sample's geometric mean, which does
#' not depend on how the bars close each sample to 1.
#'
#' Up to 60 samples are labelled with their names; more are drawn without
#' names, and with no space between bars (except between groups).
#'
#' @section Why stacked bars mislead, and what this one does:
#' No peer-reviewed paper critiques microbiome stacked bars on their own:
#' the critiques below come from the motivations of tool papers, and from
#' tutorials and forum threads.
#'
#' | Critique (source) | How this plot answers it |
#' |:------------------|:-------------------------|
#' | Segments' lengths are hard to compare across bars (Peeters et al. 2024) | Samples are ordered by their compositions, so neighbouring bars are similar and gradients read as gradients; the most abundant feature sits on the common baseline; `heatmap = TRUE` gives every feature its own aligned row |
#' | "Other" hides rare taxa, and many taxa make colour hard (Dahl et al. 2022; Peeters et al. 2024) | Hierarchical shading: one hue per group, shades for its top members, the group's own "Other" lightest; the subtitle says how much "Other" holds |
#' | Averaging samples into one bar per group distorts (tutorials and forum threads, not peer reviewed) | Never averaged: one bar per sample, and groups are panels |
#' | Closure hides absolute load (Gloor et al. 2017) | The bars are proportions only; for loads, use [qmp_scale()] and a plain [graphics::barplot()] |
#' | Stacked bars ignore phylogeny and differential abundance (Cuisiniere & Santos 2024) | Out of scope: see [plot.ecodive_ord()] and the plots of [da_test()]'s results |
#'
#' @return Invisibly, a list of what was drawn:
#'
#' * `proportions`: a matrix of the segments (rows, in stack order from the
#'   bottom, named by feature, `"Other <group>"` or `"Other"`) by the drawn
#'   samples (columns, in their drawn order). Each column sums to 1.
#' * `col`: the segments' colours, named by segment.
#' * `order`: the drawn samples' positions in `counts`, in their drawn
#'   order.
#' * `x`: the bars' positions on the x axis, as [graphics::barplot()]
#'   returns them.
#' * `clr`: with `heatmap = TRUE`, the heatmap's values (drawn features by
#'   drawn samples).
#'
#' @references
#' Dahl, E. M., Neer, E., Bowie, K. R., Leung, E. T., & Karstens, L.
#' (2022). microshades: an R package for improving color accessibility and
#' organization of microbiome data. *Microbiology Resource Announcements*,
#' 11, e00795-22. \doi{10.1128/mra.00795-22}
#'
#' Peeters, J., Bot, D. M., Rovelo Ruiz, G., & Aerts, J. (2024). Snowflake:
#' visualizing microbiome abundance tables as multivariate bipartite graphs.
#' *Frontiers in Bioinformatics*, 4, 1331043. \doi{10.3389/fbinf.2024.1331043}
#'
#' Cuisiniere, T., & Santos, M. M. (2024). StackbarExtended: a user-friendly
#' stacked bar-plot representation incorporating phylogenetic information
#' and microbiota differential abundance analysis. *F1000Research*.
#' \doi{10.12688/f1000research.151662.1}
#'
#' Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
#' (2017). Microbiome datasets are compositional: and this is not optional.
#' *Frontiers in Microbiology*, 8, 2224. \doi{10.3389/fmicb.2017.02224}
#'
#' @seealso [aitchison()], [pcoa()], [normalize()], [qmp_scale()]
#'
#' @export
#' @examples
#'     # Four body sites (ex_counts has samples in rows).
#'     composition_barplot(ex_counts)
#'
#'     # A larger table: 30 samples along a gradient, 40 features in 5 phyla.
#'     n <- 30; p <- 40
#'     grad   <- seq(0, 1, length.out = n)
#'     counts <- round(outer(1 - grad, 2^(-(1:p) / 6)) * 400 +
#'                     outer(grad, 2^(-(p:1) / 6)) * 400)
#'     dimnames(counts) <- list(paste0('S', 1:n), paste0('Genus', 1:p))
#'     counts <- counts[order((1:n * 7) %% n), ]   # shuffled
#'     phylum <- rep(c('Firmicutes', 'Bacteroidota', 'Proteobacteria',
#'                     'Actinobacteriota', 'Verrucomicrobiota'), each = 8)
#'
#'     composition_barplot(counts, feature_groups = phylum, n_feature_groups = 3)
#'
#'     # With panels by group, and the clr heatmap below.
#'     site <- ifelse(as.integer(sub('S', '', rownames(counts))) <= n / 2, 'A', 'B')
#'     res  <- composition_barplot(counts, groups = site, heatmap = TRUE, n_features = 8)
#'     res$proportions[, 1:4]
composition_barplot <- function (
    counts, groups = NULL, data = NULL, feature_groups = NULL,
    n_features = NULL, n_feature_groups = 5L, clust = 'pcoa',
    heatmap = FALSE, margin = 1L, pseudocount = NULL,
    legend = 'right', col = NULL, cpus = n_cpus(), ... ) {

  args <- list(...)

  validate_counts()
  validate_margin()
  validate_cpus()
  validate_clust()
  if (is.null(clust)) clust <- 'none'

  if (!isTRUE(heatmap) && !isFALSE(heatmap))
    stop('`heatmap` must be TRUE or FALSE.', call. = FALSE)
  legend <- cb_legend(legend)

  if (!is.null(n_features)) {
    n_features <- plot_count(n_features, 'n_features')
    if (n_features < 1) stop('`n_features` must be at least 1.', call. = FALSE)
  }
  n_feature_groups <- plot_count(n_feature_groups, 'n_feature_groups')
  if (n_feature_groups < 1) stop('`n_feature_groups` must be at least 1.', call. = FALSE)

  vals <- stored_values(counts)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'


  #--------------------------------------------------
  # Names and sizes, by `margin`.
  #--------------------------------------------------
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  n      <- dims[[margin]]
  p      <- dims[[3L - margin]]
  snames <- dn[[margin]]
  fnames <- dn[[3L - margin]]
  if (is.null(fnames)) fnames <- paste('Feature', seq_len(p))
  fnames <- as.character(fnames)


  #--------------------------------------------------
  # Each feature's mean proportion, over the samples with any counts.
  #--------------------------------------------------
  mt   <- .Call(C_comp_bar_means, counts, margin, cpus)
  keep <- which(mt$totals > 0)
  if (length(keep) == 0)
    stop('`counts` has no sample with counts above 0, so there is nothing to draw.', call. = FALSE)
  if (length(keep) < n) {
    drop <- setdiff(seq_len(n), keep)
    warning(
      sprintf('Left out %d sample%s whose counts are all zero: %s',
        length(drop), if (length(drop) == 1) '' else 's',
        name_list(if (is.null(snames)) drop else snames[drop]) ),
      call. = FALSE )
  }

  seg <- cb_select(mt$means, fnames, feature_groups, n_features, n_feature_groups, !is.null(dn[[3L - margin]]))
  K   <- length(seg$label)


  #--------------------------------------------------
  # The samples' order: by `clust`, then by group.
  #--------------------------------------------------
  need_pc <- heatmap || (clust != 'none' && length(keep) >= 3)
  pc <- if (need_pc || !is.null(pseudocount)) cb_pseudocount(counts, pseudocount) else 0

  base <- seq_along(keep)
  if (clust != 'none' && length(keep) >= 3 && !(pc == 0 && min(vals) == max(vals))) {
    d <- aitchison(counts, margin = margin, pseudocount = if (pc > 0) pc, cpus = cpus)
    if (length(keep) < n) d <- cb_dist_subset(d, keep)
    # Identical samples (all distances 0) have no PCoA axis: pcoa() warns,
    # and their input order is kept.
    base <- withCallingHandlers(
      dist_clust_order(d, clust, cpus),
      warning = function (w) {
        if (grepl('eigenvalues are > 0', conditionMessage(w))) invokeRestart('muffleWarning')
      } )
  }

  g <- cb_groups(groups, data, n, snames)
  panels <- NULL
  ord    <- base
  if (!is.null(g)) {
    codes  <- as.integer(g)[keep]
    ord    <- base[order(codes[base])]   # stable: keeps `base` within groups
    runs   <- rle(codes[ord])
    lev    <- levels(g)
    panels <- list(
      starts = cumsum(c(1L, utils::head(runs$lengths, -1L))),
      ends   = cumsum(runs$lengths),
      labels = ifelse(is.na(lev), 'NA', lev)[runs$values] )
  }
  samples <- keep[ord]
  m       <- length(samples)


  #--------------------------------------------------
  # The segments' proportions, and the heatmap's clr.
  #--------------------------------------------------
  rows <- if (heatmap) seg$heat_row else NULL
  res  <- .Call(
    C_comp_bar_segments, counts, margin, seg$seg, K, seg$other,
    samples - 1L, rows, length(seg$heat), as.double(pc), cpus )

  cnames <- if (is.null(snames)) as.character(samples) else as.character(snames[samples])
  P <- res$proportions
  dimnames(P) <- list(seg$label, cnames)

  clr <- NULL
  if (heatmap) {
    clr <- res$clr
    clr <- clr - rowMeans(clr)
    dimnames(clr) <- list(seg$label[seg$heat], cnames)
  }

  cols <- cb_colours(seg, col)


  #--------------------------------------------------
  # Draw.
  #--------------------------------------------------
  x <- cb_draw(P, cols, seg, clr, panels, clust, !is.null(g), legend, args)

  out <- list(proportions = P, col = cols, order = samples, x = x)
  if (heatmap) out$clr <- clr
  invisible(out)
}



# Checks `legend`: "right", "bottom", or NULL, FALSE or NA for none.
cb_legend <- function (legend) {
  if (is.null(legend) || isFALSE(legend) || identical(legend, NA)) return (NULL)
  if (!is.character(legend) || length(legend) != 1 || !legend %in% c('right', 'bottom'))
    stop('`legend` must be "right", "bottom", or NULL for none.', call. = FALSE)
  legend
}


# The pseudocount of the ordering and the heatmap, by the shared rule
# (log_pseudocount()), without its warning about zeros: it changes no bar.
cb_pseudocount <- function (counts, pseudocount) {
  what <- 'the Aitchison ordering'
  if (is.null(pseudocount))
    withCallingHandlers(
      log_pseudocount(counts, NULL, what, constant_ok = TRUE),
      warning = function (w) {
        if (grepl('^Zeros detected', conditionMessage(w))) invokeRestart('muffleWarning')
      } )
  else log_pseudocount(counts, pseudocount, what, constant_ok = TRUE)
}


# The distances of `d` between the samples `keep` only, as a dist.
cb_dist_subset <- function (d, keep) {
  n  <- attr(d, 'Size', exact = TRUE)
  m  <- length(keep)
  ij <- which(lower.tri(matrix(0, m, m)), arr.ind = TRUE)
  i  <- keep[ij[, 2]]   # the smaller index of each pair, column by column
  j  <- keep[ij[, 1]]
  k  <- n * (i - 1) - i * (i - 1) / 2 + j - i
  lab <- attr(d, 'Labels', exact = TRUE)
  structure(
    as.vector(d)[k], Size = m, Labels = if (!is.null(lab)) lab[keep],
    Diag = FALSE, Upper = FALSE, method = attr(d, 'method', exact = TRUE),
    class = class(d) )
}


# Each sample's group, as a factor with NA as a level (if any). Every
# distinct value is a group, numbers too. NULL without `groups`.
cb_groups <- function (groups, data, n, labels) {

  if (is.null(groups)) return (NULL)

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, dn, nrow(data), ref = '`counts`')
  }
  v <- per_sample(groups, 'groups', data, idx, labels, n, ref = '`counts`')

  f <- if (is.factor(v)) v else factor(v)
  if (anyNA(f)) f <- addNA(f, ifany = TRUE)
  f
}



# Which features make which segments. Returns a list:
#   label     each segment's name, in stack order (bottom first);
#   group     each segment's group (NA without feature groups, or for
#             "Other");
#   kind      each segment's kind: "feature", "group_other" or "other";
#   shade     each segment's position within its group (1 = darkest), and
#   n_shades  its group's number of drawn members plus 1;
#   groups    the drawn groups, most abundant first (NULL without);
#   seg       each feature's 0-based segment, or -1 (in "Other" or none);
#   other     the 0-based "Other" segment, or -1;
#   heat      the segments that are single features (the heatmap's rows,
#             in legend order), and heat_row each feature's 0-based row
#             among them, or -1.
cb_select <- function (means, fnames, feature_groups, n_features, n_feature_groups, named) {

  p    <- length(means)
  rank <- function (v) order(-v, seq_along(v))   # largest first; ties in order
  segs <- list()
  add  <- function (label, group, kind, feats, shade = NA, n_shades = NA)
    segs[[length(segs) + 1L]] <<- list(
      label = label, group = group, kind = kind, feats = feats,
      shade = shade, n_shades = n_shades )

  if (is.null(feature_groups)) {

    nf  <- if (is.null(n_features)) 12 else n_features
    pos <- rank(means)
    pos <- pos[means[pos] > 0]
    top <- utils::head(pos, nf)
    for (j in top) add(fnames[[j]], NA_character_, 'feature', j)
    rest <- setdiff(pos, top)
    if (length(rest)) add('Other', NA_character_, 'other', rest)
    drawn_groups <- NULL
  }

  else {

    fg  <- cb_feature_groups(feature_groups, fnames, p, named)
    lev <- attr(fg, 'levels_order')
    gm  <- vapply(lev, function (g) sum(means[fg == g]), double(1))
    gpos <- rank(gm)
    gpos <- gpos[gm[gpos] > 0]
    nfg  <- if (is.null(n_features)) 4 else n_features
    drawn_groups <- lev[utils::head(gpos, n_feature_groups)]

    for (g in drawn_groups) {
      mem <- which(fg == g)
      mem <- mem[rank(means[mem])]
      mem <- mem[means[mem] > 0]
      top <- utils::head(mem, nfg)
      rest <- setdiff(mem, top)
      k    <- length(top)
      for (i in seq_along(top)) add(fnames[[top[[i]]]], g, 'feature', top[[i]], i, k + 1L)
      if (length(rest)) add(paste('Other', g), g, 'group_other', rest, k + 1L, k + 1L)
    }

    rest <- which(!fg %in% drawn_groups & means > 0)
    if (length(rest)) add('Other', NA_character_, 'other', rest)
  }

  label <- make.unique(vapply(segs, `[[`, '', 'label'))
  kind  <- vapply(segs, `[[`, '', 'kind')
  K     <- length(segs)

  seg <- rep(-1L, p)
  for (k in seq_len(K)) if (kind[[k]] != 'other') seg[segs[[k]]$feats] <- k - 1L

  heat     <- which(kind == 'feature')
  heat_row <- rep(-1L, p)
  for (r in seq_along(heat)) heat_row[segs[[heat[[r]]]]$feats] <- r - 1L

  list(
    label    = label,
    group    = vapply(segs, function (s) s$group, ''),
    kind     = kind,
    shade    = vapply(segs, function (s) as.integer(s$shade), 0L),
    n_shades = vapply(segs, function (s) as.integer(s$n_shades), 0L),
    groups   = drawn_groups,
    seg      = seg,
    other    = if (any(kind == 'other')) which(kind == 'other') - 1L else -1L,
    heat     = heat,
    heat_row = heat_row )
}


# `feature_groups` as one character value per feature ("Unassigned" for
# NA), with attribute `levels_order`: the groups in the order that breaks
# ties (a factor's levels, else first appearance).
cb_feature_groups <- function (fg, fnames, p, named) {

  if (!is.atomic(fg) || !is.null(dim(fg)) || length(fg) == 0)
    stop('`feature_groups` must be a vector with one value per feature of `counts`.', call. = FALSE)

  nm <- names(fg)
  if (named && !is.null(nm) && !anyNA(nm) && !anyDuplicated(nm) && all(fnames %in% nm)) {
    fg <- fg[match(fnames, nm)]
  }
  else if (length(fg) != p) {
    stop(
      sprintf('`feature_groups` has %d value%s, but `counts` has %d feature%s.',
        length(fg), if (length(fg) == 1) '' else 's', p, if (p == 1) '' else 's'),
      if (named && !is.null(nm))
        '\nName its values by feature to match them to the features by name.',
      call. = FALSE )
  }

  lev <- if (is.factor(fg)) levels(fg) else NULL
  fg  <- as.character(fg)
  fg[is.na(fg)] <- 'Unassigned'
  lev <- unique(c(lev, fg))
  lev <- lev[lev %in% fg]
  attr(fg, 'levels_order') <- lev
  fg
}


# The segments' colours, named by segment: microshades' scheme in HCL,
# then the user's `col`.
cb_colours <- function (seg, col) {

  K    <- length(seg$label)
  cols <- rep('grey80', K)
  feat <- seg$kind != 'other'

  if (is.null(seg$groups)) {
    cols[feat] <- plot_palette(sum(feat))
  }
  else {
    G   <- length(seg$groups)
    hue <- seq(15, 375, length.out = G + 1L)
    for (k in which(feat)) {
      h <- hue[[match(seg$group[[k]], seg$groups)]]
      l <- seq(35, 80, length.out = seg$n_shades[[k]])
      cols[[k]] <- grDevices::hcl(h, c = 55, l = l[[seg$shade[[k]]]])
    }
  }
  names(cols) <- seg$label

  if (!is.null(col)) {
    if (!is.character(col) && !is.numeric(col))
      stop('`col` must be a vector of colours.', call. = FALSE)
    nm <- names(col)
    if (!is.null(nm) && all(nzchar(nm))) {
      bad <- setdiff(nm, seg$label)
      if (length(bad))
        stop('`col` names segments that are not drawn: ', name_list(bad), '.', call. = FALSE)
      cols[nm] <- col
    }
    else if (length(col) == K) {
      cols[] <- col
    }
    else {
      stop(
        sprintf('`col` must have one colour per segment (%d), or be named by segment.', K),
        call. = FALSE )
    }
  }

  cols
}



# How the samples were ordered, for the subtitle.
cb_order_name <- function (clust, grouped) {
  how <- if (identical(clust, 'pcoa')) 'Samples ordered by Aitchison PCoA axis 1'
         else if (identical(clust, 'none')) 'Samples in input order'
         else sprintf('Samples clustered by Aitchison distance (%s linkage)', clust)
  if (grouped && !identical(clust, 'none')) how <- paste(how, 'within groups')
  how
}


# The legend's entries, top down: each group's name (bold, no box), then
# its segments; segments in stack order. Returns list(text, fill, font).
cb_legend_entries <- function (seg, cols) {
  text <- character(0); fill <- character(0); font <- integer(0)
  for (k in seq_along(seg$label)) {
    g <- seg$group[[k]]
    if (!is.na(g) && (k == 1L || !identical(seg$group[[k - 1L]], g))) {
      text <- c(text, g); fill <- c(fill, NA); font <- c(font, 2L)
    }
    text <- c(text, seg$label[[k]]); fill <- c(fill, cols[[k]]); font <- c(font, 1L)
  }
  list(text = text, fill = fill, font = font)
}


# Draws the bars (and the heatmap, the legend and the group labels).
# Returns the bars' x positions.
cb_draw <- function (P, cols, seg, clr, panels, clust, grouped, legend, args) {

  K <- nrow(P)
  m <- ncol(P)
  heatmap <- !is.null(clr)

  # Bars of width 1; 0.1 apart (none above 60 samples), and one more bar
  # width before each group but the first.
  space <- rep(if (m > 60) 0 else 0.1, m)
  if (!is.null(panels)) {
    st <- panels$starts[-1L]
    space[st] <- space[st] + 1
  }
  right <- cumsum(space + 1)
  left  <- right - 1
  mids  <- left + 0.5
  xlim  <- c(min(left), max(right))

  names_on <- m <= 60
  lcex     <- 0.7

  sub <- if ('sub' %in% names(args)) args[['sub']] else {
    s <- cb_order_name(clust, grouped)
    if (any(seg$kind == 'other'))
      s <- paste0(s, "; 'Other' holds ",
                  plot_num(100 * mean(P[seg$kind == 'other', ]), 2), '% of the total')
    s
  }
  args <- plot_drop(args, 'sub')


  #--------------------------------------------------
  # Margins, in inches: room for the names, group labels and legend.
  #--------------------------------------------------
  old <- graphics::par(c('mai', 'cex', 'mfrow'))
  on.exit({
    if (heatmap) graphics::par(mfrow = old$mfrow)
    graphics::par(mai = old$mai, cex = old$cex)
  })
  if (heatmap) graphics::layout(matrix(1:2), heights = c(2, 1.4))

  csi   <- graphics::par('csi')
  namew <- if (names_on) max(graphics::strwidth(colnames(P), 'inch', cex = 0.6), 0) + 0.1 else 0
  grph  <- if (!is.null(panels)) 1.3 * csi * 0.8 else 0
  ent   <- cb_legend_entries(seg, cols)
  entw  <- max(graphics::strwidth(ent$text, 'inch', cex = lcex, font = 2), 0) + 0.45
  entw  <- entw + 0.15

  mai <- graphics::par('mai')
  mai[[3]] <- max(mai[[3]], 0.9)
  if (heatmap) {
    heatw    <- max(graphics::strwidth(rownames(clr), 'inch', cex = 0.7), 0) + 0.25
    mai[[2]] <- max(mai[[2]], heatw)
  }
  if (identical(legend, 'right')) mai[[4]] <- max(mai[[4]], entw)
  # The heatmap's key (dist_plot_key()) sits 4% to 8% of the plot's width
  # right of it, with its labels after: room for all of that.
  if (heatmap)
    mai[[4]] <- max(mai[[4]], (0.1 * (graphics::par('fin')[[1]] - mai[[2]]) + 0.55) / 1.1)

  lrows <- 0
  lcols <- 1L
  if (identical(legend, 'bottom')) {
    pw    <- graphics::par('fin')[[1]] - mai[[2]] - mai[[4]]
    lcols <- max(1L, floor(pw / entw))
    lrows <- ceiling(length(ent$text) / lcols)
  }
  legh <- if (lrows) lrows * csi * lcex * 1.2 + 0.2 else 0
  bot  <- namew + grph + 0.25

  # The bars' panel: names, group labels and the legend below it unless the
  # heatmap is.
  mai1 <- mai
  mai1[[1]] <- if (heatmap) 0.15 else max(mai[[1]], bot + legh)
  graphics::par(mai = mai1)

  plot_do(
    graphics::barplot, c(list(height = P), args),
    list(
      space = space, width = 1, col = cols, border = NA, xlim = xlim,
      ylim = c(0, 1), names.arg = colnames(P), axisnames = names_on && !heatmap,
      las = 2, cex.names = 0.6, ylab = 'Proportion',
      main = 'Composition Barplot' ))
  graphics::mtext(sub, side = 3, line = 0.4, cex = 0.8)

  if (!is.null(legend) && identical(legend, 'right'))
    cb_draw_legend(ent, lcex, 'right', 0)

  bottom_panel <- function (offset) {
    if (!is.null(panels)) {
      gm <- (left[panels$starts] + right[panels$ends]) / 2
      graphics::mtext(
        panels$labels, side = 1, at = gm, line = offset / csi + 0.3, cex = 0.8, font = 2 )
    }
    if (identical(legend, 'bottom'))
      cb_draw_legend(ent, lcex, 'bottom', namew + grph + 0.15, lcols)
  }

  if (!heatmap) {
    bottom_panel(namew)
    return (invisible(mids))
  }


  #--------------------------------------------------
  # The heatmap: one cell per bar, NA in the gaps.
  #--------------------------------------------------
  mai2 <- mai
  mai2[[1]] <- max(0.3, bot + legh)
  mai2[[3]] <- 0.1
  graphics::par(mai = mai2)

  # Cells: each bar's own [left, right], and an NA cell for each gap.
  R    <- nrow(clr)
  gap  <- c(FALSE, left[-1L] > right[-m] + 1e-9)
  bnd  <- c(left[[1]], as.vector(rbind(ifelse(gap, left, NA), right)))
  bnd  <- bnd[!is.na(bnd)]
  idx  <- as.vector(rbind(ifelse(gap, NA_integer_, -1L), seq_len(m)))
  idx  <- idx[is.na(idx) | idx > 0L]
  z    <- t(clr)[idx, rev(seq_len(R)), drop = FALSE]   # first row at the top
  fin <- is.finite(z)
  top <- if (any(fin)) max(abs(z[fin])) else 0
  if (!(top > 0)) top <- 1
  pal <- grDevices::hcl.colors(65, 'Blue-Red 3')
  brk <- seq(-top, top, length.out = length(pal) + 1L)

  graphics::image(
    x = bnd, y = seq(0.5, R + 0.5), z = z, col = pal, breaks = brk,
    xlim = xlim, xaxs = 'r', axes = FALSE, xlab = '', ylab = '' )
  rcex <- min(0.7, 0.9 * graphics::par('pin')[[2]] / R / csi)
  graphics::axis(2, at = R:1, labels = rownames(clr), las = 1, tick = FALSE, cex.axis = rcex)
  if (names_on)
    graphics::axis(1, at = mids, labels = colnames(clr), las = 2, tick = FALSE, cex.axis = 0.6)
  dist_plot_key(pal, brk, 'clr')

  bottom_panel(namew)
  invisible(mids)
}


# The legend, outside the plot region: on the right (top-aligned), or
# below it, `below` inches under the plot region's bottom, in `ncol`
# columns.
cb_draw_legend <- function (ent, cex, where, below, ncol = 1L) {
  usr <- graphics::par('usr')
  border <- ifelse(is.na(ent$fill), NA, 'grey40')
  if (where == 'right') {
    x <- graphics::grconvertX(graphics::grconvertX(usr[[2]], 'user', 'inches') + 0.1, 'inches', 'user')
    graphics::legend(
      x, usr[[4]], legend = ent$text, fill = ent$fill, border = border,
      text.font = ent$font, cex = cex, bty = 'n', xjust = 0, yjust = 1, xpd = NA )
  }
  else {
    y <- graphics::grconvertY(graphics::grconvertY(usr[[3]], 'user', 'inches') - below, 'inches', 'user')
    graphics::legend(
      mean(usr[1:2]), y, legend = ent$text, fill = ent$fill, border = border,
      text.font = ent$font, cex = cex, bty = 'n', xjust = 0.5, yjust = 1,
      ncol = ncol, xpd = NA )
  }
}
