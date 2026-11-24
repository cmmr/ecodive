# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The depth grid (np.linspace(1, max_depth, num = steps, dtype = int)) and
# the resampling design (each depth and iteration rarefied, then every
# metric computed from that table) follow q2-diversity's alpha-rarefaction
# (BSD-3, (c) 2016-2026 QIIME 2 development team). No code is ported. The
# closed forms (src/rarefaction_curve.c) are written from Hurlbert (1971),
# Heck et al. (1975), Chao & Jost (2012) and Chao et al. (2014). See
# plans/rarefaction_curves.md.



#' Rarefaction Curves
#'
#' Alpha diversity against sequencing depth: each sample's diversity when it
#' is subsampled (rarefied) to a range of depths. A curve that levels off
#' says the sample was sequenced deeply enough to see most of its features;
#' one still rising says deeper sequencing would find more. The curves also
#' show which depth to rarefy to, and how many samples that depth keeps.
#'
#' Two ways to compute a curve:
#'
#' * **Exactly, in closed form** (`'observed'`): the expected number of
#'   features in a random subsample of \eqn{d} of a sample's \eqn{N}
#'   observations, drawn without replacement (Hurlbert 1971),
#'   \deqn{E[S_d] = \sum_i \left[1 - \binom{N - x_i}{d} \Big/ \binom{N}{d}\right],}
#'   and its standard deviation (Heck et al. 1975). This is what
#'   `vegan::rarefy()` and `vegan::rarecurve()` give. There is no randomness,
#'   and `times` and `seed` are not used.
#' * **By resampling**, for any alpha metric (QIIME 2's approach): the
#'   table is rarefied `times` times at each depth, and the metric's values
#'   are averaged over the rarefactions.
#'
#' Each row also has the **expected sample coverage** of a depth-\eqn{d}
#' subsample (Chao & Jost 2012), always in closed form,
#' \deqn{C_d = 1 - \sum_i \frac{x_i}{N} \binom{N - x_i}{d} \Big/ \binom{N - 1}{d},}
#' the share of the community's individuals that belong to features seen in
#' such a subsample. At the sample's own total (\eqn{d = N}), it is the
#' estimated coverage of the whole sample, \eqn{1 - \frac{f_1}{N}
#' \frac{(N - 1) f_1}{(N - 1) f_1 + 2 f_2}}, where \eqn{f_1} and \eqn{f_2}
#' are the numbers of features seen once and twice (Chao et al. 2014, Eq. 12;
#' with \eqn{f_2 = 0}, \eqn{f_1 (f_1 - 1)} replaces \eqn{f_1^2}, as in their
#' Appendix G). These are the sample-completeness curves, and the x axis of
#' coverage-based curves; [rarefy()] uses them to rarefy to equal coverage.
#'
#' @inheritParams alpha_div
#' @inheritParams rarefied_diversity
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features) of whole
#'        numbers; see **Input Types**. A `phyloseq`, `rbiom` or
#'        `TreeSummarizedExperiment` object's tree is used for the
#'        phylogenetic metrics.
#'
#' @param metric  One or more alpha diversity metrics, such as
#'        `c('observed', 'shannon')`: any of [alpha_div()]'s, matched as it
#'        matches them (case-insensitive, with prefixes and aliases). See
#'        [list_metrics()]. Every rarefied table is used for all of them.
#'        Default: `'observed'`.
#'
#' @param depths  The depths to compute the curves at: whole numbers of at
#'        least 1 (sorted, and duplicates dropped). Default: `NULL`, a grid
#'        of `n_depths` depths from 1 to `max_depth` (see "The depths").
#'
#' @param n_depths  The number of depths in the default grid, a whole number
#'        of at least 2. Default: `20`.
#'
#' @param max_depth  The largest depth computed. Default: `NULL`, the
#'        largest library size (with `depths`, the largest of them).
#'
#' @param times  The number of rarefactions at each depth, for the metrics
#'        that are resampled, a whole number of at least 1. Default: `10`,
#'        as QIIME 2's `iterations`.
#'
#' @param analytic  Which metrics to compute in closed form. `NULL` (the
#'        default): `'observed'` in closed form and every other metric by
#'        resampling. `FALSE`: resample `'observed'` too (to compare with
#'        QIIME 2). `TRUE`: every metric in closed form, which only
#'        `'observed'` has (an error otherwise).
#'
#' @param average  How the rarefactions' values at a depth are combined, for
#'        the resampled metrics: `"mean"` (the default) or `"median"`.
#'        Case-insensitive, and a prefix is enough.
#'
#' @param seed  An integer seed for the rarefactions: rarefaction \eqn{t}
#'        at every depth uses `seed + t - 1`, as `rarefy(times = )` does (see
#'        "Nested rarefactions"). The same seed gives the same result, for
#'        any `cpus`. Default: `0`.
#'
#' @return
#' A data frame of class `c("ecodive_rarefaction_curve", "data.frame")`,
#' with one row per sample, metric and depth, ordered by sample, then
#' metric (in the order given), then depth. Columns:
#' * `sample`: the sample's name, or its position (as text) in an unnamed
#'   table.
#' * `metric`: the metric's ID, as in [list_metrics()].
#' * `depth`: the depth (an integer).
#' * `value`: the metric's expected value (closed form), or the mean or
#'   median over the `times` rarefactions.
#' * `sd`: the standard deviation of a subsample's value: Heck et al.'s
#'   (closed form; 0 at the sample's total), or over the rarefactions (with
#'   \eqn{n - 1} in the denominator; `NA` with `times = 1`).
#' * `coverage`: the expected sample coverage at that depth. It does not
#'   depend on the metric.
#'
#' Attributes:
#' * `totals`: each sample's library size, named by sample.
#' * `analytic`: the IDs of the metrics computed in closed form (possibly
#'   none).
#' * `times` and `average`: how the resampled metrics' `value` and `sd`
#'   were made.
#'
#' @section The depths:
#' With `depths = NULL`, the depths are `n_depths` points from 1 to
#' `max_depth`, as QIIME 2's `np.linspace(1, max_depth, num = n_depths,
#' dtype = int)`: \eqn{1 + i (M - 1) / (n - 1)}, \eqn{i = 0, \ldots, n - 1},
#' rounded down, with the last exactly \eqn{M}, and duplicates dropped.
#'
#' **Each sample's own library size is added to its depths** (when it is at
#' most `max_depth`), so that each curve ends at the sample's observed
#' value. A sample has no rows at depths above its library size (QIIME 2
#' gives `NaN` there), and an empty sample has none at all.
#'
#' * QIIME 2's `alpha-rarefaction`: `n_depths = 10`, `analytic = FALSE`,
#'   `average = "median"`, with `max_depth` as given to it.
#' * `vegan::rarecurve()` (every depth): `depths = seq_len(max(totals))`,
#'   where `totals` are the library sizes; it computes richness in closed
#'   form, as `rarefaction_curve()` does.
#'
#' @section Nested rarefactions:
#' Rarefaction \eqn{t} at depth \eqn{d} is exactly `rarefy(counts, d, seed
#' = seed + t - 1)`, with the same seed at every depth. [rarefy()] keeps
#' observations by Knuth's selection sampling, one random draw per
#' observation, so with the same seed a smaller depth keeps a subset of
#' what a larger one keeps: each rarefaction is one accumulation of
#' observations, and its richness never decreases with depth. The curves
#' are smooth rather than jagged. QIIME 2 draws each depth independently;
#' the averages agree (up to random error), but single rarefactions differ.
#'
#' @section Choosing a depth:
#' Rarefy to a depth where most curves have levelled off, and that keeps the
#' samples you need: `plot(x, depth = 1000)` draws a line at that depth,
#' labelled with how many samples have at least that many observations
#' (from the `totals` attribute). To compare samples at equal completeness
#' instead, see the `coverage` column (`plot(x, type = "completeness")`) and
#' `rarefy(counts, coverage = )`.
#'
#' @section Same names in other packages:
#' None, but `vegan::rarecurve()` draws the same curves for richness (and
#' `vegan::rarefy()` gives their values), and `iNEXT::iNEXT()`'s
#' interpolated sample coverage (`SC`) is the `coverage` column.
#'
#' @references
#' Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique
#' and alternative parameters. *Ecology*, 52, 577--586.
#' \doi{10.2307/1934145}
#'
#' Heck, K. L., van Belle, G., & Simberloff, D. (1975). Explicit calculation
#' of the rarefaction diversity measurement and the determination of
#' sufficient sample size. *Ecology*, 56, 1459--1461. \doi{10.2307/1934716}
#'
#' Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. *Ecology*, 93,
#' 2533--2547. \doi{10.1890/11-1952.1}
#'
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
#' R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
#' numbers: a framework for sampling and estimation in species diversity
#' studies. *Ecological Monographs*, 84, 45--67. \doi{10.1890/13-0133.1}
#'
#' Bolyen, E., et al. (2019). Reproducible, interactive, scalable and
#' extensible microbiome data science using QIIME 2. *Nature
#' Biotechnology*, 37, 852--857. \doi{10.1038/s41587-019-0209-9}
#'
#' @seealso [plot.ecodive_rarefaction_curve()], [rarefy()],
#'   [alpha_div_rarefied()], [alpha_div()]
#'
#' @export
#' @examples
#'     # Expected richness, in closed form.
#'     rc <- rarefaction_curve(ex_counts)
#'     rc
#'     plot(rc, depth = 1000)
#'
#'     # Shannon diversity by resampling, coloured by body area.
#'     area <- c(Saliva = 'Mouth', Gums = 'Mouth', Nose = 'Skin', Stool = 'Gut')
#'     rc <- rarefaction_curve(ex_counts, c('observed', 'shannon'), times = 5)
#'     rc
#'     plot(rc, 'shannon', groups = area)
#'
#'     # Sample completeness: how coverage grows with depth.
#'     plot(rc, type = 'completeness')
#'
rarefaction_curve <- function (
    counts,
    metric    = 'observed',
    depths    = NULL,
    n_depths  = 20L,
    max_depth = NULL,
    times     = 10L,
    analytic  = NULL,
    average   = 'mean',
    cutoff    = 10L,
    digits    = 3L,
    tree      = NULL,
    q         = 2,
    bias_corrected = TRUE,
    margin    = 1L,
    seed      = 0L,
    cpus      = n_cpus() ) {

  # The metrics first, so that a misspelt metric is the error reported.
  metrics <- rc_match_metrics(metric)
  ids     <- vapply(metrics, `[[`, '', 'id')

  counts  <- rarefied_counts(environment())
  times   <- rarefied_times(times)

  depths    <- rc_check_depths(depths)
  n_depths  <- rc_check_whole(n_depths, 'n_depths', 2, 'a whole number of at least 2')
  if (!is.null(max_depth))
    max_depth <- rc_check_whole(max_depth, 'max_depth', 1, 'NULL or a whole number of at least 1')
  analytic  <- rc_check_analytic(analytic, ids)
  average   <- rarefied_average(average, 'alpha')

  validate_cutoff()
  validate_digits()
  validate_q()
  validate_bias_corrected()


  # The samples, and each one's depths.
  info   <- count_samples(counts, margin)
  n      <- info$n
  labels <- if (is.null(info$labels)) as.character(seq_len(n)) else info$labels
  totals <- sample_sums(counts, margin)
  names(totals) <- labels

  grid <- rc_grid(depths, n_depths, max_depth, totals)
  top  <- attr(grid, 'max_depth')

  per_sam <- lapply(totals, function (N) {
    d <- grid[grid <= N]
    if (N >= 1 && N <= top) d <- sort(unique(c(d, N)))
    as.integer(d)
  })
  n_rows  <- lengths(per_sam)
  row_pos <- as.integer(c(0L, cumsum(n_rows)))
  sam_idx <- rep(seq_len(n), n_rows)
  depth   <- as.integer(unlist(per_sam, use.names = FALSE))
  if (is.null(depth)) depth <- integer(0)
  R       <- length(depth)


  # Closed forms: expected richness, its sd, and coverage, for every row.
  cf <- .Call(C_rarefaction_curve, counts, margin, row_pos, depth, cpus)


  # Resampled metrics.
  resampled <- setdiff(ids, analytic)
  vals      <- list()
  sds       <- list()

  if (length(resampled) && R > 0) {

    seeds <- rarefied_seeds(seed, times)
    mets  <- metrics[match(resampled, ids)]
    args  <- list(
      cutoff = cutoff, digits = digits, tree = tree, margin = margin,
      q = q, bias_corrected = bias_corrected, cpus = cpus )
    margs <- lapply(mets, function (m) args[intersect(names(args), m$params)])

    # The grid depths that some sample reaches, and the rows at each; the
    # rows at a sample's own total that aren't grid depths.
    gd      <- grid[grid <= max(totals)]
    at_rows <- lapply(gd, function (g) which(depth == g))
    at_keep <- lapply(gd, function (g) totals >= g)
    tot_row <- which(depth == totals[sam_idx] & !depth %in% gd)
    tot_sam <- sam_idx[tot_row]

    acc <- lapply(resampled, function (m) list(
      sum = numeric(R), mean = numeric(R), m2 = numeric(R),
      M = if (average == 'median') matrix(NA_real_, R, times) else NULL ))
    names(acc) <- resampled

    record <- function (m, rows, a, t) {
      st <- acc[[m]]
      st$sum[rows] <- st$sum[rows] + a
      if (average == 'mean') {
        delta <- a - st$mean[rows]
        st$mean[rows] <- st$mean[rows] + delta / t
        st$m2[rows]   <- st$m2[rows] + delta * (a - st$mean[rows])
      }
      else {
        st$M[rows, t] <- a
      }
      acc[[m]] <<- st
    }

    subset_sams <- function (m, keep) {
      if (margin == 1L) m[keep, , drop = FALSE] else m[, keep, drop = FALSE]
    }

    metric_of <- function (k, tbl) {
      as.vector(do.call(mets[[k]]$func, c(list(counts = tbl), margs[[k]])))
    }

    warns <- new.env(parent = emptyenv())
    rarefied_quietly(warns, {

      # At a sample's own total, rarefy() returns the sample as it is.
      base <- list()
      if (length(tot_row)) {
        tbl  <- subset_sams(counts, seq_len(n) %in% tot_sam)
        kept <- which(seq_len(n) %in% tot_sam)
        for (k in seq_along(mets))
          base[[k]] <- metric_of(k, tbl)[match(tot_sam, kept)]
      }

      for (t in seq_len(times)) {
        for (j in seq_along(gd)) {
          tbl <- .Call(C_rarefy, counts, gd[[j]], seeds[[t]], margin, cpus)
          tbl <- subset_sams(tbl, at_keep[[j]])
          for (k in seq_along(mets))
            record(resampled[[k]], at_rows[[j]], metric_of(k, tbl), t)
        }
        if (length(tot_row))
          for (k in seq_along(mets))
            record(resampled[[k]], tot_row, base[[k]], t)
      }
    })

    for (m in resampled) {
      st <- acc[[m]]
      if (average == 'mean') {
        vals[[m]] <- st$sum / times
        sds[[m]]  <- if (times > 1L) sqrt(st$m2 / (times - 1L)) else rep(NA_real_, R)
      }
      else {
        vals[[m]] <- .Call(C_row_medians, st$M, cpus)
        if (times > 1L) {
          mu <- rowMeans(st$M)
          s2 <- numeric(R)
          for (t in seq_len(times)) s2 <- s2 + (st$M[, t] - mu)^2
          sds[[m]] <- sqrt(s2 / (times - 1L))
        }
        else sds[[m]] <- rep(NA_real_, R)
      }
    }

    rarefied_rewarn(warns)
  }

  for (m in analytic) {
    vals[[m]] <- cf$value
    sds[[m]]  <- cf$sd
  }
  for (m in setdiff(ids, names(vals))) {   # no rows to resample
    vals[[m]] <- numeric(0)
    sds[[m]]  <- numeric(0)
  }


  # One block of rows per metric, then ordered by sample, metric, depth.
  k   <- length(ids)
  ord <- order(rep(sam_idx, k), rep(seq_len(k), each = R), rep(depth, k))
  out <- data.frame(
    sample   = rep(labels[sam_idx], k)[ord],
    metric   = rep(ids, each = R)[ord],
    depth    = rep(depth, k)[ord],
    value    = unlist(vals[ids], use.names = FALSE)[ord],
    sd       = unlist(sds[ids],  use.names = FALSE)[ord],
    coverage = rep(cf$coverage, k)[ord],
    stringsAsFactors = FALSE )
  if (R == 0) out <- out[0, , drop = FALSE]
  rownames(out) <- NULL

  attr(out, 'totals')   <- totals
  attr(out, 'analytic') <- analytic
  attr(out, 'times')    <- times
  attr(out, 'average')  <- average
  class(out) <- c('ecodive_rarefaction_curve', 'data.frame')
  out
}



# The alpha metrics named by `metric`, each matched as alpha_div() matches
# one, without duplicates (by ID). A list of match_metric() results.
rc_match_metrics <- function (metric) {
  if (!is.character(metric) || length(metric) == 0 || anyNA(metric))
    stop('`metric` must be one or more alpha diversity metric names.', call. = FALSE)
  mets <- lapply(metric, match_metric, div = 'alpha')
  ids  <- vapply(mets, `[[`, '', 'id')
  mets[!duplicated(ids)]
}


# `depths`, checked: NULL, or whole numbers of at least 1, sorted, without
# duplicates, as integers.
rc_check_depths <- function (depths) {
  if (is.null(depths)) return (NULL)
  if (!is.numeric(depths) || !is.null(dim(depths)) || length(depths) == 0 ||
      anyNA(depths) || any(depths < 1) || any(depths %% 1 != 0) ||
      any(depths > .Machine$integer.max))
    stop('`depths` must be NULL or whole numbers of at least 1.', call. = FALSE)
  sort(unique(as.integer(depths)))
}


# A single whole number `v` of at least `lo`, as an integer.
rc_check_whole <- function (v, arg, lo, must) {
  if (!is.numeric(v) || length(v) != 1 || !is.null(dim(v)) || is.na(v) ||
      v < lo || v %% 1 != 0 || v > .Machine$integer.max)
    stop('`', arg, '` must be ', must, '.', call. = FALSE)
  as.integer(v)
}


# The metrics with a closed form.
RC_CLOSED_FORMS <- 'observed'


# `analytic` (NULL, TRUE or FALSE), as the IDs among `ids` to compute in
# closed form.
rc_check_analytic <- function (analytic, ids) {
  if (is.null(analytic)) return (intersect(ids, RC_CLOSED_FORMS))
  if (identical(analytic, FALSE)) return (character(0))
  if (!identical(analytic, TRUE))
    stop('`analytic` must be NULL, TRUE or FALSE.', call. = FALSE)
  no <- setdiff(ids, RC_CLOSED_FORMS)
  if (length(no))
    stop(
      '`analytic = TRUE` needs a closed form, which only ',
      paste0('"', RC_CLOSED_FORMS, '"', collapse = ', '), ' has; not ',
      paste0('"', no, '"', collapse = ', '), '. ',
      'Use `analytic = NULL` to resample the other metrics.', call. = FALSE )
  ids
}


# The depths of the curves: `depths`, or else `n_depths` points from 1 to
# `max_depth`, as np.linspace(1, max_depth, num = n_depths, dtype = int)
# (each i * step + 1, the last exactly max_depth, then rounded down), with
# duplicates dropped. Only depths up to `max_depth` (default: the largest of
# `depths`, else the largest library size) are kept. The depth cap is the
# attribute `max_depth`.
rc_grid <- function (depths, n_depths, max_depth, totals) {

  if (is.null(max_depth))
    max_depth <- if (is.null(depths)) max(c(0, totals)) else max(depths)

  if (is.null(depths)) {
    if (max_depth < 1) {
      grid <- integer(0)
    }
    else {
      step <- (max_depth - 1) / (n_depths - 1)
      g    <- (seq_len(n_depths) - 1) * step + 1
      g[n_depths] <- max_depth
      grid <- unique(as.integer(floor(g)))
    }
  }
  else {
    grid <- depths[depths <= max_depth]
  }

  attr(grid, 'max_depth') <- max_depth
  grid
}




#' Plot rarefaction curves
#'
#' Draws the curves from [rarefaction_curve()]: one line per sample, of a
#' metric against depth (or coverage), or of coverage against depth.
#'
#' @param x  The result of [rarefaction_curve()].
#'
#' @param metric  Which of `x`'s metrics to draw, by ID or name (matched as
#'        [alpha_div()] matches them). Default: the first in `x`.
#'
#' @param groups  Colours the samples by group, with a legend, and draws
#'        each group's median curve as a thicker line: the name of a column
#'        of `data`, or a vector with one value per sample (all the samples
#'        in the `totals` attribute, matched by name when both are named, as
#'        for [permanova()]).
#'
#' @param data  A data frame of sample metadata, with one row per sample
#'        (matched by row name, else in order), when `groups` names one of
#'        its columns.
#'
#' @param type  What to draw: `"depth"` (the default), the metric against
#'        depth; `"coverage"`, the metric against expected coverage
#'        (coverage-based rarefaction curves); or `"completeness"`, expected
#'        coverage against depth (sample completeness curves, Chao et al.
#'        2014, Fig. 4). A prefix is enough.
#'
#' @param depth  A depth to mark with a vertical line (for `type = "depth"`
#'        or `"completeness"`), such as one you plan to rarefy to. It is
#'        labelled with how many of the samples have at least that many
#'        observations. Default: `NULL` (none).
#'
#' @param legend  Where to put the legend of `groups`: a keyword for
#'        [graphics::legend()], such as `"bottomright"` (the default), or
#'        `NULL` for none.
#'
#' @param col  The lines' colours: one colour for every sample, or with
#'        `groups`, one per group (in order, or named by group). Defaults to
#'        grey without groups, and to [grDevices::hcl.colors()]'s `"Dark 3"`
#'        palette with them.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `ylab`, `xlim`, `ylim` or `log`. They replace the
#'        defaults made from `x`.
#'
#' @details
#' The title is "Rarefaction Curves of" and the metric's name (see
#' [list_metrics()]), or "Sample Completeness Curves"; the subtitle says
#' how the values were made: "Expected values" for a closed form, or, for
#' instance, "Mean of 10 rarefactions". A group's median curve is the
#' median of its samples' values at each depth that at least two of them
#' share (or every depth, for a group of one); past a sample's library size,
#' the sample has no value, so the median is over fewer samples (as QIIME
#' 2's group medians).
#'
#' @return `x`, invisibly.
#'
#' @references
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
#' R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
#' numbers: a framework for sampling and estimation in species diversity
#' studies. *Ecological Monographs*, 84, 45--67. \doi{10.1890/13-0133.1}
#'
#' @export
#' @examples
#'     rc <- rarefaction_curve(ex_counts, c('observed', 'shannon'), times = 5)
#'
#'     plot(rc, depth = 1000)
#'     plot(rc, 'shannon', groups = c('A', 'A', 'B', 'C'))
#'     plot(rc, type = 'coverage')
#'     plot(rc, type = 'completeness', main = 'How complete is each sample?')
plot.ecodive_rarefaction_curve <- function (
    x, metric = NULL, groups = NULL, data = NULL, type = 'depth',
    depth = NULL, legend = 'bottomright', col = NULL, ... ) {

  pos  <- plot_legend_pos(legend)
  type <- rc_plot_type(type)

  have <- unique(x$metric)
  if (length(have) == 0)
    stop('`x` has no rows: there is nothing to draw.', call. = FALSE)
  if (is.null(metric)) {
    id <- have[[1]]
  }
  else {
    id <- match_metric(metric, div = 'alpha')$id
    if (!id %in% have)
      stop(
        '`metric = "', metric, '"` is not in `x`, which has ',
        paste0('"', have, '"', collapse = ', '), '.', call. = FALSE )
  }

  if (!is.null(depth)) {
    if (type == 'coverage')
      stop(
        '`depth` marks a depth on the x axis, which is coverage here; ',
        'use type = "depth" or "completeness".', call. = FALSE )
    if (!is.numeric(depth) || length(depth) != 1 || is.na(depth) || depth < 0)
      stop('`depth` must be a single number of at least 0, or NULL.', call. = FALSE)
  }

  totals <- attr(x, 'totals', exact = TRUE)
  labs   <- names(totals)
  n      <- length(totals)
  unnamed <- identical(labs, as.character(seq_len(n)))

  df  <- x[x$metric == id, , drop = FALSE]
  sam <- match(df$sample, labs)
  xv  <- if (type == 'coverage') df$coverage else df$depth
  yv  <- if (type == 'completeness') df$coverage else df$value

  lab <- rc_plot_labels(x, id, type)

  g    <- plot_groups(n, if (unnamed) NULL else labs, groups, data)
  if (!is.null(g$values)) g <- list(f = factor(g$values), values = NULL, title = g$title)
  look <- ord_plot_look(g, col, NULL)
  lcol <- if (is.null(g$f)) {
            if (is.null(col)) rep(PLOT_MUTED, n) else rep_len(col, n)
          } else look$col

  args <- list(...)
  plot_do(
    graphics::plot.default, c(list(NA, NA), args),
    list(
      type = 'n',
      xlim = range(xv, finite = TRUE),
      ylim = range(yv, finite = TRUE),
      xlab = lab$xlab, ylab = lab$ylab, main = lab$main, sub = lab$sub ))

  lwd <- plot_arg(args, 'lwd', 1)
  lty <- plot_arg(args, 'lty', 1)
  for (s in unique(sam)) {
    i <- which(sam == s)
    graphics::lines(xv[i], yv[i], col = lcol[[s]], lwd = lwd, lty = lty)
  }

  # Each group's median curve.
  if (!is.null(g$f)) {
    f <- g$f[sam]
    for (lv in seq_along(levels(g$f))) {
      i <- which(as.integer(f) == lv)
      if (length(i) == 0) next
      key  <- if (type == 'coverage') df$depth[i] else xv[i]
      tab  <- table(key)
      ns   <- length(unique(sam[i]))
      keep <- if (ns > 1) as.numeric(names(tab)[tab >= 2]) else as.numeric(names(tab))
      keep <- sort(keep)
      if (length(keep) == 0) next
      mx <- vapply(keep, function (k) stats::median(xv[i][key == k]), 0)
      my <- vapply(keep, function (k) stats::median(yv[i][key == k]), 0)
      graphics::lines(mx, my, col = look$pal[[lv]], lwd = 3 * lwd)
    }
    plot_legend(pos, look$legend, look$legend_col, NA, g$title, lty = 1)
  }

  # A planned depth, and the samples it keeps.
  if (!is.null(depth)) {
    graphics::abline(v = depth, lty = 2, col = PLOT_HIGHLIGHT)
    usr <- graphics::par('usr')
    plot_label(
      depth, usr[4] - 0.03 * diff(usr[3:4]),
      sprintf('%d of %d samples', sum(totals >= depth), n),
      col = PLOT_HIGHLIGHT, pos = 4, cex = 0.8 )
  }

  invisible(x)
}


# The default title, subtitle and axis labels of plot(x, metric = id, type =
# type): from the metric's display name (METRICS) and x's attributes.
rc_plot_labels <- function (x, id, type) {
  name <- METRICS$name[match(id, METRICS$id)]
  if (is.na(name)) name <- id
  analytic <- attr(x, 'analytic', exact = TRUE)
  times    <- attr(x, 'times',    exact = TRUE)
  average  <- attr(x, 'average',  exact = TRUE)
  list(
    main = if (type == 'completeness') 'Sample Completeness Curves'
           else paste('Rarefaction Curves of', name),
    sub  = if (type == 'completeness' || id %in% analytic) 'Expected values'
           else sprintf(
             '%s of %s rarefaction%s', plot_cap(average),
             format(times, scientific = FALSE), if (identical(as.integer(times), 1L)) '' else 's' ),
    xlab = if (type == 'coverage') 'Sample Coverage' else 'Depth',
    ylab = if (type == 'completeness') 'Sample Coverage' else name )
}


# `type`, matched (exactly, else by a unique prefix, ignoring case).
rc_plot_type <- function (type) {
  choices <- c('depth', 'coverage', 'completeness')
  if (is.character(type) && length(type) == 1 && !is.na(type) && nzchar(type)) {
    needle <- tolower(type)
    hit    <- if (needle %in% choices) needle else choices[startsWith(choices, needle)]
    if (length(hit) == 1) return (hit)
  }
  stop('`type` must be one of ', paste0('"', choices, '"', collapse = ', '), '.', call. = FALSE)
}
