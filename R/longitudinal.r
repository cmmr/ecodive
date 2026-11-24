# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# first_differences(), first_distances(), volatility() and nmit(): ports of
# q2-longitudinal 2026.10.0's first-differences, first-distances,
# volatility and nmit actions (q2_longitudinal/_utilities.py:
# _first_differences(), _get_group_pairs(), _nmit(), _temporal_corr(),
# _temporal_distance(); q2_longitudinal/_vega_specs/volatility/data.py for
# the volatility chart's statistics), restructured for ecodive (see
# plans/longitudinal.md). The chart's confidence intervals follow Vega's
# bootstrapCI() (vega-statistics), as src/longitudinal.c says, where the
# correlations of nmit() are computed.
#
# q2-longitudinal's notice:
#
#   Copyright (c) 2017-2026, QIIME 2 development team.
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   3. Neither the name of the copyright holder nor the names of its
#      contributors may be used to endorse or promote products derived from
#      this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' First Differences Within Subjects
#'
#' For repeated measures: how much a value changes within each subject from
#' one time point (state) to the next, or from a baseline. These paired
#' changes remove each subject's own level, so they show change over time
#' that pooled comparisons of time points can't (Bokulich et al. 2018).
#' The same as q2-longitudinal's `first-differences` action.
#'
#' @inheritParams documentation
#'
#' @param x  What to difference: a numeric vector with one value per sample
#'        (such as an [alpha_div()] result); the name of a numeric column of
#'        `data` (q2's `metric`); or a count table (any of the types listed
#'        under **Input Types**), whose features are all differenced at
#'        once, after `norm`.
#'
#' @param subjects  The subject (individual) of each sample: one value per
#'        sample, matched to the samples by name when both are named, else
#'        in the samples' order; or, with `data`, the name of one of its
#'        columns (q2's `individual_id_column`). No `NA`s.
#'
#' @param states  The time point of each sample, as a number (days, weeks,
#'        or visit numbers): given as `subjects` is (q2's `state_column`).
#'        At least two distinct values; no `NA`s.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding the columns named in `x`, `subjects` and
#'        `states`. Rows are matched to the samples by row name (the names
#'        of `x`, or the sample names of a count table or `dm`), and rows
#'        for other samples are ignored; without row names, they must be in
#'        the samples' order. When `x` names one of its columns, its rows
#'        are the samples. Default: `NULL`.
#'
#' @param baseline  `NULL` (the default) for first differences, between
#'        consecutive states. Or one of the states: then each state is
#'        compared with it instead (static differences, such as changes
#'        from before a treatment).
#'
#' @param replicates  What to do when a subject has more than one sample at
#'        a state (q2's `replicate_handling`): `'error'` (the default),
#'        `'drop'` (that subject is left out of the pairs that involve the
#'        state) or `'random'` (one of them is picked at random, with
#'        `seed`).
#'
#' @param norm  For a count table, how to normalize it first: `'none'` (the
#'        default; the values as given), `'percent'` (relative abundances,
#'        q2's relative frequencies), `'clr'`, or any other `norm` of
#'        [normalize()]. Ignored otherwise.
#'
#' @param seed  An integer seed for `replicates = 'random'`. Default: `0`.
#'
#' @return A data frame of class `c("ecodive_first_differences",
#'   "data.frame")`, one row per pair of samples, ordered by pair of states
#'   and then by subject (in the order of their factor levels, or sorted),
#'   with columns:
#'   \describe{
#'     \item{`sample`}{The later sample of the pair (q2's index): its name,
#'       or its position when the samples have no names.}
#'     \item{`subject`}{The subject.}
#'     \item{`state_from`}{The earlier state (the baseline, with one).}
#'     \item{`state`}{The later state.}
#'     \item{`difference`}{The value at `state` minus the value at
#'       `state_from`. For a count table, one column per feature instead,
#'       named by feature.}
#'   }
#'   Pairs whose difference is `NA` are left out, as q2 does. With
#'   `baseline`, attribute `baseline` holds it. `plot()` draws the
#'   differences against the state.
#'
#' @details
#' **Pairs of samples**: the states are the distinct values of `states`
#' over all the samples, sorted. Without `baseline`, each pair of
#' consecutive states \eqn{(t_k, t_{k+1})} pairs, for every subject with a
#' sample at both, those two samples; with `baseline = b`, each pair is
#' \eqn{(b, t)}. A subject that skipped a state contributes nothing on
#' either side of it: one sampled at \eqn{t_1} and \eqn{t_3} but not
#' \eqn{t_2}, when another subject was sampled at \eqn{t_2}, gives no pair
#' at all. This is q2-longitudinal's rule, kept for parity. To compare each
#' subject's consecutive samples whatever their times, give `states` as
#' each subject's visit number (1, 2, 3, ...) instead.
#'
#' To test the differences, use [stats::wilcox.test()] or
#' [stats::t.test()] on one pair of states (q2's `pairwise-differences`
#' action: `first_differences(x, ..., baseline = s1)` and the rows whose
#' `state` is `s2`), or a linear mixed model with a random intercept per
#' subject, as q2-longitudinal's tutorial does.
#'
#' **Differences from q2-longitudinal**:
#' \itemize{
#'   \item q2 returns the pairs in the arbitrary order of a Python set of
#'     subjects; ecodive sorts them, and returns the subjects and states
#'     alongside the differences.
#'   \item q2 takes one feature of a table at a time (its `metric`);
#'     ecodive differences them all, after `norm`.
#'   \item With `replicates = 'random'`, q2 draws with Python's unseeded
#'     `random.choice()` for each pair of states, and so does ecodive, but
#'     with its own seeded generator (PCG, one stream per pair of states,
#'     subject and side of the pair): the picks differ from q2's, as from
#'     one q2 run to the next.
#'   \item q2 reads the states from every row of the metadata; ecodive from
#'     the samples of `x` alone.
#'   \item `NA` in `subjects` or `states` is an error, as `groups` is
#'     elsewhere in ecodive; q2 has no rule for it.
#' }
#'
#' @references
#' Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
#' Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
#' longitudinal and paired-sample analyses of microbiome data. *mSystems*,
#' 3, e00219-18. \doi{10.1128/mSystems.00219-18}
#'
#' @seealso [first_distances()], [volatility()], [nmit()].
#'
#' @family longitudinal
#' @export
#' @examples
#'     # Four subjects, each sampled on days 0, 7 and 14.
#'     meta <- data.frame(
#'       subject = rep(c('a', 'b', 'c', 'd'), each = 3),
#'       day     = rep(c(0, 7, 14), 4),
#'       row.names = sprintf('s%02d', 1:12) )
#'     shannon <- setNames(c(2.1, 2.4, 2.2, 1.8, 2.0, 2.3, 2.5, 2.2, 2.6, 1.9, 2.1, 2.0), rownames(meta))
#'
#'     first_differences(shannon, 'subject', 'day', meta)
#'
#'     # Changes from day 0.
#'     fd <- first_differences(shannon, 'subject', 'day', meta, baseline = 0)
#'     fd
#'     plot(fd)
#'
#'     # Every feature of a count table, as relative abundances.
#'     counts <- matrix(rpois(12 * 5, 20), 12, dimnames = list(rownames(meta), letters[1:5]))
#'     first_differences(counts, 'subject', 'day', meta, norm = 'percent')
#'
first_differences <- function (
    x, subjects, states, data = NULL, baseline = NULL, replicates = 'error',
    margin = 1L, norm = 'none', pseudocount = NULL, seed = 0L ) {

  validate_seed()
  replicates <- long_choice(replicates, 'replicates', c('error', 'drop', 'random'))

  inp <- long_x(x, data, margin, norm, pseudocount, table = TRUE)
  md  <- long_meta(inp, subjects, states, missing(subjects), missing(states))
  pr  <- long_pairs(md$subjects, md$states, baseline, replicates, seed)

  vals <- inp$values[pr$to, , drop = FALSE] - inp$values[pr$from, , drop = FALSE]
  keep <- if (inp$table) rep(TRUE, nrow(vals)) else !is.na(vals[, 1])
  if (!any(keep)) long_empty()

  res <- long_frame(inp$labels, md$subjects, md$states, pr, keep)
  if (inp$table) {
    feat <- inp$features
    for (j in seq_along(feat)) res[[feat[[j]]]] <- unname(vals[keep, j])
  }
  else {
    res$difference <- unname(vals[keep, 1])
  }

  if (!is.null(baseline)) attr(res, 'baseline') <- pr$baseline
  class(res) <- c('ecodive_first_differences', 'data.frame')
  res
}



#' First Distances Within Subjects
#'
#' For repeated measures: how far each subject's community moves from one
#' time point (state) to the next, or from a baseline: the distance between
#' the subject's two samples, from a distance matrix. The same as
#' q2-longitudinal's `first-distances` action. Typical uses are the rate of
#' change of each subject's community (consecutive states), and how far it
#' is from its baseline (with `baseline`, such as recovery after
#' antibiotics).
#'
#' @inheritParams first_differences
#'
#' @param dm  A distance matrix over samples: a `dist` object, such as from
#'        [beta_div()] or [bray()], or a symmetric numeric matrix. Its
#'        sample names are matched to `data`'s row names and to the names
#'        of `subjects` and `states`.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding the columns named in `subjects` and `states`.
#'        Rows are matched to the samples in `dm` by row name, and rows for
#'        other samples are ignored; without row names, they must be in the
#'        same order as `dm`. Default: `NULL`.
#'
#' @return A data frame of class `c("ecodive_first_distances",
#'   "data.frame")`, one row per pair of samples, ordered by pair of states
#'   and then by subject, with columns `sample` (the later sample), `subject`,
#'   `state_from`, `state` and `distance`. With `baseline`, attribute
#'   `baseline` holds it. Its title (see [ecodive_display]) is "First
#'   Distances of" and the title of `dm`. `plot()` draws the distances
#'   against the state.
#'
#' @details
#' The pairs of samples are those of [first_differences()] (see its
#' **Pairs of samples**). Distances are read from the packed triangle of
#' `dm`, without forming the full matrix.
#'
#' **Differences from q2-longitudinal**: those of [first_differences()],
#' and q2 reads the states from every row of its metadata, even samples
#' that are not in the distance matrix (whose pairs it then skips); ecodive
#' takes them from the samples of `dm` alone.
#'
#' @references
#' Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
#' Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
#' longitudinal and paired-sample analyses of microbiome data. *mSystems*,
#' 3, e00219-18. \doi{10.1128/mSystems.00219-18}
#'
#' @family longitudinal
#' @export
#' @examples
#'     meta <- data.frame(
#'       subject = rep(c('a', 'b', 'c', 'd'), each = 3),
#'       day     = rep(c(0, 7, 14), 4),
#'       row.names = sprintf('s%02d', 1:12) )
#'     counts <- matrix(rpois(12 * 5, 20), 12, dimnames = list(rownames(meta), letters[1:5]))
#'
#'     fd <- first_distances(bray(counts), 'subject', 'day', meta)
#'     fd
#'     plot(fd)
#'
#'     # Distance from day 0.
#'     first_distances(bray(counts), 'subject', 'day', meta, baseline = 0)
#'
first_distances <- function (
    dm, subjects, states, data = NULL, baseline = NULL, replicates = 'error',
    seed = 0L ) {

  validate_seed()
  replicates <- long_choice(replicates, 'replicates', c('error', 'drop', 'random'))

  dv <- dist_vector(dm)
  n  <- dv$n
  if (n < 2)
    stop('`dm` must have at least 2 samples.', call. = FALSE)

  inp <- c(long_data(data, dv$labels, n, '`dm`'), list(labels = dv$labels, n = n, ref = '`dm`'))
  md  <- long_meta(inp, subjects, states, missing(subjects), missing(states))
  pr  <- long_pairs(md$subjects, md$states, baseline, replicates, seed)

  # Each pair's place in the packed lower triangle (columns j < rows i).
  i   <- as.double(pmax(pr$from, pr$to))
  j   <- as.double(pmin(pr$from, pr$to))
  pos <- (j - 1) * n - j * (j - 1) / 2 + (i - j)

  res <- long_frame(dv$labels, md$subjects, md$states, pr, rep(TRUE, length(pos)))
  res$distance <- as.vector(dv$dm[pos])

  if (!is.null(baseline)) attr(res, 'baseline') <- pr$baseline
  class(res) <- c('ecodive_first_distances', 'data.frame')
  set_display(res, title = title_of('First Distances', dv$title))
}



#' Volatility of a Variable Over Time
#'
#' The numbers behind q2-longitudinal's volatility chart (a control chart of
#' a variable over time), and the chart itself with `plot()`: each group's
#' mean at each state with a bootstrap confidence interval, global control
#' limits (the mean \eqn{\pm} 2 and 3 standard deviations over all the
#' samples), and each subject's trajectory.
#'
#' @inheritParams first_differences
#'
#' @param x  The variable: a numeric vector with one value per sample (such
#'        as an [alpha_div()] result), or the name of a numeric column of
#'        `data`. `NA`s are left out of the means and limits. For the
#'        features of a count table, call `volatility()` for each.
#'
#' @param subjects  Optional: the subject of each sample, for the
#'        trajectories (q2's `individual_id_column`), given as `x` is. No
#'        `NA`s. `NULL` for none.
#'
#' @param groups  Optional: the group of each sample (q2's
#'        `default_group_column`), given as `x` is: a mean line per group.
#'        No `NA`s. Default: `NULL`, one group of every sample.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding the columns named in `x`, `subjects`, `states`
#'        and `groups`. When `x` names one of its columns, its rows are the
#'        samples; otherwise they are matched to the names of `x` by row
#'        name, or else taken in order. Default: `NULL`.
#'
#' @param conf_level  The confidence level of the intervals: a number
#'        strictly between 0 and 1. Default: `0.95`, as q2's chart.
#'
#' @param n_resamples  The number of bootstrap resamples for each interval.
#'        Default: `1000`, as q2's chart.
#'
#' @param seed  An integer seed for the bootstrap resamples. Default: `0`.
#'
#' @return A list of class `ecodive_volatility`:
#'   \describe{
#'     \item{`means`}{A data frame, one row per group and state (in the
#'       groups' order, then by state): `group`, `state`, `mean` (of the
#'       values that are not `NA`), `lower` and `upper` (the bootstrap
#'       interval of the mean) and `n` (the number of values).}
#'     \item{`global`}{A named numeric vector over all the samples' values:
#'       `mean`, `sd` (the sample standard deviation, with \eqn{n - 1}; 0
#'       for a single value, as Vega's), `lower_3sd`, `lower_2sd`,
#'       `upper_2sd` and `upper_3sd` (the mean \eqn{\pm} 2 and 3 SDs).}
#'     \item{`trajectories`}{With `subjects`, a data frame of each sample
#'       (`sample`, `subject`, `group`, `state`, `value`), ordered by subject
#'       and state; else `NULL`.}
#'   }
#'   Attributes `conf_level` and `n_resamples` describe the intervals.
#'
#' @details
#' **The interval** is Vega's `ci0` and `ci1` aggregates, which q2's chart
#' draws (Vega's `bootstrapCI()`): `n_resamples` times, the mean of as many
#' values as the group has, drawn from them with replacement; then the
#' \eqn{\alpha/2} and \eqn{1 - \alpha/2} quantiles of those means (R's
#' default [stats::quantile()] type 7), with \eqn{\alpha = 1 -}
#' `conf_level`. Vega draws with JavaScript's unseeded `Math.random()`, so
#' its intervals change each time the chart is drawn; ecodive's come from
#' its own seeded generator (PCG, one stream per group and state), so the
#' same `seed` gives the same intervals, which differ from Vega's as one
#' draw does from another. A group of one value has the interval of width 0.
#'
#' **Differences from q2-longitudinal**: its chart lets the viewer pick the
#' variable and the grouping interactively; here, one of each per call.
#' Its tooltips count every sample of a group at a state (Vega's `count`),
#' those with a missing value too; `n` counts the values used.
#'
#' **Same names in other packages**: `TTR::volatility()` (of financial
#' prices) is unrelated; whichever package was attached last masks the
#' other's, so write `ecodive::volatility()` to be sure.
#'
#' @references
#' Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
#' Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
#' longitudinal and paired-sample analyses of microbiome data. *mSystems*,
#' 3, e00219-18. \doi{10.1128/mSystems.00219-18}
#'
#' Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
#' Bootstrap*. Chapman & Hall.
#'
#' Satyanarayan, A., Russell, R., Hoffswell, J., & Heer, J. (2016). Reactive
#' Vega: a streaming dataflow architecture for declarative interactive
#' visualization. *IEEE Transactions on Visualization and Computer
#' Graphics*, 22, 659--668. \doi{10.1109/TVCG.2015.2467091}
#'
#' @family longitudinal
#' @export
#' @examples
#'     meta <- data.frame(
#'       subject = rep(c('a', 'b', 'c', 'd'), each = 3),
#'       day     = rep(c(0, 7, 14), 4),
#'       diet    = rep(c('x', 'y'), each = 6),
#'       shannon = c(2.1, 2.4, 2.2, 1.8, 2.0, 2.3, 2.5, 2.2, 2.6, 1.9, 2.1, 2.0) )
#'
#'     v <- volatility('shannon', 'subject', 'day', 'diet', data = meta)
#'     v
#'     plot(v)
#'
volatility <- function (
    x, subjects, states, groups = NULL, data = NULL, conf_level = 0.95,
    n_resamples = 1000L, seed = 0L ) {

  validate_seed()

  if (!is.numeric(conf_level) || length(conf_level) != 1 || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop('`conf_level` must be a number strictly between 0 and 1.', call. = FALSE)
  if (!is.numeric(n_resamples) || length(n_resamples) != 1 || !is.finite(n_resamples) ||
      n_resamples < 1 || n_resamples %% 1 != 0 || n_resamples > 1e7)
    stop('`n_resamples` must be a whole number from 1 to 1e7.', call. = FALSE)
  n_resamples <- as.integer(n_resamples)

  inp <- long_x(x, data, table = FALSE)
  if (missing(subjects)) subjects <- NULL
  md  <- long_meta(inp, subjects, states, FALSE, missing(states), subjects_ok = TRUE)
  n   <- inp$n
  y   <- inp$values[, 1]
  st  <- md$states

  g <- if (is.null(groups)) NULL else long_var(groups, 'groups', inp)
  grp <- if (is.null(g)) factor(rep('all', n)) else long_factor(g$val)


  #--------------------------------------------------------
  # Each group at each state: mean, n and interval.
  #--------------------------------------------------------
  st_u  <- sort(unique(st))
  cell  <- (as.integer(grp) - 1L) * length(st_u) + match(st, st_u)
  cells <- sort(unique(cell))
  ok    <- !is.na(y)

  # The values of each cell, in sample order.
  o    <- order(cell, seq_len(n), method = 'radix')
  o    <- o[ok[o]]
  vals <- as.double(y[o])
  cnt  <- tabulate(match(cell[o], cells), length(cells))
  pos  <- c(0L, cumsum(cnt))

  mean_v <- vapply(seq_along(cells), function (k) {
    if (cnt[[k]] == 0) NA_real_ else sum(vals[seq_len(cnt[[k]]) + pos[[k]]]) / cnt[[k]]
  }, numeric(1))

  ci <- .Call(C_volatility_ci, vals, as.integer(pos), c(n_resamples, seed), 1 - conf_level)

  gi <- (cells - 1L) %/% length(st_u) + 1L
  si <- (cells - 1L) %%  length(st_u) + 1L

  means <- data.frame(
    group = factor(levels(grp)[gi], levels = levels(grp)),
    state = st_u[si],
    mean  = mean_v,
    lower = ci[1, ],
    upper = ci[2, ],
    n     = cnt )


  #--------------------------------------------------------
  # Global control limits (Vega's mean and stdev).
  #--------------------------------------------------------
  yy <- as.double(y[ok])
  gm <- if (length(yy)) sum(yy) / length(yy) else NA_real_
  gs <- if (length(yy) > 1) sqrt(sum((yy - gm)^2) / (length(yy) - 1)) else if (length(yy)) 0 else NA_real_
  global <- c(
    mean = gm, sd = gs,
    lower_3sd = gm - 3 * gs, lower_2sd = gm - 2 * gs,
    upper_2sd = gm + 2 * gs, upper_3sd = gm + 3 * gs )


  #--------------------------------------------------------
  # Each subject's trajectory.
  #--------------------------------------------------------
  traj <- NULL
  if (!is.null(md$subjects)) {
    lev <- long_levels(md$subjects)
    sc  <- match(as.character(md$subjects), lev)
    o   <- order(sc, st, seq_len(n), method = 'radix')
    traj <- data.frame(
      sample  = if (is.null(inp$labels)) o else inp$labels[o],
      subject = md$subjects[o],
      group   = grp[o],
      state   = st[o],
      value   = as.vector(y[o]),
      stringsAsFactors = FALSE )
    rownames(traj) <- NULL
  }

  res <- structure(
    list(means = means, global = global, trajectories = traj),
    class = 'ecodive_volatility',
    conf_level  = conf_level,
    n_resamples = n_resamples )

  set_display(
    res, title = 'Volatility', ylab = inp$name, xlab = md$states_col,
    legend = if (!is.null(g)) g$col )
}



#' Nonparametric Microbial Interdependence Test (NMIT) Distances
#'
#' Distances between subjects in how their features vary together over
#' time (Zhang et al. 2017), as q2-longitudinal's `nmit` action computes
#' them: for each subject, the correlation of every pair of features over
#' that subject's samples (its time points); then, for each pair of
#' subjects, a matrix norm of the difference of their correlation
#' matrices. NMIT's test is a PERMANOVA of these distances against a
#' grouping of the subjects: run [permanova()] on the result.
#'
#' @inheritParams documentation
#'
#' @param counts  A count table: samples \eqn{\times} features, or any of
#'        the types listed under **Input Types**. Every sample of a subject
#'        enters its correlations, in any order.
#'
#' @param subjects  The subject (individual) of each sample: one value per
#'        sample, matched to the samples by name when both are named, else
#'        in the samples' order; or, with `data`, the name of one of its
#'        columns (q2's `individual_id_column`). No `NA`s.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding the `subjects` column. Rows are matched to the
#'        samples in `counts` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the samples' order.
#'        Default: `NULL`.
#'
#' @param method  The correlation over time: `'kendall'` (Kendall's
#'        \eqn{\tau_b}, the default, as q2's), `'pearson'` or `'spearman'`.
#'
#' @param matrix_norm  The norm of the difference of two subjects'
#'        correlation matrices: `'frobenius'` (the default, as q2's) or
#'        `'nuclear'` (the sum of its singular values).
#'
#' @param norm  How to normalize `counts` first: `'percent'` (relative
#'        abundances, the default, as q2 expects), `'none'`, `'clr'`, or any
#'        other `norm` of [normalize()].
#'
#' @param clust  Whether and how to order the subjects by clustering, as
#'        for [beta_div()]: `NULL` (the default) or `FALSE` to keep their
#'        order, or one of `stats::hclust()`'s linkage methods (such as
#'        `"average"`), or `"pcoa"`. The `order` attribute then gives each
#'        subject's position before clustering.
#'
#' @return A `dist` object over the **subjects** (class
#'   `c("ecodive_dist", "dist")`), labelled by subject and in their order
#'   (a factor's levels, or sorted), with `method` `"nmit"` and title "NMIT
#'   Distance". It can be given as `dm` to [permanova()], [pcoa()] or
#'   [stats::hclust()].
#'
#' @details
#' **Correlations** are computed as pandas' `DataFrame.corr()` (which q2
#' calls) computes them. A pair involving a feature that doesn't vary within
#' the subject has no correlation, and counts as 0 (q2's `fillna(0)`).
#' On the diagonal, pandas puts 1 for Kendall's \eqn{\tau_b} whatever the
#' feature, but for Pearson and Spearman a constant feature's own
#' correlation is also undefined, so 0: subjects differ there too, and that
#' adds to their distance.
#'
#' **Subjects with fewer than 3 samples** have correlations of \eqn{\pm 1}
#' or 0 only (one sample: none at all), which say little about how their
#' features vary together. q2 computes them all the same; so does
#' `nmit()`, with a warning naming those subjects. Leave them out first
#' if you can.
#'
#' **Distances**: the Frobenius norm of the difference \eqn{\Delta} is
#' \eqn{\sqrt{2\sum_{a<b} \Delta_{ab}^2 + \sum_a \Delta_{aa}^2}}, computed
#' from the stored lower triangles. The nuclear norm is the sum of the
#' absolute eigenvalues of \eqn{\Delta} (symmetric, so these are its
#' singular values): an eigendecomposition of a \eqn{D \times D} matrix for
#' each pair of subjects, \eqn{O(S^2 D^3)} in all for \eqn{S} subjects and
#' \eqn{D} features, which for hundreds of features takes minutes.
#'
#' **Memory**: every subject's correlations are kept, \eqn{S D (D + 1) / 2}
#' numbers (50 subjects and 500 features: 50 MB). `nmit()` warns above
#' 500 MB; as q2's documentation advises, filter out rare features first.
#'
#' **The test**: give [permanova()] the distances and one group per
#' subject, matched by subject name, for instance from the sample metadata:
#' \preformatted{
#'   d    <- nmit(counts, 'subject', meta)
#'   subj <- unique(meta[, c('subject', 'group')])
#'   rownames(subj) <- subj$subject
#'   permanova(d, 'group', subj)
#' }
#' Its p-values equal those of other PERMANOVAs of the same distances, such
#' as QIIME 2's `beta-group-significance`, given the same permutations.
#'
#' **Differences from q2-longitudinal**: q2 labels each subject by the ID
#' of its first sample, in the order the subjects first appear; ecodive
#' labels them by subject. Distances are the same to rounding (about
#' \eqn{10^{-15}}).
#'
#' @references
#' Zhang, Y., Han, S. W., Cox, L. M., & Li, H. (2017). A multivariate
#' distance-based analytic framework for microbial interdependence
#' association test in longitudinal study. *Genetic Epidemiology*, 41,
#' 769--778. \doi{10.1002/gepi.22065}
#'
#' Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
#' Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
#' longitudinal and paired-sample analyses of microbiome data. *mSystems*,
#' 3, e00219-18. \doi{10.1128/mSystems.00219-18}
#'
#' @inheritSection documentation Input Types
#'
#' @family longitudinal
#' @export
#' @examples
#'     # Six subjects, each sampled 5 times; two diets.
#'     meta <- data.frame(
#'       subject = rep(sprintf('m%d', 1:6), each = 5),
#'       diet    = rep(c('x', 'y'), each = 15) )
#'     counts <- matrix(rpois(30 * 6, 30), 30, dimnames = list(NULL, letters[1:6]))
#'
#'     d <- nmit(counts, 'subject', meta)
#'     d
#'
#'     # NMIT's test: a PERMANOVA of the subjects.
#'     subj <- unique(meta)
#'     rownames(subj) <- subj$subject
#'     permanova(d, 'diet', subj, permutations = 99)
#'
nmit <- function (
    counts, subjects, data = NULL, method = 'kendall',
    matrix_norm = 'frobenius', margin = 1L, norm = 'percent',
    pseudocount = NULL, clust = NULL, cpus = n_cpus() ) {

  method      <- long_choice(method, 'method', c('kendall', 'pearson', 'spearman'))
  matrix_norm <- long_choice(matrix_norm, 'matrix_norm', c('frobenius', 'nuclear'))

  validate_counts()
  validate_margin()
  validate_norm()
  validate_pseudocount()
  validate_clust()
  validate_cpus()

  if (!all(is.finite(stored_values(counts))))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)

  smp <- count_samples(counts, margin)
  inp <- c(long_data(data, smp$labels, smp$n, '`counts`'),
           list(labels = smp$labels, n = smp$n, ref = '`counts`'))

  if (missing(subjects)) subjects <- NULL
  sv  <- long_var(subjects, 'subjects', inp)
  lev <- long_levels(sv$val)
  sc  <- match(as.character(sv$val), lev)
  S   <- length(lev)
  if (S < 2)
    stop('`subjects` must have at least 2 subjects.', call. = FALSE)

  dims <- switch(
    get_matrix_package(counts),
    base = dim(counts), slam = c(counts$nrow, counts$ncol), Matrix = counts@Dim )
  D <- as.double(dims[[3L - margin]])
  if (D < 1)
    stop('`counts` has no features.', call. = FALSE)

  ns  <- tabulate(sc, S)
  few <- which(ns < 3)
  if (length(few))
    warning(
      sprintf('%d subject%s %s fewer than 3 samples, so %s correlations are only +/-1 or 0: %s.',
        length(few), if (length(few) == 1) '' else 's', if (length(few) == 1) 'has' else 'have',
        if (length(few) == 1) 'its' else 'their',
        paste(sprintf("'%s' (%d)", lev[few], ns[few]), collapse = ', ')),
      call. = FALSE )

  mb <- as.double(S) * D * (D + 1) / 2 * 8 / 1e6
  if (mb > 500)
    warning(sprintf(
      'nmit() keeps every subject\'s correlations: %.0f MB for %d subjects and %d features. Consider filtering out rare features first.',
      mb, S, as.integer(D) ), call. = FALSE)

  res <- .Call(
    C_nmit, counts, margin, norm, as.double(pseudocount), sc - 1L, S,
    match(method, c('pearson', 'spearman', 'kendall')),
    match(matrix_norm, c('frobenius', 'nuclear')), cpus )

  attributes(res) <- list(
    Size = S, Labels = lev, Diag = FALSE, Upper = FALSE,
    class = c('ecodive_dist', 'dist'), method = 'nmit' )
  res <- set_display(res, title = 'NMIT Distance')
  dist_clust(res, clust, cpus)
}



#------------------------------------------------------------------------------
# Plots
#------------------------------------------------------------------------------

#' Plot Longitudinal Results
#'
#' Charts of [first_differences()], [first_distances()] and
#' [volatility()] results.
#'
#' `plot()` of a [first_differences()] or [first_distances()] result draws
#' the values against the (later) state: a box of them at each state, and
#' each subject's values joined by a faint line; a dashed line at 0 for
#' differences.
#'
#' `plot()` of a [volatility()] result draws q2-longitudinal's control
#' chart: the state on the x axis; each group's means, joined by a line,
#' with their confidence intervals as error bars, one colour per group; each
#' subject's trajectory as a thin, faint line in its group's colour (the
#' "spaghetti"); the global mean as a dashed line and the \eqn{\pm} 2 and 3
#' SD control limits as dotted lines.
#'
#' @param x  A result of [first_differences()], [first_distances()] or
#'        [volatility()].
#'
#' @param feature  For a [first_differences()] result of a count table: the
#'        feature to draw (a column name). Default: the first.
#'
#' @param spaghetti  Draw each subject's trajectory. Default: `TRUE`.
#'
#' @param limits  Draw the global mean and control limits. Default: `TRUE`.
#'
#' @param legend  Where to put the legend of groups (a keyword of
#'        [graphics::legend()], such as `"topright"`), or `NULL` for none.
#'        Only drawn for two or more groups.
#'
#' @param col  The groups' colours. Default: [grDevices::hcl.colors()]'s
#'        "Dark 3" palette.
#'
#' @param ...  Passed to [graphics::plot.default()], such as `main`, `xlab`,
#'        `ylab` or `ylim`.
#'
#' @return `x`, invisibly.
#'
#' @name plot_longitudinal
#' @examples
#'     meta <- data.frame(
#'       subject = rep(c('a', 'b', 'c', 'd'), each = 3),
#'       day     = rep(c(0, 7, 14), 4),
#'       diet    = rep(c('x', 'y'), each = 6),
#'       shannon = c(2.1, 2.4, 2.2, 1.8, 2.0, 2.3, 2.5, 2.2, 2.6, 1.9, 2.1, 2.0) )
#'
#'     plot(volatility('shannon', 'subject', 'day', 'diet', data = meta))
#'     plot(first_differences('shannon', 'subject', 'day', meta))
NULL


#' @rdname plot_longitudinal
#' @export
plot.ecodive_first_differences <- function (x, feature = NULL, ...) {

  keys <- c('sample', 'subject', 'state_from', 'state')
  vcol <- setdiff(names(x), keys)
  if (length(vcol) == 0)
    stop('`x` has no values to draw.', call. = FALSE)
  if (is.null(feature)) {
    col <- vcol[[1]]
  }
  else {
    if (!is.character(feature) || length(feature) != 1 || !feature %in% vcol)
      stop('`feature` must be the name of one of the value columns of `x`: ',
           name_list(vcol), '.', call. = FALSE)
    col <- feature
  }

  b <- attr(x, 'baseline', exact = TRUE)
  long_plot_pairs(
    x, x[[col]], zero = TRUE,
    main = 'First Differences',
    sub  = if (!is.null(b)) paste('Differences from state', format(b)) else 'Consecutive states',
    ylab = if (col == 'difference') 'Difference' else paste('Difference in', col),
    ... )
}


#' @rdname plot_longitudinal
#' @export
plot.ecodive_first_distances <- function (x, ...) {
  b <- attr(x, 'baseline', exact = TRUE)
  long_plot_pairs(
    x, x$distance, zero = FALSE,
    main = plot_title(x, 'First Distances'),
    sub  = if (!is.null(b)) paste('Distances from state', format(b)) else 'Consecutive states',
    ylab = 'Distance', ... )
}


# Values v of a first_differences() or first_distances() result x against
# its states: boxes per state and a faint line per subject.
long_plot_pairs <- function (x, v, zero, main, sub, ylab, ...) {

  if (nrow(x) == 0)
    stop('`x` has no rows: there is nothing to draw.', call. = FALSE)

  st   <- x$state
  st_u <- sort(unique(st))
  w    <- if (length(st_u) > 1) 0.5 * min(diff(st_u)) else 1
  ok   <- is.finite(v)
  yl   <- range(c(v[ok], if (zero) 0))
  if (!all(is.finite(yl))) yl <- c(-1, 1)

  args <- list(...)
  plot_do(
    graphics::plot.default, c(list(NA, NA), args),
    list(
      type = 'n', xlim = range(st_u) + c(-w, w), ylim = yl,
      xlab = 'State', ylab = ylab, main = main, sub = sub ))

  if (zero) graphics::abline(h = 0, lty = 2, col = PLOT_MUTED)

  line_col <- grDevices::adjustcolor(PLOT_MUTED, alpha.f = 0.4)
  subj     <- as.character(x$subject)
  for (s in unique(subj)) {
    i <- which(subj == s & ok)
    i <- i[order(st[i], method = 'radix')]
    if (length(i) > 1) graphics::lines(st[i], v[i], col = line_col)
  }

  f <- factor(st[ok], levels = st_u)
  graphics::boxplot(
    split(v[ok], f), at = st_u, boxwex = w, add = TRUE, axes = FALSE,
    outline = TRUE, col = NA, border = 'black', pch = 20, cex = 0.6 )

  invisible(x)
}


#' @rdname plot_longitudinal
#' @export
plot.ecodive_volatility <- function (
    x, spaghetti = TRUE, limits = TRUE, legend = 'topright', col = NULL, ... ) {

  pos  <- plot_legend_pos(legend)
  m    <- x$means
  gl   <- x$global
  tr   <- if (isTRUE(spaghetti)) x$trajectories
  lev  <- levels(m$group)
  k    <- length(lev)
  pal  <- if (is.null(col)) plot_palette(k) else rep_len(col, k)

  ys <- c(m$mean, m$lower, m$upper, if (!is.null(tr)) tr$value,
          if (isTRUE(limits)) gl[c('lower_3sd', 'upper_3sd')])
  yl <- range(ys[is.finite(ys)])
  if (!all(is.finite(yl))) yl <- c(-1, 1)
  xs <- c(m$state, if (!is.null(tr)) tr$state)

  xlab <- get_display(x, 'xlab')
  ylab <- get_display(x, 'ylab')

  args <- list(...)
  plot_do(
    graphics::plot.default, c(list(NA, NA), args),
    list(
      type = 'n', xlim = range(xs), ylim = yl,
      xlab = if (is.null(xlab)) 'State' else xlab,
      ylab = if (is.null(ylab)) 'Value' else ylab,
      main = plot_title(x, 'Volatility') ))

  if (isTRUE(limits) && is.finite(gl[['mean']])) {
    graphics::abline(h = gl[['mean']], lty = 2, col = PLOT_MUTED)
    graphics::abline(
      h = gl[c('lower_3sd', 'lower_2sd', 'upper_2sd', 'upper_3sd')],
      lty = 3, col = PLOT_MUTED )
  }

  if (!is.null(tr)) {
    gi  <- match(as.character(tr$group), lev)
    key <- paste(as.character(tr$subject), gi, sep = '\r')
    for (s in unique(key)) {
      i <- which(key == s & is.finite(tr$value))
      if (length(i) > 1)
        graphics::lines(
          tr$state[i], tr$value[i],
          col = grDevices::adjustcolor(pal[[gi[i[[1]]]]], alpha.f = 0.3) )
    }
  }

  for (g in seq_len(k)) {
    i <- which(as.integer(m$group) == g & is.finite(m$mean))
    if (length(i) == 0) next
    graphics::lines(m$state[i], m$mean[i], col = pal[[g]], lwd = 2)
    graphics::points(m$state[i], m$mean[i], col = pal[[g]], pch = 19)
    j <- i[is.finite(m$lower[i]) & m$upper[i] > m$lower[i]]
    if (length(j))
      graphics::arrows(
        m$state[j], m$lower[j], m$state[j], m$upper[j],
        angle = 90, code = 3, length = 0.04, col = pal[[g]] )
  }

  if (k > 1)
    plot_legend(pos, lev, pal, 19, get_display(x, 'legend'), lty = 1)

  invisible(x)
}


#' @export
print.ecodive_volatility <- function (x, n = NULL, width = NULL, digits = 3L, ...) {
  gl <- x$global
  cat(plot_title(x, 'Volatility'), '\n', sep = '')
  cat(sprintf(
    'Global mean %s, SD %s; limits %s to %s (2 SD), %s to %s (3 SD)\n',
    plot_num(gl[['mean']], digits), plot_num(gl[['sd']], digits),
    plot_num(gl[['lower_2sd']], digits), plot_num(gl[['upper_2sd']], digits),
    plot_num(gl[['lower_3sd']], digits), plot_num(gl[['upper_3sd']], digits) ))
  cat(sprintf(
    'Means with %s%% bootstrap intervals (%d resamples):\n',
    format(100 * attr(x, 'conf_level', exact = TRUE)),
    as.integer(attr(x, 'n_resamples', exact = TRUE)) ))
  head <- sprintf('A data.frame: %d %s %d', nrow(x$means), print_sym()$times, length(x$means))
  cat(print_table_lines(x$means, head, n = n, width = width, digits = digits), sep = '\n')
  tr <- x$trajectories
  if (!is.null(tr))
    cat(sprintf('Trajectories of %d subjects ($trajectories)\n', length(unique(tr$subject))))
  invisible(x)
}



#------------------------------------------------------------------------------
# Helpers
#------------------------------------------------------------------------------

# `v` (a string) matched to `choices`: exactly, or by a unique prefix,
# ignoring case.
long_choice <- function (v, arg, choices) {
  if (is.character(v) && length(v) == 1 && !is.na(v) && nzchar(v)) {
    needle <- tolower(v)
    hit    <- if (needle %in% choices) needle else choices[startsWith(choices, needle)]
    if (length(hit) == 1) return (hit)
  }
  stop('`', arg, '` must be one of ', paste0("'", choices, "'", collapse = ', '), '.', call. = FALSE)
}


# `data` checked and matched to the n samples (named `labels`): list(data,
# idx), or NULLs.
long_data <- function (data, labels, n, ref) {
  if (is.null(data)) return (list(data = NULL, idx = NULL))
  if (!is.data.frame(data) && !is.matrix(data))
    stop('`data` must be a data.frame.', call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
  list(data = data, idx = match_rows(labels, n, rn, nrow(data), ref = ref))
}


# The samples of `x` for first_differences() and volatility(): a numeric
# vector, the name of a numeric column of `data`, or (with `table`) a count
# table. Returns list(values (n x p), table, features, name, labels, n,
# data, idx, ref).
long_x <- function (x, data, margin = 1L, norm = 'none', pseudocount = NULL, table = FALSE) {

  # A column of `data`: its rows are the samples.
  if (is.character(x) && length(x) == 1 && !is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!x %in% names(data))
      stop('`x` must be the name of a column of `data`, but ', name_list(x), ' is not one.', call. = FALSE)
    v <- data[[x]]
    if (!is.numeric(v) || is.factor(v))
      stop('`data` column ', name_list(x), ' must be numeric.', call. = FALSE)
    n      <- nrow(data)
    labels <- if (.row_names_info(data) > 0) rownames(data) else NULL
    return (list(
      values = matrix(as.double(v), ncol = 1), table = FALSE, name = x,
      labels = labels, n = n, data = data, idx = seq_len(n), ref = '`data`' ))
  }

  # A vector, one value per sample.
  if (is.atomic(x) && is.null(dim(x))) {
    if (!is.numeric(x) || is.factor(x))
      stop(
        '`x` must be a numeric vector with one value per sample',
        if (!is.null(data)) ', or the name of a numeric column of `data`',
        if (table) ', or a count table', '.', call. = FALSE )
    n      <- length(x)
    labels <- names(x)
    if (!is.null(labels) && (anyNA(labels) || anyDuplicated(labels) || !all(nzchar(labels))))
      labels <- NULL
    if (n == 0)
      stop('`x` has no values.', call. = FALSE)
    return (c(
      list(values = matrix(as.double(x), ncol = 1), table = FALSE, name = NULL,
           labels = labels, n = n, ref = '`x`'),
      long_data(data, labels, n, '`x`') ))
  }

  if (!table)
    stop(
      '`x` must be a numeric vector with one value per sample, or the name of ',
      'a numeric column of `data`.', call. = FALSE )

  # A count table.
  counts <- x
  validate_counts()
  validate_margin()
  smp <- count_samples(counts, margin)
  tab <- normalize(counts, norm, margin, pseudocount, cpus = 1L)
  tab <- as.matrix(tab)
  if (margin == 2L) tab <- t(tab)
  storage.mode(tab) <- 'double'
  feat <- colnames(tab)
  if (is.null(feat)) feat <- as.character(seq_len(ncol(tab)))

  c(list(values = unname(tab), table = TRUE, features = feat, name = NULL,
         labels = smp$labels, n = smp$n, ref = '`x`'),
    long_data(data, smp$labels, smp$n, '`x`') )
}


# One per-sample variable (`subjects`, `states`, `groups`) of the samples
# in `inp` (from long_x() or long_data()): list(val, col (its column name
# in `data`, or NULL)). No NAs.
long_var <- function (v, arg, inp) {

  if (is.null(v))
    stop('`', arg, '` is required: one value per sample, or a column name of `data`.', call. = FALSE)

  val <- per_sample(v, arg, inp$data, inp$idx, inp$labels, inp$n, ref = inp$ref)
  col <- if (!is.null(inp$data) && is.character(v) && length(v) == 1 && v %in% names(inp$data)) v

  if (anyNA(val))
    stop(
      if (is.null(col)) paste0('`', arg, '` must not contain NAs.')
      else paste0('Missing values (NA) in `data` column ', name_list(col), '.'),
      call. = FALSE )

  list(val = val, col = col)
}


# `subjects` and `states` of the samples in `inp`: list(subjects, states,
# states_col). With `subjects_ok`, `subjects` may be NULL.
long_meta <- function (inp, subjects, states, no_subjects, no_states, subjects_ok = FALSE) {

  if (no_states) states <- NULL
  if (no_subjects && !subjects_ok) subjects <- NULL

  if (is.character(subjects) && is.character(states) && length(subjects) == 1 &&
      length(states) == 1 && identical(subjects, states))
    stop('`subjects` and `states` must be different columns.', call. = FALSE)

  sv <- if (!is.null(subjects) || !subjects_ok) long_var(subjects, 'subjects', inp)
  tv <- long_var(states, 'states', inp)

  st <- tv$val
  if (!is.numeric(st) || is.factor(st))
    stop('`states` must be numeric, such as days or visit numbers.', call. = FALSE)
  if (!all(is.finite(st)))
    stop('`states` must be finite numbers.', call. = FALSE)
  st <- as.double(st)

  if (length(unique(st)) < 2)
    stop('`states` must have at least two distinct values (it has only ',
         format(st[[1]]), ').', call. = FALSE)

  list(subjects = sv$val, states = st, states_col = tv$col)
}


# The distinct subjects, in order: a factor's levels (those used), or sorted.
long_levels <- function (v) {
  if (is.factor(v)) return (levels(droplevels(v)))
  as.character(sort(unique(v), method = 'radix'))
}

# `v` as a factor, in long_levels()'s order.
long_factor <- function (v) {
  factor(as.character(v), levels = long_levels(v))
}


#------------------------------------------------------------------------------
# The pairs of samples (q2's _first_differences() and _get_group_pairs()).
#
#   subj, st     Each sample's subject and state.
#   baseline     NULL (consecutive states) or one of the states.
#   replicates   'error', 'drop' or 'random'.
#
# Returns list(from, to (sample indices), k (the pair of states), baseline):
# the pairs in order of pair of states, then subject.
#------------------------------------------------------------------------------
long_pairs <- function (subj, st, baseline, replicates, seed) {

  st_u <- sort(unique(st))

  if (!is.null(baseline)) {
    if (!is.numeric(baseline) || length(baseline) != 1 || !is.finite(baseline))
      stop('`baseline` must be NULL or one of the states (a number).', call. = FALSE)
    if (!baseline %in% st_u)
      stop(sprintf('`baseline` must be one of the states: %s is not in %s.',
                   format(baseline), paste(vapply(st_u, format, ''), collapse = ', ')), call. = FALSE)
    baseline <- st_u[match(baseline, st_u)]
    from_st  <- rep(baseline, length(st_u) - 1L)
    to_st    <- st_u[st_u != baseline]
  }
  else {
    from_st <- st_u[-length(st_u)]
    to_st   <- st_u[-1L]
  }

  lev <- long_levels(subj)
  S   <- length(lev)
  sc  <- match(as.character(subj), lev)
  ti  <- match(st, st_u)
  n_t <- length(st_u)

  # Each (subject, state) cell's samples, in input order.
  cell  <- (ti - 1L) * S + sc
  cnt   <- tabulate(cell, S * n_t)
  first <- match(seq_len(S * n_t), cell)

  rep_cells <- which(cnt > 1L)
  if (length(rep_cells) && replicates == 'error') {
    c1 <- rep_cells[[1]]
    s1 <- (c1 - 1L) %% S + 1L
    t1 <- (c1 - 1L) %/% S + 1L
    stop(sprintf(
      paste0('Detected replicate samples for subject %s at state %s (%d samples). ',
             'Remove the replicate samples, or set `replicates` to "drop" or "random" ',
             'to select how replicates are handled.'),
      lev[[s1]], format(st_u[[t1]]), cnt[[c1]]), call. = FALSE)
  }

  members <- if (length(rep_cells)) split(seq_along(cell), factor(cell, levels = rep_cells))

  pick <- function (k, s, side, t) {
    c1 <- (t - 1L) * S + s
    if (cnt[[c1]] == 1L) return (first[[c1]])
    m <- members[[as.character(c1)]]
    stream <- ((k - 1) * S + (s - 1)) * 2 + side
    u <- .Call(C_rng_dist, 1L, 1L, 0, seed, stream)
    m[[min(length(m), floor(u * length(m)) + 1L)]]
  }

  from <- to <- kk <- integer(0)
  for (k in seq_along(to_st)) {
    ta  <- match(from_st[[k]], st_u)
    tb  <- match(to_st[[k]], st_u)
    ca  <- cnt[(ta - 1L) * S + seq_len(S)]
    cb  <- cnt[(tb - 1L) * S + seq_len(S)]
    ok  <- if (replicates == 'drop') ca == 1L & cb == 1L else ca >= 1L & cb >= 1L
    for (s in which(ok)) {
      from <- c(from, pick(k, s, 0L, ta))
      to   <- c(to,   pick(k, s, 1L, tb))
      kk   <- c(kk, k)
    }
  }

  if (length(from) == 0) long_empty()

  list(from = from, to = to, k = kk, baseline = baseline)
}


long_empty <- function () {
  stop(
    'No pairs of samples: no subject has samples at two states that are ',
    'compared (or they were dropped as replicates, or their values are NA). ',
    'Check `subjects`, `states`, `baseline` and `replicates`.', call. = FALSE )
}


# The key columns of a first_differences() or first_distances() result:
# the pairs `pr` whose `keep` is TRUE.
long_frame <- function (labels, subj, st, pr, keep) {
  to   <- pr$to[keep]
  from <- pr$from[keep]
  res  <- data.frame(
    sample     = if (is.null(labels)) to else labels[to],
    subject    = subj[to],
    state_from = st[from],
    state      = st[to],
    stringsAsFactors = FALSE )
  rownames(res) <- NULL
  res
}


# The uniforms that volatility() draws for its cell `cell` (0-based) of
# `n` values: n_resamples x n of them, resample after resample, for tests
# of its intervals against a direct implementation of Vega's bootstrapCI().
volatility_details <- function (n, n_resamples, seed, cell) {
  .Call(C_rng_dist, as.integer(n * n_resamples), 1L, 0, as.integer(seed), cell)
}
