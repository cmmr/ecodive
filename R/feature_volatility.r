# feature_volatility(): q2-longitudinal's feature-volatility pipeline.
# Ported from q2-longitudinal 2026.10.0 (feature_volatility(),
# plot_feature_volatility(), _importance_filtering() and
# _summarize_feature_stats()) and q2-sample-classifier 2026.10.0
# (regress_samples(): split_table(), fit_regressor(), the accuracy of the
# held-out samples' predictions), with its recursive feature elimination
# (_rfecv_feature_selection()). The random forest is random_forest()'s
# (src/random_forest.c: scikit-learn 1.9.1), the recursive feature
# elimination is scikit-learn's RFECV (rf_rfe() in R/random_forest.r), and
# the train/test split is scikit-learn's train_test_split() (ShuffleSplit)
# with PCG draws (scikit-learn: BSD-3, Copyright (c) 2007-2026 The
# scikit-learn developers; see inst/COPYRIGHTS).
#
# BSD 3-Clause License
#
# Copyright (c) 2017-2026, QIIME 2 development team.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#' Feature Volatility
#'
#' Which features change most predictably over time, and how: q2-longitudinal's
#' `feature-volatility` pipeline. A random-forest regression of each sample's
#' state (time point) on its features ranks the features by importance; the
#' most important ones are summarized (their mean change between states,
#' global mean, spread and so on), and `plot()` draws the importances beside
#' the volatility chart of a chosen feature.
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#' @inheritSection documentation Input Types
#'
#' @param counts  A count table: samples \eqn{\times} features, or any of the
#'        table types listed under **Input Types** (sparse tables stay
#'        sparse).
#'
#' @param subjects  Optional: the subject (individual) of each sample, for
#'        the chart's trajectories (q2's `individual_id_column`): one value
#'        per sample, or a column name of `data`. No `NA`s. `NULL` for none.
#'
#' @param states  The state (time point) of each sample, the forest's
#'        target: numeric, such as days or visit numbers (q2's
#'        `state_column`), one value per sample or a column name of `data`.
#'        No `NA`s; at least two distinct values.
#'
#' @param groups  Optional: the group of each sample, a mean line per group
#'        in the chart, given as `states` is. No `NA`s. Default: `NULL`, one
#'        group of every sample.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding the columns named in `subjects`, `states` and
#'        `groups`. Rows are matched to the samples in `counts` by row name,
#'        and rows for other samples are ignored; without row names, they
#'        must be in the samples' order. Default: `NULL`.
#'
#' @param n_features  How many of the most important features to keep
#'        (q2's `feature_count`), of those recursive feature elimination
#'        selected: a whole number, or `Inf` for all of them. Default: `100`,
#'        as q2.
#'
#' @param test_size  The fraction of the samples held out of the forest, to
#'        measure its accuracy on (q2's `test_size`): \eqn{\lceil}
#'        `test_size` \eqn{\times n \rceil} samples, drawn at random. `0`
#'        trains on every sample (accuracy then only out of bag). Default:
#'        `0.2`, as q2.
#'
#' @param rfe  Whether to select features first by recursive feature
#'        elimination with cross-validation (scikit-learn's `RFECV`), as q2's
#'        pipeline always does (its `optimize_feature_selection`). Default:
#'        `FALSE`: one forest on every feature ranks them all. `TRUE` is
#'        q2's pipeline, at the cost of about \eqn{(}`rfe_folds`
#'        \eqn{+ 1) / } `rfe_step` forests instead of one.
#'
#' @param rfe_step  How many features recursive feature elimination drops at
#'        each round: a fraction of the features (below 1) or a count.
#'        Default: `0.05`, as q2 (5%, at least 1).
#'
#' @param rfe_folds  The number of cross-validation folds that choose how
#'        many features to keep. Default: `5`, as q2.
#'
#' @param importance  What ranks the features: `'impurity'` (the forest's
#'        mean decrease in impurity, scikit-learn's `feature_importances_`,
#'        as q2; the default) or `'permutation'` (out-of-bag permutation
#'        importance, slower). See [random_forest()].
#'
#' @param norm  How the forest's table is normalized: `'percent'` (relative
#'        abundances, the default), `'none'` (the values as given, as q2
#'        does), or any other `norm` of [normalize()]. The statistics and
#'        the chart always use relative abundances, as q2's do.
#'
#' @param n_trees  The number of trees. Default: `100`, as q2.
#'
#' @param conf_level  The confidence level of the chart's intervals: a
#'        number strictly between 0 and 1. Default: `0.95`, as q2's chart.
#'
#' @param n_resamples  The number of bootstrap resamples for each interval.
#'        Default: `1000`, as q2's chart.
#'
#' @param seed  An integer seed for the split, the forest and the bootstrap
#'        intervals. Default: `0`.
#'
#' @param ...  More settings of the forest, passed to [random_forest()]:
#'        `max_features`, `min_samples_leaf`, `min_samples_split`,
#'        `max_depth`, `criterion`, `max_samples`, `n_folds` (default `0`:
#'        no cross-validation, as q2's fit; the held-out samples measure the
#'        accuracy) and `keep_forest`.
#'
#' @return A list of class `ecodive_feature_volatility`:
#'
#' \describe{
#'   \item{`features`}{A data frame of the `n_features` most important
#'     features, most important first (q2's `feature_metadata.tsv`):
#'     `feature`, `importance`, `cumulative_avg_decrease` and
#'     `cumulative_avg_increase` (the sums of the negative and of the
#'     positive differences between the mean relative abundances at
#'     consecutive states), `net_avg_change` (their sum: the last state's
#'     mean minus the first's), and over all the samples `global_variance`,
#'     `global_mean`, `global_median`, `global_sd` and `global_cv` (SD /
#'     mean).}
#'   \item{`means`}{A data frame, one row per feature, group and state:
#'     `feature`, `group`, `state`, `mean` (of the relative abundances),
#'     `lower` and `upper` (the bootstrap interval of the mean, as
#'     [volatility()]'s) and `n`.}
#'   \item{`abundances`}{The relative abundances of those features: a
#'     samples \eqn{\times} features matrix.}
#'   \item{`samples`}{A data frame, one row per sample, in the input order:
#'     `sample`, `subject` (with `subjects`), `group`, `state` and `set`
#'     (`"train"` or `"test"`).}
#'   \item{`predictions`}{A data frame of the held-out samples (in the
#'     input order): `sample`, `observed` (the state) and `predicted` (by
#'     the forest). `NULL` with `test_size = 0`.}
#'   \item{`performance`}{A named numeric vector: of the held-out samples'
#'     predictions, `mse`, `r_squared` (\eqn{1 - SS_{res}/SS_{tot}}),
#'     `slope`, `intercept` and `r` (of the predicted on the observed
#'     states, as q2's scatter plot); and `oob_r_squared`, the forest's out
#'     of bag. The held-out values are `NA` with `test_size = 0`.}
#'   \item{`rfe`}{With `rfe = TRUE`, recursive feature elimination's
#'     results: `ranking`, a named integer vector over every feature (1 for
#'     the selected ones, then 2, 3, ... for those dropped later to earlier),
#'     and `scores`, a data frame of each number of features tried
#'     (`n_features`), the mean and SD (divisor n, as scikit-learn's) of the
#'     folds' \eqn{R^2}, and each fold's (`fold_1`, ...). `NULL` with
#'     `rfe = FALSE` (the default).}
#'   \item{`forest`}{The [random_forest()] result, trained on the training
#'     samples (with `rfe`, on the selected features only): every feature's
#'     importances (0 for the features not selected), out-of-bag
#'     predictions, and (with `keep_forest`, the default) the forest for
#'     `predict()`.}
#' }
#'
#' Attributes `conf_level` and `n_resamples` describe the intervals, and
#' `test_size` the held-out samples. Its title (see [ecodive_display]) is
#' "Feature Volatility".
#'
#' @details
#' **The pipeline**, as q2's:
#' 1. Hold out a random `test_size` of the samples (scikit-learn's
#'    `train_test_split()`, not stratified).
#' 2. With `rfe`, select features by recursive feature elimination on the
#'    training samples (see below).
#' 3. Grow a random forest of `n_trees` regression trees of the state on
#'    the (selected) features of the training samples, and predict the
#'    held-out samples' states, to see how well the features track time.
#' 4. Rank the features by the forest's importance and keep the
#'    `n_features` most important.
#' 5. Summarize those features' relative abundances (of the whole table,
#'    every sample): q2's `_summarize_feature_stats()`. The mean at each
#'    state is over every sample at that state, whatever its group.
#' 6. For the chart, each feature's mean and bootstrap interval for each
#'    group at each state, as [volatility()] computes them.
#'
#' **Recursive feature elimination** is scikit-learn's `RFECV` as q2 runs
#' it. The training samples are split, in order and without shuffling (they
#' are already in random order), into `rfe_folds` folds. In each, a forest
#' is grown on the other folds' samples and scored (\eqn{R^2}) on the fold;
#' the `rfe_step` least important features (impurity importance; ties in
#' the table's feature order) are dropped, and so on down to one feature.
#' The number of features with the largest total score over the folds (the
#' fewest, on ties) is kept: the same elimination on all the training
#' samples, stopped there, selects them. Every forest has the settings of
#' the final one. This grows about \eqn{(}`rfe_folds` \eqn{+ 1) / }
#' `rfe_step` forests (126 by default), so it takes about that many times
#' as long as `rfe = FALSE`.
#'
#' **Differences from q2-longitudinal**:
#' * **Random draws** (the split, the forest and the intervals) are
#'   ecodive's own seeded PCG streams, not numpy's or JavaScript's, so they
#'   differ from q2's as one `random_state` does from another. Given q2's
#'   split and `random_state`, the selected features, the forest and its
#'   importances are the same as q2's (see [random_forest()]).
#' * **The forest's table** is relative abundances by default; q2 grows it
#'   on the table it is given (usually counts): use `norm = 'none'` for that.
#' * **Features absent from every training sample** are kept, with
#'   importance 0 (so recursive feature elimination drops them first, but
#'   they count towards the fraction `rfe_step`); q2 drops them. Ties in
#'   importance are ranked in the table's feature order.
#' * **`global_cv`** is the standard deviation over the mean, as q2's
#'   "Global CV (%)" column, which despite its name is not multiplied by
#'   100.
#' * q2's chart lets the viewer pick the feature and the grouping
#'   interactively; here `plot()` draws one feature per call, and the
#'   grouping is set by `groups`.
#'
#' @references
#' Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
#' Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
#' longitudinal and paired-sample analyses of microbiome data. *mSystems*,
#' 3, e00219-18. \doi{10.1128/mSystems.00219-18}
#'
#' Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley,
#' G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning
#' tools for microbiome classification and regression. *Journal of Open
#' Source Software*, 3, 934. \doi{10.21105/joss.00934}
#'
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5-32.
#' \doi{10.1023/A:1010933404324}
#'
#' @seealso [plot.ecodive_feature_volatility()], [random_forest()],
#'   [volatility()]
#'
#' @family longitudinal
#' @export
#' @examples
#'     # 8 subjects sampled on days 0, 7, 14 and 21; feature 1 rises.
#'     meta <- data.frame(
#'       subject = rep(sprintf('s%d', 1:8), each = 4),
#'       day     = rep(c(0, 7, 14, 21), 8),
#'       diet    = rep(c('x', 'y'), each = 16) )
#'     counts <- outer(seq_len(32), 1:6, function (i, j) (i * j) %% 11 + 5)
#'     counts[, 1] <- counts[, 1] + 3 * meta$day
#'     colnames(counts) <- sprintf('F%d', 1:6)
#'
#'     fv <- feature_volatility(counts, 'subject', 'day', 'diet', data = meta,
#'                              n_trees = 20, cpus = 1)
#'     fv
#'     plot(fv)
feature_volatility <- function (
    counts, subjects, states, groups = NULL, data = NULL,
    n_features = 100L, test_size = 0.2, rfe = FALSE, rfe_step = 0.05,
    rfe_folds = 5L, importance = 'impurity',
    margin = 1L, norm = 'percent', pseudocount = NULL, n_trees = 100L,
    conf_level = 0.95, n_resamples = 1000L, seed = 0L, cpus = n_cpus(),
    ... ) {

  if (missing(subjects)) subjects <- NULL
  if (missing(states))   states   <- NULL
  fv_run(environment(), list(...))
}


# The same, with the internals the tests need:
#   numpy_random_state   NULL, or an integer: split and grow as
#                        q2-sample-classifier does with that random_state
#                        (numpy's RandomState permutation, scikit-learn's
#                        tree seeds) instead of PCG.
#   sparse               As random_forest_details()'s.
# Returns the result with `split`: list(train, test), sample positions in
# the order drawn.
feature_volatility_details <- function (
    counts, subjects, states, groups = NULL, data = NULL,
    n_features = 100L, test_size = 0.2, rfe = FALSE, rfe_step = 0.05,
    rfe_folds = 5L, importance = 'impurity',
    margin = 1L, norm = 'percent', pseudocount = NULL, n_trees = 100L,
    conf_level = 0.95, n_resamples = 1000L, seed = 0L, cpus = n_cpus(),
    ..., numpy_random_state = NULL, sparse = NULL ) {

  if (missing(subjects)) subjects <- NULL
  if (missing(states))   states   <- NULL
  fv_run(environment(), list(...), numpy_random_state = numpy_random_state,
         sparse = sparse, details = TRUE)
}


# The forest settings that `...` may give.
FV_FOREST_ARGS <- c(
  'max_features', 'min_samples_leaf', 'min_samples_split', 'max_depth',
  'criterion', 'max_samples', 'n_folds', 'keep_forest' )


fv_run <- function (env, dots, numpy_random_state = NULL, sparse = NULL, details = FALSE) {

  # The forest's own arguments, before validation changes them.
  rf <- mget(c('counts', 'margin', 'norm', 'pseudocount', 'n_trees', 'importance',
               'seed', 'cpus'), envir = env)

  with(env, {
    validate_counts()
    validate_margin()
    validate_seed()
    validate_cpus()
  })
  counts <- env$counts; margin <- env$margin; seed <- env$seed; cpus <- env$cpus

  if (length(dots)) {
    nm <- names(dots)
    if (is.null(nm) || any(!nzchar(nm)))
      stop('Arguments in `...` must be named: they are settings of the forest.', call. = FALSE)
    bad <- setdiff(nm, FV_FOREST_ARGS)
    if (length(bad))
      stop('Unknown argument', if (length(bad) > 1) 's', ' ', name_list(bad),
           '. `...` takes these settings of random_forest(): ',
           paste0('`', FV_FOREST_ARGS, '`', collapse = ', '), '.', call. = FALSE)
    if (anyDuplicated(nm))
      stop('Arguments in `...` must not be repeated.', call. = FALSE)
  }

  n_features <- env$n_features
  if (!is.numeric(n_features) || length(n_features) != 1 || is.na(n_features) ||
      n_features < 1 || (is.finite(n_features) && n_features != floor(n_features)))
    stop('`n_features` must be a whole number, at least 1, or `Inf`.', call. = FALSE)

  test_size <- env$test_size
  if (!is.numeric(test_size) || length(test_size) != 1 || is.na(test_size) ||
      test_size < 0 || test_size >= 1)
    stop('`test_size` must be a number from 0 up to (not including) 1.', call. = FALSE)

  rfe <- env$rfe
  if (!is.logical(rfe) || length(rfe) != 1 || is.na(rfe))
    stop('`rfe` must be TRUE or FALSE.', call. = FALSE)
  rfe_step <- env$rfe_step
  if (!is.numeric(rfe_step) || length(rfe_step) != 1 || !is.finite(rfe_step) || rfe_step <= 0 ||
      (rfe_step >= 1 && rfe_step != floor(rfe_step)))
    stop('`rfe_step` must be a fraction below 1 or a whole number of features.', call. = FALSE)
  rfe_folds <- env$rfe_folds
  if (!is.numeric(rfe_folds) || length(rfe_folds) != 1 || !is.finite(rfe_folds) ||
      rfe_folds < 2 || rfe_folds != floor(rfe_folds) || rfe_folds > 63)
    stop('`rfe_folds` must be a whole number from 2 to 63.', call. = FALSE)

  imp_type <- env$importance
  if (!is.character(imp_type) || length(imp_type) != 1 || !imp_type %in% c('impurity', 'permutation'))
    stop("`importance` must be 'impurity' or 'permutation'.", call. = FALSE)

  conf_level <- env$conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1 || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop('`conf_level` must be a number strictly between 0 and 1.', call. = FALSE)
  n_resamples <- env$n_resamples
  if (!is.numeric(n_resamples) || length(n_resamples) != 1 || !is.finite(n_resamples) ||
      n_resamples < 1 || n_resamples %% 1 != 0 || n_resamples > 1e7)
    stop('`n_resamples` must be a whole number from 1 to 1e7.', call. = FALSE)
  n_resamples <- as.integer(n_resamples)


  #--------------------------------------------------------
  # The samples' subjects, states and groups.
  #--------------------------------------------------------
  smp <- count_samples(counts, margin)
  n   <- smp$n
  if (n < 2) stop('`counts` must have at least two samples.', call. = FALSE)
  inp <- c(list(labels = smp$labels, n = n, ref = '`counts`'),
           long_data(env$data, smp$labels, n, '`counts`'))
  md  <- long_meta(inp, env$subjects, env$states, FALSE, is.null(env$states), subjects_ok = TRUE)
  st  <- md$states
  g   <- if (is.null(env$groups)) NULL else long_var(env$groups, 'groups', inp)
  grp <- if (is.null(g)) factor(rep('all', n)) else long_factor(g$val)
  labels <- if (is.null(smp$labels)) as.character(seq_len(n)) else smp$labels


  #--------------------------------------------------------
  # The split (scikit-learn's ShuffleSplit): the first
  # ceil(test_size * n) of a random permutation are held out.
  #--------------------------------------------------------
  n_test  <- as.integer(ceiling(test_size * n))
  n_train <- n - n_test
  if (n_train < 2)
    stop(sprintf('`test_size` (%s) leaves fewer than two of the %d samples to train on.',
                 format(test_size), n), call. = FALSE)
  if (n_test == 0) {
    train <- seq_len(n)
    test  <- integer(0)
  }
  else {
    perm <- if (!is.null(numpy_random_state))
              as.integer(.Call(C_rf_rng_test, 5L, as.double(numpy_random_state), as.double(n), n)) + 1L
            else .Call(C_rf_shuffle, seq_len(n), seed, 2^62 + 2^48)
    test  <- perm[seq_len(n_test)]
    train <- perm[n_test + seq_len(n_train)]
  }


  #--------------------------------------------------------
  # The forest, on the training samples.
  #--------------------------------------------------------
  fenv <- new.env()
  fvals <- c(rf, list(
    y = unname(st), data = NULL, max_features = NULL, min_samples_leaf = 1L,
    min_samples_split = 2L, max_depth = NULL, criterion = NULL,
    class_weight = NULL, max_samples = NULL, n_folds = 0L, folds = NULL,
    keep_forest = TRUE ))
  fvals[names(dots)] <- dots
  for (k in names(fvals)) assign(k, fvals[[k]], envir = fenv)
  fit <- rf_run(
    fenv, numpy_random_state = numpy_random_state, sparse = sparse, details = details, train = train,
    rfe = if (rfe) list(step = rfe_step, n_folds = as.integer(rfe_folds)) )

  ho <- fit$holdout
  fit$holdout <- NULL
  rfe_res <- fit$rfe
  fit$rfe <- NULL
  feats <- count_features(counts, margin)
  if (!is.null(rfe_res)) names(rfe_res$ranking) <- feats


  #--------------------------------------------------------
  # The most important features.
  #--------------------------------------------------------
  imp   <- fit$importance
  if (!is.null(rfe_res))
    imp <- imp[rfe_res$ranking[match(imp$feature, feats)] == 1L, , drop = FALSE]
  top   <- utils::head(imp, min(nrow(imp), n_features))
  j     <- match(top$feature, feats)
  rel   <- fv_relative(counts, margin, j, n, cpus)
  dimnames(rel) <- list(smp$labels, top$feature)

  stats <- fv_stats(rel, st)
  features <- data.frame(
    feature    = top$feature,
    importance = top[[imp_type]],
    stats,
    stringsAsFactors = FALSE )
  rownames(features) <- NULL


  #--------------------------------------------------------
  # The chart's means and intervals: cell k of feature f
  # draws from PCG stream 2^52 + (f - 1) * cells + k.
  #--------------------------------------------------------
  vc <- long_vol_cells(rel, st, grp, n_resamples, seed, 1 - conf_level, stream_hi = 2^20)
  means <- vc$means
  means$variable <- top$feature[means$variable]
  names(means)[[1]] <- 'feature'
  rownames(means) <- NULL


  #--------------------------------------------------------
  # Samples, held-out predictions and performance.
  #--------------------------------------------------------
  set <- rep('train', n)
  set[test] <- 'test'
  samples <- data.frame(sample = labels, stringsAsFactors = FALSE)
  if (!is.null(md$subjects)) samples$subject <- md$subjects
  samples$group <- grp
  samples$state <- st
  samples$set   <- set

  oob_r2 <- unname(fit$performance[['oob_r_squared']])
  if (length(ho$index)) {
    obs  <- st[ho$index]
    pred <- ho$value
    predictions <- data.frame(
      sample = labels[ho$index], observed = obs, predicted = pred,
      stringsAsFactors = FALSE )
    lr <- rf_linregress(obs, pred)
    performance <- c(
      mse = mean((pred - obs)^2), r_squared = rf_r2(obs, pred),
      slope = lr[[1]], intercept = lr[[2]], r = lr[[3]], oob_r_squared = oob_r2 )
  }
  else {
    predictions <- NULL
    performance <- c(
      mse = NA_real_, r_squared = NA_real_, slope = NA_real_, intercept = NA_real_,
      r = NA_real_, oob_r_squared = oob_r2 )
  }

  res <- list(
    features = features, means = means, abundances = rel, samples = samples,
    predictions = predictions, performance = performance, rfe = rfe_res, forest = fit )
  if (details) res$split <- list(train = train, test = test)

  res <- structure(
    res,
    class       = 'ecodive_feature_volatility',
    method      = 'feature_volatility',
    conf_level  = conf_level,
    n_resamples = n_resamples,
    test_size   = test_size )

  set_display(
    res, title = 'Feature Volatility', xlab = md$states_col,
    legend = if (!is.null(g)) g$col )
}


# The feature names of `counts` (as random_forest() names them).
count_features <- function (counts, margin) {
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  f <- 3L - margin
  if (length(dn[[f]])) dn[[f]] else sprintf('F%d', seq_len(dims[[f]]))
}


# Relative abundances of features `j` (q2's relative frequencies): a dense
# n x length(j) matrix, without densifying the rest of a sparse table.
fv_relative <- function (counts, margin, j, n, cpus) {

  if (get_matrix_package(counts) == 'base') {
    x <- .Call(C_normalize, counts, margin, NORM_PERCENT, numeric(0), TRUE, cpus)
    if (margin == 2L) x <- t(x)
    x <- x[, j, drop = FALSE]
    storage.mode(x) <- 'double'
    dimnames(x) <- NULL
    return (x)
  }

  tr <- .Call(C_normalize, counts, margin, NORM_PERCENT, numeric(0), FALSE, cpus)
  ri <- if (margin == 2L) tr$j else tr$i
  rj <- if (margin == 2L) tr$i else tr$j
  k  <- match(rj + 1L, j)
  ok <- !is.na(k)
  x  <- matrix(0, n, length(j))
  x[cbind(ri[ok] + 1L, k[ok])] <- as.double(tr$x[ok])
  x
}


# q2-longitudinal's _summarize_feature_stats(): the columns of `rel`
# (samples x features) summarized over states `st`. A data frame, one row
# per feature.
fv_stats <- function (rel, st) {

  st_u <- sort(unique(st))
  si   <- match(st, st_u)
  n    <- nrow(rel)

  one <- function (v) {
    sm  <- vapply(seq_along(st_u), function (s) { w <- v[si == s]; sum(w) / length(w) }, 0)
    d   <- diff(sm)
    dec <- sum(d[d < 0])
    inc <- sum(d[d > 0])
    gm  <- sum(v) / n
    gv  <- sum((gm - v)^2) / (n - 1)
    gs  <- sqrt(gv)
    c(dec, inc, dec + inc, gv, gm, stats::median(v), gs, gs / gm)
  }

  m <- if (ncol(rel)) t(vapply(seq_len(ncol(rel)), function (k) one(rel[, k]), numeric(8)))
       else matrix(numeric(0), 0, 8)
  colnames(m) <- c(
    'cumulative_avg_decrease', 'cumulative_avg_increase', 'net_avg_change',
    'global_variance', 'global_mean', 'global_median', 'global_sd', 'global_cv' )
  as.data.frame(m)
}



#' Plot a Feature Volatility Result
#'
#' The importances of the most important features as bars, beside the
#' volatility chart of one of them (as [volatility()]'s `plot()`): its
#' relative abundance over time, each group's mean with a bootstrap
#' interval, each subject's trajectory, and the global mean and control
#' limits.
#'
#' @inheritParams plot_longitudinal
#'
#' @param x  A [feature_volatility()] result.
#'
#' @param feature  The feature to chart: its name, or its rank (`1`, the
#'        default, is the most important).
#'
#' @param n_bars  How many of the most important features to draw as bars.
#'        Default: `20`.
#'
#' @param ...  Passed to [graphics::plot.default()] for the chart, such as
#'        `main`, `xlab`, `ylab` or `ylim`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     meta <- data.frame(
#'       subject = rep(sprintf('s%d', 1:8), each = 4),
#'       day     = rep(c(0, 7, 14, 21), 8) )
#'     counts <- outer(seq_len(32), 1:6, function (i, j) (i * j) %% 11 + 5)
#'     counts[, 1] <- counts[, 1] + 3 * meta$day
#'     colnames(counts) <- sprintf('F%d', 1:6)
#'     fv <- feature_volatility(counts, 'subject', 'day', data = meta,
#'                              n_trees = 20, cpus = 1)
#'
#'     plot(fv)
#'     plot(fv, feature = fv$features$feature[[1]], spaghetti = FALSE, limits = FALSE)
plot.ecodive_feature_volatility <- function (
    x, feature = 1L, n_bars = 20L, spaghetti = TRUE, limits = TRUE,
    legend = 'topright', col = NULL, ... ) {

  ft <- x$features
  if (nrow(ft) == 0) stop('`x` has no features to draw.', call. = FALSE) # nocov
  if (is.character(feature) && length(feature) == 1 && feature %in% ft$feature) {
    f <- match(feature, ft$feature)
  }
  else if (is.numeric(feature) && length(feature) == 1 && !is.na(feature) &&
           feature >= 1 && feature <= nrow(ft) && feature == floor(feature)) {
    f <- as.integer(feature)
  }
  else {
    stop('`feature` must be the name of one of the features of `x` or its ',
         'rank, from 1 to ', nrow(ft), '.', call. = FALSE)
  }
  nb  <- plot_count(n_bars, 'n_bars')
  if (nb < 1) stop('`n_bars` must be at least 1.', call. = FALSE)
  pos <- plot_legend_pos(legend)
  fname <- ft$feature[[f]]

  # The chosen feature's volatility result.
  m <- x$means[x$means$feature == fname, -1L]
  rownames(m) <- NULL
  gm <- ft$global_mean[[f]]; gs <- ft$global_sd[[f]]
  global <- c(
    mean = gm, sd = gs, lower_3sd = gm - 3 * gs, lower_2sd = gm - 2 * gs,
    upper_2sd = gm + 2 * gs, upper_3sd = gm + 3 * gs )
  smp <- x$samples
  tr  <- NULL
  if (!is.null(smp$subject)) {
    sc <- match(as.character(smp$subject), long_levels(smp$subject))
    o  <- order(sc, smp$state, seq_len(nrow(smp)), method = 'radix')
    tr <- data.frame(
      sample = smp$sample[o], subject = smp$subject[o], group = smp$group[o],
      state = smp$state[o], value = unname(x$abundances[o, f]),
      stringsAsFactors = FALSE )
  }
  vol <- structure(
    list(means = m, global = global, trajectories = tr),
    class = 'ecodive_volatility',
    conf_level = attr(x, 'conf_level', exact = TRUE),
    n_resamples = attr(x, 'n_resamples', exact = TRUE) )
  vol <- set_display(
    vol, title = fname, ylab = 'Relative abundance',
    xlab = get_display(x, 'xlab'), legend = get_display(x, 'legend') )

  old <- graphics::par(c('mfrow', 'mai'))
  on.exit({
    graphics::par(mfrow = old$mfrow)
    graphics::par(mai = old$mai)
  })
  graphics::layout(matrix(1:2, 1), widths = c(1, 2))

  # The importance bars, most important at the top.
  rows <- rev(seq_len(min(nb, nrow(ft))))
  kk   <- length(rows)
  v    <- ft$importance[rows]
  linch <- max(graphics::strwidth(ft$feature[rows], 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai  <- old$mai
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)
  lo <- min(0, v, na.rm = TRUE)
  hi <- max(0, v, na.rm = TRUE)
  if (!(hi > lo)) hi <- lo + 1
  graphics::plot.default(
    0, 0, type = 'n', yaxt = 'n', xlim = c(lo, hi), ylim = c(0.5, kk + 0.5),
    main = plot_title(x, 'Feature Volatility'), xlab = 'Importance', ylab = '' )
  y    <- seq_len(kk)
  bcol <- ifelse(rows == f, PLOT_HIGHLIGHT, PLOT_MUTED)
  graphics::rect(0, y - 0.35, v, y + 0.35, col = bcol, border = NA)
  graphics::axis(2, at = y, labels = ft$feature[rows], tick = FALSE, las = 1, cex.axis = 0.8)
  graphics::abline(v = 0, col = PLOT_MUTED)

  # The chart.
  graphics::par(mai = old$mai)
  plot.ecodive_volatility(vol, spaghetti = spaghetti, limits = limits, legend = legend,
                          col = col, ...)

  invisible(x)
}


#' @export
print.ecodive_feature_volatility <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  fo   <- x$forest
  pf   <- x$performance
  smp  <- x$samples
  n_te <- sum(smp$set == 'test')
  cat(plot_title(x, 'Feature Volatility'), '\n', sep = '')
  cat(sprintf(
    'Random forest of %d trees, trained on %d of %d samples',
    as.integer(attr(fo, 'n_trees', exact = TRUE)), nrow(smp) - n_te, nrow(smp) ))
  if (n_te > 0) {
    cat(sprintf(
      '; on the %d held out: MSE %s, R^2 %s, r %s',
      n_te, plot_num(pf[['mse']], digits),
      plot_num(pf[['r_squared']], digits), plot_num(pf[['r']], digits) ))
  }
  cat(sprintf('; OOB R^2 %s\n', plot_num(pf[['oob_r_squared']], digits)))
  ft  <- x$features
  nf  <- nrow(fo$importance)
  rf  <- x$rfe
  if (!is.null(rf)) {
    ns <- sum(rf$ranking == 1L)
    sc <- rf$scores
    cat(sprintf(
      'Recursive feature elimination kept %d of %d features (CV R^2 %s)\n',
      ns, nf, plot_num(sc$mean[match(ns, sc$n_features)], digits) ))
    nf <- ns
  }
  cat(sprintf(
    'The %d most important of %d %sfeatures, with relative abundances over time:\n',
    nrow(ft), nf, if (!is.null(rf)) 'selected ' else '' ))
  head <- sprintf('A data.frame: %d %s %d', nrow(ft), print_sym()$times, length(ft))
  cat(print_table_lines(ft, head, n = n, width = width, digits = digits), sep = '\n')
  cat(sprintf(
    'Means by state with %s%% bootstrap intervals in $means; trajectories in $abundances\n',
    format(100 * attr(x, 'conf_level', exact = TRUE)) ))
  invisible(x)
}
