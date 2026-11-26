# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# random_forest(): random forests for classification and regression. The
# trees, forest, bootstrap, out-of-bag predictions and impurity importances
# are ported from scikit-learn 1.9.1 (BSD-3) in src/random_forest.c; the
# workflow (cross-validation, mean importances over folds, the accuracy
# summary and ROC averages) follows q2-sample-classifier 2026.10.0 (BSD-3:
# Copyright (c) 2017-2026, QIIME 2 development team; see inst/COPYRIGHTS),
# and the folds and roc_curve() follow scikit-learn's StratifiedKFold,
# KFold and sklearn/metrics/_ranking.py.


#' Random forests
#'
#' Grows a random forest (Breiman 2001) of classification or regression
#' trees on a feature table: to predict each sample's class (disease or
#' health, body site) or a number (age, a clinical score), to rank the
#' features by importance, and to estimate how well that works, by
#' cross-validation and out of bag. The trees are an exact port of
#' scikit-learn's, and the workflow is q2-sample-classifier's.
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#' @inheritSection ml_section Machine learning
#' @inheritSection documentation Input Types
#'
#' @param counts  A numeric matrix of feature values (samples \eqn{\times}
#'        features), or any of the table types listed under "Input Types"
#'        (sparse tables stay sparse). Typically counts.
#'
#' @param y   What to predict: a vector with one value per sample, or the
#'   name of a column of `data`. A factor, character or logical `y` makes a
#'   classification forest (the classes are its values, or a factor's
#'   levels); a numeric `y`, a regression forest. Samples whose `y` is `NA`
#'   are left out, with a warning.
#'
#' @param data   A data frame of sample metadata, with one row per sample
#'   (matched to `counts` by row names when both have them), when `y` names
#'   one of its columns. Default: `NULL`
#'
#' @param norm   How to normalize each sample before growing the trees:
#'   `'percent'` (relative abundances, the default), `'none'` (the values
#'   as given, for any numeric table), or one of [normalize()]'s other
#'   choices. Trees don't change under monotone transforms of a feature, so
#'   only transforms across features (percent, clr, chord, ...) matter.
#'   Default: `'percent'`
#'
#' @param n_trees   The number of trees in each forest. Default: `500`
#'
#' @param max_features   How many features each split chooses among: a
#'   count, a fraction of the features (below 1), `'sqrt'`, `'log2'`, or
#'   `'all'`; `NULL` for scikit-learn's defaults (`'sqrt'` for
#'   classification, `'all'` for regression). A vector of these makes a
#'   grid (see Details). Default: `NULL`
#'
#' @param min_samples_leaf   The fewest (distinct) samples a leaf may hold.
#'   A vector of values makes a grid (see Details). Default: `1`
#'
#' @param min_samples_split   The fewest (distinct) samples a node needs to
#'   be split; at least twice `min_samples_leaf` is used. Default: `2`
#'
#' @param max_depth   The deepest a tree may grow, or `NULL` for no limit.
#'   Default: `NULL`
#'
#' @param criterion   The impurity to split on: `'gini'` (the default for
#'   classification) or `'entropy'`; `'squared_error'` for regression.
#'   Default: `NULL`
#'
#' @param class_weight   Weights for the classes (classification only):
#'   `NULL` (equal), `'balanced'` (\eqn{n / (k n_c)} for a class of
#'   \eqn{n_c} of the \eqn{n} training samples, \eqn{k} classes),
#'   `'balanced_subsample'` (the same, computed on each tree's bootstrap
#'   sample), or a named numeric vector of weights (classes not named get
#'   1). As in scikit-learn 1.9, weights other than
#'   `'balanced_subsample'` drive a weighted bootstrap (samples are drawn in
#'   proportion to their class weight). Default: `NULL`
#'
#' @param max_samples   The size of each tree's bootstrap sample: `NULL`
#'   (as many draws as training samples), a count, or a fraction (below 1)
#'   of the training samples (of their total class weight, with
#'   `class_weight`). Default: `NULL`
#'
#' @param importance   Which feature importances to compute: `'impurity'`
#'   (mean decrease in impurity, free), `'permutation'` (out-of-bag
#'   permutation importance), or `'both'`. Default: `'impurity'`
#'
#' @param n_folds   The number of cross-validation folds; `0` or `1` for
#'   none (out-of-bag estimates only). Default: `5`
#'
#' @param folds   A fold for each sample instead (any labels; each label is
#'   a fold), which overrides `n_folds` and `seed` for the folds. Matched to
#'   the samples like `y`. Default: `NULL`
#'
#' @param keep_forest   Keep the full-data forest in the result, for
#'   [predict()][predict.ecodive_random_forest()]. Default: `TRUE`
#'
#' @param seed   An integer seed for the bootstraps, the trees' feature
#'   draws, the folds and the permutations. The same `seed` gives the same
#'   forest for any `cpus`. Default: `0`
#'
#' @details
#'
#' ## The trees
#'
#' Each tree is grown on a bootstrap sample of the training samples
#' (drawn with replacement; the number of draws of each sample is its
#' weight), splitting each node at the threshold, over a random subset of
#' `max_features` features, that most reduces the weighted impurity: Gini
#' or entropy for classes, the squared error for numbers. Trees grow until
#' their leaves are pure or reach `min_samples_leaf`, `min_samples_split`
#' or `max_depth`. A tree predicts the (weighted) class fractions or mean
#' of the leaf a sample falls in, and the forest averages its trees: class
#' probabilities (soft voting, as scikit-learn), whose largest gives the
#' predicted class (the first class on ties), or the mean of the trees'
#' predictions.
#'
#' This is scikit-learn 1.9.1's `RandomForestClassifier` and
#' `RandomForestRegressor`, ported to C: given the same random draws, the
#' trees are identical (the tests check every node against scikit-learn's).
#' The draws themselves come from ecodive's PCG generator, one stream per
#' tree, so a forest is the same on any number of threads, and differs
#' from scikit-learn's as a different seed would.
#'
#' ## Out of bag and cross-validated
#'
#' About a third of the samples are left out of each tree's bootstrap.
#' Predicting each sample from only the trees that never saw it gives
#' **out-of-bag** (OOB) predictions and accuracy (or \eqn{R^2}) for free.
#' **Cross-validation** grows a separate forest on each fold's training
#' samples and predicts its held-out samples (q2-sample-classifier's
#' `classify-samples-ncv` and `regress-samples-ncv`): with `n_folds` folds,
#' stratified by class for classification (as scikit-learn's
#' `StratifiedKFold`), or plain (`KFold`) for regression. Both are honest
#' estimates of how well the forest predicts new samples; cross-validation
#' also gives ROC curves, a confusion matrix and fold-to-fold variation. A
#' full-data forest is also grown: it gives the OOB results, the reported
#' importances and [predict()][predict.ecodive_random_forest()].
#'
#' ## Choosing settings
#'
#' A vector of `max_features` or `min_samples_leaf` values makes a grid
#' (all their combinations). Within each training set (each fold's, and
#' the full data), a forest is grown for every grid point and the one with
#' the lowest out-of-bag error is kept (the first on ties). The held-out
#' samples never choose a setting, so the cross-validated estimates stay
#' honest. The `tuning` element lists every grid point's OOB error.
#'
#' ## Importance
#'
#' `'impurity'` importance (mean decrease in impurity, scikit-learn's
#' `feature_importances_`) is each feature's share of the trees' total
#' impurity decrease, averaged over trees and scaled to sum to 1. It is
#' biased toward features with many distinct values (Strobl et al. 2007).
#' `'permutation'` importance (Breiman 2001) permutes each feature among
#' each tree's out-of-bag samples and records how much the tree's
#' out-of-bag error (misclassification rate, or mean squared error) grows;
#' it is the mean over all trees (a tree that doesn't split on the feature
#' adds 0), with its standard deviation over trees. Importances measure
#' association with the prediction, not effects; correlated features share
#' it.
#'
#' Features should not be filtered by their association with `y` before
#' cross-validation: that leaks the held-out samples into the model and
#' inflates the estimates. Relative abundances are compositional, so a
#' feature's importance can reflect changes in others; see "Machine
#' learning" below on transforms (Yerke et al. 2024 found relative
#' abundances did best with random forests).
#'
#' ## Recipes
#'
#' * **q2-sample-classifier's settings:** `n_trees = 100`, `n_folds = 5`,
#'   `norm = 'none'` (q2 passes the table it's given).
#' * **randomForest's scaled permutation importance** (its default
#'   `importance(fit)`): `importance = 'permutation'`, then
#'   `permutation / (permutation_sd / sqrt(n_trees))`.
#' * **A microbiota maturity index** (Subramanian et al. 2014): a
#'   regression forest of age on healthy children's samples, then
#'   [predict()][predict.ecodive_random_forest()] for the others; their
#'   relative maturity is the predicted minus a spline of predicted age on
#'   age in the healthy children.
#'
#' @section Differences from scikit-learn, q2-sample-classifier, randomForest and ranger:
#'
#' * **Random draws:** PCG streams, not numpy's MT19937, so forests differ
#'   from scikit-learn's as with a different `random_state`. (A parity mode
#'   that reproduces scikit-learn's draws is used by the tests.)
#' * **Defaults:** 500 trees (randomForest's and ranger's; scikit-learn's
#'   and q2's 100), with scikit-learn's `max_features` (all features for
#'   regression, where randomForest uses a third) and leaf sizes (1 for
#'   regression, where randomForest and ranger use 5). Tables are relative
#'   abundances by default (`norm = 'percent'`).
#' * **Soft voting:** classes are predicted from the mean of the trees'
#'   class probabilities, as scikit-learn; randomForest counts votes.
#' * **Importances:** impurity importances sum to 1, as scikit-learn's;
#'   randomForest's `MeanDecreaseGini` and ranger's `'impurity'` are not
#'   scaled. Permutation importance is unscaled (randomForest's
#'   `scale = FALSE`, ranger's).
#' * **Never out of bag:** a sample that every bootstrap drew gets `NA` as
#'   its OOB prediction, with a warning, and is left out of the OOB scores;
#'   scikit-learn predicts 0 for it.
#' * **Folds:** as scikit-learn's `StratifiedKFold` (classification) and
#'   `KFold` (regression), shuffled by PCG instead of numpy. All features
#'   are kept in every fold (q2's DictVectorizer drops features absent from
#'   a fold's training samples).
#' * **Float32:** values are rounded to single precision before growing and
#'   predicting, as scikit-learn does; values that differ only beyond that
#'   (a relative 6e-8) are ties.
#' * **Not ported:** q2's recursive feature elimination
#'   (`optimize_feature_selection`) and randomized parameter search, and
#'   scikit-learn's `bootstrap = FALSE`, `max_leaf_nodes`,
#'   `min_impurity_decrease`, `min_weight_fraction_leaf`, `ccp_alpha`,
#'   `monotonic_cst` and missing values (`NA` in `counts` is an error).
#'
#' @return A list of class `ecodive_random_forest`:
#'
#' \describe{
#'   \item{`importance`}{A data frame, one row per feature, sorted by the
#'     requested importance (largest first): `feature`, `impurity` (the
#'     full-data forest's), `impurity_cv_mean` and `impurity_cv_sd` (over
#'     the cross-validation folds' forests; q2 reports the mean), and, with
#'     `importance = 'permutation'` or `'both'`, `permutation` and
#'     `permutation_sd`.}
#'   \item{`predictions`}{A data frame, one row per sample: `sample`,
#'     `observed`, `predicted` (cross-validated: by the forest of the fold
#'     that held the sample out), `fold`, `oob_predicted`, and for
#'     classification one cross-validated probability per class
#'     (`prob_<class>`).}
#'   \item{`performance`}{A named numeric vector. Classification:
#'     `accuracy` (cross-validated), `accuracy_fold_mean`,
#'     `accuracy_fold_sd` (over folds, `sd` with divisor n, as q2),
#'     `balanced_accuracy`, `baseline_accuracy` (the largest class's
#'     share), `accuracy_ratio`, `auc` (two classes) or `auc_macro` and
#'     `auc_micro`, and `oob_accuracy`. Regression: `mse`, `r_squared`,
#'     `mse_fold_mean`, `mse_fold_sd`, `slope`, `intercept` and `r` (of the
#'     predicted on the observed values, as q2), and `oob_r_squared`.
#'     Without cross-validation, only the OOB values are not `NA`.}
#'   \item{`confusion`}{Classification: the observed (rows) by predicted
#'     (columns) counts, cross-validated (out of bag without
#'     cross-validation).}
#'   \item{`roc`}{Classification: a data frame `class`, `fpr`, `tpr`,
#'     `threshold`: scikit-learn's `roc_curve()` of each class against the
#'     rest (only the second class with two classes), and for more than two
#'     classes q2's `'micro'` and `'macro'` averages; from the
#'     cross-validated probabilities (out of bag without
#'     cross-validation).}
#'   \item{`tuning`}{With a grid: a data frame of `training_set` (`'full'`
#'     or `'fold <k>'`), `max_features`, `min_samples_leaf`, `oob_error`
#'     and `chosen`.}
#'   \item{`forest`}{With `keep_forest`: the full-data forest, for
#'     [predict()][predict.ecodive_random_forest()].}
#' }
#'
#' Attributes: `method` (`'random_forest'`), `type` (`'classification'` or
#' `'regression'`), `n_trees`, `max_features` and `min_samples_leaf` (of
#' the full-data forest: as chosen, with a grid), `criterion`,
#' `class_weight` and `n_folds`.
#'
#' @references
#'
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5-32.
#' \doi{10.1023/A:1010933404324}
#'
#' Breiman, L., Friedman, J. H., Olshen, R. A., & Stone, C. J. (1984).
#' *Classification and Regression Trees*. Wadsworth.
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830.
#'
#' Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley,
#' G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning
#' tools for microbiome classification and regression. *Journal of Open
#' Source Software*, 3, 934. \doi{10.21105/joss.00934}
#'
#' Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias in
#' random forest variable importance measures. *BMC Bioinformatics*, 8, 25.
#' \doi{10.1186/1471-2105-8-25}
#'
#' Subramanian, S., et al. (2014). Persistent gut microbiota immaturity in
#' malnourished Bangladeshi children. *Nature*, 510, 417-421.
#' \doi{10.1038/nature13421}
#'
#' Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
#' normalizations outperform compositional data transformations in machine
#' learning applications. *Microbiome*, 12, 45.
#' \doi{10.1186/s40168-023-01747-z}
#'
#' Fawcett, T. (2006). An introduction to ROC analysis. *Pattern
#' Recognition Letters*, 27, 861-874. \doi{10.1016/j.patrec.2005.10.010}
#'
#' @seealso [predict.ecodive_random_forest()], [plot.ecodive_random_forest()],
#'   [plsda()], [lda()]
#'
#' @export
#' @examples
#'     # Which body site? 20 samples from 4 sites, 5-fold cross-validation.
#'     counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
#'     site   <- rep(rownames(ex_counts), each = 5)
#'     fit    <- random_forest(counts, site, n_trees = 50, cpus = 1)
#'     fit
#'
#'     # Predict new samples.
#'     predict(fit, ex_counts, cpus = 1)
random_forest <- function (
    counts, y, data = NULL, margin = 1L, norm = 'percent',
    pseudocount = NULL, n_trees = 500L, max_features = NULL,
    min_samples_leaf = 1L, min_samples_split = 2L, max_depth = NULL,
    criterion = NULL, class_weight = NULL, max_samples = NULL,
    importance = 'impurity', n_folds = 5L, folds = NULL,
    keep_forest = TRUE, seed = 0L, cpus = n_cpus() ) {

  if (missing(y)) y <- NULL
  rf_run(environment())
}


# The same, with the internals the tests need:
#   numpy_random_state   NULL, or an integer: draw as scikit-learn does with
#                        that random_state (numpy's MT19937 and our_rand_r())
#                        instead of PCG. Every training set uses the same
#                        random_state, as q2-sample-classifier's refits do.
#   sparse               NULL (sparse for sparse tables, unless norm is clr
#                        or log), TRUE or FALSE.
# Returns the result with `forests`: each training set's chosen forest (the
# full data's first), with its inbag counts and node impurities.
random_forest_details <- function (
    counts, y, data = NULL, margin = 1L, norm = 'percent',
    pseudocount = NULL, n_trees = 500L, max_features = NULL,
    min_samples_leaf = 1L, min_samples_split = 2L, max_depth = NULL,
    criterion = NULL, class_weight = NULL, max_samples = NULL,
    importance = 'impurity', n_folds = 5L, folds = NULL,
    keep_forest = TRUE, seed = 0L, cpus = n_cpus(),
    numpy_random_state = NULL, sparse = NULL ) {

  if (missing(y)) y <- NULL
  rf_run(environment(), numpy_random_state = numpy_random_state,
         sparse = sparse, details = TRUE)
}



# train: NULL, or the 1-based positions of the samples to train on, in that
# order (feature_volatility()'s split); the others' predictions by the
# full-data forest are then returned as `holdout`, list(index, value).
rf_run <- function (env, numpy_random_state = NULL, sparse = NULL, details = FALSE,
                    train = NULL) {

  with(env, {
    validate_counts()
    validate_norm()
    validate_pseudocount()
    validate_margin()
    validate_cpus()
    validate_seed()
  })
  counts <- env$counts; margin <- env$margin; norm <- env$norm
  pseudocount <- env$pseudocount; cpus <- env$cpus; seed <- env$seed

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  pkg <- get_matrix_package(counts)
  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }
  n <- as.integer(dims[[1]])
  p <- as.integer(dims[[2]])
  if (n < 2) stop('`counts` must have at least two samples.', call. = FALSE)
  if (p < 1) stop('`counts` must have at least one feature.', call. = FALSE)
  feats <- if (is.null(dn[[2]])) sprintf('F%d', seq_len(p)) else dn[[2]]
  sams  <- dn[[1]]


  # y and folds, matched to the samples.
  data <- env$data
  idx  <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(sams, n, rn, nrow(data), ref = '`counts`')
  }
  if (is.null(env$y))
    stop('`y` is required: one value per sample, or a column name of `data`.', call. = FALSE)
  yv <- per_sample(env$y, 'y', data, idx, sams, n, ref = '`counts`')
  fv <- per_sample(env$folds, 'folds', data, idx, sams, n, ref = '`counts`')

  # The training samples: those with a `y`, or (feature_volatility()) the
  # given ones, in the given order.
  if (is.null(train)) {
    keep <- !is.na(yv)
    if (!all(keep)) {
      warning(sprintf(
        'Left out %d sample%s whose `y` is NA.', sum(!keep), if (sum(!keep) == 1) '' else 's'),
        call. = FALSE )
      if (sum(keep) < 2) stop('Fewer than two samples have a `y`.', call. = FALSE)
    }
    keep <- which(keep)
  }
  else {
    keep <- as.integer(train)
  }
  rows <- keep - 1L
  yv   <- yv[keep]
  if (!is.null(fv)) fv <- fv[keep]
  sam_out <- if (is.null(sams)) as.character(keep) else sams[keep]
  nk <- length(rows)

  if (is.factor(yv) || is.character(yv) || is.logical(yv)) {
    type <- 'classification'
    if (is.factor(yv)) {
      if (anyNA(levels(yv))) yv <- factor(as.character(yv))
      unused <- setdiff(levels(yv), unique(as.character(yv)))
      if (length(unused))
        warning('Dropped unused `y` level', if (length(unused) > 1) 's', ': ',
                name_list(unused), '.', call. = FALSE)
      yv <- droplevels(yv)
    }
    else if (is.logical(yv)) yv <- factor(yv, levels = c(FALSE, TRUE))[, drop = TRUE]
    else yv <- factor(yv, levels = sort(unique(yv), method = 'radix'))
    classes <- levels(yv)
    k <- length(classes)
    if (k < 2) stop('`y` must have at least two classes.', call. = FALSE)
    ycode <- as.double(as.integer(yv) - 1L)
  }
  else if (is.numeric(yv)) {
    type <- 'regression'
    if (!all(is.finite(yv))) stop('`y` must be finite.', call. = FALSE)
    classes <- NULL
    k <- 0L
    ycode <- as.double(yv)
  }
  else {
    stop('`y` must be a factor, character, logical or numeric vector.', call. = FALSE)
  }
  cls <- type == 'classification'


  # Settings.
  set <- rf_settings(env, p, cls, classes)
  imp_type <- match.arg(env$importance, c('impurity', 'permutation', 'both'))
  keep_forest <- env$keep_forest
  if (!is.logical(keep_forest) || length(keep_forest) != 1 || is.na(keep_forest))
    stop('`keep_forest` must be TRUE or FALSE.', call. = FALSE)


  # The folds.
  if (!is.null(fv)) {
    if (anyNA(fv)) stop('`folds` must not contain NAs.', call. = FALSE)
    ff <- if (is.factor(fv)) fv else factor(fv, levels = sort(unique(fv), method = 'radix'))
    if (any(tabulate(as.integer(ff), nlevels(ff)) == 0))
      stop('Every fold in `folds` must have samples.', call. = FALSE)
    fold <- as.integer(ff)
    n_folds <- nlevels(ff)
    if (n_folds < 2) stop('`folds` must have at least two folds.', call. = FALSE)
  }
  else {
    n_folds <- env$n_folds
    if (!is.numeric(n_folds) || length(n_folds) != 1 || is.na(n_folds) ||
        n_folds < 0 || n_folds != floor(n_folds))
      stop('`n_folds` must be a whole number, at least 0.', call. = FALSE)
    n_folds <- as.integer(n_folds)
    if (n_folds == 1L) n_folds <- 0L
    if (n_folds > nk)
      stop(sprintf('`n_folds` (%d) must be at most the number of samples (%d).', n_folds, nk),
           call. = FALSE)
    fold <- if (n_folds >= 2) rf_folds(yv, cls, n_folds, seed) else NULL
  }
  if (!is.null(fold)) for (f in seq_len(n_folds)) {
    tr <- fold != f
    if (cls && length(unique(yv[tr])) < 2)
      stop(sprintf('The training samples of fold %d have only one class.', f), call. = FALSE)
    if (sum(tr) < 2)
      stop(sprintf('Fold %d leaves fewer than two training samples.', f), call. = FALSE) # nocov
  }


  # The normalized table: dense samples x features, or CSC.
  nm_norm <- names(which(c(none = NORM_NONE, percent = NORM_PERCENT, chord = NORM_CHORD,
    hellinger = NORM_HELLINGER, binary = NORM_BINARY, clr = NORM_CLR, rclr = NORM_RCLR,
    log = NORM_LOG) == norm))
  clr <- identical(norm, NORM_CLR) || identical(norm, NORM_LOG)
  use_sparse <- if (is.null(sparse)) !clr && pkg != 'base' && !inherits(counts, 'dgeMatrix')
                else isTRUE(sparse)
  if (clr && use_sparse) use_sparse <- FALSE
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  x <- rf_table(counts, margin, norm, pseudocount, use_sparse, n, p, cpus)


  # numpy mode's tree seeds.
  np_seeds <- NULL
  if (!is.null(numpy_random_state))
    np_seeds <- as.integer(.Call(C_rf_rng_test, 3L, as.double(numpy_random_state), 0, set$n_trees))


  # Grow the forests: the full data (set 0), then each fold's training set.
  grid  <- set$grid
  sets  <- c(list(rows), if (n_folds >= 2) lapply(seq_len(n_folds), function (f) rows[fold != f]))
  pos   <- c(list(seq_len(nk)), if (n_folds >= 2) lapply(seq_len(n_folds), function (f) which(fold != f)))
  tuning  <- NULL
  forests <- vector('list', length(sets))
  chosen  <- integer(length(sets))

  for (s in seq_along(sets)) {
    best <- NULL; best_err <- Inf
    for (g in seq_len(nrow(grid))) {
      fr  <- rf_grow(x, sets[[s]], ycode[pos[[s]]], k, set, grid[g, ], seed,
                     (s - 1) * 2^32, np_seeds, cpus, classes)
      err <- NA_real_
      if (nrow(grid) > 1 || s == 1) {
        oob <- rf_oob(fr, x, sets[[s]], ycode[pos[[s]]], k, cpus)
        fr$oob <- oob
        err <- oob$error
        if (nrow(grid) > 1)
          tuning <- rbind(tuning, data.frame(
            training_set = if (s == 1) 'full' else paste('fold', s - 1),
            max_features = grid$max_features[[g]], min_samples_leaf = grid$min_samples_leaf[[g]],
            oob_error = err, chosen = FALSE, stringsAsFactors = FALSE ))
      }
      if (is.null(best) || (!is.na(err) && err < best_err)) {
        best <- fr; best_err <- if (is.na(err)) Inf else err; chosen[[s]] <- g
      }
    }
    forests[[s]] <- best
  }
  if (!is.null(tuning))
    tuning$chosen <- unlist(lapply(seq_along(sets), function (s) seq_len(nrow(grid)) == chosen[[s]]))

  full <- forests[[1]]
  oob  <- full$oob
  if (any(oob$n_trees == 0))
    warning(sprintf(
      '%d sample%s never out of bag: %s OOB prediction%s NA. Use more trees.',
      sum(oob$n_trees == 0), if (sum(oob$n_trees == 0) == 1) ' was' else 's were',
      if (sum(oob$n_trees == 0) == 1) 'its' else 'their', if (sum(oob$n_trees == 0) == 1) ' is' else 's are'),
      call. = FALSE )


  # Cross-validated predictions and fold importances.
  kk   <- max(k, 1L)
  cvv  <- matrix(NA_real_, kk, nk)
  imp_cv <- NULL
  if (n_folds >= 2) {
    imp_cv <- matrix(NA_real_, p, n_folds)
    for (f in seq_len(n_folds)) {
      te <- which(fold == f)
      pr <- .Call(C_rf_predict, forests[[f + 1]]$forest, x, rows[te], FALSE, cpus)
      cvv[, te] <- pr$value
      imp_cv[, f] <- .Call(C_rf_importance, forests[[f + 1]]$forest, p)
    }
  }


  # Importances.
  imp_full <- .Call(C_rf_importance, full$forest, p)
  imp <- data.frame(
    feature          = feats,
    impurity         = imp_full,
    impurity_cv_mean = if (is.null(imp_cv)) NA_real_ else rowMeans(imp_cv),
    impurity_cv_sd   = if (is.null(imp_cv)) NA_real_ else
                         sqrt(rowMeans((imp_cv - rowMeans(imp_cv))^2)),
    stringsAsFactors = FALSE )
  if (imp_type != 'impurity') {
    pm <- .Call(C_rf_permutation, full$forest, x, rows, ycode, k, seed, cpus)
    imp$permutation    <- pm$mean
    imp$permutation_sd <- pm$sd
  }
  ord <- order(-(if (imp_type == 'impurity') imp$impurity else imp$permutation), seq_len(p))
  imp <- imp[ord, , drop = FALSE]
  rownames(imp) <- NULL


  # Predictions and performance.
  res <- rf_summaries(cls, classes, yv, ycode, cvv, fold, n_folds, oob, sam_out)

  res <- c(list(importance = imp), res)
  if (!is.null(tuning)) res$tuning <- tuning
  if (keep_forest) {
    fo <- full$forest
    fo[c('impurity', 'n_node_samples', 'weighted_n_node_samples', 'inbag')] <- list(NULL)
    res$forest <- list(
      trees = fo, classes = classes, features = feats, norm = nm_norm,
      pseudocount = pseudocount, type = type )
  }
  if (!is.null(train)) {
    te <- setdiff(seq_len(n), keep)
    res$holdout <- list(
      index = te,
      value = if (length(te)) as.vector(.Call(C_rf_predict, full$forest, x, te - 1L, FALSE, cpus)$value)
              else numeric(0) )
  }
  if (details) {
    res$forests <- lapply(forests, function (fr) fr$forest)
    res$folds   <- fold
    res$oob_value <- oob$value
  }

  g1 <- grid[chosen[[1]], ]
  structure(
    res,
    class            = 'ecodive_random_forest',
    method           = 'random_forest',
    type             = type,
    n_trees          = set$n_trees,
    max_features     = g1$max_features,
    min_samples_leaf = g1$min_samples_leaf,
    criterion        = set$criterion,
    class_weight     = set$class_weight_label,
    n_folds          = n_folds )
}



# Checks and resolves the forest's settings. Returns list(n_trees, grid
# (data frame of max_features and min_samples_leaf counts), max_features
# (the values as given, per grid row), min_samples_split, max_depth (NA:
# none), criterion, crit_code, class_weight, class_weight_label,
# balanced_sub, max_samples).
rf_settings <- function (env, p, cls, classes) {

  whole <- function (v, arg, lo) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < lo || v != floor(v) || !is.finite(v))
      stop('`', arg, '` must be a whole number, at least ', lo, '.', call. = FALSE)
    as.integer(v)
  }

  n_trees <- whole(env$n_trees, 'n_trees', 1)
  min_samples_split <- whole(env$min_samples_split, 'min_samples_split', 2)
  max_depth <- if (is.null(env$max_depth)) NA_integer_ else whole(env$max_depth, 'max_depth', 1)

  # max_features: counts, fractions, or keywords.
  mf <- env$max_features
  if (is.null(mf)) mf <- if (cls) 'sqrt' else 'all'
  if (!(is.character(mf) || is.numeric(mf)) || length(mf) == 0 || anyNA(mf))
    stop("`max_features` must be counts, fractions, 'sqrt', 'log2' or 'all'.", call. = FALSE)
  mf_count <- vapply(as.list(mf), function (v) {
    if (is.character(v)) {
      if (v %in% c('sqrt', 'log2', 'all'))
        return (switch(v,
          sqrt = max(1L, as.integer(floor(sqrt(p)))),
          log2 = max(1L, as.integer(floor(log2(p)))),
          all  = p ))
      v <- suppressWarnings(as.numeric(v))
      if (is.na(v))
        stop("`max_features` must be counts, fractions, 'sqrt', 'log2' or 'all'.", call. = FALSE)
    }
    if (!is.finite(v) || v <= 0)
      stop('`max_features` must be positive.', call. = FALSE)
    if (v < 1) return (max(1L, as.integer(floor(v * p))))
    if (v != floor(v))
      stop('`max_features` must be a whole number (a count), or a fraction below 1.', call. = FALSE)
    if (v > p)
      stop(sprintf('`max_features` (%d) must be at most the number of features (%d).', as.integer(v), p),
           call. = FALSE)
    as.integer(v)
  }, 1L)

  msl <- env$min_samples_leaf
  if (!is.numeric(msl) || length(msl) == 0 || anyNA(msl) || any(msl < 1) ||
      any(msl != floor(msl)) || any(!is.finite(msl)))
    stop('`min_samples_leaf` must be whole numbers, at least 1.', call. = FALSE)
  msl <- as.integer(msl)

  grid <- expand.grid(
    max_features = unique(mf_count), min_samples_leaf = unique(msl),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE )

  # criterion
  crit <- env$criterion
  if (is.null(crit)) crit <- if (cls) 'gini' else 'squared_error'
  if (!is.character(crit) || length(crit) != 1 || is.na(crit))
    stop('`criterion` must be one string.', call. = FALSE)
  if (cls) {
    if (!crit %in% c('gini', 'entropy'))
      stop("`criterion` must be 'gini' or 'entropy' for classification.", call. = FALSE)
  }
  else {
    if (crit == 'mse') crit <- 'squared_error'
    if (crit != 'squared_error')
      stop("`criterion` must be 'squared_error' for regression.", call. = FALSE)
  }
  crit_code <- switch(crit, gini = 0L, entropy = 1L, squared_error = 2L)

  # class_weight
  cw <- env$class_weight
  cw_label <- NULL
  balanced_sub <- FALSE
  if (!is.null(cw)) {
    if (!cls) stop('`class_weight` is for classification only.', call. = FALSE)
    if (is.character(cw) && length(cw) == 1 && cw %in% c('balanced', 'balanced_subsample')) {
      cw_label <- cw
      balanced_sub <- cw == 'balanced_subsample'
    }
    else if (is.numeric(cw) && !is.null(names(cw)) && length(cw) > 0 &&
             all(is.finite(cw)) && all(cw > 0)) {
      bad <- setdiff(names(cw), classes)
      if (length(bad))
        stop('`class_weight` names classes that `y` lacks: ', name_list(bad), '.', call. = FALSE)
      cw_label <- cw
    }
    else stop("`class_weight` must be 'balanced', 'balanced_subsample', ",
              'or a named vector of positive weights.', call. = FALSE)
  }

  # max_samples
  ms <- env$max_samples
  if (!is.null(ms)) {
    if (!is.numeric(ms) || length(ms) != 1 || is.na(ms) || !is.finite(ms) || ms <= 0 ||
        (ms >= 1 && ms != floor(ms)))
      stop('`max_samples` must be a count, or a fraction below 1.', call. = FALSE)
    ms <- if (ms >= 1) c(1, ms) else c(2, ms)
  }

  list(
    n_trees = n_trees, grid = grid, min_samples_split = min_samples_split,
    max_depth = max_depth, criterion = crit, crit_code = crit_code,
    class_weight = cw, class_weight_label = cw_label, balanced_sub = balanced_sub,
    max_samples = ms )
}



# The normalized table, as C_rf_grow() reads it: a dense samples x
# features double matrix, or list(i (0-based rows), p (column pointers),
# x, n) for a CSC table.
rf_table <- function (counts, margin, norm, pseudocount, sparse, n, p, cpus) {

  if (!sparse) {
    x <- .Call(C_normalize, counts, margin, norm, as.double(pseudocount), TRUE, cpus)
    if (margin == 2L) x <- t(x)
    storage.mode(x) <- 'double'
    dimnames(x) <- NULL
    if (!all(is.finite(x)))
      stop('The normalized `counts` contain NA, NaN, or infinite values.', call. = FALSE) # nocov
    return (x)
  }

  tr <- .Call(C_normalize, counts, margin, norm, as.double(pseudocount), FALSE, cpus)
  if (margin == 2L) {
    o <- order(tr$i, tr$j, method = 'radix')
    ri <- tr$j[o]; rj <- tr$i[o]; rx <- tr$x[o]
  }
  else { ri <- tr$i; rj <- tr$j; rx <- tr$x }
  if (!all(is.finite(rx)))
    stop('The normalized `counts` contain NA, NaN, or infinite values.', call. = FALSE) # nocov
  list(i = as.integer(ri), p = c(0L, cumsum(tabulate(rj + 1L, nbins = p))),
       x = as.double(rx), n = as.integer(n))
}



# StratifiedKFold (classification) or KFold (regression), shuffled with
# PCG streams 2^62 + class (or 2^62). Returns each sample's fold, 1-based.
rf_folds <- function (yv, cls, n_folds, seed) {

  n <- length(yv)
  if (!cls) {
    idx  <- .Call(C_rf_shuffle, seq_len(n), seed, 2^62)
    size <- rep(n %/% n_folds, n_folds)
    size[seq_len(n %% n_folds)] <- size[seq_len(n %% n_folds)] + 1L
    fold <- integer(n)
    fold[idx] <- rep(seq_len(n_folds), size)
    return (fold)
  }

  # Classes encoded in order of first appearance.
  enc <- match(as.character(yv), unique(as.character(yv)))
  kc  <- max(enc)
  cnt <- tabulate(enc, kc)
  if (all(n_folds > cnt))
    stop(sprintf('`n_folds` (%d) is larger than every class.', n_folds), call. = FALSE)
  if (n_folds > min(cnt))
    warning(sprintf(
      'The smallest class in `y` has only %d sample%s, fewer than `n_folds` (%d).',
      min(cnt), if (min(cnt) == 1) '' else 's', n_folds), call. = FALSE)
  y_order <- sort(enc)
  alloc <- t(vapply(seq_len(n_folds), function (i)
    tabulate(y_order[seq(i, n, by = n_folds)], kc), integer(kc)))
  if (kc == 1) alloc <- matrix(alloc, ncol = 1) # nocov
  fold <- integer(n)
  for (c in seq_len(kc)) {
    ff <- rep(seq_len(n_folds), alloc[, c])
    ff <- .Call(C_rf_shuffle, ff, seed, 2^62 + c - 1)
    fold[enc == c] <- ff
  }
  fold
}



# Grows one forest. Returns list(forest, rows).
rf_grow <- function (x, rows, y, k, set, g, seed, stream_base, np_seeds, cpus, classes) {

  cw <- NULL
  if (!is.null(set$class_weight) && !set$balanced_sub) {
    if (is.character(set$class_weight)) {
      # 'balanced': n / (k * n_c) over the classes present.
      cnt <- tabulate(as.integer(y) + 1L, k)
      kp  <- sum(cnt > 0)
      w   <- length(y) / (kp * cnt)
      cw  <- w[as.integer(y) + 1L]
    }
    else {
      w <- rep(1, k)
      w[match(names(set$class_weight), classes)] <- set$class_weight
      cw <- w[as.integer(y) + 1L]
    }
  }

  params <- c(
    set$n_trees, g$max_features, g$min_samples_leaf, set$min_samples_split,
    set$max_depth, set$crit_code, as.integer(set$balanced_sub) )

  forest <- .Call(
    C_rf_grow, x, as.integer(rows), as.double(y), as.integer(k),
    as.integer(params), if (is.null(set$max_samples)) NULL else as.double(set$max_samples),
    if (is.null(cw)) NULL else as.double(cw), seed, as.double(stream_base),
    np_seeds, cpus )

  list(forest = forest, rows = rows)
}



# Out-of-bag predictions of a forest's training rows. Returns list(value (k
# x n), n_trees, pred (class code or value; NA where never OOB), error
# (misclassification rate or MSE over the rows with OOB predictions),
# score (accuracy or R^2)).
rf_oob <- function (fr, x, rows, y, k, cpus) {
  pr <- .Call(C_rf_predict, fr$forest, x, as.integer(rows), TRUE, cpus)
  ok <- pr$n_trees > 0
  v  <- pr$value
  v[, !ok] <- NA_real_
  if (k > 0) {
    pred <- rep(NA_real_, length(ok))
    pred[ok] <- max.col(t(v[, ok, drop = FALSE]), ties.method = 'first') - 1
    err   <- if (any(ok)) mean(pred[ok] != y[ok]) else NA_real_
    score <- 1 - err
  }
  else {
    pred  <- as.vector(v)
    err   <- if (any(ok)) mean((pred[ok] - y[ok])^2) else NA_real_
    score <- if (any(ok)) rf_r2(y[ok], pred[ok]) else NA_real_
  }
  list(value = v, n_trees = pr$n_trees, pred = pred, error = err, score = score)
}


# sklearn's r2_score(): 1 - SS_res / SS_tot.
rf_r2 <- function (obs, pred) {
  ss_tot <- sum((obs - mean(obs))^2)
  ss_res <- sum((obs - pred)^2)
  if (ss_tot == 0) return (if (ss_res == 0) 1 else 0)
  1 - ss_res / ss_tot
}



# The predictions data frame, performance, confusion and ROC.
rf_summaries <- function (cls, classes, yv, ycode, cvv, fold, n_folds, oob, sam_out) {

  nk  <- length(yv)
  has_cv <- n_folds >= 2

  if (cls) {
    k <- length(classes)
    cv_code <- if (has_cv) max.col(t(cvv), ties.method = 'first') else rep(NA_integer_, nk)
    pred_cv <- factor(classes[cv_code], levels = classes)
    oob_pred <- factor(classes[oob$pred + 1], levels = classes)
    pr <- data.frame(
      sample = sam_out, observed = yv, predicted = pred_cv,
      fold = if (has_cv) fold else NA_integer_, oob_predicted = oob_pred,
      stringsAsFactors = FALSE )
    probs <- t(cvv)
    colnames(probs) <- paste0('prob_', classes)
    pr <- cbind(pr, as.data.frame(probs, stringsAsFactors = FALSE))

    # What the confusion matrix and ROC use: CV, or else OOB.
    use_p    <- if (has_cv) t(cvv) else t(oob$value)
    use_pred <- if (has_cv) pred_cv else oob_pred
    ok <- !is.na(use_pred)

    acc_fold <- if (has_cv) vapply(seq_len(n_folds), function (f)
      mean(pred_cv[fold == f] == yv[fold == f]), 0) else NA_real_
    acc  <- if (has_cv) mean(pred_cv == yv) else NA_real_
    bacc <- if (has_cv) mean(vapply(classes, function (c) mean(pred_cv[yv == c] == c), 0)) else NA_real_
    base <- max(tabulate(as.integer(yv), k)) / nk

    roc  <- rf_roc(ycode[ok], use_p[ok, , drop = FALSE], classes)
    aucs <- attr(roc, 'auc')
    attr(roc, 'auc') <- NULL

    perf <- c(
      accuracy           = acc,
      accuracy_fold_mean = if (has_cv) mean(acc_fold) else NA_real_,
      accuracy_fold_sd   = if (has_cv) sqrt(mean((acc_fold - mean(acc_fold))^2)) else NA_real_,
      balanced_accuracy  = bacc,
      baseline_accuracy  = base,
      accuracy_ratio     = acc / base )
    if (k == 2) perf <- c(perf, auc = if (has_cv) aucs[[classes[[2]]]] else NA_real_)
    else perf <- c(perf, auc_macro = if (has_cv) aucs[['macro']] else NA_real_,
                         auc_micro = if (has_cv) aucs[['micro']] else NA_real_)
    perf <- c(perf, oob_accuracy = oob$score)

    conf <- table(
      observed = factor(yv[ok], levels = classes),
      predicted = factor(use_pred[ok], levels = classes) )
    conf <- matrix(as.integer(conf), k, k, dimnames = list(observed = classes, predicted = classes))

    return (list(predictions = pr, performance = perf, confusion = conf, roc = roc))
  }

  pred_cv <- if (has_cv) as.vector(cvv) else rep(NA_real_, nk)
  pr <- data.frame(
    sample = sam_out, observed = yv, predicted = pred_cv,
    fold = if (has_cv) fold else NA_integer_, oob_predicted = oob$pred,
    stringsAsFactors = FALSE )
  mse_fold <- if (has_cv) vapply(seq_len(n_folds), function (f)
    mean((pred_cv[fold == f] - yv[fold == f])^2), 0) else NA_real_
  lr <- if (has_cv) rf_linregress(yv, pred_cv) else c(NA_real_, NA_real_, NA_real_)
  perf <- c(
    mse           = if (has_cv) mean((pred_cv - yv)^2) else NA_real_,
    r_squared     = if (has_cv) rf_r2(yv, pred_cv) else NA_real_,
    mse_fold_mean = if (has_cv) mean(mse_fold) else NA_real_,
    mse_fold_sd   = if (has_cv) sqrt(mean((mse_fold - mean(mse_fold))^2)) else NA_real_,
    slope         = lr[[1]],
    intercept     = lr[[2]],
    r             = lr[[3]],
    oob_r_squared = oob$score )
  list(predictions = pr, performance = perf)
}


# scipy.stats.linregress(observed, predicted): slope, intercept, r.
rf_linregress <- function (x, y) {
  mx <- mean(x); my <- mean(y)
  ssxm <- mean((x - mx)^2); ssym <- mean((y - my)^2); ssxym <- mean((x - mx) * (y - my))
  slope <- if (ssxm == 0) NA_real_ else ssxym / ssxm
  r <- if (ssxm == 0 || ssym == 0) 0 else max(-1, min(1, ssxym / sqrt(ssxm * ssym)))
  c(slope, my - slope * mx, r)
}



# scikit-learn's roc_curve(y_true, y_score) (drop_intermediate = TRUE):
# returns list(fpr, tpr, threshold).
rf_roc_curve <- function (truth, score) {
  o  <- order(score, decreasing = TRUE, method = 'radix')
  s  <- score[o]
  t  <- truth[o]
  n  <- length(s)
  ix <- c(which(diff(s) != 0), n)
  tps <- cumsum(t)[ix]
  fps <- ix - tps
  thr <- s[ix]
  if (length(fps) > 2) {
    keep <- c(TRUE, diff(fps, differences = 2) != 0 | diff(tps, differences = 2) != 0, TRUE)
    fps <- fps[keep]; tps <- tps[keep]; thr <- thr[keep]
  }
  tps <- c(0, tps); fps <- c(0, fps); thr <- c(Inf, thr)
  fpr <- if (fps[[length(fps)]] <= 0) rep(NA_real_, length(fps)) else fps / fps[[length(fps)]]
  tpr <- if (tps[[length(tps)]] <= 0) rep(NA_real_, length(tps)) else tps / tps[[length(tps)]]
  list(fpr = fpr, tpr = tpr, threshold = thr)
}


# The trapezoid rule, as sklearn's auc().
rf_auc <- function (x, y) {
  if (length(x) < 2 || anyNA(x) || anyNA(y)) return (NA_real_)
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}


# numpy.interp(x, xp, fp) for increasing xp (with ties).
rf_interp <- function (x, xp, fp) {
  n <- length(xp)
  j <- findInterval(x, xp)
  out <- numeric(length(x))
  for (a in seq_along(x)) {
    if (x[[a]] >= xp[[n]]) out[[a]] <- fp[[n]]
    else if (j[[a]] == 0) out[[a]] <- fp[[1]]
    else {
      jj <- j[[a]]
      slope <- (fp[[jj + 1]] - fp[[jj]]) / (xp[[jj + 1]] - xp[[jj]])
      out[[a]] <- slope * (x[[a]] - xp[[jj]]) + fp[[jj]]
    }
  }
  out
}


# ROC curves of class probabilities: each class against the rest (the
# second class only, with two), and with more than two classes q2's micro
# and macro averages. Returns a data frame (class, fpr, tpr, threshold) with
# attribute `auc` (named by class, 'micro', 'macro').
rf_roc <- function (ycode, probs, classes) {
  k <- length(classes)
  if (nrow(probs) == 0) {
    df <- data.frame(class = character(0), fpr = numeric(0), tpr = numeric(0),
                     threshold = numeric(0), stringsAsFactors = FALSE)
    attr(df, 'auc') <- stats::setNames(rep(NA_real_, k), classes)
    return (df)
  }
  bin <- outer(ycode, seq_len(k) - 1, `==`) + 0
  curves <- lapply(seq_len(k), function (c) rf_roc_curve(bin[, c], probs[, c]))
  names(curves) <- classes
  aucs <- vapply(curves, function (r) rf_auc(r$fpr, r$tpr), 0)

  if (k > 2) {
    micro <- rf_roc_curve(as.vector(t(bin)), as.vector(t(probs)))
    all_fpr <- sort(unique(unlist(lapply(curves, `[[`, 'fpr'))))
    mean_tpr <- numeric(length(all_fpr))
    for (c in seq_len(k)) mean_tpr <- mean_tpr + rf_interp(all_fpr, curves[[c]]$fpr, curves[[c]]$tpr)
    mean_tpr <- mean_tpr / k
    macro <- list(fpr = all_fpr, tpr = mean_tpr, threshold = rep(NA_real_, length(all_fpr)))
    curves <- c(curves, list(micro = micro, macro = macro))
    aucs <- c(aucs, micro = rf_auc(micro$fpr, micro$tpr), macro = rf_auc(all_fpr, mean_tpr))
  }
  else curves <- curves[2]

  df <- do.call(rbind, lapply(names(curves), function (nm) data.frame(
    class = nm, fpr = curves[[nm]]$fpr, tpr = curves[[nm]]$tpr,
    threshold = curves[[nm]]$threshold, stringsAsFactors = FALSE )))
  attr(df, 'auc') <- aucs
  df
}



#' Predict with a random forest
#'
#' Predicts the class (or class probabilities), or the value, of new
#' samples from the full-data forest of [random_forest()].
#'
#' Features are matched by name. Features of the forest missing from
#' `counts` are taken as 0, with a warning saying how many (count tables
#' often lack rare features); features the forest doesn't know are
#' dropped. The table is then normalized as the forest's was (same `norm`
#' and `pseudocount`), over the forest's features only.
#'
#' @param object   A [random_forest()] result, grown with `keep_forest =
#'   TRUE`.
#' @param counts   The new samples' feature table, of any type
#'   [random_forest()] takes.
#' @param type   For classification: `'class'` (the predicted classes, a
#'   factor) or `'prob'` (a samples x classes matrix of probabilities).
#'   Ignored for regression, which returns numbers. Default: `'class'`
#' @param margin   As for [random_forest()]: `1` if samples are rows of
#'   `counts`, `2` if columns. Default: `1`
#' @param cpus   The number of CPU threads. Default: `n_cpus()`
#' @param ...   Not used.
#'
#' @return A factor of classes, a matrix of class probabilities, or a
#'   numeric vector, named by sample.
#'
#' @seealso [random_forest()]
#'
#' @export
#' @examples
#'     counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
#'     site   <- rep(rownames(ex_counts), each = 5)
#'     fit    <- random_forest(counts, site, n_trees = 50, n_folds = 0, cpus = 1)
#'     predict(fit, ex_counts, type = 'prob', cpus = 1)
predict.ecodive_random_forest <- function (
    object, counts, type = c('class', 'prob'), margin = 1L, cpus = n_cpus(), ... ) {

  if (is.null(object$forest))
    stop('The forest was not kept: grow it with `keep_forest = TRUE`.', call. = FALSE)
  if (missing(counts)) stop('`counts` is required.', call. = FALSE)
  type <- match.arg(type)

  validate_counts()
  validate_margin()
  validate_cpus()

  fo  <- object$forest
  pkg <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }
  n <- as.integer(dims[[1]])
  feats <- fo$features
  p <- length(feats)

  if (is.null(dn[[2]])) {
    if (dims[[2]] != p)
      stop(sprintf('`counts` has no feature names, and %d features where the forest has %d.',
                   dims[[2]], p), call. = FALSE)
    col <- seq_len(p)
  }
  else col <- match(dn[[2]], feats)
  n_missing <- p - length(unique(col[!is.na(col)]))
  if (n_missing == p)
    stop('`counts` has none of the forest\'s features.', call. = FALSE)
  if (n_missing > 0)
    warning(sprintf('%d of the forest\'s %d features %s not in `counts`; taken as 0.',
                    n_missing, p, if (n_missing == 1) 'is' else 'are'), call. = FALSE)

  # The forest's features only (as triplets), then normalized as before.
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  tr <- .Call(C_normalize, counts, margin, NORM_NONE, NA_real_, FALSE, cpus)
  ri <- if (margin == 2L) tr$j else tr$i
  rj <- if (margin == 2L) tr$i else tr$j
  cj <- col[rj + 1L]
  ok <- !is.na(cj)
  stm <- structure(
    list(i = as.integer(ri[ok] + 1L), j = as.integer(cj[ok]), v = as.double(tr$x[ok]),
         nrow = n, ncol = p, dimnames = list(dn[[1]], feats)),
    class = 'simple_triplet_matrix' )
  norm <- switch(fo$norm,
    none = NORM_NONE, percent = NORM_PERCENT, chord = NORM_CHORD,
    hellinger = NORM_HELLINGER, binary = NORM_BINARY, clr = NORM_CLR,
    rclr = NORM_RCLR, log = NORM_LOG )
  clr <- fo$norm %in% c('clr', 'log')
  x <- rf_table(stm, 1L, norm, fo$pseudocount, !clr, n, p, cpus)

  pr <- .Call(C_rf_predict, fo$trees, x, seq_len(n) - 1L, FALSE, cpus)
  nms <- dn[[1]]

  if (fo$type == 'regression') {
    out <- as.vector(pr$value)
    names(out) <- nms
    return (out)
  }
  probs <- t(pr$value)
  dimnames(probs) <- list(nms, fo$classes)
  if (type == 'prob') return (probs)
  out <- factor(fo$classes[max.col(probs, ties.method = 'first')], levels = fo$classes)
  names(out) <- nms
  out
}



#' @export
print.ecodive_random_forest <- function (x, n = 10L, ...) {

  type <- attr(x, 'type', exact = TRUE)
  cls  <- type == 'classification'
  pr   <- x$predictions
  perf <- x$performance
  nf   <- attr(x, 'n_folds', exact = TRUE)
  imp  <- x$importance
  fmt  <- function (v) if (is.na(v)) 'NA' else formatC(v, digits = 3, format = 'fg', flag = '#')

  cat(sprintf('# Random forest %s: %d samples x %d features\n',
              type, nrow(pr), nrow(imp)))
  cw <- attr(x, 'class_weight', exact = TRUE)
  cat(sprintf('# %d trees; max_features = %d, min_samples_leaf = %d, criterion = %s%s\n',
              attr(x, 'n_trees', exact = TRUE), attr(x, 'max_features', exact = TRUE),
              attr(x, 'min_samples_leaf', exact = TRUE), attr(x, 'criterion', exact = TRUE),
              if (is.null(cw)) '' else if (is.character(cw)) paste0(', class_weight = ', cw)
              else ', class weights given'))
  if (cls) cat(sprintf('# Classes: %s\n', paste(levels(pr$observed), collapse = ', ')))

  if (cls) {
    if (nf >= 2) {
      cat(sprintf('# %d-fold CV accuracy: %s (folds %s +/- %s; baseline %s, ratio %s)\n',
                  nf, fmt(perf[['accuracy']]), fmt(perf[['accuracy_fold_mean']]),
                  fmt(perf[['accuracy_fold_sd']]), fmt(perf[['baseline_accuracy']]),
                  fmt(perf[['accuracy_ratio']])))
      if ('auc' %in% names(perf)) cat(sprintf('# CV AUC: %s\n', fmt(perf[['auc']])))
      else cat(sprintf('# CV AUC: macro %s, micro %s\n', fmt(perf[['auc_macro']]), fmt(perf[['auc_micro']])))
    }
    cat(sprintf('# Out-of-bag accuracy: %s\n', fmt(perf[['oob_accuracy']])))
  }
  else {
    if (nf >= 2)
      cat(sprintf('# %d-fold CV R^2: %s, MSE: %s (folds %s +/- %s)\n',
                  nf, fmt(perf[['r_squared']]), fmt(perf[['mse']]),
                  fmt(perf[['mse_fold_mean']]), fmt(perf[['mse_fold_sd']])))
    cat(sprintf('# Out-of-bag R^2: %s\n', fmt(perf[['oob_r_squared']])))
  }

  k <- min(n, nrow(imp))
  if (k > 0) {
    col <- if ('permutation' %in% names(imp) && !identical(imp$permutation, imp$impurity) &&
               is.unsorted(-imp$impurity)) 'permutation' else 'impurity'
    cat(sprintf('# Top %d feature%s by %s importance:\n', k, if (k == 1) '' else 's', col))
    top <- imp[seq_len(k), c('feature', intersect(c('impurity', 'permutation'), names(imp))), drop = FALSE]
    for (j in seq_along(top)[-1]) top[[j]] <- signif(top[[j]], 3)
    print(top, row.names = FALSE)
  }
  invisible(x)
}



#' Plot a random forest
#'
#' Five views of a [random_forest()] result, drawn from the result alone:
#'
#' * `'importance'` (the default): horizontal bars of the `n_features`
#'   most important features, largest at the top. With impurity importance,
#'   the full-data forest's values, with whiskers of the cross-validation
#'   folds' standard deviation; with permutation importance, the mean over
#'   trees, with whiskers of its standard error (SD / sqrt(trees)).
#' * `'roc'` (classification): the ROC curve of each class against the
#'   rest (one curve for two classes), and for more classes the micro and
#'   macro averages, with AUCs in the legend, as q2-sample-classifier draws
#'   them. Cross-validated (out of bag without cross-validation).
#' * `'confusion'` (classification): a heatmap of the confusion matrix,
#'   each observed class's row as proportions (shading) and counts
#'   (written).
#' * `'predictions'` (regression): cross-validated (or out-of-bag)
#'   predictions against the observed values, with the 1:1 line.
#' * `'tuning'` (with a grid): each grid point's out-of-bag error, one line
#'   per training set.
#'
#' @param x   A [random_forest()] result.
#' @param type   Which view: `'importance'`, `'roc'`, `'confusion'`,
#'   `'predictions'` or `'tuning'`. Default: `'importance'`
#' @param n_features   For `'importance'`: how many features to show.
#'   Default: `20`
#' @param importance   For `'importance'`: `'impurity'` or `'permutation'`;
#'   `NULL` for the one the result was sorted by. Default: `NULL`
#' @param ...   Passed on to the base plotting function (e.g. `main`,
#'   `xlab`, `col`).
#'
#' @return `x`, invisibly.
#'
#' @seealso [random_forest()]
#'
#' @export
#' @examples
#'     counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
#'     site   <- rep(rownames(ex_counts), each = 5)
#'     fit    <- random_forest(counts, site, n_trees = 50, cpus = 1)
#'     plot(fit)
#'     plot(fit, 'roc')
#'     plot(fit, 'confusion')
plot.ecodive_random_forest <- function (
    x, type = c('importance', 'roc', 'confusion', 'predictions', 'tuning'),
    n_features = 20, importance = NULL, ... ) {

  type <- match.arg(type)
  args <- list(...)
  cls  <- attr(x, 'type', exact = TRUE) == 'classification'
  cv   <- attr(x, 'n_folds', exact = TRUE) >= 2
  how  <- if (cv) 'cross-validated' else 'out of bag'

  if (type %in% c('roc', 'confusion') && !cls)
    stop("`type = '", type, "'` is for classification forests.", call. = FALSE)
  if (type == 'predictions' && cls)
    stop("`type = 'predictions'` is for regression forests; use 'confusion' or 'roc'.", call. = FALSE)
  if (type == 'tuning' && is.null(x$tuning))
    stop("`type = 'tuning'` needs a grid of `max_features` or `min_samples_leaf` values.", call. = FALSE)

  switch(type,
    importance  = rf_plot_importance(x, n_features, importance, args),
    roc         = rf_plot_roc(x, how, args),
    confusion   = rf_plot_confusion(x, how, args),
    predictions = rf_plot_predictions(x, cv, how, args),
    tuning      = rf_plot_tuning(x, args) )

  invisible(x)
}


rf_plot_importance <- function (x, n_features, importance, args) {

  n   <- plot_count(n_features, 'n_features')
  if (n < 1) stop('`n_features` must be at least 1.', call. = FALSE)
  imp <- x$importance
  has_perm <- 'permutation' %in% names(imp)
  if (is.null(importance))
    importance <- if (has_perm && !is.unsorted(-imp$permutation)) 'permutation' else 'impurity'
  importance <- match.arg(importance, c('impurity', 'permutation'))
  if (importance == 'permutation' && !has_perm)
    stop("The forest has no permutation importances: grow it with `importance = 'permutation'`.",
         call. = FALSE)

  if (importance == 'impurity') {
    val <- imp$impurity
    err <- imp$impurity_cv_sd
    sub <- 'Mean decrease in impurity (whiskers: SD over CV folds)'
  }
  else {
    val <- imp$permutation
    err <- imp$permutation_sd / sqrt(attr(x, 'n_trees', exact = TRUE))
    sub <- 'Out-of-bag permutation importance (whiskers: SE over trees)'
  }
  rows <- utils::head(order(-val, seq_along(val)), n)
  rows <- rev(rows)                 # bottom to top
  kk   <- length(rows)
  v    <- val[rows]
  e    <- err[rows]
  e[!is.finite(e)] <- 0

  opar <- graphics::par('mai')
  on.exit(graphics::par(mai = opar))
  linch <- max(graphics::strwidth(imp$feature[rows], 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)

  lo <- min(0, v - e, na.rm = TRUE)
  hi <- max(0, v + e, na.rm = TRUE)
  if (!(hi > lo)) hi <- lo + 1
  y  <- seq_len(kk)

  plot_do(
    graphics::plot.default, c(list(0, 0), args),
    list(
      type = 'n', yaxt = 'n', xlim = c(lo, hi), ylim = c(0.5, max(kk, 1) + 0.5),
      main = 'Random Forest Feature Importance', xlab = 'Importance', ylab = '',
      sub = sub ))
  col <- plot_arg(args, 'col', PLOT_HIGHLIGHT)
  graphics::rect(0, y - 0.35, v, y + 0.35, col = col, border = NA)
  ok <- e > 0
  if (any(ok)) graphics::segments(v[ok] - e[ok], y[ok], v[ok] + e[ok], y[ok], col = PLOT_MUTED)
  graphics::axis(2, at = y, labels = imp$feature[rows], tick = FALSE, las = 1, cex.axis = 0.8)
  graphics::abline(v = 0, col = PLOT_MUTED)
}


rf_plot_roc <- function (x, how, args) {
  roc  <- x$roc
  perf <- x$performance
  cl   <- levels(x$predictions$observed)
  k    <- length(cl)
  curves <- unique(roc$class)
  aucs <- vapply(curves, function (c) {
    r <- roc[roc$class == c, ]
    rf_auc(r$fpr, r$tpr)
  }, 0)

  plot_do(
    graphics::plot.default, c(list(c(0, 1), c(0, 1)), args),
    list(
      type = 'n', xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i',
      main = sprintf('Random Forest ROC (%s)', how),
      xlab = 'False positive rate', ylab = 'True positive rate' ))
  graphics::abline(0, 1, col = PLOT_MUTED, lty = 2)

  cols <- if (k == 2) PLOT_HIGHLIGHT else c(plot_palette(k), 'navy', 'lightblue3')
  ltys <- if (k == 2) 1 else c(rep(1, k), 3, 3)
  names(cols) <- names(ltys) <- if (k == 2) curves else c(cl, 'micro', 'macro')
  for (c in curves) {
    r <- roc[roc$class == c, ]
    graphics::lines(r$fpr, r$tpr, col = cols[[c]], lty = ltys[[c]], lwd = 2)
  }
  lab <- sprintf('%s (AUC = %s)',
    ifelse(curves %in% c('micro', 'macro'), paste0(curves, '-average'), curves),
    vapply(aucs, plot_num, ''))
  graphics::legend('bottomright', legend = lab, col = cols[curves], lty = ltys[curves],
                   lwd = 2, bty = 'n', cex = 0.8)
}


rf_plot_confusion <- function (x, how, args) {
  cm   <- x$confusion
  k    <- nrow(cm)
  rs   <- rowSums(cm)
  prop <- cm / ifelse(rs == 0, 1, rs)
  perf <- x$performance
  acc  <- if (how == 'cross-validated') perf[['accuracy']] else perf[['oob_accuracy']]
  sub  <- sprintf('Accuracy %s (%s); baseline %s; ratio %s', plot_num(acc), how,
    plot_num(perf[['baseline_accuracy']]), plot_num(acc / perf[['baseline_accuracy']]))

  pal <- grDevices::colorRampPalette(c('white', PLOT_HIGHLIGHT))(101)
  # image() draws column 1 at the left and row 1 at the bottom: observed
  # classes from the top.
  z <- t(prop[rev(seq_len(k)), , drop = FALSE])
  plot_do(
    graphics::image, c(list(x = seq_len(k), y = seq_len(k), z = z), args),
    list(
      col = pal, zlim = c(0, 1), axes = FALSE,
      main = 'Random Forest Confusion Matrix', xlab = 'Predicted', ylab = 'Observed',
      sub = sub ))
  graphics::axis(1, at = seq_len(k), labels = colnames(cm), tick = FALSE, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(k), labels = rev(rownames(cm)), tick = FALSE, las = 1, cex.axis = 0.8)
  graphics::box()
  cz <- t(cm[rev(seq_len(k)), , drop = FALSE])
  graphics::text(rep(seq_len(k), k), rep(seq_len(k), each = k), labels = as.vector(cz),
                 col = ifelse(as.vector(z) > 0.6, 'white', 'black'), cex = 0.9)
}


rf_plot_predictions <- function (x, cv, how, args) {
  pr   <- x$predictions
  pred <- if (cv) pr$predicted else pr$oob_predicted
  perf <- x$performance
  sub  <- if (cv) sprintf('R^2 = %s, MSE = %s (%s)', plot_num(perf[['r_squared']]),
                          plot_num(perf[['mse']]), how)
          else sprintf('R^2 = %s (%s)', plot_num(perf[['oob_r_squared']]), how)
  rng  <- range(c(pr$observed, pred), na.rm = TRUE)
  plot_do(
    graphics::plot.default, c(list(pr$observed, pred), args),
    list(
      xlim = rng, ylim = rng, pch = 19, col = PLOT_HIGHLIGHT,
      main = sprintf('Random Forest Predictions (%s)', how),
      xlab = 'Observed', ylab = 'Predicted', sub = sub ))
  graphics::abline(0, 1, col = PLOT_MUTED, lty = 2)
}


rf_plot_tuning <- function (x, args) {
  tu  <- x$tuning
  lab <- sprintf('%d / %d', tu$max_features, tu$min_samples_leaf)
  pts <- unique(lab)
  sets <- unique(tu$training_set)
  cols <- c(PLOT_HIGHLIGHT, rep(PLOT_MUTED, length(sets) - 1))
  xi <- match(lab, pts)
  plot_do(
    graphics::plot.default, c(list(xi, tu$oob_error), args),
    list(
      type = 'n', xaxt = 'n', xlim = c(0.5, length(pts) + 0.5),
      main = 'Random Forest Tuning', xlab = 'max_features / min_samples_leaf',
      ylab = 'Out-of-bag error',
      sub = 'Red: the full data; grey: each CV fold\'s training samples; filled: chosen' ))
  graphics::axis(1, at = seq_along(pts), labels = pts, cex.axis = 0.8)
  for (s in rev(seq_along(sets))) {
    r <- tu$training_set == sets[[s]]
    graphics::lines(xi[r], tu$oob_error[r], col = cols[[s]], type = 'b',
                    pch = ifelse(tu$chosen[r], 19, 1))
  }
}
