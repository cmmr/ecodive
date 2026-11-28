# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# coda_lasso() and coda_lasso_pairwise(): penalized log-contrast regression.
# coda_lasso() ports zeroSum 2.0.7's solver, lambda path and CV driver (MIT:
# Copyright (c) 2019 Thorsten Rehberg; see src/coda_lasso.c and
# inst/COPYRIGHTS). coda_lasso_pairwise() follows coda4microbiome 0.2.4's
# coda_glmnet() (MIT: Copyright (c) 2022 coda4microbiome authors): its zero
# handling (here ecodive's shared rule instead), covariate offset and
# performance summary, on ecodive's own coordinate descent
# (src/coda_lasso_pairwise.c) in place of glmnet's.


#' Compositional lasso: penalized log-contrast regression
#'
#' Predicts an outcome (a measurement, or one of two classes) from a
#' feature table and picks out the features that matter, in a way that
#' respects compositionality: the **log-contrast model** (Aitchison &
#' Bacon-Shone 1984; Lin et al. 2014), fitted with an elastic-net penalty
#' and the penalty chosen by cross-validation. `coda_lasso()` constrains the
#' features' coefficients to sum to zero, as zeroSum does
#' (Altenbuchinger et al. 2017); `coda_lasso_pairwise()` puts the penalty on
#' all pairwise log ratios, as coda4microbiome's `coda_glmnet()` does (Calle
#' et al. 2023).
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#' @inheritSection ml_section Machine learning
#' @inheritSection documentation Input Types
#'
#' @param counts  A numeric matrix of feature values (samples
#'        \eqn{\times} features), or any of the table types listed under
#'        "Input Types". Typically counts; any positive scale (counts,
#'        proportions, the output of [replace_zeros()]) gives the same fit.
#'
#' @param y   What to predict: a vector with one value per sample, or the
#'   name of a column of `data`. Numeric `y` is a measurement (`family =
#'   'gaussian'`); a factor, logical or character `y` with two values is a
#'   class (`'binomial'`; the second level, or `TRUE`, is the "event" whose
#'   log odds are modelled). Samples whose `y` is `NA` are left out, with a
#'   warning.
#'
#' @param data   A data frame of sample metadata, with one row per sample
#'   (matched to `counts` by row names when both have them), when `y`,
#'   `given` or `folds` name its columns. Default: `NULL`
#'
#' @param given   Covariates to adjust for (age, sex, BMI): a one-sided
#'   formula over the columns of `data` (such as `~ age + sex`), the names
#'   of columns of `data`, or a data frame with one row per sample.
#'   `coda_lasso()` fits them alongside the features, unpenalized and
#'   outside the zero-sum constraint (zeroSum's `penalty.factor = 0` and
#'   `zeroSum.weights = 0`); `coda_lasso_pairwise()` first regresses `y` on
#'   them alone and passes that fit as an offset (coda4microbiome's
#'   `covar`). Samples with `NA` in `given` are left out, with a warning;
#'   covariates that are linear combinations of others are dropped, with a
#'   warning. Default: `NULL`
#'
#' @param family   `'gaussian'` (least squares) or `'binomial'` (logistic
#'   regression), or `NULL` to choose from `y` (see `y`). Default: `NULL`
#'
#' @param l1_ratio   The elastic-net mixing parameter, between 0 (ridge)
#'   and 1 (lasso): the penalty is \eqn{\lambda[(1 - r)/2 \|\beta\|_2^2 +
#'   r \|\beta\|_1]}. glmnet and zeroSum call it `alpha`. Default: `1` for
#'   `coda_lasso()` (zeroSum's), `0.9` for `coda_lasso_pairwise()`
#'   (coda4microbiome's).
#'
#' @param lambda   The penalties to fit, or `NULL` for a path of
#'   `n_lambdas` values from the smallest that makes every coefficient zero
#'   down to `lambda_min_ratio` times it, evenly spaced on the log scale.
#'   Default: `NULL`
#'
#' @param n_lambdas   The length of the default path. Default: `100`
#'
#' @param lambda_min_ratio   The smallest penalty of the default path, as a
#'   fraction of the largest; `NULL` for 0.01 when there are fewer samples
#'   than coefficients, else 1e-4 (zeroSum's and glmnet's rule). Default:
#'   `NULL`
#'
#' @param standardize   Penalize each coefficient in proportion to the
#'   standard deviation of its column (log abundance, or log ratio), as if
#'   the columns were scaled to unit variance; the constraint and the
#'   reported coefficients stay on the original scale. Default: `FALSE` for
#'   `coda_lasso()` (zeroSum's), `TRUE` for `coda_lasso_pairwise()`
#'   (coda4microbiome's, through glmnet's default).
#'
#' @param n_folds   The number of cross-validation folds, which choose
#'   \eqn{\lambda}; `0` or `1` for none (then `lambda` must be a single
#'   value). Default: `10`
#'
#' @param folds   A fold for each sample instead (any labels; each label is
#'   a fold, as glmnet's and zeroSum's `foldid`), matched to the samples
#'   like `y`. Overrides `n_folds` and `seed`. Default: `NULL`
#'
#' @param select   Which \eqn{\lambda} the result reports: `'min'`, the one
#'   with the best cross-validated error, or `'1se'`, a larger (sparser)
#'   one within a standard error of it (see Details). Default: `'min'` for
#'   `coda_lasso()` (zeroSum's), `'1se'` for `coda_lasso_pairwise()`
#'   (coda4microbiome's).
#'
#' @param tol   The convergence tolerance: zeroSum's `precision` (the
#'   change in the cost function between sweeps) for `coda_lasso()`;
#'   glmnet's `thresh` (the largest change of any coefficient's
#'   contribution, relative to the null deviance) for
#'   `coda_lasso_pairwise()`. Default: `1e-8` and `1e-7`.
#'
#' @param max_iters   The most iterations: zeroSum's outer (active-set)
#'   steps for each \eqn{\lambda} in `coda_lasso()`; glmnet's `maxit`, the
#'   most passes over the coefficients in all, for `coda_lasso_pairwise()`.
#'   Default: `100` and `1e5`.
#'
#' @param seed   An integer seed for the folds. Default: `0`
#'
#' @details
#'
#' ## The log-contrast model
#'
#' For samples \eqn{i} and features \eqn{j}, with \eqn{g} the identity
#' (gaussian) or the logit (binomial),
#' \deqn{g(E[y_i]) = \beta_0 + \sum_j \beta_j \log x_{ij},
#'   \quad \sum_j \beta_j = 0.}
#' Because the coefficients sum to zero, multiplying a sample by any
#' constant changes nothing: the fit is the same on counts, proportions or
#' rarefied counts, and \eqn{\sum_j \beta_j \log x_{ij} = \sum_j \beta_j
#' \mathrm{clr}(x_i)_j}. The features with positive coefficients form the
#' numerator, and those with negative ones the denominator, of a single
#' weighted log ratio (the "signature"; Calle et al. 2023).
#'
#' `coda_lasso()` minimizes
#' \deqn{-\frac{1}{n}\ell(\beta_0, \beta) + \lambda \sum_j \left[\frac{1 -
#'   r}{2} v_j^2 \beta_j^2 + r v_j |\beta_j|\right]}
#' (\eqn{\ell} the log-likelihood, \eqn{r} = `l1_ratio`, \eqn{v_j} the
#' column's standard deviation with `standardize`, else 1) subject to the
#' constraint, by zeroSum's coordinate descent over pairs of coefficients
#' (raising one while lowering another keeps the sum at zero). For
#' gaussian `y` the ridge term is divided by the standard deviation of
#' `y`, as zeroSum and glmnet do.
#'
#' `coda_lasso_pairwise()` fits the same model through the log ratios of
#' every pair of features, \eqn{\log(x_j / x_k)}, \eqn{j < k}, as
#' coda4microbiome does: \eqn{\beta_j = \sum_k \theta_{jk} - \sum_k
#' \theta_{kj}} sums to zero by construction, and the penalty falls on the
#' \eqn{\theta}. It never forms the \eqn{n \times p(p-1)/2} table of log
#' ratios. With a pure lasso (`l1_ratio = 1`), no standardization and no
#' covariates, its fit at \eqn{2\lambda} is `coda_lasso()`'s at
#' \eqn{\lambda}, since the cheapest \eqn{\theta} giving a \eqn{\beta} has
#' \eqn{\|\theta\|_1 = \|\beta\|_1 / 2}. With the defaults (standardized,
#' `l1_ratio = 0.9`) they differ.
#'
#' ## Choosing lambda
#'
#' Each fold's samples are predicted from a fit to the others along the
#' whole path. `coda_lasso()` follows zeroSum: the error is the mean
#' squared error (gaussian) or the binomial deviance per sample, averaged
#' over folds, with standard error \eqn{\mathrm{sd}/\sqrt{K}} over the
#' \eqn{K} folds; `'1se'` is the \eqn{\lambda} (no smaller than the best)
#' whose error is closest to the best error plus its standard error.
#' `coda_lasso_pairwise()` follows `cv.glmnet()` as coda4microbiome calls
#' it: the mean squared error (gaussian) or the AUC (binomial; the
#' binomial deviance, with a warning, when there are fewer than 10 samples
#' per fold, as glmnet does), averaged over folds weighted by their size;
#' `'1se'` is the largest \eqn{\lambda} whose error is within a standard
#' error of the best. The default path stops early, as glmnet's does, once
#' the fit explains 99.9% of the deviance or stops improving.
#'
#' Folds are drawn at random with `seed` (for binomial `y`, within each
#' class, so that each fold has both). Penalized fits come with no
#' p-values: a feature's selection is not a test of it.
#'
#' ## Zeros and covariates
#'
#' Zeros have no logarithm: replace them first (see [replace_zeros()]),
#' or a `pseudocount` is added (see "Pseudocount" below).
#'
#' @section Differences from zeroSum and coda4microbiome:
#'
#' * **zeroSum:** the largest penalty of the default path includes every
#'   pair of features (zeroSum's `lambdaMax` skips the first feature), and
#'   with `given` it is computed from the residuals of `y` on the
#'   covariates; the path is always fitted to its end (zeroSum's `cvStop =
#'   0.1` stops it once the CV error rises); folds are drawn with ecodive's
#'   random numbers (and stratified by class); `max_iters` caps zeroSum's
#'   fixed 100 outer steps. Coefficients agree with zeroSum's to its
#'   convergence tolerance; both reach the same cost to about 1e-10.
#' * **coda4microbiome:** zeros get ecodive's shared pseudocount rule
#'   (coda4microbiome adds 1, or half the smallest value below 1, to every
#'   value); `coefficient` is \eqn{\beta} as fitted, not rescaled to
#'   \eqn{\sum|\beta| = 2}; the apparent \eqn{R^2} and AUC are computed
#'   whether or not features are selected (coda4microbiome reports 0 and
#'   0.5 when none are); the reported CV AUC is the AUC even when the
#'   deviance chose \eqn{\lambda}. glmnet's solution is matched to its
#'   tolerance, not its iterations.
#'
#' @section Same names in other packages:
#' None known.
#'
#' @return A data frame of class `ecodive_coda_lasso` (or
#' `ecodive_coda_lasso_pairwise`), one row per feature in table order:
#' `feature`, `coefficient` (\eqn{\beta_j} at the chosen \eqn{\lambda}) and
#' `selected` (`coefficient != 0`). Its attributes:
#'
#' \describe{
#'   \item{`method`}{`'coda_lasso'` or `'coda_lasso_pairwise'`.}
#'   \item{`intercept`, `lambda`, `select`, `family`, `l1_ratio`,
#'     `pseudocount`}{The fit's intercept, the chosen \eqn{\lambda} and the
#'     rule that chose it, and the settings that qualify the
#'     coefficients.}
#'   \item{`covariates`}{With `given`: the covariates' coefficients
#'     (`coda_lasso()`), or those of the offset model
#'     (`coda_lasso_pairwise()`, with its intercept).}
#'   \item{`folds`}{Each sample's fold.}
#'   \item{`path`}{A data frame: `lambda`, `cv_error`, `cv_se` and
#'     `n_selected`, one row per \eqn{\lambda}; `cv_measure` names the
#'     error (`'mse'`, `'deviance'` or `'auc'`).}
#'   \item{`coef_path`}{The features' coefficients along the path
#'     (features \eqn{\times} \eqn{\lambda}).}
#'   \item{`fitted`}{A data frame: `sample`, `link` (the linear predictor,
#'     the "signature score") and, for binomial `y`, `probability`;
#'     `cv_link` holds each sample's cross-validated linear predictor.}
#'   \item{`performance`}{Gaussian: `r_squared` (the squared correlation of
#'     the fitted and observed values), `cv_mse` and `cv_mse_se`. Binomial:
#'     `auc` (apparent), `cv_auc` and `cv_auc_se` (the folds' AUCs averaged
#'     as glmnet does, at the chosen \eqn{\lambda}).}
#'   \item{`y`}{The outcome as fitted (a factor for binomial `y`).}
#'   \item{`pairs`}{`coda_lasso_pairwise()` only: the non-zero
#'     \eqn{\theta_{jk}}, a data frame of `numerator`, `denominator` and
#'     `coefficient`.}
#' }
#'
#' @references
#'
#' Aitchison, J., & Bacon-Shone, J. (1984). Log contrast models for
#' experiments with mixtures. *Biometrika*, 71, 323-330.
#' \doi{10.1093/biomet/71.2.323}
#'
#' Lin, W., Shi, P., Feng, R., & Li, H. (2014). Variable selection in
#' regression with compositional covariates. *Biometrika*, 101, 785-797.
#' \doi{10.1093/biomet/asu031}
#'
#' Altenbuchinger, M., Rehberg, T., Zacharias, H. U., et al. (2017).
#' Reference point insensitive molecular data analysis. *Bioinformatics*,
#' 33, 219-226. \doi{10.1093/bioinformatics/btw598}
#'
#' Calle, M. L., Pujolassos, M., & Susin, A. (2023). coda4microbiome:
#' compositional data analysis for microbiome cross-sectional and
#' longitudinal studies. *BMC Bioinformatics*, 24, 82.
#' \doi{10.1186/s12859-023-05205-3}
#'
#' Susin, A., Wang, Y., Lê Cao, K.-A., & Calle, M. L. (2020). Variable
#' selection in microbiome compositional data analysis. *NAR Genomics and
#' Bioinformatics*, 2, lqaa029. \doi{10.1093/nargab/lqaa029}
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths
#' for generalized linear models via coordinate descent. *Journal of
#' Statistical Software*, 33(1), 1-22. \doi{10.18637/jss.v033.i01}
#'
#' Tibshirani, R., Bien, J., Friedman, J., et al. (2012). Strong rules for
#' discarding predictors in lasso-type problems. *Journal of the Royal
#' Statistical Society B*, 74, 245-266.
#' \doi{10.1111/j.1467-9868.2011.01004.x}
#'
#' @seealso [predict.ecodive_coda_lasso()], [plot.ecodive_coda_lasso()],
#'   [log_ratio()], [random_forest()]
#'
#' @export
#' @examples
#'     # Which features predict a measurement? 30 samples, 8 features.
#'     set.seed(1)
#'     counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
#'     y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
#'     fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
#'     fit
#'
#'     # The same with all pairwise log ratios.
#'     coda_lasso_pairwise(counts, y, n_folds = 5, cpus = 1)
coda_lasso <- function (
    counts, y, data = NULL, given = NULL, family = NULL, margin = 1L,
    pseudocount = NULL, l1_ratio = 1, lambda = NULL, n_lambdas = 100L,
    lambda_min_ratio = NULL, standardize = FALSE, n_folds = 10L, folds = NULL,
    select = c('min', '1se'), tol = 1e-8, max_iters = 100L, seed = 0L,
    cpus = n_cpus() ) {

  if (missing(y)) y <- NULL
  cl_run(environment(), pairwise = FALSE)
}


#' @rdname coda_lasso
#' @export
coda_lasso_pairwise <- function (
    counts, y, data = NULL, given = NULL, family = NULL, margin = 1L,
    pseudocount = NULL, l1_ratio = 0.9, lambda = NULL, n_lambdas = 100L,
    lambda_min_ratio = NULL, standardize = TRUE, n_folds = 10L, folds = NULL,
    select = c('1se', 'min'), tol = 1e-7, max_iters = 1e5L, seed = 0L,
    cpus = n_cpus() ) {

  if (missing(y)) y <- NULL
  cl_run(environment(), pairwise = TRUE)
}



#------------------------------------------------------------------------------
# The shared front end. `zero_sum = FALSE` (coda_lasso() only; internal, for
# the tests' glmnet cross-check) drops the constraint.
#------------------------------------------------------------------------------
cl_run <- function (env, pairwise, zero_sum = TRUE) {

  with(env, {
    validate_counts()
    validate_margin()
    validate_cpus()
    validate_seed()
  })
  counts <- env$counts; margin <- env$margin; cpus <- env$cpus; seed <- env$seed
  method <- if (pairwise) 'coda_lasso_pairwise' else 'coda_lasso'

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
  if (n < 3) stop('`counts` must have at least three samples.', call. = FALSE)
  if (p < 2) stop('`counts` must have at least two features.', call. = FALSE)
  feats <- if (is.null(dn[[2]])) sprintf('F%d', seq_len(p)) else dn[[2]]
  sams  <- dn[[1]]


  # Settings.
  l1_ratio <- env$l1_ratio
  if (!is.numeric(l1_ratio) || length(l1_ratio) != 1 || is.na(l1_ratio) ||
      l1_ratio < 0 || l1_ratio > 1)
    stop('`l1_ratio` must be a number from 0 to 1.', call. = FALSE)
  l1_ratio <- as.double(l1_ratio)
  standardize <- env$standardize
  if (!is.logical(standardize) || length(standardize) != 1 || is.na(standardize))
    stop('`standardize` must be TRUE or FALSE.', call. = FALSE)
  select <- match.arg(env$select, if (pairwise) c('1se', 'min') else c('min', '1se'))
  tol <- env$tol
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0)
    stop('`tol` must be a positive number.', call. = FALSE)
  max_iters <- cl_count(env$max_iters, 'max_iters', 1)
  n_lambdas <- cl_count(env$n_lambdas, 'n_lambdas', 1)
  lmr <- env$lambda_min_ratio
  if (!is.null(lmr) && (!is.numeric(lmr) || length(lmr) != 1 || !is.finite(lmr) ||
                        lmr <= 0 || lmr >= 1))
    stop('`lambda_min_ratio` must be a number between 0 and 1.', call. = FALSE)
  lambda <- env$lambda
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) < 1 || !all(is.finite(lambda)) || any(lambda < 0))
      stop('`lambda` must be non-negative numbers.', call. = FALSE)
    lambda <- sort(as.double(lambda), decreasing = TRUE)
  }


  # y, folds and given, matched to the samples.
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
  gm <- cl_given(env$given, data, idx, sams, n)

  keep <- !is.na(yv)
  if (!all(keep))
    warning(sprintf('Left out %d sample%s whose `y` is NA.', sum(!keep),
                    if (sum(!keep) == 1) '' else 's'), call. = FALSE)
  if (!is.null(gm)) {
    gok <- stats::complete.cases(gm$mm)
    if (!all(gok[keep]))
      warning(sprintf('Left out %d sample%s with NA in `given`.', sum(!gok[keep]),
                      if (sum(!gok[keep]) == 1) '' else 's'), call. = FALSE)
    keep <- keep & gok
  }
  keep <- which(keep)
  if (length(keep) < 3) stop('Fewer than three samples are left to fit.', call. = FALSE)
  yv <- yv[keep]
  if (!is.null(fv)) {
    fv <- fv[keep]
    if (anyNA(fv)) stop('`folds` must not contain NAs.', call. = FALSE)
  }
  sam_out <- if (is.null(sams)) as.character(keep) else sams[keep]
  nk <- length(keep)


  # The family, and y as fitted.
  family <- env$family
  if (!is.null(family)) family <- match.arg(family, c('gaussian', 'binomial'))
  if (is.factor(yv) || is.character(yv) || is.logical(yv)) {
    if (is.factor(yv)) yf <- droplevels(yv)
    else if (is.logical(yv)) yf <- factor(yv, levels = c(FALSE, TRUE))[, drop = TRUE]
    else yf <- factor(yv, levels = sort(unique(yv), method = 'radix'))
    if (nlevels(yf) > 2)
      stop('`y` has ', nlevels(yf), ' classes; only two are supported (binomial).', call. = FALSE)
    if (nlevels(yf) < 2) stop('`y` must have two classes.', call. = FALSE)
    if (identical(family, 'gaussian'))
      stop("`family = 'gaussian'` needs a numeric `y`.", call. = FALSE)
    family <- 'binomial'
    ycode  <- as.double(as.integer(yf) - 1L)
    y_out  <- yf
  }
  else if (is.numeric(yv)) {
    if (!all(is.finite(yv))) stop('`y` must be finite.', call. = FALSE)
    if (identical(family, 'binomial')) {
      if (!all(yv %in% c(0, 1)))
        stop("`family = 'binomial'` needs `y` of two classes, or 0 and 1.", call. = FALSE)
      if (length(unique(yv)) < 2) stop('`y` must have two classes.', call. = FALSE)
      y_out <- factor(yv, levels = c(0, 1))
    }
    else {
      family <- 'gaussian'
      y_out  <- as.double(yv)
    }
    ycode <- as.double(yv)
  }
  else {
    stop('`y` must be numeric, or a factor, character or logical vector.', call. = FALSE)
  }
  binom <- family == 'binomial'


  # Folds.
  if (!is.null(fv)) {
    ff <- if (is.factor(fv)) fv else factor(fv, levels = sort(unique(fv), method = 'radix'))
    if (any(tabulate(as.integer(ff), nlevels(ff)) == 0))
      stop('Every fold in `folds` must have samples.', call. = FALSE)
    if (nlevels(ff) < 2) stop('`folds` must have at least two folds.', call. = FALSE)
    fold <- as.integer(ff)
    K    <- nlevels(ff)
  }
  else {
    K <- cl_count(env$n_folds, 'n_folds', 0)
    if (K == 1L) K <- 0L
    if (K > nk) stop(sprintf('`n_folds` (%d) is larger than the number of samples (%d).', K, nk),
                     call. = FALSE)
    fold <- if (K >= 2) rf_folds(if (binom) y_out else ycode, binom, K, seed) else integer(nk)
  }
  if (K >= 2) {
    for (f in seq_len(K)) {
      tr <- fold != f
      if (sum(tr) < 2)
        stop('Each fold must leave at least two samples to train on.', call. = FALSE)
      if (binom && length(unique(ycode[tr])) < 2)
        stop('The training samples of fold ', f, ' all have one class; use fewer folds.',
             call. = FALSE)
    }
  }
  else if (is.null(lambda) || length(lambda) > 1) {
    stop('Choosing `lambda` needs cross-validation (`n_folds` of 2 or more), ',
         'or give a single `lambda`.', call. = FALSE)
  }


  # The log table.
  pc <- log_pseudocount(counts, env$pseudocount, 'log-contrast regression')
  # log(x + pc) of the values as given (as zeroSum and coda4microbiome):
  # with `standardize`, the columns' sds depend on the samples' totals.
  L <- .Call(C_normalize, counts, margin, NORM_NONE, 0, TRUE, cpus)
  if (margin == 2L) L <- t(L)
  storage.mode(L) <- 'double'
  dimnames(L) <- NULL
  L <- log(L[keep, , drop = FALSE] + pc)
  if (!all(is.finite(L)))
    stop('The log `counts` contain NA, NaN, or infinite values.', call. = FALSE) # nocov

  Z <- NULL
  if (!is.null(gm)) Z <- cl_given_matrix(gm, keep)


  fit <- if (pairwise)
    cl_pairwise(L, ycode, Z, fold, K, binom, l1_ratio, lambda, n_lambdas, lmr,
                standardize, select, tol, max_iters, cpus)
  else
    cl_zerosum(L, ycode, Z, fold, K, binom, l1_ratio, lambda, n_lambdas, lmr,
               standardize, select, tol, max_iters, cpus, zero_sum)


  # The result at the chosen lambda.
  s    <- fit$sel
  beta <- fit$coef_path[, s]
  link <- fit$link
  res  <- data.frame(feature = feats, coefficient = beta, selected = beta != 0,
                     stringsAsFactors = FALSE)

  cvl <- if (K >= 2) fit$cv_link[, s] else rep(NA_real_, nk)
  fitted <- data.frame(sample = sam_out, link = link, stringsAsFactors = FALSE)
  if (binom) fitted$probability <- stats::plogis(link)

  if (binom) {
    perf <- c(auc = cl_auc(link, ycode), cv_auc = NA_real_, cv_auc_se = NA_real_)
    if (K >= 2) {
      a <- cl_fold_auc(fit$cv_link[, s, drop = FALSE], ycode, fold, K)
      perf[['cv_auc']] <- a$m; perf[['cv_auc_se']] <- a$se
    }
  }
  else {
    r2 <- if (stats::sd(link) > 0) stats::cor(link, ycode)^2 else 0
    perf <- c(r_squared = r2, cv_mse = NA_real_, cv_mse_se = NA_real_)
    if (K >= 2) {
      if (identical(fit$measure, 'mse')) {
        perf[['cv_mse']] <- fit$cv_error[s]; perf[['cv_mse_se']] <- fit$cv_se[s]
      }
    }
  }

  path <- data.frame(
    lambda = fit$lambda, cv_error = fit$cv_error, cv_se = fit$cv_se,
    n_selected = colSums(fit$coef_path != 0) )
  coef_path <- fit$coef_path
  dimnames(coef_path) <- list(feats, NULL)

  attr(res, 'method')      <- method
  attr(res, 'intercept')   <- fit$intercept[s]
  attr(res, 'lambda')      <- fit$lambda[s]
  attr(res, 'select')      <- select
  attr(res, 'family')      <- family
  attr(res, 'l1_ratio')    <- l1_ratio
  attr(res, 'pseudocount') <- pc
  if (!is.null(Z)) attr(res, 'covariates') <- fit$covariates
  attr(res, 'folds')       <- if (K >= 2) fold else NULL
  attr(res, 'path')        <- path
  attr(res, 'cv_measure')  <- if (K >= 2) fit$measure else NULL
  attr(res, 'coef_path')   <- coef_path
  attr(res, 'fitted')      <- fitted
  attr(res, 'cv_link')     <- if (K >= 2) cvl else NULL
  attr(res, 'performance') <- perf
  attr(res, 'y')           <- y_out
  if (pairwise) {
    pr <- fit$pairs
    attr(res, 'pairs') <- data.frame(
      numerator = feats[pr$j], denominator = feats[pr$k], coefficient = pr$t,
      stringsAsFactors = FALSE )
  }
  if (!is.null(gm)) attr(res, 'given_model') <- gm[c('terms', 'xlevels', 'cols')]
  class(res) <- c(paste0('ecodive_', method), 'data.frame')
  res
}


# A whole number >= `lo`.
cl_count <- function (x, arg, lo) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x %% 1 != 0 || x < lo)
    stop('`', arg, '` must be a whole number', if (lo > 0) sprintf(' of at least %d', lo), '.',
         call. = FALSE)
  as.integer(x)
}



#------------------------------------------------------------------------------
# `given`: the covariates' model frame terms, in x's sample order (NAs kept).
# Returns NULL or list(mm (n x k, NAs kept), terms, xlevels, cols).
#------------------------------------------------------------------------------
cl_given <- function (given, data, idx, labels, n) {

  if (is.null(given)) return (NULL)
  cols <- NULL

  if (inherits(given, 'formula')) {
    if (is.null(data))
      stop('`given` is a formula, but `data` is not given.', call. = FALSE)
    if (length(given) != 2)
      stop('`given` must be a one-sided formula, such as ~ age + sex.', call. = FALSE)
    df <- data[idx, , drop = FALSE]
    tt <- stats::terms(given, data = df)
  }
  else if (is.character(given) && is.null(dim(given))) {
    if (is.null(data))
      stop('`given` names columns, but `data` is not given.', call. = FALSE)
    miss <- setdiff(given, names(data))
    if (length(miss))
      stop('`given` names ', name_list(miss), ', not column', if (length(miss) > 1) 's',
           ' of `data`.', call. = FALSE)
    df   <- data[idx, given, drop = FALSE]
    tt   <- stats::terms(~ ., data = df)
    cols <- given
  }
  else {
    if (is.matrix(given)) given <- as.data.frame(given)
    if (!is.data.frame(given))
      stop('`given` must be a formula, column names of `data`, or a data frame ',
           'with one row per sample.', call. = FALSE)
    rn <- if (.row_names_info(given) > 0) rownames(given) else NULL
    if (!is.null(rn) && (anyNA(rn) || anyDuplicated(rn))) rn <- NULL
    df <- given[match_rows(labels, n, rn, nrow(given), ref = '`counts`', what = '`given`'),
                , drop = FALSE]
    tt   <- stats::terms(~ ., data = df)
    cols <- names(given)
  }

  mf <- stats::model.frame(tt, df, na.action = stats::na.pass)
  xlevels <- stats::.getXlevels(tt, mf)
  ok <- stats::complete.cases(mf)
  mm <- matrix(NA_real_, n, 0)
  if (any(ok)) {
    mok <- stats::model.matrix(tt, mf[ok, , drop = FALSE])
    mok <- mok[, colnames(mok) != '(Intercept)', drop = FALSE]
    mm  <- matrix(NA_real_, n, ncol(mok), dimnames = list(NULL, colnames(mok)))
    mm[ok, ] <- mok
  }
  if (ncol(mm) == 0) stop('`given` has no covariates.', call. = FALSE)
  list(mm = mm, terms = stats::delete.response(tt), xlevels = xlevels, cols = cols)
}


# The covariates of the kept samples, without aliased columns (as lm()).
cl_given_matrix <- function (gm, keep) {
  Z <- gm$mm[keep, , drop = FALSE]
  q <- qr(cbind(1, Z))
  if (q$rank < ncol(Z) + 1L) {
    piv  <- q$pivot[seq_len(q$rank)]
    kept <- sort(piv[piv > 1L]) - 1L
    drop <- setdiff(seq_len(ncol(Z)), kept)
    warning('Dropped `given` covariate', if (length(drop) > 1) 's', ' ',
            name_list(colnames(Z)[drop]), ': linear combination',
            if (length(drop) > 1) 's', ' of the others.', call. = FALSE)
    Z <- Z[, kept, drop = FALSE]
  }
  if (ncol(Z) == 0) return (NULL) # nocov
  attr(Z, 'kept') <- colnames(Z)
  Z
}


# AUC (Mann-Whitney, ties count 1/2) of `score` for 0/1 `y`.
cl_auc <- function (score, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return (NA_real_)
  r <- rank(score)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


# Per-fold errors (K x L) averaged as cv.glmnet: weighted by fold size, with
# SE sqrt(weighted variance / (K - 1)).
cl_fold_summary <- function (fm, w) {
  K  <- nrow(fm)
  m  <- apply(fm, 2, function (e) stats::weighted.mean(e, w, na.rm = TRUE))
  se <- sqrt(apply(sweep(fm, 2, m)^2, 2, function (e) stats::weighted.mean(e, w, na.rm = TRUE)) /
             (K - 1))
  list(m = m, se = se)
}

cl_fold_auc <- function (link, y, fold, K) {
  fm <- t(vapply(seq_len(K), function (f) {
    i <- fold == f
    apply(link[i, , drop = FALSE], 2, cl_auc, y = y[i])
  }, numeric(ncol(link))))
  if (ncol(link) == 1) fm <- matrix(fm, ncol = 1)
  cl_fold_summary(fm, tabulate(fold, K))
}



#------------------------------------------------------------------------------
# coda_lasso(): zeroSum.
#------------------------------------------------------------------------------
cl_zerosum <- function (L, y, Z, fold, K, binom, alpha, lambda, n_lambdas, lmr,
                        standardize, select, tol, max_iters, cpus, zero_sum = TRUE) {

  n <- nrow(L); p <- ncol(L); q <- if (is.null(Z)) 0L else ncol(Z)
  X <- if (q) cbind(L, Z) else L
  u <- c(rep(1, p), rep(0, q))
  v <- c(rep(1, p), rep(0, q))
  if (!zero_sum) u[] <- 1

  constant <- !binom && all(y == y[[1]])
  if (constant)
    warning('`y` is constant, so every coefficient is 0.', call. = FALSE)

  if (is.null(lambda)) {
    lmax <- if (constant) 1 else cl_zs_lambda_max(L, y, Z, binom, alpha, standardize, zero_sum)
    if (!(lmax > 0) || !is.finite(lmax)) lmax <- 1   # nocov
    if (is.null(lmr)) lmr <- if (n < p + q) 0.01 else 1e-4
    lambda <- exp(seq(log(lmax), log(lmax * lmr), length.out = n_lambdas))
  }
  nl <- length(lambda)

  if (constant) {
    out <- list(beta = matrix(0, p + q, nl), intercept = rep(y[[1]], nl),
                cv_error = rep(0, nl), cv_sd = rep(0, nl),
                cv_link = matrix(y[[1]], n, nl), n_fit = nl)
    if (K < 2) { out$cv_error[] <- NA; out$cv_sd[] <- NA; out$cv_link[] <- NA }
  }
  else {
    out <- .Call(
      C_coda_lasso, X, as.double(y), as.integer(if (K >= 2) fold else rep(1L, n)),
      as.integer(K), as.integer(binom), alpha, as.double(lambda), as.double(u),
      as.double(v), as.integer(standardize), as.double(tol), as.integer(max_iters),
      1L, as.integer(zero_sum), as.integer(cpus) )
    if (out$n_fit < nl) {
      if (out$n_fit < 1) stop('The fit failed at the first `lambda`.', call. = FALSE) # nocov
      warning(sprintf('The fit diverged at lambda %d of %d; the path stops there.',
                      out$n_fit + 1L, nl), call. = FALSE)
      keep <- seq_len(out$n_fit)
      lambda <- lambda[keep]; nl <- length(lambda)
      out$beta <- out$beta[, keep, drop = FALSE]; out$intercept <- out$intercept[keep]
      out$cv_error <- out$cv_error[keep]; out$cv_sd <- out$cv_sd[keep]
      out$cv_link <- out$cv_link[, keep, drop = FALSE]
    }
  }

  # zeroSum's choice: the smallest CV error (the first, on ties); '1se' the
  # lambda no smaller whose error is closest to that plus its SD.
  if (K >= 2) {
    imin <- which.min(out$cv_error)
    sel  <- if (select == 'min') imin else
      which.min(abs(out$cv_error[seq_len(imin)] - (out$cv_error[imin] + out$cv_sd[imin])))
  }
  else sel <- 1L

  B   <- out$beta
  cov <- NULL
  if (q) {
    cov <- B[p + seq_len(q), sel]
    names(cov) <- colnames(Z)
  }
  link <- drop(X %*% B[, sel]) + out$intercept[sel]
  list(lambda = lambda, coef_path = B[seq_len(p), , drop = FALSE],
       intercept = out$intercept, covariates = cov, link = link,
       cv_error = if (K >= 2) out$cv_error else rep(NA_real_, nl),
       cv_se = if (K >= 2) out$cv_sd else rep(NA_real_, nl),
       cv_link = out$cv_link, sel = sel,
       measure = if (binom) 'deviance' else 'mse')
}


# zeroSum's lambdaMax: the largest |u_k g_s / u_s - g_k| / (alpha (v_k +
# v_s u_k/u_s)) over pairs of features, g the gradient at the null model
# (here, the model of the covariates alone).
cl_zs_lambda_max <- function (L, y, Z, binom, alpha, standardize, zero_sum) {
  n <- nrow(L)
  D <- if (is.null(Z)) matrix(1, n, 1) else cbind(1, Z)
  if (binom) {
    mu <- if (is.null(Z)) rep(mean(y), n) else
      suppressWarnings(stats::glm.fit(D, y, family = stats::binomial())$fitted.values)
    res <- (y - mu) / n
  }
  else {
    res <- if (is.null(Z)) y - mean(y) else stats::lm.fit(D, y)$residuals
    res <- res / n
  }
  a <- if (alpha == 0) 0.01 else alpha
  g <- drop(crossprod(L, res))
  v <- if (standardize) sqrt(colMeans(sweep(L, 2, colMeans(L))^2)) else rep(1, ncol(L))
  if (!zero_sum) return (max(abs(g) / (a * v)))
  if (!standardize) return ((max(g) - min(g)) / (2 * a))
  m <- 0
  for (j in seq_len(ncol(L) - 1L)) {
    k <- (j + 1L):ncol(L)
    m <- max(m, abs(g[j] - g[k]) / (a * (v[j] + v[k])))
  }
  m
}



#------------------------------------------------------------------------------
# coda_lasso_pairwise(): coda4microbiome's model on glmnet's conventions.
#------------------------------------------------------------------------------
cl_pairwise <- function (L, y, Z, fold, K, binom, alpha, lambda, n_lambdas, lmr,
                         standardize, select, tol, max_iters, cpus) {

  n <- nrow(L); p <- ncol(L)

  # The covariates' offset (coda4microbiome's model0), on all the samples.
  off <- rep(0, n)
  cov <- NULL
  if (!is.null(Z)) {
    D <- cbind(1, Z)
    m0 <- if (binom) suppressWarnings(stats::glm.fit(D, y, family = stats::binomial()))
          else stats::lm.fit(D, y)
    cov <- m0$coefficients
    names(cov) <- c('(Intercept)', colnames(Z))
    off <- drop(D %*% cov)
  }

  constant <- !binom && all(y == y[[1]])
  if (constant) warning('`y` is constant, so every coefficient is 0.', call. = FALSE)

  early <- is.null(lambda)
  if (early) {
    lmax <- if (constant) 1 else cl_pw_lambda_max(L, y, off, binom, alpha, standardize)
    if (!(lmax > 0) || !is.finite(lmax)) lmax <- 1   # nocov
    if (is.null(lmr)) lmr <- if (n < p * (p - 1) / 2) 0.01 else 1e-4
    lambda <- exp(seq(log(lmax), log(lmax * lmr), length.out = n_lambdas))
  }
  nl <- length(lambda)

  if (constant) {
    out <- list(beta = matrix(0, p, nl), intercept = rep(y[[1]], nl),
                cv_link = matrix(y[[1]], n, nl), n_done = rep(nl, K + 1L),
                maxit = rep(0L, K + 1L), pairs = list(l = integer(0), j = integer(0),
                                                       k = integer(0), t = numeric(0)))
  }
  else {
    out <- .Call(
      C_coda_lasso_pairwise, L, as.double(y), as.double(off),
      as.integer(if (K >= 2) fold else rep(1L, n)), as.integer(K), as.integer(binom),
      alpha, as.double(lambda), as.integer(standardize), as.double(tol),
      as.integer(max_iters), as.integer(early), as.integer(cpus) )
  }
  if (any(out$maxit > 0))
    warning('`max_iters` was reached before convergence; the path stops there.', call. = FALSE)
  nf <- out$n_done[[K + 1L]]
  if (nf < 1) stop('`max_iters` was reached at the first `lambda`.', call. = FALSE)
  keep <- seq_len(nf)
  lambda <- lambda[keep]

  # Folds that stopped early predict from their last fit beyond (done in C);
  # those cut short by max_iters stay NA there.
  cvm <- cvse <- rep(NA_real_, nf)
  measure <- NULL
  if (K >= 2) {
    link <- out$cv_link[, keep, drop = FALSE]
    w <- tabulate(fold, K)
    measure <- if (!binom) 'mse' else if (n / K < 10) 'deviance' else 'auc'
    if (binom && measure == 'deviance')
      warning('Fewer than 10 samples per fold: lambda is chosen by the binomial ',
              'deviance instead of the AUC, as glmnet does.', call. = FALSE)
    fm <- switch(measure,
      mse = t(vapply(seq_len(K), function (f) {
        i <- fold == f
        colMeans((y[i] - link[i, , drop = FALSE])^2) }, numeric(nf))),
      deviance = t(vapply(seq_len(K), function (f) {
        i <- fold == f
        pp <- pmin(pmax(stats::plogis(link[i, , drop = FALSE]), 1e-5), 1 - 1e-5)
        colMeans(-2 * (y[i] * log(pp) + (1 - y[i]) * log(1 - pp))) }, numeric(nf))),
      auc = t(vapply(seq_len(K), function (f) {
        i <- fold == f
        apply(link[i, , drop = FALSE], 2, cl_auc, y = y[i]) }, numeric(nf))) )
    if (nf == 1) fm <- matrix(fm, ncol = 1)
    s <- cl_fold_summary(fm, w)
    cvm <- s$m; cvse <- s$se
    if (measure == 'auc') {
      imin <- which.max(cvm)
      sel  <- if (select == 'min') imin else min(which(cvm >= cvm[imin] - cvse[imin]))
    }
    else {
      imin <- which.min(cvm)
      sel  <- if (select == 'min') imin else min(which(cvm <= cvm[imin] + cvse[imin]))
    }
  }
  else sel <- 1L

  B <- out$beta[, keep, drop = FALSE]
  pr <- out$pairs
  ps <- pr$l == sel
  pairs <- list(j = pr$j[ps], k = pr$k[ps], t = pr$t[ps])
  link <- drop(L %*% B[, sel]) + out$intercept[sel] + off
  list(lambda = lambda, coef_path = B, intercept = out$intercept[keep],
       covariates = cov, link = link, cv_error = cvm, cv_se = cvse,
       cv_link = if (K >= 2) out$cv_link[, keep, drop = FALSE] else NULL,
       sel = sel, measure = measure, pairs = pairs)
}


# glmnet's largest lambda: max over pairs of |<x_c, r>| / (n a_c alpha), r
# the null model's residuals (with the offset), a_c the pair's sd with
# standardize; alpha below 1e-3 counts as 1e-3.
cl_pw_lambda_max <- function (L, y, off, binom, alpha, standardize) {
  n <- nrow(L)
  if (binom) {
    b <- suppressWarnings(stats::glm.fit(matrix(1, n, 1), y, offset = off,
                                         family = stats::binomial())$coefficients)
    r <- y - stats::plogis(off + b)
  }
  else {
    yy <- y - off
    r  <- yy - mean(yy)
  }
  a  <- max(alpha, 1e-3)
  Lc <- sweep(L, 2, colMeans(L))
  g  <- drop(crossprod(Lc, r)) / n
  if (!standardize) return ((max(g) - min(g)) / a)
  dv <- colMeans(Lc^2)
  m  <- 0
  for (j in seq_len(ncol(L) - 1L)) {
    k  <- (j + 1L):ncol(L)
    cj <- drop(crossprod(Lc[, k, drop = FALSE], Lc[, j])) / n
    sd <- sqrt(pmax(dv[j] + dv[k] - 2 * cj, 0))
    ok <- sd > 0
    if (any(ok)) m <- max(m, abs(g[j] - g[k])[ok] / sd[ok])
  }
  m / a
}



#------------------------------------------------------------------------------
# predict()
#------------------------------------------------------------------------------

#' Predictions from a compositional lasso
#'
#' The fitted model's linear predictor (the signature score) or
#' probabilities for new samples.
#'
#' @param object   The result of [coda_lasso()] or [coda_lasso_pairwise()].
#'
#' @param counts   A feature table of the new samples, with (at least) the
#'   fit's selected features, matched by name (by position when the table
#'   has no feature names and as many features as the fit).
#'
#' @param data   Metadata of the new samples, for a fit with a `given`
#'   formula or column names. Default: `NULL`
#'
#' @param given   For a fit whose `given` was a data frame: the new
#'   samples' covariates, as a data frame with the same columns. Default:
#'   `NULL`
#'
#' @param margin   `1` if the new samples are rows, `2` if columns.
#'   Default: `1`
#'
#' @param type   `'link'`, the linear predictor (the log odds for a
#'   binomial fit), or `'response'` (the probability of the second class
#'   for a binomial fit; the same as `'link'` for gaussian). Default:
#'   `'link'`
#'
#' @param ...   Not used.
#'
#' @return A numeric vector, one value per new sample.
#'
#' @details The fit's `pseudocount` is added before taking logs (with
#'   none, zeros in the selected features are an error). Only the selected
#'   features enter, so the others needn't be in `counts`; because the
#'   coefficients sum to zero, each sample's total doesn't matter.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
#'     y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
#'     fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
#'     predict(fit, counts[1:3, ])
predict.ecodive_coda_lasso <- function (
    object, counts, data = NULL, given = NULL, margin = 1L,
    type = c('link', 'response'), ... ) {

  type <- match.arg(type)
  if (missing(counts)) stop('`counts` (the new samples) is required.', call. = FALSE)
  with(environment(), { validate_counts(); validate_margin() })

  feats <- object$feature
  beta  <- object$coefficient
  sel   <- which(beta != 0)
  pc    <- attr(object, 'pseudocount', exact = TRUE)

  dn <- switch(get_matrix_package(counts),
               base = dimnames(counts), slam = counts$dimnames, Matrix = counts@Dimnames)
  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) dn <- rev(dn)
  x <- .Call(C_normalize, counts, margin, NORM_NONE, 0, TRUE, 1L)
  if (margin == 2L) x <- t(x)
  dimnames(x) <- NULL
  nn <- nrow(x)
  cn <- dn[[2]]
  if (is.null(cn)) {
    if (ncol(x) != length(feats))
      stop('`counts` has no feature names and ', ncol(x), ' features; the fit has ',
           length(feats), '.', call. = FALSE)
    cols <- sel
  }
  else {
    cols <- match(feats[sel], cn)
    if (anyNA(cols))
      stop('`counts` lacks selected feature', if (sum(is.na(cols)) > 1) 's', ' ',
           name_list(feats[sel][is.na(cols)]), '.', call. = FALSE)
  }
  eta <- rep(attr(object, 'intercept', exact = TRUE), nn)
  if (length(sel)) {
    xs <- x[, cols, drop = FALSE] + pc
    if (any(xs <= 0))
      stop('`counts` has zeros in selected features, and the fit used no pseudocount.',
           call. = FALSE)
    eta <- eta + drop(log(xs) %*% beta[sel])
  }

  gm <- attr(object, 'given_model', exact = TRUE)
  if (!is.null(gm)) {
    df <- if (!is.null(given)) as.data.frame(given) else data
    if (is.null(df))
      stop('The fit adjusted for `given`; supply the new samples\' ',
           if (is.null(gm$cols)) '`data`' else '`given` (or `data`)', '.', call. = FALSE)
    df <- as.data.frame(df)
    if (nrow(df) != nn)
      stop('`data` (or `given`) must have one row per new sample.', call. = FALSE)
    mf <- stats::model.frame(gm$terms, df, na.action = stats::na.pass, xlev = gm$xlevels)
    mm <- stats::model.matrix(gm$terms, mf)
    cov <- attr(object, 'covariates', exact = TRUE)
    eta <- eta + drop(mm[, names(cov), drop = FALSE] %*% cov)
  }

  if (type == 'response' && attr(object, 'family', exact = TRUE) == 'binomial')
    eta <- stats::plogis(eta)
  names(eta) <- dn[[1]]
  eta
}

#' @rdname predict.ecodive_coda_lasso
#' @export
predict.ecodive_coda_lasso_pairwise <- predict.ecodive_coda_lasso



#------------------------------------------------------------------------------
# print()
#------------------------------------------------------------------------------

# The registry's display names (README: facts about a method by `method`).
CODA_LASSO_NAMES <- c(
  coda_lasso          = 'Compositional Lasso',
  coda_lasso_pairwise = 'Compositional Lasso (pairwise log ratios)' )

cl_name <- function (x) {
  m <- attr(x, 'method', exact = TRUE)
  if (is.null(m) || !m %in% names(CODA_LASSO_NAMES)) return (NULL)
  CODA_LASSO_NAMES[[m]]
}


#' @export
print.ecodive_coda_lasso <- function (x, n = 10L, ...) {
  nm <- cl_name(x)
  if (is.null(nm)) return (print_table_result(x, ...))  # rows subset with `[`
  sym <- print_sym()
  df  <- as.data.frame(unclass(x)[c('feature', 'coefficient', 'selected')],
                       stringsAsFactors = FALSE)
  fam  <- attr(x, 'family', exact = TRUE)
  perf <- attr(x, 'performance', exact = TRUE)
  cat(sprintf('%s, %s y: a data.frame %d %s %d\n', nm, fam, nrow(df), sym$times, 3L))
  cat(sprintf('lambda = %s (%s), l1_ratio = %s: %d of %d features selected\n',
              format(attr(x, 'lambda', exact = TRUE), digits = 3),
              attr(x, 'select', exact = TRUE),
              format(attr(x, 'l1_ratio', exact = TRUE), digits = 3),
              sum(df$selected), nrow(df)))
  cat(cl_perf_text(perf, fam), '\n', sep = '')

  side <- function (i, label) {
    if (!length(i)) return (cat(label, ': none\n', sep = ''))
    i <- i[order(-abs(df$coefficient[i]))]
    more <- length(i) - n
    i <- utils::head(i, n)
    s <- paste0(df$feature[i], ' (', formatC(df$coefficient[i], digits = 3, format = 'g'), ')')
    cat(strwrap(paste0(label, ': ', paste(s, collapse = ', '),
                       if (more > 0) sprintf(', and %d more', more)),
                exdent = 2), sep = '\n')
  }
  side(which(df$coefficient > 0), 'Numerator (coefficient > 0)')
  side(which(df$coefficient < 0), 'Denominator (coefficient < 0)')
  cat(print_wrap(paste('Attributes:', paste(print_attrs(x, c('names', 'row.names', 'class')),
                                           collapse = ', ')),
                 getOption('width', 80L)), sep = '\n')
  invisible(x)
}

#' @export
print.ecodive_coda_lasso_pairwise <- print.ecodive_coda_lasso


cl_perf_text <- function (perf, fam) {
  f <- function (v) formatC(v, digits = 3, format = 'f')
  if (fam == 'binomial')
    sprintf('Apparent AUC %s; cross-validated AUC %s +/- %s', f(perf[['auc']]),
            f(perf[['cv_auc']]), f(perf[['cv_auc_se']]))
  else
    sprintf('Apparent R^2 %s; cross-validated MSE %s +/- %s', f(perf[['r_squared']]),
            formatC(perf[['cv_mse']], digits = 3, format = 'g'),
            formatC(perf[['cv_mse_se']], digits = 3, format = 'g'))
}



#------------------------------------------------------------------------------
# plot()
#------------------------------------------------------------------------------

#' Plot a compositional lasso
#'
#' `type = 'signature'` draws coda4microbiome's two standard figures side
#' by side: the selected features' coefficients as bars (numerator,
#' positive, in red; denominator, negative, in blue), and the signature
#' score (the linear predictor) by outcome: a box per class for binomial
#' `y`, against the observed value for gaussian `y`. `type = 'path'` draws
#' the coefficients along the \eqn{\lambda} path, and `type = 'cv'` the
#' cross-validated error with its standard errors; both mark the chosen
#' \eqn{\lambda}.
#'
#' coda4microbiome rescales the coefficients to \eqn{\sum|\beta| = 2}
#' before drawing them; that changes neither their signs nor their order.
#'
#' @param x   The result of [coda_lasso()] or [coda_lasso_pairwise()].
#' @param type   `'signature'`, `'path'` or `'cv'`. Default: `'signature'`
#' @param n_features   The most features to show in the signature's bars
#'   (the largest \eqn{|\beta|}). Default: `30`
#' @param ...   Passed to the plotting functions (`main`, `xlab`, `col`,
#'   ...).
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
#'     y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
#'     fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
#'     plot(fit)
#'     plot(fit, type = 'path')
#'     plot(fit, type = 'cv')
plot.ecodive_coda_lasso <- function (x, type = c('signature', 'path', 'cv'),
                                     n_features = 30, ...) {
  type <- match.arg(type)
  args <- list(...)
  if (!is.numeric(n_features) || length(n_features) != 1 || is.na(n_features) ||
      n_features < 1)
    stop('`n_features` must be a positive number.', call. = FALSE)
  main <- cl_name(x)
  if (is.null(main)) stop('`x` is not a whole coda_lasso() result.', call. = FALSE)
  switch(type,
    signature = cl_plot_signature(x, as.integer(n_features), main, args),
    path      = cl_plot_path(x, main, args),
    cv        = cl_plot_cv(x, main, args) )
  invisible(x)
}

#' @rdname plot.ecodive_coda_lasso
#' @export
plot.ecodive_coda_lasso_pairwise <- plot.ecodive_coda_lasso


cl_plot_signature <- function (x, n_features, main, args) {

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(1, 1.2))

  b   <- x$coefficient
  nms <- x$feature
  sel <- which(b != 0)
  fam <- attr(x, 'family', exact = TRUE)
  fit <- attr(x, 'fitted', exact = TRUE)
  y   <- attr(x, 'y', exact = TRUE)

  if (!length(sel)) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, 'No features are selected\nat the chosen lambda.')
    graphics::plot.new()
    return (invisible())
  }

  shown <- sel[order(-abs(b[sel]))]
  sub   <- NULL
  if (length(shown) > n_features) {
    sub   <- sprintf('%d of %d selected features shown', n_features, length(sel))
    shown <- shown[seq_len(n_features)]
  }
  shown <- shown[order(b[shown])]      # most negative at the bottom
  cols  <- ifelse(b[shown] > 0, PLOT_HIGHLIGHT, 'royalblue3')
  lab_w <- max(graphics::strwidth(nms[shown], units = 'inches'))
  graphics::par(mai = c(graphics::par('mai')[1], min(lab_w + 0.2, 2.5),
                        graphics::par('mai')[3], 0.2))
  rng <- range(c(0, b[shown]))
  rng <- rng + c(-0.15, 0.15) * diff(rng)
  mp <- plot_do(
    graphics::barplot, c(list(b[shown]), args[setdiff(names(args), c('xlab', 'ylab'))]),
    list(horiz = TRUE, names.arg = nms[shown], las = 1, col = cols, border = NA,
         xlim = rng, main = main, sub = sub, xlab = 'Coefficient') )
  graphics::abline(v = 0, col = 'grey40')
  graphics::text(b[shown], mp, formatC(b[shown], digits = 2, format = 'f'),
                 pos = ifelse(b[shown] > 0, 4, 2), cex = 0.75, xpd = TRUE)

  # The signature score by outcome.
  graphics::par(mai = op$mai)
  perf <- attr(x, 'performance', exact = TRUE)
  ok   <- !is.na(fit$link)
  if (fam == 'binomial') {
    f  <- factor(y)
    bx <- as.integer(f)
    xy <- log_ratio_spread(bx[ok], fit$link[ok])
    graphics::plot.default(
      NA_real_, NA_real_, xlim = c(0.5, nlevels(f) + 0.5), ylim = range(fit$link[ok]),
      xaxt = 'n', las = 1, xlab = '', ylab = 'Signature score (log odds)',
      sub = cl_perf_text(perf, fam), main = 'Signature score by class' )
    graphics::boxplot(split(fit$link[ok], f[ok]), add = TRUE, at = seq_len(nlevels(f)),
                      outline = FALSE, border = 'grey60', col = NA, axes = FALSE, boxwex = 0.6)
    graphics::axis(1, at = seq_len(nlevels(f)), labels = levels(f))
    graphics::points(xy$x, xy$y, col = PLOT_MUTED)
  }
  else {
    graphics::plot.default(
      fit$link[ok], y[ok], las = 1, xlab = 'Signature score', ylab = 'Observed y',
      col = PLOT_MUTED, sub = cl_perf_text(perf, fam), main = 'Observed by signature score' )
  }
}


cl_plot_path <- function (x, main, args) {
  cp  <- attr(x, 'coef_path', exact = TRUE)
  lam <- attr(x, 'path', exact = TRUE)$lambda
  if (length(lam) < 2) stop('The fit has a single lambda: no path to plot.', call. = FALSE)
  nz  <- rowSums(cp != 0) > 0
  col <- ifelse(x$coefficient > 0, PLOT_HIGHLIGHT,
                ifelse(x$coefficient < 0, 'royalblue3', PLOT_MUTED))
  plot_do(
    graphics::matplot, c(list(log(lam), t(cp[nz, , drop = FALSE])), args),
    list(type = 'l', lty = 1, col = col[nz], las = 1, main = main,
         xlab = 'log(lambda)', ylab = 'Coefficient',
         xlim = rev(range(log(lam)))) )
  graphics::abline(h = 0, col = 'grey70')
  graphics::abline(v = log(attr(x, 'lambda', exact = TRUE)), lty = 2, col = 'grey40')
}


cl_plot_cv <- function (x, main, args) {
  path <- attr(x, 'path', exact = TRUE)
  if (all(is.na(path$cv_error)))
    stop('The fit has no cross-validation to plot.', call. = FALSE)
  ll <- log(path$lambda)
  lo <- path$cv_error - path$cv_se
  hi <- path$cv_error + path$cv_se
  ylab <- switch(attr(x, 'cv_measure', exact = TRUE),
                 mse = 'Mean squared error', deviance = 'Binomial deviance', auc = 'AUC')
  plot_do(
    graphics::plot.default, c(list(ll, path$cv_error), args),
    list(pch = 20, col = PLOT_HIGHLIGHT, las = 1, main = main,
         xlab = 'log(lambda)', ylab = ylab,
         ylim = range(c(lo, hi), na.rm = TRUE), xlim = rev(range(ll))) )
  graphics::segments(ll, lo, ll, hi, col = 'grey60')
  graphics::abline(v = log(attr(x, 'lambda', exact = TRUE)), lty = 2, col = 'grey40')
  graphics::axis(3, at = ll, labels = path$n_selected, tick = FALSE, cex.axis = 0.7,
                 line = -0.8)
}
