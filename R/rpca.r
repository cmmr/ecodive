# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Ported from gemelli 0.0.13 (gemelli/rpca.py), which carries this notice:
#
#   BSD 3-Clause License
#
#   Copyright (c) 2019, Cameron Martino
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   * Neither the name of the copyright holder nor the names of its
#     contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' Robust Aitchison PCA (RPCA)
#'
#' Ordinates samples straight from a count table, treating zeros as missing
#' values rather than replacing them with pseudocounts. Robust Aitchison PCA
#' (Martino et al. 2019) is the method of DEICODE and of QIIME 2's `gemelli`
#' plugin: a robust centered log-ratio transform of the non-zero counts,
#' completion of the missing entries with a low-rank matrix (OptSpace), and a
#' PCA of the completed matrix. It returns sample scores and feature loadings
#' for a biplot.
#'
#' @inheritParams documentation
#'
#' @param counts   A numeric matrix of count data (samples \eqn{\times}
#'        features), or any of the other objects listed in "Input Types".
#'        Values must be finite and non-negative; they need not be integers
#'        (proportions give the same result as counts).
#'
#' @param n_axes   The rank of the low-rank model, which is also the number of
#'        axes returned: a whole number from 1 to 10 (gemelli suggests
#'        fewer than 10). Default: `3`
#'
#' @param max_iters   The largest number of OptSpace gradient descent
#'        iterations; fewer are run if `tol` is reached first. More is not
#'        better: see "Iterations". Default: `5`
#'
#' @param tol   Stop early once the root mean squared error of the fit to the
#'        (rescaled) observed values falls below `tol`. Default: `1e-5`
#'
#' @param min_sample_count   Drop samples whose total count (after the
#'        feature filters) is not above this. Default: `0`
#'
#' @param min_feature_count   Drop features whose total count is not above
#'        this. Default: `0`
#'
#' @param min_feature_frequency   Drop features that are non-zero in no more
#'        than this percentage of samples (0 to 100). Default: `0`
#'
#' @param distances   If `TRUE`, also return gemelli's distances between
#'        samples, as a `dist` object in attribute `dm`. See "Distances"
#'        below. Default: `FALSE`
#'
#' @return A numeric matrix with one row per kept sample and one column per
#'   axis, named `"RPCA1"`, `"RPCA2"`, etc., as described in "Ordination
#'   result" below. Samples removed by the filters (see "Method") have no
#'   row. Rows are named after the samples, or numbered by their position in
#'   `counts` if it has no sample names, so that the kept samples can always
#'   be identified. Its attributes are:
#'
#'   * `method`: `"rpca"`.
#'   * `eig`: The variance of each axis: \eqn{s^2 / (n - 1)} for \eqn{n}
#'     kept samples, as for [pca()]. See "Relation to gemelli's results".
#'   * `prop`: The proportion of the centered completed matrix's total
#'     variance that each axis explains, `eig / trace`. This is gemelli's
#'     `proportion_explained`.
#'   * `trace`: The centered completed matrix's total variance: the sum of
#'     all `n_axes` axes' variances. (The completed matrix has rank `n_axes`.)
#'   * `loadings`: A matrix with one row per kept feature and one column per
#'     axis: the unit-length feature loadings, for a biplot.
#'   * `cost`: The root mean squared error of the low-rank fit to the
#'     observed values, after OptSpace's initialization and after each
#'     iteration (on OptSpace's rescaled data, as gemelli measures it).
#'   * `n_iters`: The number of OptSpace iterations run.
#'   * `dm`: Only with `distances = TRUE`. A `dist` object of distances
#'     between the kept samples, of class `c("ecodive_dist", "dist")`, as
#'     [beta_div()] returns.
#'
#' @details
#'
#' **Method**
#'
#' 1. **Filtering.** Features are kept if their total count is above
#'    `min_feature_count` and they are non-zero in more than
#'    `min_feature_frequency` percent of samples; then samples are kept if
#'    their total over those features is above `min_sample_count`; then
#'    features that are zero in every remaining sample are dropped. These
#'    are gemelli's filters with its defaults: at `0`, they only drop
#'    samples and features that are entirely zero. (gemelli's documentation
#'    suggests a `min_sample_count` of at least 500.) Samples with a single
#'    non-zero feature are also dropped, with a warning: see "Differences
#'    from gemelli". The result has rows for the kept samples only.
#'
#' 2. **Robust clr.** For each sample, the log of each non-zero count minus
#'    the mean log over its non-zero counts (as [normalize()] with
#'    `norm = 'rclr'`). Zeros are missing values: they are not part of the
#'    mean and are not given a value.
#'
#' 3. **Matrix completion.** OptSpace (Keshavan, Montanari & Oh 2010) fits a
#'    rank-`n_axes` matrix \eqn{USV^T} to the observed values only: it starts from
#'    their top `n_axes` singular vectors and runs up to `max_iters` steps of
#'    gradient descent on the Grassmann manifold, re-fitting the
#'    \eqn{k \times k} matrix \eqn{S} after each step.
#'
#' 4. **PCA.** The completed matrix \eqn{X = USV^T} is centered (columns,
#'    then rows) and decomposed as \eqn{X = u\,\mathrm{diag}(s)\,v^T}. The
#'    sample scores are \eqn{u\,\mathrm{diag}(s)} and the feature loadings
#'    \eqn{v}. So `ord %*% t(attr(ord, 'loadings'))` is the centered
#'    completed matrix, and [pcoa()] of the Euclidean distances between its
#'    rows would give the same scores.
#'
#' `rpca()` takes a count table, not distances, so like [pca()] it is called
#' directly rather than through [dist_ord()].
#'
#'
#' **Relation to gemelli's results**
#'
#' gemelli reports the singular values \eqn{s} as its "eigenvalues", the
#' unscaled \eqn{u} as its sample loadings, and \eqn{s^2 / \sum s^2} as the
#' proportion explained. Here `eig` is a variance, \eqn{s^2 / (n - 1)}, on the
#' same scale as [pca()]'s, and the scores are \eqn{u\,\mathrm{diag}(s)}, as
#' for [pca()] and [pcoa()]. So, for an ordination `ord` of \eqn{n} samples,
#' gemelli's eigenvalues are `sqrt(attr(ord, 'eig') * (n - 1))` and its
#' sample loadings are `sweep(ord, 2, sqrt(attr(ord, 'eig') * (n - 1)), '/')`.
#' `prop` and the feature loadings are the same as gemelli's.
#'
#'
#' **Speed and memory**
#'
#' gemelli works with dense \eqn{n \times p} arrays. `rpca()` never forms
#' one: every step works from the non-zero counts and the
#' \eqn{n \times k} and \eqn{p \times k} factors, in parallel C. The initial
#' singular vectors come from block Lanczos iteration using sparse products
#' (the same solver as [pca()] and [pcoa()]); the descent's cost and
#' gradients are sums over the non-zero counts; the line search's cost along
#' the search direction is a quartic polynomial whose coefficients take one
#' pass; and the final PCA reduces to \eqn{k \times k} problems.
#'
#' Each descent step takes time proportional to the number of non-zero
#' counts times \eqn{k^2}, plus \eqn{k^4/4} per sample for re-fitting
#' \eqn{S}. The initial SVD takes a number of steps that depends on the gap
#' between the `n_axes`-th and next singular values, each a pass over the
#' non-zero counts plus work on a basis of \eqn{\min(n, p)}-long vectors.
#' Memory beyond the input is about 24 bytes per non-zero count, plus about
#' \eqn{80 (n + p) k} bytes for the factors, plus about 1.5 kB per sample or
#' per feature, whichever there are fewer of, for the Lanczos basis.
#'
#' Threads (`cpus`) share every pass over the non-zero counts; small tables
#' use one thread. Results are identical from run to run and for any number
#' of `cpus`; R's random number generator is neither used nor changed. Each
#' axis is flipped so that its largest score in absolute value is positive,
#' with the loadings flipped to match.
#'
#'
#' **Biplots**
#'
#' The scores are \eqn{u\,\mathrm{diag}(s)} and the loadings are the unit
#' vectors \eqn{v}: a distance (form) biplot, in which distances between
#' samples are preserved and each feature's arrow shows the direction in
#' which its robust clr increases. (gemelli, and QIIME 2's Emperor, plot
#' \eqn{u} and \eqn{v}.) The loadings are much shorter than the scores;
#' `plot(ord)` scales them to fit among the samples and draws them as
#' arrows (see [plot.ecodive_ord()]):
#'
#' ```r
#' ord <- rpca(ex_counts, n_axes = 2)
#' plot(ord)
#' ```
#'
#'
#' **Iterations**
#'
#' OptSpace minimizes the error on the observed values only, so with many
#' iterations it can overfit them, and the completed values (and the
#' eigenvalues) can grow far beyond the data. For example, on a 20
#' \eqn{\times} 10 table with half its counts zero, the first singular value
#' was 15.5 after 5 iterations, 30.2 after 50, and 41.7 after 400; gemelli's
#' own OptSpace behaves the same. The default of 5 is gemelli's and is
#' recommended; `cost` shows how the fit changes.
#'
#'
#' **Distances**
#'
#' With `distances = TRUE`, `rpca()` also returns the distances that gemelli
#' returns with its ordination (which it describes as the Aitchison distance
#' of the sample loadings): Euclidean distances between the rows of
#' OptSpace's sample factor \eqn{U}, which are normalized so that each of its
#' `n_axes` columns has a squared length of about \eqn{n}, before the singular
#' values in \eqn{S} are applied. All `n_axes` axes count equally in them. They
#' are computed in parallel straight into a `dist` object, so they can go to
#' [pcoa()], [nmds()] or [permanova()], and nothing of size
#' \eqn{n \times n} is made without `distances = TRUE`.
#'
#' They differ from two other distances:
#'
#' * `dist(ord)`: Euclidean distances between the returned scores, i.e.
#'   between the samples' rows of the centered completed matrix, where each
#'   axis counts in proportion to its variance.
#' * [robust_aitchison()]: Euclidean distances between the samples' robust
#'   clr values with no completion: missing values count as 0, so each pair
#'   of samples is compared on the features observed in either, and no
#'   low-rank model is fitted. For a table with no zeros and `n_axes` at least
#'   the rank of its clr matrix, the completed matrix is the clr matrix, so
#'   `dist(ord)` is then the same as `robust_aitchison()` and [aitchison()].
#'
#'
#' **Differences from gemelli**
#'
#' `rpca()` is a port of `gemelli.rpca.rpca()` (gemelli 0.0.13), and gives
#' its results to within rounding error and the initial SVD's tolerance, with
#' these exceptions. On the package's test tables, results agree with
#' gemelli's to about \eqn{10^{-13}} (relative). Ill-conditioned problems
#' agree less closely: when the `n_axes`-th and next singular values of the rclr
#' matrix are nearly tied, the initial subspace (in gemelli as here) is only
#' determined to about \eqn{10^{-15}} divided by their relative gap, and
#' results can then differ from gemelli's by \eqn{10^{-10}} or more.
#'
#' * gemelli's OptSpace treats any robust clr value that is exactly 0 as
#'   missing, which depends on rounding. Here a value is missing exactly
#'   when its count is 0.
#' * gemelli keeps samples with a single non-zero feature. Their robust clr
#'   is 0, which gemelli treats as missing, so they have no observed values
#'   and their position only reflects the other samples. `rpca()` drops them
#'   with a warning.
#' * gemelli requires \eqn{2 \le k \le \min(n, p) - 1}. `rpca()` accepts
#'   `n_axes = 1`, and reduces `n_axes` above \eqn{\min(n, p) - 1} with a warning.
#'   It also limits `n_axes` to 10, because re-fitting \eqn{S} takes time in
#'   proportion to \eqn{k^4}; gemelli has no limit but suggests
#'   \eqn{k < 10}. gemelli adds an axis of zeros when \eqn{k = 2};
#'   `rpca()` doesn't.
#' * gemelli's `rpca()` takes `max_iterations` but doesn't pass it on, so it
#'   always runs 5 iterations; here `max_iters` is used. `tol` is OptSpace's
#'   own.
#' * `eig` and the scores are scaled as described in "Relation to gemelli's
#'   results".
#' * Axes whose variance is at rounding-error level (below
#'   `sqrt(.Machine$double.eps)` times the larger of the first axis's and
#'   the rclr values' sum of squares, i.e. variation below about
#'   \eqn{10^{-4}} of the data's in amplitude) are dropped with a warning;
#'   if no axis is left (e.g. all samples have the same composition),
#'   `rpca()` stops with an error. gemelli returns such axes.
#' * The line search's trial costs are computed from a polynomial rather
#'   than by re-summing, which differs only by rounding; in a borderline
#'   case that can change which step length is accepted.
#'
#' @references
#' Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
#' Knight, R., & Zengler, K. (2019). A novel sparse compositional technique
#' reveals microbial perturbations. *mSystems*, 4(1), e00016-19.
#' \doi{10.1128/mSystems.00016-19}
#'
#' Keshavan, R. H., Montanari, A., & Oh, S. (2010). Matrix completion from a
#' few entries. *IEEE Transactions on Information Theory*, 56(6), 2980-2998.
#' \doi{10.1109/TIT.2010.2046205}
#'
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman & Hall, London.
#'
#' Martino, C. (2019). *gemelli* (Python package, version 0.0.13). BSD
#' 3-Clause License. \url{https://pypi.org/project/gemelli/}
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [pca()], which with `norm = 'rclr'` keeps zeros as 0 rather than
#'   completing them; [rpca_phylo()], the same over the nodes of a
#'   phylogenetic tree; [pcoa()]; [robust_aitchison()].
#'
#' @export
#' @examples
#'     ord <- rpca(iris[, 1:4])
#'     plot(ord, 'Species', iris)
#'     ord
#'     str(ord)
#'
rpca <- function (
    counts,
    margin                = 1L,
    n_axes                = 3L,
    max_iters             = 5L,
    tol                   = 1e-5,
    min_sample_count      = 0,
    min_feature_count     = 0,
    min_feature_frequency = 0,
    distances             = FALSE,
    cpus                  = n_cpus() ) {

  validate_counts()
  validate_margin()
  validate_cpus()

  k         <- ord_arg(n_axes, '`n_axes` must be a whole number from 1 to 10.', whole = TRUE, lo = 1, hi = 10)
  max_iters <- ord_arg(max_iters, '`max_iters` must be a whole number from 1 to 1000000.', whole = TRUE, lo = 1, hi = 1e6)
  tol       <- ord_arg(tol,       '`tol` must be a number of at least 0.')

  min_sample_count      <- ord_arg(min_sample_count,      '`min_sample_count` must be a number of at least 0.')
  min_feature_count     <- ord_arg(min_feature_count,     '`min_feature_count` must be a number of at least 0.')
  min_feature_frequency <- ord_arg(min_feature_frequency, '`min_feature_frequency` must be a number from 0 to 100.', hi = 100)

  if (!isTRUE(distances) && !isFALSE(distances))
    stop('`distances` must be TRUE or FALSE.', call. = FALSE)


  # Sample and feature names, with samples first. Unnamed ones are numbered,
  # so that the kept samples and features can be told apart.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  sam_names <- if (is.null(dn[[1]])) as.character(seq_len(dims[[1]])) else as.character(dn[[1]])
  fea_names <- if (is.null(dn[[2]])) as.character(seq_len(dims[[2]])) else as.character(dn[[2]])


  res <- .Call(
    C_rpca, counts, margin, k, max_iters, tol,
    min_sample_count, min_feature_count, min_feature_frequency,
    distances, FALSE, 1e4, cpus )


  if (length(res$single) > 0)
    warning(
      call. = FALSE,
      length(res$single), ' sample(s) with a single non-zero feature ',
      'removed (their robust clr is 0, which carries no information): ',
      rpca_list(sam_names[res$single]) )

  rpca_finish(
    res, k, 'rpca', sam_names[res$samples], fea_names[res$features],
    'features', distances )
}



# The ordination from C_rpca()'s or C_rpca_phylo()'s list `res`: warnings,
# dropping rounding-error axes, and the result. `k` is the rank asked for,
# `samples` and `features` the kept samples' and features' names, `what` the
# word for features in a warning, and `...` further attributes.
rpca_finish <- function (res, k, method, samples, features, what, distances, ...) {

  if (!res$lanczos_converged) # nocov start
    warning(
      call. = FALSE,
      "OptSpace's initial SVD did not converge; results may differ from gemelli's." )
  # nocov end

  if (res$k < k)
    warning(
      call. = FALSE,
      '`n_axes` reduced to ', res$k, ': at most one less than the number of ',
      'samples and of ', what, ' (', length(samples), ' and ',
      length(features), ' after filtering).' )

  n     <- length(samples)
  s2    <- res$eig    # s^2, largest first
  total <- sum(s2)

  # Axes with rounding-error variance, measured against the variance in
  # the data (the rclr values' sum of squares) or on the first axis.
  keep <- s2 > sqrt(.Machine$double.eps) * max(res$ss, s2[[1]])
  if (!any(keep))
    stop(
      call. = FALSE,
      'There is no variation to ordinate: after the robust clr transform, ',
      'all samples are the same (to within rounding error).' )
  if (!all(keep))
    warning(sprintf('only %d of the first %d eigenvalues are > 0', sum(keep), length(keep)), call. = FALSE)

  s2       <- s2[keep]
  loadings <- res$loadings[, keep, drop = FALSE]
  dimnames(loadings) <- list(features, paste0('RPCA', seq_len(sum(keep))))

  dm <- NULL
  if (distances)
    dm <- structure(
      .Data   = res$dist,
      Size    = n,
      Labels  = samples,
      Diag    = FALSE,
      Upper   = FALSE,
      method  = method,
      display = list(title =
        if (method == 'rpca_phylo') 'Phylogenetic Robust Aitchison Distance'
        else 'Robust Aitchison Distance (RPCA)' ),
      class   = c('ecodive_dist', 'dist') )

  denom <- n - 1L

  ord_result(
    res$scores[, keep, drop = FALSE], method, samples, prefix = 'RPCA',
    eig      = s2 / denom,
    prop     = s2 / total,
    trace    = total / denom,
    loadings = loadings,
    cost     = res$cost,
    n_iters  = res$n_iters,
    dm       = dm,
    ... )
}



# "a, b, c" or "a, b, c, ... (n more)".
rpca_list <- function (x, max = 5L) {
  if (length(x) <= max) return (paste(x, collapse = ', '))
  paste0(paste(x[seq_len(max)], collapse = ', '), ', ... (', length(x) - max, ' more)')
}
