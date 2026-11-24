# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# balance() is one coordinate of ilr() with a sequential binary partition:
# its basis row comes from sbp_basis() (R/log_ratio_transforms.r), which
# follows sbp_basis() from scikit-bio 0.7.4 (skbio/stats/composition/
# _base.py), Copyright (c) 2013--, scikit-bio development team, under the
# BSD-3 license given in src/log_ratio.c and inst/COPYRIGHTS.



#' Log Ratio of Summed Abundances
#'
#' One value per sample: the log of the summed abundance of the `numerator`
#' features over that of the `denominator` features, as Qurro plots it.
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of a zero replacement
#'        are all accepted: each sample's scale cancels out of the ratio.
#'
#' @param numerator,denominator  The features on the top and the bottom of
#'        the ratio, each given as feature names, (1-based) feature indices,
#'        or a logical vector with one value per feature. Both are required,
#'        each must select at least one feature, and they must not share
#'        any. The two may use different forms (names for one, indices for
#'        the other).
#'
#' @param pseudocount  Added to every selected value before summing, so
#'        \eqn{r c} is added to a numerator of \eqn{r} features. `0` (the
#'        default) adds nothing; `NULL` applies the rule of the other
#'        log-ratio functions (half the smallest non-zero value of the whole
#'        table, with a warning, when it has zeros). See **Zeros**.
#'
#' @return A numeric vector with one value per sample, in the order of
#'   `counts` and named after the samples (unnamed if they have no names),
#'   whatever the `margin`. Larger values mean relatively more of the
#'   numerator. Attributes:
#'
#'   * `method`: `'log_ratio'`.
#'   * `numerator`, `denominator`: the selected features' names (their
#'     indices when the features have no names), in the order of `counts`.
#'   * `pseudocount`: the pseudocount used (`0` for none).
#'
#'   Its class, `ecodive_log_ratio`, gives it a short `print()` (a header
#'   naming the features, then the values) and a `plot()` by group or
#'   against a numeric variable ([plot.ecodive_log_ratio()]). Everything
#'   else treats it as a plain numeric vector (`is.numeric()` is `TRUE`);
#'   subsetting (`x[1:5]`), `sort()`, and arithmetic or maths (`-x`, `x -
#'   mean(x)`, `abs(x)`) return a plain named vector without the
#'   attributes, which no longer describe the values.
#'
#'   The vector drops straight into metadata (`data$ratio <- log_ratio(...)`)
#'   and into any test or fit that takes a per-sample variable, such as
#'   [envfit()] or [mantel()].
#'
#' @details
#' For a sample \eqn{x}, numerator features \eqn{N} and denominator features
#' \eqn{D}, the log ratio is
#' \deqn{\log\frac{\sum_{i \in N} (x_i + c)}{\sum_{j \in D} (x_j + c)}}
#' with \eqn{c} the pseudocount (natural logs). With one feature on each
#' side, this is the log ratio of two features, and equals the matching
#' [balance()] divided by \eqn{\sqrt{1/2}}.
#'
#' Ratios between features can be interpreted without knowing each sample's
#' total load: if the ratio of A to B goes up, it went up whatever happened
#' to the total (the "reference frames" of Morton et al. 2019). Qurro
#' (Fedarko et al. 2020) is a tool for choosing the numerator and
#' denominator interactively from feature rankings; this function computes
#' the same log ratio (`log(Num_Sum / Denom_Sum)` in Qurro's `qarcoal`) for
#' use in R.
#'
#' **Log ratio or balance?** Sums are robust to zeros in single features and
#' easy to explain ("the log of the ratio of total Firmicutes to total
#' Bacteroidetes"), but a sum is dominated by its most abundant features.
#' [balance()] uses geometric means instead: every feature counts equally,
#' it is a coordinate of the isometric log-ratio geometry the CoDa theory is
#' built on, and its scale makes balances of different sizes comparable;
#' but a zero in any one of its features needs replacing first.
#'
#' @section Zeros:
#'   Zeros in single features are harmless as long as each sum is positive,
#'   so no pseudocount is needed in almost every case, and none is used by
#'   default (unlike the other log-ratio functions, whose default is
#'   `NULL`). A sample whose numerator or denominator sums to 0 has no
#'   defined ratio: its value is `NA`, with one warning naming such samples.
#'   (Qurro drops those samples; ecodive keeps their place, so the result
#'   lines up with `counts`.) A positive `pseudocount` gives them a value.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Morton, J. T., Marotz, C., Washburne, A., et al. (2019). Establishing
#' microbial composition measurement standards with reference frames.
#' *Nature Communications*, 10, 2719.
#'
#' Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
#' 'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
#' Bioinformatics*, 2, lqaa023.
#'
#' @seealso [balance()] for the ratio of geometric means; [alr()] for every
#'   feature against one; [vlr()], [rho_p()] and [phi_s()] for finding
#'   features that move together. [log_ratio_features()] to take the
#'   numerator and denominator from the two ends of a ranking, and
#'   [plot.ecodive_log_ratio()] to plot the result.
#'
#' @export
#' @examples
#'     # Oral streptococci and Haemophilus against skin taxa
#'     r <- log_ratio(
#'       ex_counts,
#'       numerator   = c('Streptococcus', 'Haemophilus'),
#'       denominator = c('Propionibacterium', 'Staphylococcus', 'Corynebacterium') )
#'     r
#'
#'     # The same, by index; zero sums give NA
#'     log_ratio(ex_counts, c(1, 4), c(3, 5, 6))
#'     suppressWarnings(log_ratio(ex_counts, 2, 6))
#'
log_ratio <- function (
    counts, numerator, denominator, margin = 1L, pseudocount = 0 ) {

  ratio_args_given(missing(numerator), missing(denominator))
  info <- log_ratio_input('a log ratio', zero_ok = TRUE, dense_out = FALSE)
  sel  <- ratio_parts(numerator, denominator, info)

  res <- .Call(
    C_log_ratio, counts, margin, sel$num - 1L, sel$den - 1L, pseudocount )

  ratio_na_warning(res, info,
    'a numerator or denominator that sums to 0, so %s NA.',
    c('its log ratio is', 'their log ratios are'),
    'A positive `pseudocount` gives %s a value.' )

  ratio_result(res, 'log_ratio', sel, info, pseudocount)
}



#' Isometric Log-Ratio Balance
#'
#' One value per sample: the balance between two groups of features, the
#' scaled log ratio of their geometric means. This is an isometric
#' log-ratio balance (Egozcue & Pawlowsky-Glahn 2005), not a balanced
#' design.
#'
#' @inheritParams log_ratio
#'
#' @param pseudocount  Added to every value of `counts` before taking logs.
#'        `NULL` (the default) adds nothing when `counts` has no zeros, and
#'        otherwise uses half the smallest non-zero value of the whole
#'        table, with a warning, as [ilr()] does. `0` adds nothing, and
#'        gives `NA` for samples with a zero in a selected feature. See
#'        **Pseudocount**.
#'
#' @return A numeric vector with one value per sample, in the order of
#'   `counts` and named after the samples (unnamed if they have no names),
#'   whatever the `margin`. Larger values mean relatively more of the
#'   numerator. Attributes:
#'
#'   * `method`: `'balance'`.
#'   * `numerator`, `denominator`: the selected features' names (their
#'     indices when the features have no names), in the order of `counts`.
#'   * `pseudocount`: the pseudocount used (`0` for none).
#'
#'   Its class, `ecodive_log_ratio`, gives it a short `print()` (a header
#'   naming the features, then the values) and a `plot()` by group or
#'   against a numeric variable ([plot.ecodive_log_ratio()]). Everything
#'   else treats it as a plain numeric vector (`is.numeric()` is `TRUE`);
#'   subsetting (`x[1:5]`), `sort()`, and arithmetic or maths (`-x`, `x -
#'   mean(x)`, `abs(x)`) return a plain named vector without the
#'   attributes, which no longer describe the values.
#'
#' @details
#' For a sample \eqn{x}, numerator features \eqn{N} (\eqn{r} of them) and
#' denominator features \eqn{D} (\eqn{s} of them),
#' \deqn{\mathrm{balance}(x) = \sqrt{\frac{rs}{r + s}}\,\log\frac{g(x_N +
#' c)}{g(x_D + c)}}
#' where \eqn{g} is the geometric mean and \eqn{c} the pseudocount (natural
#' logs). This is the inner product of \eqn{\mathrm{clr}(x)} with the basis
#' vector that has \eqn{\sqrt{s / (r (r + s))}} on the numerator features
#' and \eqn{-\sqrt{r / (s (r + s))}} on the denominator ones: one row of
#' scikit-bio's `sbp_basis()`. So `balance()` equals, to the last bit, the
#' coordinate of [ilr()] whose row of a sequential binary partition has `1`
#' on `numerator` and `-1` on `denominator` (with the same `pseudocount`),
#' without building the full basis. Features in neither group don't affect
#' the value.
#'
#' The factor \eqn{\sqrt{rs/(r+s)}} puts balances of different sizes on the
#' same scale: each is a coordinate of an orthonormal basis, so their
#' variances can be compared directly. Swapping `numerator` and
#' `denominator` flips the sign. Like every log ratio, a balance is
#' unchanged when a sample is multiplied by a constant (with no
#' pseudocount). See [log_ratio()] for when to prefer the ratio of sums.
#'
#' Balances are what selbal (Rivera-Pinto et al. 2018) selects, gneiss and
#' PhILR build from trees, and coda4microbiome generalizes to log-contrasts.
#'
#' **Same names in other packages**: `compositions::balance()` computes
#' the same balances, with the same sign, from a formula (`~ (A*B)/C`);
#' whichever package was attached last masks the other's. Write
#' `ecodive::balance()` to be sure.
#'
#' @section Pseudocount:
#'   A zero in any numerator or denominator feature makes its geometric mean
#'   0, so zeros need replacing first. With `pseudocount = NULL` (the
#'   default), no pseudocount is used when `counts` has no zeros; otherwise
#'   half the smallest non-zero value of the **whole table** is added to
#'   every value, with a warning, exactly as [ilr()] does, so a balance
#'   agrees with the matching `ilr()` coordinate. To avoid a pseudocount,
#'   replace the zeros first with [replace_zeros()]: its result has none.
#'   With `pseudocount = 0`, samples with a zero in a selected feature get
#'   `NA`, with one warning naming them.
#'
#'   The choice of pseudocount is not neutral: it can distort the results,
#'   especially for sparse tables. See [aitchison()] for references.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
#' balances in compositional data analysis. *Mathematical Geology*, 37,
#' 795–828.
#'
#' Rivera-Pinto, J., Egozcue, J. J., Pawlowsky-Glahn, V., et al. (2018).
#' Balances: a new perspective for microbiome analysis. *mSystems*, 3,
#' e00053-18.
#'
#' @seealso [log_ratio()] for the ratio of sums; [ilr()] for all the
#'   balances of a sequential binary partition at once; [replace_zeros()];
#'   [vlr()], [rho_p()] and [phi_s()] for finding features that move
#'   together.
#'
#' @export
#' @examples
#'     # Replace zeros first, then take the balance
#'     x <- replace_zeros(ex_counts, 'mult')
#'     b <- balance(
#'       x,
#'       numerator   = c('Streptococcus', 'Haemophilus'),
#'       denominator = c('Propionibacterium', 'Staphylococcus', 'Corynebacterium') )
#'     b
#'
#'     # The first coordinate of ilr() with a partition that starts with it
#'     sbp <- rbind(
#'       c(1,  0, -1,  1, -1, -1),
#'       c(1,  0,  0, -1,  0,  0),
#'       c(0,  0,  1,  0, -1, -1),
#'       c(0,  0,  0,  0,  1, -1),
#'       c(1, -1,  1,  1,  1,  1) )
#'     all.equal(unname(b), unname(ilr(x, basis = sbp)[, 1]))
#'
balance <- function (
    counts, numerator, denominator, margin = 1L, pseudocount = NULL ) {

  ratio_args_given(missing(numerator), missing(denominator))
  info <- log_ratio_input('a balance', zero_ok = TRUE, dense_out = FALSE)
  sel  <- ratio_parts(numerator, denominator, info)

  # The balance's row of sbp_basis(), as ilr_basis() stores it.
  sbp <- matrix(0, 1L, info$D)
  sbp[sel$num] <-  1
  sbp[sel$den] <- -1
  row  <- sbp_basis(sbp)
  feat <- which(row != 0)

  res <- .Call(
    C_balance, counts, margin, feat - 1L, row[feat], pseudocount )

  ratio_na_warning(res, info,
    'a zero in a numerator or denominator feature, so %s NA.',
    c('its balance is', 'their balances are'),
    'Use `pseudocount = NULL`, or replace_zeros() first, to give %s a value.' )

  ratio_result(res, 'balance', sel, info, pseudocount)
}



# Both sides of a ratio are required.
ratio_args_given <- function (no_num, no_den) {
  if (no_num || no_den)
    stop(
      '`', if (no_num) 'numerator' else 'denominator', '` is required: ',
      'the features on the ', if (no_num) 'top' else 'bottom', ' of the ratio, ',
      'as names, indices, or a logical vector.', call. = FALSE )
}



# The numerator's and denominator's features, checked: each a non-empty set
# (feature_set()) and the two disjoint.
ratio_parts <- function (numerator, denominator, info) {

  num <- feature_set(numerator,   info$features, info$D, 'numerator')
  den <- feature_set(denominator, info$features, info$D, 'denominator')

  shared <- intersect(num, den)
  if (length(shared))
    stop(sprintf(
      '`numerator` and `denominator` must not share features, but they share %s: %s.',
      if (length(shared) == 1L) '1 feature' else paste(length(shared), 'features'),
      name_list(if (is.null(info$features)) shared else info$features[shared]) ),
      call. = FALSE )

  list(num = num, den = den)
}



# The 1-based indices, increasing, of a set of features given by name,
# index or a logical vector with one value per feature. At least one, and
# none twice. (feature_index() is the one-feature version, for alr().)
# `src` names the features' source in messages.
feature_set <- function (x, features, D, arg, src = '`counts`') {

  msg <- sprintf(
    '`%s` must select features: their names, their 1-based indices, or a logical vector with one value per feature.',
    arg )

  if (is.factor(x)) x <- as.character(x)
  if (is.null(x) || !is.atomic(x) || !(is.character(x) || is.logical(x) || is.numeric(x)))
    stop(msg, call. = FALSE)
  if (length(x) == 0L)
    stop(sprintf('`%s` selects no features.', arg), call. = FALSE)
  if (anyNA(x))
    stop(sprintf('`%s` must not contain NA.', arg), call. = FALSE)

  if (is.logical(x)) {
    if (length(x) != D)
      stop(sprintf(
        'A logical `%s` needs one value per feature: %s has %d features, not %d.',
        arg, src, D, length(x) ), call. = FALSE)
    if (!any(x))
      stop(sprintf('`%s` selects no features: it is all FALSE.', arg), call. = FALSE)
    return (which(x))
  }

  if (is.character(x)) {
    if (is.null(features))
      stop(sprintf('%s has no feature names; give `%s` as indices.', plot_cap(src), arg), call. = FALSE)
    i <- match(x, features)
    if (anyNA(i)) {
      bad <- unique(x[is.na(i)])
      stop(sprintf(
        '`%s` names %s not in %s: %s.', arg,
        if (length(bad) == 1L) 'a feature' else paste(length(bad), 'features'),
        src, name_list(bad) ), call. = FALSE )
    }
  }
  else {
    if (any(!is.finite(x)) || any(x %% 1 != 0))
      stop(msg, call. = FALSE)
    bad <- unique(x[x < 1 | x > D])
    if (length(bad))
      stop(sprintf(
        '`%s` has %s outside 1 to %d (%s has %d features): %s.', arg,
        if (length(bad) == 1L) 'an index' else 'indices', D, src, D,
        paste(utils::head(bad, 5), collapse = ', ') ), call. = FALSE )
    i <- as.integer(x)
  }

  if (anyDuplicated(i)) {
    dup <- unique(x[duplicated(i)])
    stop(sprintf(
      '`%s` selects %s more than once: %s.', arg,
      if (length(dup) == 1L) 'a feature' else 'features',
      if (is.character(dup)) name_list(dup) else paste(utils::head(dup, 5), collapse = ', ') ),
      call. = FALSE )
  }

  sort(i)
}



# One warning naming the samples whose value is NA. `why` and `fix` each
# have one %s, filled with `what` (singular, plural) and 'it' / 'them'.
ratio_na_warning <- function (res, info, why, what, fix) {
  na <- which(is.na(res))
  if (!length(na)) return (invisible())
  one <- length(na) == 1L
  ids <- if (is.null(info$samples)) na else info$samples[na]
  warning(sprintf(
    '%d %s (%s) %s %s %s', length(na), if (one) 'sample' else 'samples',
    name_list(ids), if (one) 'has' else 'have',
    sprintf(why, what[[2L - one]]), sprintf(fix, if (one) 'it' else 'them') ),
    call. = FALSE )
}



# The per-sample vector with its names and attributes.
ratio_result <- function (res, method, sel, info, pseudocount) {
  names(res) <- info$samples
  lab <- function (i) if (is.null(info$features)) i else info$features[i]
  attr(res, 'method')      <- method
  attr(res, 'numerator')   <- lab(sel$num)
  attr(res, 'denominator') <- lab(sel$den)
  attr(res, 'pseudocount') <- pseudocount
  class(res) <- 'ecodive_log_ratio'
  res
}



#' Features for a Log Ratio from a Ranking
#'
#' The `n` highest- and `n` lowest-ranked features of a differential
#' abundance result or of an ordination's loadings, as the `numerator` and
#' `denominator` of [log_ratio()] or [balance()]: Qurro's "autoselection"
#' (Fedarko et al. 2020), after the reference frames of Morton et al.
#' (2019).
#'
#' @param x  What ranks the features:
#'
#'   * a differential abundance result with a signed effect: from
#'     [dirmult_ttest()], [ancombc()], [ancombc2()], [da_ttest()],
#'     [da_wilcoxon()], [da_lm()], or [linda()] (any method whose `effect`
#'     column [list_da_tests()] names); its effect column ranks them;
#'   * an ordination with feature loadings, such as [pca()], [rpca()],
#'     [rpca_phylo()] or [rda()]: the loadings on axis `axes` rank them;
#'   * a numeric vector of effects, named after the features.
#'
#' @param n  How many features on each side: a whole number of at least 1,
#'        with `2 * n` no more than the number of ranked features. Default
#'        `5`.
#'
#' @param covariate  For a model result ([ancombc()], [ancombc2()],
#'        [da_lm()], [linda()]): which covariate's effects rank the
#'        features, as named in the result's `covariate` column. Defaults to
#'        the one their plots show (see [plot.ecodive_ancombc()]).
#'
#' @param axes  For an ordination: which axis's loadings rank the features,
#'        one axis number. Default `1`.
#'
#' @return A list of two vectors, `numerator` (the `n` features with the
#'   largest effects, largest first) and `denominator` (the `n` with the
#'   smallest, smallest first): feature names, or indices into the
#'   features when they have none. So `do.call(log_ratio, c(list(counts),
#'   sel))` is the log ratio of the top over the bottom of the ranking.
#'
#' @details
#' Features whose effect isn't finite (an estimate that couldn't be made,
#' `NaN`) are left out first. The rest are sorted by their effect,
#' ascending, with ties left in the result's order (Qurro's stable sort);
#' the last `n` are the numerator and the first `n` the denominator.
#' [ancom()] has no signed effect, and so no reference frame: its `W`
#' counts significant log ratios without a direction.
#'
#' **Phylogenetic RPCA.** The loadings of [rpca_phylo()] include the tree's
#' internal nodes (`"n12"`), which aren't columns of the count table. Only
#' the tips are ranked, with a message. To look at log ratios of clades,
#' run [rpca_phylo()] with `node_counts = TRUE` and rank its loadings as a
#' named vector: `log_ratio_features(attr(ord, 'loadings')[, 1], n)` with
#' `log_ratio(attr(ord, 'node_counts'), ...)`.
#'
#' **Circularity.** Choosing features by how much they differ between two
#' groups, then testing their log ratio between the same groups, finds the
#' difference by construction: the p-value means nothing. Use the log
#' ratio to look at the samples, or select its features on one subset of
#' the samples and test it on another.
#'
#' @references
#' Morton, J. T., Marotz, C., Washburne, A., et al. (2019). Establishing
#' microbial composition measurement standards with reference frames.
#' *Nature Communications*, 10, 2719. \doi{10.1038/s41467-019-10656-5}
#'
#' Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
#' 'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
#' Bioinformatics*, 2, lqaa023. \doi{10.1093/nargab/lqaa023}
#'
#' @seealso [log_ratio()], [balance()], and their plot
#'   ([plot.ecodive_log_ratio()]); the rank plots of DA results and
#'   ordinations (`type = "rank"`, see [plot.ecodive_ancombc()]).
#'
#' @export
#' @examples
#'     x <- rbind(
#'       c(1, 10, 16, 0, 0, 12, 0, 7), c(0, 10, 45, 2, 1, 25, 1, 9),
#'       c(4,  8, 22, 1, 5, 14, 1, 5), c(1, 11, 22, 1, 4, 14, 1, 3),
#'       c(4,  7, 36, 0, 0, 19, 0, 9), c(2, 30, 26, 0, 0, 12, 1, 0),
#'       c(1, 36, 10, 1, 5, 15, 2, 9), c(1, 25, 11, 2, 2,  8, 0, 5),
#'       c(1, 36, 23, 0, 2, 11, 1, 2), c(3, 60, 33, 3, 4, 14, 1, 6),
#'       c(3, 57, 30, 3, 1,  3, 0, 4), c(0, 65, 23, 0, 2,  2, 0, 8),
#'       c(2, 46, 19, 0, 1,  2, 0, 1), c(0, 41, 26, 0, 2,  0, 0, 3),
#'       c(3, 46, 20, 0, 0,  3, 0, 0) )
#'     dimnames(x) <- list(paste0('S', 1:15), paste0('F', 1:8))
#'     meta <- data.frame(
#'       status = rep(c('mild', 'moderate', 'severe'), each = 5),
#'       row.names = rownames(x) )
#'     fit <- ancombc2(x, ~ status, meta, groups = 'status')
#'
#'     # The two features that rose most with moderate status, over the two
#'     # that fell most.
#'     sel <- log_ratio_features(fit, n = 2)
#'     sel
#'
#'     # Where they sit in the ranking, and their log ratio by status.
#'     plot(fit, type = 'rank', numerator = sel$numerator, denominator = sel$denominator)
#'     plot(do.call(log_ratio, c(list(x), sel)), 'status', meta)
#'
#'     # From an ordination's loadings.
#'     log_ratio_features(pca(iris[, 1:4]), n = 1)
log_ratio_features <- function (x, n = 5L, covariate = NULL, axes = 1L) {

  src <- rank_source(x, covariate, axes, has_axes = !missing(axes))
  eff <- src$effect
  ok  <- is.finite(eff)

  # rpca_phylo(): only the tips are columns of the count table.
  tree <- if (inherits(x, 'ecodive_ord') && identical(attr(x, 'method', exact = TRUE), 'rpca_phylo'))
    attr(x, 'tree', exact = TRUE)
  if (!is.null(tree)) {
    tip <- src$feature %in% tree$tip.label
    if (any(!tip & ok))
      message(sprintf(
        'Ranking the %d tips only: the %d internal nodes of rpca_phylo() aren\'t features of the count table.',
        sum(tip & ok), sum(!tip & ok) ))
    ok <- ok & tip
  }

  m <- sum(ok)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != floor(n) || n < 1)
    stop('`n` must be a whole number of at least 1.', call. = FALSE)
  if (2 * n > m)
    stop(sprintf(
      '`n` must be at most %d: %d %s ranked, and `2 * n` must not exceed that.',
      m %/% 2L, m, if (m == 1L) 'feature is' else 'features are' ), call. = FALSE)
  n <- as.integer(n)

  i <- which(ok)
  o <- i[order(eff[i])]
  m <- length(o)
  list(
    numerator   = src$feature[rev(o[(m - n + 1L):m])],
    denominator = src$feature[o[seq_len(n)]] )
}
