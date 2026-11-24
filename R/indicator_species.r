# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Written from the papers (Dufrene & Legendre 1997; De Caceres & Legendre
# 2009; De Caceres, Legendre & Moretti 2010; Tichy & Chytry 2006) and from
# the outputs and documentation of indicspecies::multipatt() and
# strassoc() (GPL; no code from them). See plans/indicator_species.md.



#' Indicator Species Analysis
#'
#' Finds the features (species, taxa) that characterise a group of samples,
#' or a combination of groups: found mostly in its samples (specificity)
#' and in most of them (fidelity). For each feature, the association with
#' every group and combination of groups is measured by the indicator value
#' (IndVal, Dufrêne & Legendre 1997) or a correlation index (De Cáceres &
#' Legendre 2009), the best one is kept, and a permutation test gives it a
#' p-value: De Cáceres, Legendre & Moretti's (2010) multilevel pattern
#' analysis. Gives the same results as `indicspecies::multipatt()`, but in
#' parallel C on sparse tables.
#'
#' @inheritParams documentation
#'
#' @param counts  A feature table, such as taxa abundances: a numeric matrix
#'        (samples in rows, or in columns with `margin = 2`), a sparse
#'        matrix from \pkg{Matrix} or \pkg{slam}, or a `phyloseq`, `rbiom`,
#'        or `SummarizedExperiment` object. Sparse tables are never made
#'        dense. Values must be non-negative. They are used as given: see
#'        **Normalization**.
#'
#' @param groups  Each sample's group: a vector or factor with one value per
#'        sample, matched to the samples in `counts` by name when both are
#'        named, else in the same order. Or, when `data` is given, the name
#'        of one of its columns (or a vector with one value per row of
#'        `data`). There must be at least two groups. A factor's levels set
#'        the groups' order (unused levels are dropped); other values are
#'        sorted.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the columns named in `groups` or `strata`. Rows are matched
#'        to the samples in `counts` by row name, and rows for other samples
#'        are ignored; without row names, they must be in the same order as
#'        the samples. Default: `NULL`.
#'
#' @param index  The association index: `'indval_g'` (the default), the
#'        group-equalized indicator value; `'indval'`, the indicator value
#'        without equalization; `'r'`, the point-biserial correlation (the
#'        phi coefficient for presence-absence data); or `'r_g'`, its
#'        group-equalized form. `indicspecies`' names (`'IndVal.g'`,
#'        `'IndVal'`, `'r'`, `'r.g'`) are accepted too. See **Indices**.
#'
#' @param min_order,max_order  The smallest and largest number of groups in
#'        a combination: `max_order = 1` considers single groups only, as
#'        Dufrêne & Legendre (1997) did. Defaults: `1` and `NULL`, all the
#'        groups (for `'r'` and `'r_g'`, all but one: see **Combinations of
#'        groups**).
#'
#' @param combinations  The combinations of groups to consider, when not
#'        all of them make sense: a character vector of names, such as
#'        `c('a', 'b', 'a+b')` (groups joined with `"+"`), or the
#'        combinations' `combination_id`s (see **Value**). Those outside
#'        `min_order` to `max_order` are dropped. Default: `NULL`, all of
#'        them.
#'
#' @param negative  For `'r'` and `'r_g'` only: whether a negative
#'        association (a feature avoiding the combination's samples) can be
#'        the best one, by the absolute value of the index. Default:
#'        `FALSE`.
#'
#' @param alpha  The significance level of the adjusted p-values: `reject`
#'        is `p_adjusted <= alpha`. Strictly between 0 and 1. Default:
#'        `0.05`.
#'
#' @param p_adjust  How to adjust the p-values for testing every feature:
#'        any method of [stats::p.adjust()] (`'holm'`, the default, as
#'        ecodive's DA methods; `'BH'`, `'hochberg'`, `'hommel'`,
#'        `'bonferroni'`, `'BY'`, `'fdr'`, or `'none'`, as
#'        `indicspecies::multipatt()`). Features without a p-value are left
#'        out of the family.
#'
#' @param permutations  The number of permutations. `0` skips the test, and
#'        gives `NA` p-values. Default: `999`. Alternatively, a matrix with
#'        one permutation of `1:n` per row, as from
#'        `permute::shuffleSet()`, to use instead of random ones: sample `i`
#'        then takes the group of sample `permutations[k, i]`, as in
#'        `indicspecies::multipatt(permutations = )`.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`.
#'
#' @return A data frame with one row per feature, in the order of `counts`,
#'   of class `c('ecodive_indicator_species', 'data.frame')`, with columns:
#'
#'   * `feature`: The feature's name (its column number, as text, when
#'     `counts` has no feature names).
#'   * `combination`: The group, or combination of groups, the feature is
#'     most associated with, as `"a"` or `"a+b"` (groups in their order).
#'     `NA` when no index is defined: for a feature absent from every
#'     sample, or (correlations) with the same value in every sample.
#'   * `combination_id`: The combination's number among all
#'     \eqn{2^K - 1} combinations of the \eqn{K} groups, numbered as
#'     `indicspecies` numbers them: the single groups first, then the pairs,
#'     and so on, each size in the groups' order (with three groups: 1 =
#'     a, 2 = b, 3 = c, 4 = a+b, 5 = a+c, 6 = b+c, 7 = a+b+c). These are
#'     the ids `combinations` takes.
#'   * `n_groups`: How many groups the combination has.
#'   * `stat`: The index for that combination (see **Indices**).
#'   * `a`, `b`: For the indicator values, its two parts: the specificity
#'     \eqn{A} and the fidelity \eqn{B}, with `stat` \eqn{= \sqrt{A B}}.
#'     `NA` for the correlation indices.
#'   * `p_value`: The permutation p-value (see **Permutation test**). `NA`
#'     when the combination is all the groups (there are no other samples to
#'     compare with), with `permutations = 0`, and where `combination` is
#'     `NA`.
#'   * `p_adjusted`: `p_value` adjusted with `p_adjust`.
#'   * `reject`: `p_adjusted <= alpha` (`FALSE` where `NA`).
#'
#'   Attributes: `method` (`"indicator_species"`), `index`, `groups` (the
#'   groups, in order), `combinations` (the combinations considered: their
#'   `combination_id`s, named by combination), `alpha`, `p_adjust`,
#'   `n_perm` (the number of permutations), `permutation` (`"free"`,
#'   `"within strata"` or `"supplied matrix"`), and `means` (each group's
#'   mean of each feature: a groups \eqn{\times} features matrix, for
#'   `plot()`).
#'
#' @details
#'
#' **Indices**
#'
#' For a feature with value \eqn{a_i} in sample \eqn{i}, and a combination
#' \eqn{C} of groups with \eqn{n_C} samples in all:
#'
#' * **Specificity** \eqn{A}, the probability that a sample belongs to
#'   \eqn{C} given that it has the feature (De Cáceres' "positive predictive
#'   value"). `'indval'` estimates it from the abundances,
#'   \eqn{A = \sum_{i \in C} a_i / \sum_i a_i}; `'indval_g'` from the
#'   groups' mean abundances \eqn{\bar a_k}, so that every group weighs
#'   the same whatever its size, \eqn{A_g = \sum_{k \in C} \bar a_k /
#'   \sum_k \bar a_k}.
#' * **Fidelity** \eqn{B}, the probability of finding the feature in a
#'   sample of \eqn{C} (its "sensitivity"): the proportion of \eqn{C}'s
#'   samples that have it, \eqn{B = \#\{i \in C : a_i > 0\} / n_C}.
#' * The **indicator value** is \eqn{\sqrt{A B}}, between 0 and 1: 1 for a
#'   feature found in every sample of \eqn{C} and nowhere else. (Dufrêne &
#'   Legendre's IndVal is \eqn{A B}, which their equalized \eqn{A} makes
#'   `'indval_g'` here; `indicspecies` and ecodive report the square root.)
#' * `'r'` is the Pearson correlation between the feature's values and
#'   membership of \eqn{C} (1 in \eqn{C}, 0 elsewhere): the point-biserial
#'   correlation, or the phi coefficient for presence-absence data.
#'   `'r_g'` is the same correlation with each sample weighted by one over
#'   its group's size, so that every group weighs the same (Tichý & Chytrý
#'   2006). They range from -1 to 1.
#'
#' The group-equalized indices (`'indval_g'`, `'r_g'`) suit groups of
#' unequal sizes, and are the ones the `indicspecies` vignette recommends.
#' The indicator values measure how well a feature predicts the
#' combination (for identifying it in the field, or for monitoring); the
#' correlations measure a feature's preference for it over the other
#' samples, and can be negative (De Cáceres & Legendre 2009).
#'
#'
#' **Combinations of groups**
#'
#' A feature found in two of three habitats is a poor indicator of either
#' one, but a good one of the pair. So every combination of groups is
#' considered (De Cáceres et al. 2010), within `min_order` to `max_order`
#' groups and `combinations`, and the one with the largest index is the
#' feature's `combination`. With `negative = TRUE`, the largest absolute
#' value: a correlation is negative for a combination when it is positive
#' for the other groups (the complement), by the same amount, so the
#' negative association with the complement ties the positive one, and the
#' tie goes to the combination of fewer groups. Values within a relative
#' \eqn{10^{-12}} of the largest are ties (for rounding), and ties go to
#' the combination of fewest groups, then to the first in the groups'
#' order. The combination of all the groups can be a feature's best for an
#' indicator value (a feature as common everywhere as anywhere), and is
#' then not tested; the correlations are undefined for it (there are no
#' other samples), so it is never considered for them. With many groups,
#' restrict the combinations (`max_order`): \eqn{K} groups have \eqn{2^K -
#' 1}, and more than \eqn{2^{20}} are refused.
#'
#'
#' **Permutation test**
#'
#' The samples' groups are shuffled (only within `strata`, when given), and
#' the index computed again for every combination. As the best combination
#' was chosen as the largest of many, the test compares it with the largest
#' value over all the combinations considered in each permutation (the
#' largest absolute value, with `negative = TRUE`), as `multipatt()`
#' does. The p-value is \eqn{(1 + \#\{perm \ge obs\}) / (1 +
#' permutations)}, comparing the values exactly, as `multipatt()` does. A
#' permutation that leaves the best combination's samples in it gives
#' exactly the observed value, and counts.
#'
#' Each feature is tested on its own, so the p-values need adjusting for
#' the number of features tested (`p_adjust`): `multipatt()` reports them
#' unadjusted, and its vignette warns about statements such as "the number
#' of indicator species is X" (De Cáceres & Legendre 2009; Legendre &
#' Legendre 2012).
#'
#'
#' **Normalization**
#'
#' The specificity uses abundances, and the fidelity only presence. So the
#' results depend on how `counts` was normalized: give relative abundances
#' (`normalize(counts, 'percent')`) unless the samples' totals are
#' comparable, and remember that sequencing depth (and rarefaction) change
#' which samples have a rare feature, and so its fidelity.
#'
#'
#' **Speed and memory**
#'
#' Each feature costs one pass over the samples that have it, plus a few
#' operations per combination, for the observed data and for each
#' permutation; nothing of the size of features \eqn{\times} combinations
#' is kept. Permutations run in parallel, in blocks of a fixed size, and
#' only integer counts are combined, so results are identical for any
#' `cpus`.
#'
#'
#' **Differences from indicspecies::multipatt()**
#'
#' * One row per feature, with the combination named (`"a+b"`) and its
#'   `combination_id`, rather than `multipatt()`'s membership columns
#'   (`s.a`, `s.b`, ...) and `index`. `multipatt()`'s `index` is the
#'   combination's position among those considered, so it differs from
#'   `combination_id` when `min.order` or `restcomb` leaves out the first
#'   ones.
#' * `combinations` takes names as well as ids; `restcomb` takes ids.
#'   `duleg = TRUE` is `max_order = 1`.
#' * `p_adjusted` and `reject` are additions; `p_adjust = 'none'` gives
#'   `multipatt()`'s p-values.
#' * Ties (and values within a relative \eqn{10^{-5}} of the largest) are
#'   broken at random by `multipatt()`, so the same data can give a
#'   different combination from one run to the next; ecodive breaks them
#'   by the rule above, the same way every time.
#' * When the best combination is all the groups, the p-value is `NA` here,
#'   as in `multipatt()` without `restcomb`; with `restcomb`,
#'   `multipatt()` tests it and gives 1.
#' * Random permutations come from ecodive's own generator (PCG, with
#'   `seed`), so p-values differ from `multipatt()`'s by sampling error; a
#'   supplied `permutations` matrix gives the same p-values. `strata` is
#'   `permute::how(blocks = )`.
#' * `strassoc()`'s other indices, its bootstrap confidence intervals, and
#'   `indicators()` (combinations of species as indicators) are not
#'   provided.
#'
#'
#' @section Same names in other packages:
#'   None: `labdsv::indval()` computes `index = 'indval_g'` with
#'   `max_order = 1`, and `indicspecies::multipatt()` the whole analysis.
#'   [da_test()] runs this function as `"indicator_species"`, and also
#'   as `"indval"` or `"multipatt"`.
#'
#' @references
#' De Cáceres, M., & Legendre, P. (2009). Associations between species and
#' groups of sites: indices and statistical inference. *Ecology*, 90,
#' 3566-3574. \doi{10.1890/08-1823.1}
#'
#' De Cáceres, M., Legendre, P., & Moretti, M. (2010). Improving indicator
#' species analysis by combining groups of sites. *Oikos*, 119, 1674-1684.
#' \doi{10.1111/j.1600-0706.2010.18334.x}
#'
#' Dufrêne, M., & Legendre, P. (1997). Species assemblages and indicator
#' species: the need for a flexible asymmetrical approach. *Ecological
#' Monographs*, 67, 345-366.
#' \doi{10.1890/0012-9615(1997)067[0345:SAAIST]2.0.CO;2}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English
#' ed.). Elsevier.
#'
#' Tichý, L., & Chytrý, M. (2006). Statistical determination of diagnostic
#' species for site groups of unequal size. *Journal of Vegetation
#' Science*, 17, 809-818. \doi{10.1111/j.1654-1103.2006.tb02504.x}
#'
#' @seealso [da_test()], [simper()], [lefse()],
#'   [plot.ecodive_indicator_species()]
#'
#' @export
#' @examples
#'     # Three habitats; f1 lives in 'a', f2 in 'b' and 'c', f3 everywhere.
#'     set.seed(1)
#'     grp <- rep(c('a', 'b', 'c'), c(6, 8, 10))
#'     x <- matrix(rpois(24 * 5, 2), 24, 5, dimnames = list(NULL, paste0('f', 1:5)))
#'     x[grp != 'a', 1] <- 0
#'     x[grp == 'a', 2] <- 0
#'     x[, 3] <- x[, 3] + 1
#'     res <- indicator_species(x, grp, cpus = 1)
#'     res
#'
#'     # The original analysis of Dufrene & Legendre: single groups only.
#'     indicator_species(x, grp, max_order = 1, cpus = 1)
#'
#'     # Correlations, allowing negative associations.
#'     indicator_species(x, grp, index = 'r_g', negative = TRUE, cpus = 1)
#'
#'     plot(res)
indicator_species <- function (
    counts, groups, data = NULL, margin = 1L,
    index = c('indval_g', 'indval', 'r', 'r_g'),
    min_order = 1L, max_order = NULL, combinations = NULL,
    negative = FALSE, alpha = 0.05, p_adjust = 'holm',
    permutations = 999L, strata = NULL, seed = 0L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  indsp_run(environment())$result
}



# indicator_species(), and every combination's value (for tests): returns
# list(result, strength (features x combinations considered, named)).
indicator_species_details <- function (
    counts, groups, data = NULL, margin = 1L,
    index = c('indval_g', 'indval', 'r', 'r_g'),
    min_order = 1L, max_order = NULL, combinations = NULL,
    negative = FALSE, alpha = 0.05, p_adjust = 'holm',
    permutations = 999L, strata = NULL, seed = 0L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  indsp_run(environment(), keep_str = TRUE)
}



#------------------------------------------------------------------------------
# The checks, the C call, and the result.
#------------------------------------------------------------------------------
indsp_run <- function (env, keep_str = FALSE) {

  validate_counts(env)
  validate_margin(env)
  validate_seed(env)
  validate_cpus(env)

  counts <- env$counts
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  margin <- env$margin

  index <- indsp_index(env$index)

  negative <- env$negative
  if (!isTRUE(negative) && !isFALSE(negative))
    stop('`negative` must be TRUE or FALSE.', call. = FALSE)
  if (negative && !index %in% c('r', 'r_g'))
    stop(
      "`negative = TRUE` is only for `index = 'r'` or `'r_g'`: the ",
      'indicator values are never negative.', call. = FALSE )

  alpha <- env$alpha
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1)
    stop('`alpha` must be a single number strictly between 0 and 1.', call. = FALSE)

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)


  # Dimensions and names, with samples first.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  n <- as.integer(dims[[1]])
  D <- as.integer(dims[[2]])
  if (D < 1L) stop('`counts` must have at least one feature.', call. = FALSE)
  features <- if (length(dn[[2]])) as.character(dn[[2]]) else as.character(seq_len(D))

  pm <- perm_matrix(env$permutations, n)


  #--------------------------------------------------------
  # Each sample's group (and stratum), in counts' order.
  #--------------------------------------------------------
  des <- test_group_values(env$groups, env$data, n, dn[[1]], env$strata, ref = '`counts`')
  val <- des$groups
  lev <- if (is.factor(val)) levels(val) else sort(unique(val), method = 'radix')
  gi  <- match(val, lev)
  lev <- as.character(lev)
  K   <- length(lev)

  if (K < 2L)
    stop('`groups` must have at least two different values.', call. = FALSE)
  if (K > 30L)
    stop('`groups` has ', K, ' groups; at most 30 are allowed.', call. = FALSE)

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  cmb <- indsp_combos(K, lev, index, env$min_order, env$max_order, env$combinations)


  #--------------------------------------------------------
  # The indices, best combinations, and permutation counts.
  #--------------------------------------------------------
  res <- .Call(
    C_indicator_species, counts, margin, as.integer(gi - 1L), K,
    as.integer(cmb$mask), match(index, INDSP_INDICES) - 1L, negative,
    pm$perms, des$strata, env$seed, pm$n_perm, keep_str, env$cpus )

  best <- res$best
  p    <- if (pm$n_perm > 0L) (res$count + 1) / (pm$n_perm + 1) else rep(NA_real_, D)
  p[is.na(res$count)] <- NA_real_
  padj <- stats::p.adjust(p, p_adjust)
  ind  <- index %in% c('indval', 'indval_g')

  out <- data.frame(
    feature        = features,
    combination    = cmb$name[best],
    combination_id = cmb$id[best],
    n_groups       = cmb$n_groups[best],
    stat           = res$stat,
    a              = if (ind) res$a else rep(NA_real_, D),
    b              = if (ind) res$b else rep(NA_real_, D),
    p_value        = p,
    p_adjusted     = padj,
    reject         = !is.na(padj) & padj <= alpha,
    stringsAsFactors = FALSE )
  out$stat[is.na(best)] <- NA_real_

  means <- t(res$means)
  dimnames(means) <- list(lev, features)

  attr(out, 'method')       <- 'indicator_species'
  attr(out, 'index')        <- index
  attr(out, 'groups')       <- lev
  attr(out, 'combinations') <- stats::setNames(cmb$id, cmb$name)
  attr(out, 'alpha')        <- as.double(alpha)
  attr(out, 'p_adjust')     <- p_adjust
  attr(out, 'n_perm')       <- pm$n_perm
  attr(out, 'permutation')  <- if (!is.null(pm$perms)) 'supplied matrix'
                               else if (!is.null(des$strata)) 'within strata'
                               else 'free'
  attr(out, 'means')        <- means
  class(out) <- c('ecodive_indicator_species', 'data.frame')

  strength <- NULL
  if (keep_str) {
    strength <- res$str
    dimnames(strength) <- list(features, cmb$name)
  }

  list(result = out, strength = strength)
}


INDSP_INDICES <- c('indval', 'indval_g', 'r', 'r_g')   # The C code's order.


# The association index, from `index`: ecodive's ids, or indicspecies'
# names ('IndVal.g', 'r.g'), ignoring case and anything but letters.
indsp_index <- function (index) {

  if (identical(index, c('indval_g', 'indval', 'r', 'r_g'))) return ('indval_g')

  msg <- "`index` must be one of 'indval_g', 'indval', 'r', or 'r_g'."
  if (!is.character(index) || length(index) != 1L || is.na(index))
    stop(msg, call. = FALSE)

  i <- match(gsub('[^a-z]', '', tolower(index)), c('indval', 'indvalg', 'r', 'rg'))
  if (is.na(i)) stop(msg, call. = FALSE)
  INDSP_INDICES[[i]]
}



#------------------------------------------------------------------------------
# The combinations of groups to consider, in indicspecies' order (by size,
# then lexicographically by the groups' positions). Returns list(mask (bit
# k - 1: group k), id (the position among all 2^K - 1), name ("a+b"),
# n_groups).
#------------------------------------------------------------------------------
indsp_combos <- function (K, lev, index, min_order, max_order, combinations) {

  whole <- function (v, arg) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v %% 1 != 0 || v < 1)
      stop('`', arg, '` must be a whole number of at least 1.', call. = FALSE)
    as.integer(v)
  }
  min_order <- whole(min_order, 'min_order')
  max_order <- if (is.null(max_order)) K else whole(max_order, 'max_order')

  if (max_order > K)
    stop('`max_order` must be at most the number of groups (', K, ').', call. = FALSE)
  if (min_order > max_order)
    stop('`min_order` must not be larger than `max_order` (', max_order, ').', call. = FALSE)

  # The correlations are undefined for the combination of all the groups.
  corr <- index %in% c('r', 'r_g')
  if (corr) max_order <- min(max_order, K - 1L)
  if (min_order > max_order)
    stop(
      "With `index = '", index, "'`, the combination of all the groups ",
      "can't be tested (there are no other samples): `min_order` must be ",
      'less than the number of groups (', K, ').', call. = FALSE )

  if (is.null(combinations)) {

    n_comb <- sum(choose(K, min_order:max_order))
    if (n_comb > 2^20)
      stop(
        K, ' groups have ', format(n_comb, big.mark = ',', scientific = FALSE),
        ' combinations of ', min_order, ' to ', max_order, ' groups, more ',
        'than 2^20: set `max_order` or `combinations`.', call. = FALSE )

    parts <- lapply(min_order:max_order, function (size) {
      m <- utils::combn(K, size)
      list(
        mask = as.integer(colSums(matrix(2^(m - 1), nrow = size))),
        id   = as.integer(indsp_offset(K, size) + seq_len(ncol(m))),
        name = apply(matrix(lev[m], nrow = size), 2L, paste, collapse = '+'),
        size = rep(size, ncol(m)) )
    })
    return (list(
      mask     = unlist(lapply(parts, `[[`, 'mask')),
      id       = unlist(lapply(parts, `[[`, 'id')),
      name     = unlist(lapply(parts, `[[`, 'name')),
      n_groups = unlist(lapply(parts, `[[`, 'size')) ))
  }


  # Named or numbered combinations.
  n_all <- 2^K - 1
  if (is.numeric(combinations)) {
    ok <- !is.na(combinations) & combinations %% 1 == 0 & combinations >= 1 & combinations <= n_all
    if (!length(combinations) || !all(ok))
      stop(
        '`combinations` must be names of groups or combinations, or ',
        'whole numbers from 1 to ', n_all, ' (the combination ids).', call. = FALSE )
    sets <- lapply(combinations, indsp_unrank, K = K)
  }
  else if (is.character(combinations) || is.factor(combinations)) {
    combinations <- as.character(combinations)
    if (!length(combinations) || anyNA(combinations))
      stop('`combinations` must not be empty or contain NAs.', call. = FALSE)
    sets <- lapply(combinations, function (s) {
      set <- indsp_parse(s, lev, seq_len(K))
      if (is.null(set))
        stop(
          '`combinations`: "', s, '" is not a group or a combination of ',
          'groups joined with "+". The groups are ',
          paste0('"', utils::head(lev, 10), '"', collapse = ', '),
          if (K > 10) ', ...', '.', call. = FALSE )
      sort(set)
    })
  }
  else {
    stop(
      '`combinations` must be NULL, names of groups or combinations ',
      '(such as "a+b"), or combination ids.', call. = FALSE )
  }

  size <- lengths(sets)
  keep <- size >= min_order & size <= max_order
  sets <- sets[keep]
  if (!length(sets))
    stop(
      'None of `combinations` has ', min_order, ' to ', max_order, ' groups',
      if (corr) " (the combination of all the groups can't be tested with correlations)",
      '.', call. = FALSE )

  id   <- vapply(sets, indsp_rank, double(1), K = K)
  o    <- order(id)
  o    <- o[!duplicated(id[o])]
  sets <- sets[o]

  list(
    mask     = vapply(sets, function (s) as.integer(sum(2^(s - 1))), integer(1)),
    id       = as.integer(id[o]),
    name     = vapply(sets, function (s) paste(lev[s], collapse = '+'), ''),
    n_groups = lengths(sets) )
}


# How many combinations have fewer than `size` of K groups.
indsp_offset <- function (K, size) sum(choose(K, seq_len(size - 1L)))


# The id of a combination (sorted group positions) among all 2^K - 1.
indsp_rank <- function (set, K) {
  size <- length(set)
  r    <- 0
  prev <- 0L
  for (j in seq_len(size)) {
    v <- prev + 1L
    while (v < set[[j]]) { r <- r + choose(K - v, size - j); v <- v + 1L }
    prev <- set[[j]]
  }
  indsp_offset(K, size) + r + 1
}


# The combination (sorted group positions) with id `id`.
indsp_unrank <- function (id, K) {
  size <- 1L
  while (id > indsp_offset(K, size + 1L)) size <- size + 1L
  r    <- id - indsp_offset(K, size) - 1
  set  <- integer(size)
  prev <- 0L
  for (j in seq_len(size)) {
    v <- prev + 1L
    while (choose(K - v, size - j) <= r) { r <- r - choose(K - v, size - j); v <- v + 1L }
    set[[j]] <- v
    prev <- v
  }
  set
}


# The groups named by `s`, group names joined with '+' (group names may
# themselves contain '+'), from those in `avail`; NULL if none.
indsp_parse <- function (s, lev, avail) {
  hit <- avail[lev[avail] == s]
  if (length(hit)) return (hit[[1]])
  for (k in avail) {
    pre <- paste0(lev[[k]], '+')
    if (startsWith(s, pre) && nchar(s) > nchar(pre)) {
      rest <- indsp_parse(substring(s, nchar(pre) + 1L), lev, setdiff(avail, k))
      if (!is.null(rest)) return (c(k, rest))
    }
  }
  NULL
}


# "Indicator species (IndVal.g)": the title of prints and plots.
indsp_label <- function (x) {
  index <- attr(x, 'index', exact = TRUE)
  nm <- c(indval = 'IndVal', indval_g = 'IndVal.g', r = 'r', r_g = 'r.g')
  if (is.character(index) && length(index) == 1L && index %in% names(nm))
    paste0('Indicator species (', nm[[index]], ')')
  else 'Indicator species'
}



#' Plot an Indicator Species Analysis
#'
#' The features found by [indicator_species()] (`reject`): their index as
#' bars, or their groups' mean abundances as a heatmap.
#'
#' @param x  The result of [indicator_species()].
#'
#' @param type  `'bar'` (the default): one horizontal bar per feature, its
#'        length the index (`stat`), grouped and coloured by `combination`,
#'        the largest first within each. `'heatmap'`: the features (rows)
#'        by the groups (columns), each tile shaded by the group's mean
#'        abundance as a proportion of the feature's largest group mean,
#'        with the groups of the feature's `combination` outlined.
#'
#' @param n_features  How many features to show, those with the largest
#'        (absolute) index. Default: `30`.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `xlab` or `sub`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     grp <- rep(c('a', 'b', 'c'), c(6, 8, 10))
#'     x <- matrix(rpois(24 * 6, 2), 24, 6, dimnames = list(NULL, paste0('f', 1:6)))
#'     x[grp != 'a', 1:2] <- 0
#'     x[grp == 'a', 3] <- 0
#'     res <- indicator_species(x, grp, p_adjust = 'none', cpus = 1)
#'     plot(res)
#'     plot(res, 'heatmap')
plot.ecodive_indicator_species <- function (x, type = c('bar', 'heatmap'), n_features = 30, ...) {

  type  <- match.arg(type)
  n     <- plot_count(n_features, 'n_features')
  grp   <- attr(x, 'groups', exact = TRUE)
  alpha <- attr(x, 'alpha', exact = TRUE)
  means <- attr(x, 'means', exact = TRUE)

  rows  <- which(x$reject %in% TRUE & is.finite(x$stat))
  n_rep <- length(rows)
  rows  <- utils::head(rows[order(-abs(x$stat[rows]))], n)

  # Grouped by combination (in id order), largest first within each, from
  # the top down.
  rows <- rows[order(x$combination_id[rows], -abs(x$stat[rows]))]
  kk   <- length(rows)
  cmbs <- unique(x$combination[rows])
  col  <- plot_palette(length(cmbs))[match(x$combination[rows], cmbs)]

  sub <- sprintf('%d of %d features at alpha = %s', n_rep, nrow(x), format(alpha))
  if (kk < n_rep) sub <- sprintf('%s; top %d shown', sub, kk)
  main <- indsp_label(x)

  opar <- graphics::par('mai')
  on.exit(graphics::par(mai = opar))
  linch <- max(graphics::strwidth(x$feature[rows], 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)

  y <- rev(seq_len(kk))   # The first row at the top.

  if (type == 'bar') {

    val  <- x$stat[rows]
    xlim <- c(min(0, val, na.rm = TRUE), max(0, val, 1e-9, na.rm = TRUE))

    plot_do(
      graphics::plot.default, c(list(0, 0), list(...)),
      list(
        type = 'n', yaxt = 'n', xlim = xlim, ylim = c(0.5, max(kk, 1) + 0.5),
        main = main, xlab = base::sub('^Indicator species \\((.*)\\)$', '\\1', main),
        ylab = '', sub = sub ))

    if (kk) {
      graphics::rect(0, y - 0.4, val, y + 0.4, col = col, border = NA)
      graphics::axis(2, at = y, labels = x$feature[rows], tick = FALSE, las = 1, cex.axis = 0.8)
      graphics::legend(
        'bottomright', legend = cmbs, fill = plot_palette(length(cmbs)),
        border = NA, bty = 'n', cex = 0.8 )
    }
    else {
      graphics::text(mean(xlim), 1, 'No features found', col = PLOT_MUTED)
    }
    graphics::abline(v = 0, col = 'grey70')
    return (invisible(x))
  }


  # Heatmap.
  k  <- length(grp)
  pal <- grDevices::hcl.colors(65, 'Blues 3', rev = TRUE)

  plot_do(
    graphics::plot.default, c(list(0, 0), list(...)),
    list(
      type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0.5, k + 0.5),
      ylim = c(0.5, max(kk, 1) + 0.5), main = main, xlab = '', ylab = '',
      sub = sub, bty = 'n' ))

  if (kk) {
    m   <- means[, x$feature[rows], drop = FALSE]
    top <- apply(m, 2L, max)
    top[!(top > 0)] <- 1
    m   <- sweep(m, 2L, top, '/')
    for (j in seq_len(kk)) {
      fill <- pal[1L + round(m[, j] * (length(pal) - 1L))]
      graphics::rect(seq_len(k) - 0.5, y[[j]] - 0.5, seq_len(k) + 0.5, y[[j]] + 0.5,
                     col = fill, border = 'white')
      inC <- match(strsplit_combination(x$combination[rows[[j]]], grp), seq_len(k))
      inC <- inC[!is.na(inC)]
      if (length(inC))
        graphics::rect(inC - 0.45, y[[j]] - 0.45, inC + 0.45, y[[j]] + 0.45,
                       border = PLOT_HIGHLIGHT, lwd = 2)
    }
    graphics::axis(2, at = y, labels = x$feature[rows], tick = FALSE, las = 1, cex.axis = 0.8)
    graphics::axis(1, at = seq_len(k), labels = grp, tick = FALSE, cex.axis = 0.8)
  }
  else {
    graphics::text((k + 1) / 2, 1, 'No features found', col = PLOT_MUTED)
  }

  invisible(x)
}


# The groups (positions in `grp`) of a combination's name.
strsplit_combination <- function (name, grp) {
  if (is.na(name)) return (integer(0))
  set <- indsp_parse(name, grp, seq_along(grp))
  if (is.null(set)) integer(0) else set
}
