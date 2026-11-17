# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Analysis of Similarities (ANOSIM)
#'
#' Tests whether samples in the same group are more similar to each other
#' than to samples in other groups, from the ranks of the distances between
#' them. Gives the same \eqn{R} statistic and p-value as `vegan::anosim()`.
#' Computes them in parallel C, from one ranking of the distances, without
#' the \eqn{n \times n} matrices that `anosim()` builds.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values.
#'
#' @param groups  Each sample's group: a vector or factor with one value per
#'        sample, matched to the samples in `dm` by name when both are named,
#'        else in the same order as `dm`. Or, when `data` is given, the name of
#'        one of its columns (or a vector with one value per row of `data`).
#'        Values are treated as categories, even if they are numbers.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `groups` (and `strata`) as columns. Rows are matched
#'        to the samples in `dm` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as `dm`.
#'        Default: `NULL`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param strata  Permute samples only within these groups (blocks), for
#'        designs such as repeated measures. The name of a column of `data`,
#'        or a vector with one value per row of `data` (or, without `data`,
#'        per value of `groups`). Default: `NULL` (free permutation).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-value. Default: `0`
#'
#' @return A list of class `c("ecodive_anosim", "ecodive_test")`, the same
#'   shape as the results of ecodive's other permutation tests, with:
#'
#'   * `statistic`: The ANOSIM statistic, `c(R = ...)`: `vegan::anosim()`'s
#'     `statistic`.
#'   * `p_value`: Permutation p-value, `anosim()`'s `signif`. `NA` without
#'     permutations.
#'   * `z`: How many standard deviations \eqn{R} lies above the mean of the
#'     permuted statistics: \eqn{(R - mean(R_{perm})) / sd(R_{perm})}, as
#'     `summary(vegan::permustats())`. `NA` with fewer than two permutations,
#'     or when they are all equal.
#'   * `n_perm`: The number of permutations.
#'   * `ranks`: A data frame with one row for the pairs of samples in
#'     different groups (`"Between"`) and one for each group, with the number
#'     of `pairs` and their `mean_rank`: the averages behind
#'     `summary(anosim(...))`'s quantiles. A group of one sample has no pairs
#'     and an `NA` mean rank.
#'   * `permutation`: How samples were permuted: `"free"`, `"within strata"`,
#'     or `"supplied matrix"`.
#'
#'   Attribute `r_perm` is a vector of the \eqn{R} statistics of each
#'   permutation, `anosim()`'s `perm`.
#'
#' @details
#'
#' **Method**
#'
#' ANOSIM (Clarke 1993) ranks all \eqn{M = n(n-1)/2} distances, from 1 for
#' the smallest to \eqn{M} for the largest, giving tied distances the
#' average of their ranks. With \eqn{\bar{r}_B} the mean rank of the
#' distances between samples in different groups, and \eqn{\bar{r}_W} that
#' of the distances within groups,
#' \deqn{R = \frac{\bar{r}_B - \bar{r}_W}{M / 2}}
#' \eqn{R} is near 0 when the grouping explains nothing, and 1 when every
#' distance within groups is smaller than every distance between them.
#' Negative values mean samples are more alike between groups than within
#' them.
#'
#' For each permutation, the samples' group labels are shuffled and
#' \eqn{R} is recomputed. The p-value is the proportion of permutations,
#' counting the observed one, with an \eqn{R} at least as large:
#' \eqn{(1 + \#\{R_{perm} \ge R\}) / (1 + permutations)}. As in `anosim()`,
#' permuted statistics within `sqrt(.Machine$double.eps)` of the observed one
#' count as ties.
#'
#' ANOSIM tests for differences in location, but is also sensitive to
#' differences in spread (dispersion) between groups (Warton et al. 2012);
#' see also [permanova()].
#'
#'
#' **Speed and memory**
#'
#' The distances are ranked once, with a radix sort. Shuffling the labels
#' doesn't change how many pairs are within groups, nor the sum of all the
#' ranks, so each permutation needs only the sum of the ranks within groups:
#' one pass over the ranks, shared by many permutations at once, with no
#' \eqn{n \times n} matrix. Permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. P-values are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `anosim()`'s, so
#' p-values agree with it only up to sampling error, unless the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from vegan::anosim()**
#'
#' * `anosim()` here takes distances, not a community table: compute them
#'   first with [beta_div()] or a metric function such as [bray()].
#' * The fields are named as in ecodive's other tests: `statistic` is named
#'   `R`, `signif` is `p_value`, and the permuted statistics are attribute
#'   `r_perm`. The full vectors of ranks and their classes (`dis.rank` and
#'   `class.vec`), each as long as `dm`, are summarized as the mean rank of
#'   each class, in field `ranks`.
#' * `strata` permutes freely within strata. For other restricted designs,
#'   pass a `permutations` matrix from `permute::shuffleSet()`.
#' * When there are fewer possible permutations than requested, `vegan`
#'   (via `permute`) runs each of them once instead. `anosim()` always draws
#'   random permutations.
#' * Negative distances are an error, where `vegan` warns.
#'
#' @references
#' Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in
#' community structure. *Australian Journal of Ecology*, 18(1), 117-143.
#' \doi{10.1111/j.1442-9993.1993.tb00438.x}
#'
#' Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based
#' multivariate analyses confound location and dispersion effects.
#' *Methods in Ecology and Evolution*, 3(1), 89-101.
#' \doi{10.1111/j.2041-210X.2011.00127.x}
#'
#' @seealso [permanova()], [beta_div()]
#'
#' @export
#' @examples
#'     # 30 samples in three treatment groups. The 'high' group has more
#'     # of the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       treatment = rep(c('ctrl', 'low', 'high'), each = 10),
#'       site      = rep(c('A', 'B'), 15) )
#'     high <- meta$treatment == 'high'
#'     counts[high, 1:5] <- counts[high, 1:5] + 5
#'
#'     dm <- bray(counts)
#'
#'     # Are samples more alike within treatment groups than between them?
#'     anosim(dm, 'treatment', meta)
#'
#'     # The same, with the groups as a vector.
#'     anosim(dm, meta$treatment)
#'
#'     # Permuting treatment labels only within each site.
#'     anosim(dm, 'treatment', meta, strata = 'site')
#'
#'     # Mean rank of the distances between and within groups.
#'     anosim(dm, 'treatment', meta)$ranks
#'
anosim <- function (
    dm, groups, data = NULL, permutations = 999L, strata = NULL,
    seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  dv <- dist_vector(dm)
  n  <- dv$n

  if (n < 2)
    stop('`dm` must have at least 2 samples.', call. = FALSE)

  if (missing(groups)) groups <- NULL

  des <- test_group_values(groups, data, n, dv$labels, strata)
  grp <- list(groups = droplevels(factor(des$groups)), strata = des$strata)
  pm  <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(grp$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  lev   <- levels(grp$groups)
  codes <- as.integer(grp$groups) - 1L
  sizes <- tabulate(codes + 1L, length(lev))
  n_w   <- sum(as.double(sizes) * (sizes - 1) / 2)  # pairs within groups

  if (length(lev) < 2)
    stop('`groups` must have at least 2 different values.', call. = FALSE)
  if (n_w == 0)
    stop('Every group has only one sample: ANOSIM needs groups with replicates.', call. = FALSE)


  #--------------------------------------------------------
  # R for the identity and each permutation, and the
  # observed sums of (doubled) ranks within each group.
  #--------------------------------------------------------
  res <- .Call(
    C_anosim, dv$dm, as.integer(n), codes, pm$perms, grp$strata,
    as.integer(c(pm$n_perm, seed, cpus)) )

  M     <- as.double(n) * (n - 1) / 2
  n_b   <- M - n_w
  r_obs <- res$r[[1]]
  r_prm <- res$r[-1]


  #--------------------------------------------------------
  # p-value and z-score.
  #--------------------------------------------------------
  eps <- sqrt(.Machine$double.eps)
  p   <- if (pm$n_perm > 0) (1 + sum(r_prm >= r_obs - eps)) / (1 + pm$n_perm) else NA_real_

  z <- perm_z(r_obs, r_prm)


  #--------------------------------------------------------
  # The mean rank of each class, and the result.
  #--------------------------------------------------------
  gp <- res$group_pairs
  ranks <- data.frame(
    pairs     = c(n_b, gp),
    mean_rank = c(
      (M * (M + 1) - res$sw[[1]]) / 2 / n_b,
      ifelse(gp > 0, res$group_sw / 2 / gp, NA_real_) ),
    row.names = c('Between', lev) )

  permutation <-
    if (!is.null(pm$perms))        'supplied matrix'
    else if (!is.null(grp$strata)) 'within strata'
    else                           'free'

  new_ecodive_test(
    class       = 'ecodive_anosim',
    title       = 'Analysis of similarities (ANOSIM)',
    statistic   = c(R = r_obs),
    p_value     = p,
    z           = z,
    n_perm      = pm$n_perm,
    ranks       = ranks,
    perm        = r_prm,
    permutation = permutation )
}


# Lines added to the printout by print.ecodive_test(): the mean ranks.
#' @export
test_details.ecodive_anosim <- function (x, digits) {
  rk <- x$ranks
  rk$mean_rank <- signif(rk$mean_rank, digits + 2L)
  c('Mean rank of the distances, between and within groups:',
    utils::capture.output(print(rk)) )
}
