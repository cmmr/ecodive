# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The SIMPER breakdown follows ecopy's simper() (ecopy 0.1.2.2, MIT,
# (c) 2017 Nathan Lemoine); see src/simper.c for its license text.



#' Similarity Percentages (SIMPER)
#'
#' Breaks down the mean Bray-Curtis dissimilarity between each pair of
#' groups into the contribution of each taxon, to rank the taxa that
#' differ most between the groups (Clarke 1993). Permutation tests give
#' each taxon's contribution a p-value. Gives the same results as
#' `vegan::simper()` and `summary()` of it, but runs in parallel C on
#' sparse tables. Its p-values count ties with the observed contribution,
#' as Clarke (1993) does, so they are never smaller than vegan's, which
#' are too small for small groups (see **Permutations**).
#'
#' @inheritParams documentation
#'
#' @param counts  A feature table, such as taxa abundances. Accepts the same
#'        objects as [beta_div()]: a numeric matrix (samples in rows, or in
#'        columns with `margin = 2`), a sparse matrix from \pkg{Matrix} or
#'        \pkg{slam}, or a `phyloseq`, `rbiom`, or `SummarizedExperiment`
#'        object. Sparse tables are never made dense. Values must be
#'        non-negative, and every sample must have a non-zero total.
#'
#' @param groups  Each sample's group: a vector or factor with one value per
#'        sample, matched to the samples in `counts` by name when both are
#'        named, else in the same order. Or, when `data` is given, the name
#'        of one of its columns (or a vector with one value per row of
#'        `data`). Every distinct value is a group, numbers included; there
#'        must be at least two. A factor's unused levels are dropped.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the columns named in `groups` or `strata`. Rows are matched
#'        to the samples in `counts` by row name, and rows for other samples
#'        are ignored; without row names, they must be in the same order as
#'        the samples. Default: `NULL`.
#'
#' @param norm  Normalize `counts` first: `'none'` (as `vegan::simper()`,
#'        which uses the values as given), `'percent'`, `'chord'`,
#'        `'hellinger'`, or `'binary'`, as in [beta_div()]. `'clr'` and
#'        `'rclr'` are not allowed: they give negative values, which
#'        Bray-Curtis can't use. Default: `'none'`.
#'
#' @param p_adjust  How to adjust each pair's p-values for testing many
#'        taxa: one of `stats::p.adjust.methods` (`'holm'`, `'hochberg'`,
#'        `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`, `'fdr'`, or `'none'`),
#'        matched as `stats::p.adjust()` matches them. The taxa of each pair
#'        are adjusted together, leaving out those with a `p_value` of `NA`.
#'        Default: `'holm'`.
#'
#' @param permutations  The number of permutations. `0` skips the tests,
#'        and gives `NA` p-values. Default: `999`. Alternatively, a matrix
#'        with one permutation of `1:n` per row, as from
#'        `permute::shuffleSet()`, to use instead of random ones: it is cut
#'        down to each pair's samples (see **Permutations**).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A data frame with one row per pair of groups and taxon, the pairs
#'   in the order of the groups' factor levels (sorted values, for other
#'   vectors), and the taxa of each pair in decreasing order of `average`.
#'   Its columns are:
#'
#'   * `pair`: The pair of groups, as `"a-b"`.
#'   * `taxon`: The taxon's name (its column number, as text, when `counts`
#'     has no taxa names).
#'   * `average`: The taxon's mean contribution to the Bray-Curtis
#'     dissimilarity between the two groups' samples. A pair's `average`s
#'     sum to its mean between-group Bray-Curtis dissimilarity.
#'   * `sd`: The standard deviation of the taxon's contribution over the
#'     pairs of samples (`NA` when the groups have one sample each).
#'   * `ratio`: `average / sd`. Taxa with a high ratio contribute
#'     consistently, not just a lot in a few pairs of samples.
#'   * `mean_a`, `mean_b`: The taxon's mean abundance (after `norm`) in the
#'     first and second group of the pair.
#'   * `cumsum`: The cumulative proportion of the pair's mean dissimilarity
#'     that this taxon and those above it account for.
#'   * `p_value`: The permutation p-value: the proportion of permutations
#'     giving the taxon a larger `average` (see **Permutations**). `NA` with
#'     `permutations = 0`, and for taxa whose `average` is 0 (those absent
#'     from both groups, or with the same abundance in every sample).
#'   * `p_adjusted`: `p_value` adjusted for multiple testing with
#'     `p_adjust`, over the taxa of the pair.
#'
#'   The data frame has class `"ecodive_simper"`, only to mark what made it.
#'   Attributes: `overall`, each pair's mean between-group Bray-Curtis
#'   dissimilarity (named by pair); `n_perm`, the number of permutations;
#'   `p_adjust`, the adjustment method; and `permutation`, `"free"`,
#'   `"within strata"` or `"supplied matrix"`.
#'
#' @details
#'
#' **Method**
#'
#' For samples \eqn{i} and \eqn{j} with totals \eqn{T_i} and \eqn{T_j}, the
#' Bray-Curtis dissimilarity \eqn{\sum_k |x_{ik} - x_{jk}| / (T_i + T_j)}
#' is a sum over the taxa, so taxon \eqn{k} contributes \eqn{c_{ijk} =
#' |x_{ik} - x_{jk}| / (T_i + T_j)}. For two groups \eqn{a} and \eqn{b}, a
#' taxon's `average` is the mean of \eqn{c_{ijk}} over the \eqn{n_a n_b}
#' pairs of a sample from \eqn{a} and one from \eqn{b}, and `sd` their
#' standard deviation (with denominator \eqn{n_a n_b - 1}). SIMPER is
#' defined for Bray-Curtis only.
#'
#'
#' **Permutations**
#'
#' A taxon's p-value (Clarke 1993) comes from shuffling the group labels of
#' the pair's samples (only within `strata`, when given), and recomputing
#' its `average`. The p-value is \eqn{(1 + \#\{average_{perm} \ge
#' average\}) / (1 + permutations)}, Clarke's (1993) "at least as large"
#' rule, with permuted averages within a relative \eqn{10^{-7}} of the
#' observed one counted as ties (for rounding). Permutations that only swap
#' samples within groups reproduce the observed `average` exactly, so they
#' count. `vegan::simper()` counts only permuted averages that are strictly
#' larger, which makes its p-values too small when groups are small: with
#' two samples per group, a third of all permutations reproduce the
#' observed `average`, and a taxon whose observed `average` is the largest
#' possible gets \eqn{p \approx 1 / (1 + permutations)} from vegan, where
#' its p-value is about 1/3. So ecodive's p-values are never smaller than
#' vegan's on the same permutations, and are equal for taxa that no
#' permutation ties. Each pair has random
#' permutations of its own, drawn from `seed` and the pair's position, as
#' the pairwise tests of [permanova()] and the others have. So `p_value`s
#' are identical for any number of `cpus`, but differ by sampling error
#' from `simper()` run by hand on the pair's samples with the same `seed`.
#' A supplied `permutations` matrix is cut down to each pair's samples by
#' following each permutation's cycles, as in the pairwise tests (see
#' [permanova()]'s **Pairwise comparisons**). Give `vegan::simper()` those
#' permutations, on the pair's samples, for the same permutations.
#'
#' The adjusted p-values (`p_adjusted`) are an addition: `vegan::simper()`
#' doesn't adjust. A study of a few groups and hundreds of taxa tests many
#' hypotheses, so look at `p_adjusted`, not `p_value`.
#'
#'
#' **Interpretation: location and dispersion**
#'
#' Warton, Wright & Wang (2012) show that SIMPER confounds differences in
#' mean abundance (location) with differences in variability (dispersion):
#' the most variable taxa get the largest contributions even when their
#' mean abundance is the same in both groups, so a high `average` doesn't
#' mean a taxon differs between the groups. Check the groups' dispersions
#' with [permdisp()], compare `mean_a` with `mean_b`, and test the
#' differences of interest taxon by taxon with models of abundance (such as
#' negative binomial GLMs), rather than reading the ranking as evidence.
#' The p-values test a change in a taxon's contribution, which a change in
#' its variance alone can also give.
#'
#'
#' **Speed and memory**
#'
#' A taxon's contribution is 0 for pairs of samples that both lack it, so
#' each taxon costs about \eqn{n/2} operations per sample that has it, and
#' sparse tables are never made dense. The observed sums run in parallel
#' over taxa, and the permutations in parallel over permutations; sums are
#' always in the same order, so results are identical for any `cpus`.
#'
#'
#' **Differences from vegan::simper()**
#'
#' * The result is one data frame (as `summary()` of vegan's result, with
#'   its pairs stacked), with the pair and taxon as columns, and vegan's
#'   `ava`, `avb`, and `cusum` named `mean_a`, `mean_b`, and `cumsum`.
#' * Pairs are named `"a-b"` and come in the order of the groups' factor
#'   levels; vegan names them `"a_b"`, in the order the groups first appear.
#' * vegan permutes all the samples' groups for every pair (a supplied
#'   matrix covers all of them), so its p-values for a pair also depend on
#'   the other groups' samples. ecodive permutes only the pair's samples, as
#'   its pairwise tests do.
#' * p-values count permutations whose `average` ties the observed one, as
#'   Clarke (1993) does; vegan counts only larger ones (see
#'   **Permutations**).
#' * `p_adjusted`, `norm`, and `strata` given directly.
#'
#'
#' @section Same names in other packages:
#'   `vegan::simper()` has the same name: whichever package was attached
#'   last masks the other. Call `ecodive::simper()` or `vegan::simper()` to
#'   be sure.
#'
#' @references
#' Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in
#' community structure. *Australian Journal of Ecology*, 18, 117-143.
#' \doi{10.1111/j.1442-9993.1993.tb00438.x}
#'
#' Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based
#' multivariate analyses confound location and dispersion effects.
#' *Methods in Ecology and Evolution*, 3, 89-101.
#' \doi{10.1111/j.2041-210X.2011.00127.x}
#'
#' @seealso [permdisp()], [permanova()], [bray()]
#'
#' @export
#' @examples
#'     # Oral (Saliva, Gums) versus other (Nose, Stool) samples.
#'     grp <- c('oral', 'oral', 'other', 'other')
#'     # Four samples are far too few for permutation p-values.
#'     res <- simper(ex_counts, grp, permutations = 0, cpus = 1)
#'     res
#'
#'     # The averages sum to the mean between-group Bray-Curtis dissimilarity.
#'     sum(res$average)
#'     mean(as.matrix(bray(ex_counts))[1:2, 3:4])
simper <- function (
    counts, groups, data = NULL, margin = 1L, norm = 'none',
    p_adjust = 'holm', permutations = 999L, strata = NULL, seed = 0L,
    cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL

  if (is.character(norm) && length(norm) == 1 && tolower(norm) %in% c('clr', 'rclr'))
    stop(
      "`norm = '", norm, "'` gives negative values, which Bray-Curtis ",
      "can't use. Use 'none', 'percent', 'chord', 'hellinger', or 'binary'.",
      call. = FALSE )

  validate_counts()
  validate_norm()
  validate_margin()
  validate_seed()
  validate_cpus()

  p_adjust <- pairwise_method(TRUE, p_adjust, TRUE)


  # Dimensions and names, with samples first.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  n <- as.integer(dims[[1]])
  p <- as.integer(dims[[2]])
  taxa <- if (is.null(dn[[2]])) as.character(seq_len(p)) else as.character(dn[[2]])

  pm <- perm_matrix(permutations, n)


  #--------------------------------------------------------
  # Each sample's group (and stratum), in counts' order.
  #--------------------------------------------------------
  des <- test_group_values(groups, data, n, dn[[1]], strata, ref = '`counts`')
  grp <- factor(des$groups, exclude = NULL)  # a factor's NA level is a group
  lev <- levels(grp)
  gi  <- as.integer(grp)

  if (length(lev) < 2)
    stop('`groups` must have at least two different values.', call. = FALSE)

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  #--------------------------------------------------------
  # Each pair of groups: its samples, labels, permutations.
  #--------------------------------------------------------
  combos <- utils::combn(length(lev), 2)
  labs   <- paste(lev[combos[1, ]], lev[combos[2, ]], sep = '-')

  pairs <- lapply(seq_len(ncol(combos)), function (k) {
    idx <- which(gi %in% combos[, k])
    list(
      idx    = idx - 1L,
      lab    = as.integer(gi[idx] == combos[2, k]),
      perms  = if (is.null(pm$perms)) NULL
               else t(induced_perms(permutations, idx)) - 1L,
      strata = if (is.null(des$strata)) NULL
               else as.integer(factor(des$strata[idx])) - 1L,
      seed   = pair_seed(seed, k),
      n_perm = pm$n_perm )
  })

  res <- .Call(C_simper, counts, margin, norm, gi - 1L, length(lev), pairs, cpus)

  if (any(res$totals == 0)) {
    empty <- if (is.null(dn[[1]])) which(res$totals == 0) else dn[[1]][res$totals == 0]
    stop(
      'Samples with no counts: ', paste(utils::head(empty, 10), collapse = ', '),
      if (length(empty) > 10) ', ...', '. Bray-Curtis is undefined between ',
      'two empty samples; remove them first.', call. = FALSE )
  }


  #--------------------------------------------------------
  # One block of rows per pair, by decreasing average.
  #--------------------------------------------------------
  overall <- numeric(length(labs))

  blocks <- lapply(seq_along(labs), function (k) {

    r    <- res$pairs[[k]]
    np   <- prod(tabulate(gi[gi %in% combos[, k]], length(lev))[combos[, k]])
    avg  <- r$sum / np
    sd   <- if (np > 1) sqrt(r$ss / (np - 1)) else rep(NA_real_, p)
    pv   <- if (pm$n_perm > 0) (r$count + 1) / (pm$n_perm + 1) else rep(NA_real_, p)
    pv[avg == 0] <- NA_real_

    ord  <- order(avg, decreasing = TRUE)
    tot  <- sum(avg)
    overall[[k]] <<- tot

    data.frame(
      pair       = labs[[k]],
      taxon      = taxa[ord],
      average    = avg[ord],
      sd         = sd[ord],
      ratio      = avg[ord] / sd[ord],
      mean_a     = res$means[ord, combos[1, k]],
      mean_b     = res$means[ord, combos[2, k]],
      cumsum     = cumsum(avg[ord]) / tot,
      p_value    = pv[ord],
      p_adjusted = stats::p.adjust(pv[ord], p_adjust),
      stringsAsFactors = FALSE )
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL

  names(overall) <- labs
  attr(out, 'overall')     <- overall
  attr(out, 'n_perm')      <- pm$n_perm
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'permutation') <- if (!is.null(pm$perms)) 'supplied matrix'
                              else if (!is.null(des$strata)) 'within strata'
                              else 'free'
  class(out) <- c('ecodive_simper', 'data.frame')
  out
}
