# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Baselga's partition of beta diversity into turnover and nestedness (and
# of abundance-based dissimilarity into balanced variation and abundance
# gradients): the eight component metrics, beta_partition() and
# beta_partition_multi().
#
# The pairwise components follow abdiv 0.2.0's R/beta_components.R
# (jaccard_turnover(), sorenson_turnover(), bray_curtis_balanced(),
# ruzicka_balanced() and the nestedness / gradient functions, computed as
# total minus turnover), ported to C in src/beta_div.c:
#
#   Copyright (c) 2020 Kyle Bittinger
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
#   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
#   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# The multiple-site values (beta_partition_multi()) are written from the
# papers (Baselga 2010, 2012, 2017) and checked against betapart 1.6.1's
# outputs only: betapart is GPL, and its source was not read.


# Integer IDs for C code (src/beta_div.c).
BDIV_SORENSEN_TURNOVER   <- 25L
BDIV_SORENSEN_NESTEDNESS <- 26L
BDIV_JACCARD_TURNOVER    <- 27L
BDIV_JACCARD_NESTEDNESS  <- 28L
BDIV_BRAY_BALANCED       <- 29L
BDIV_BRAY_GRADIENT       <- 30L
BDIV_SOERGEL_BALANCED    <- 31L
BDIV_SOERGEL_GRADIENT    <- 32L


# The families of beta_partition() and beta_partition_multi(), keyed by
# `index`: the C code, the display name of the total dissimilarity, the
# names of the three parts, and the metric IDs of the three matrices.
PARTITIONS <- data.frame(
  index      = c('sorensen', 'jaccard', 'bray', 'ruzicka'),
  code       = 1:4,
  name       = c('Sorensen Dissimilarity', 'Jaccard Dissimilarity',
                 'Bray-Curtis Dissimilarity', 'Ruzicka Dissimilarity'),
  incidence  = c(TRUE, TRUE, FALSE, FALSE),
  part1      = c('turnover', 'turnover', 'balanced', 'balanced'),
  part2      = c('nestedness', 'nestedness', 'gradient', 'gradient'),
  label1     = c('Turnover', 'Turnover', 'Balanced variation', 'Balanced variation'),
  label2     = c('Nestedness', 'Nestedness', 'Abundance gradient', 'Abundance gradient'),
  metric1    = c('sorensen_turnover', 'jaccard_turnover', 'bray_balanced', 'soergel_balanced'),
  metric2    = c('sorensen_nestedness', 'jaccard_nestedness', 'bray_gradient', 'soergel_gradient'),
  metric3    = c('sorensen', 'jaccard', 'bray', 'soergel'),
  stringsAsFactors = FALSE )


# The PARTITIONS row of `index`: 'sorensen', 'jaccard', 'bray' or
# 'ruzicka' ('soergel' is taken as 'ruzicka'), in any case.
part_index <- function (index) {
  ok <- is.character(index) && length(index) == 1 && !is.na(index)
  if (ok) {
    index <- tolower(index)
    if (identical(index, 'soergel')) index <- 'ruzicka'
    ok <- index %in% PARTITIONS$index
  }
  if (!ok)
    stop(
      '`index` must be one of "sorensen", "jaccard", "bray" or "ruzicka".',
      call. = FALSE )
  as.list(PARTITIONS[match(index, PARTITIONS$index), ])
}


# Checks `norm` for an abundance partition: the norms that give negative
# values are refused, as min() and max() of abundances need values of at
# least 0. With `incidence`, only 'none' is accepted. Other values are left
# to validate_norm().
part_norm <- function (norm, incidence) {
  if (!is.character(norm) || length(norm) != 1 || is.na(norm))
    stop('`norm` must be a single string.', call. = FALSE)
  m <- pmatch(tolower(norm), c(
    'none', 'percent', 'chord', 'hellinger', 'binary', 'clr', 'rclr', 'log' ))
  if (is.na(m)) return (invisible())
  if (incidence && m != 1L)
    stop(
      "`norm` is only used by the abundance indices ('bray', 'ruzicka'); ",
      "the incidence indices use presence and absence.", call. = FALSE )
  if (m >= 6L)
    stop(
      "`norm = '", norm, "'` gives negative values, which the abundance ",
      "partition can't use. Use 'none', 'percent', 'chord', 'hellinger' ",
      "or 'binary'.", call. = FALSE )
  invisible()
}


# Refuses negative values in `counts` (already validated).
part_nonneg <- function (counts) {
  if (isTRUE(any(stored_values(counts) < 0)))
    stop(
      '`counts` contains negative values; the partition of beta diversity ',
      'needs abundances of at least 0.', call. = FALSE )
  invisible()
}


# One component metric, once validate_args() has checked the arguments in
# the caller's frame: runs kernel `code` and labels the result with
# `metric`.
part_metric <- function (code, metric, env = parent.frame()) {
  part_nonneg(env$counts)
  dist_result(
    .Call(C_beta_div, code, env$counts, env$margin, env$norm, env$pairs, env$cpus, 0, NULL),
    metric, env = env )
}




#' Turnover and nestedness of presence/absence dissimilarity
#'
#' Baselga's partition of the Sorensen and Jaccard dissimilarities between
#' two samples into **turnover** (features replaced by others) and
#' **nestedness-resultant** dissimilarity (one sample holding a subset of
#' the other's features).
#'
#' @inherit documentation
#'
#' @param counts   A numeric matrix of count data (samples \eqn{\times}
#'        features). Counts are converted to presence/absence. Values must
#'        not be negative.
#'
#' @family Partitioned metrics
#' @seealso [beta_partition()] for the three matrices of a family at once,
#'   [beta_partition_multi()] for whole sets of samples, [beta_div()].
#'
#' @details
#' With \eqn{a} features in both samples, and \eqn{b} and \eqn{c} only in
#' the first and only in the second (Baselga 2010, 2012):
#'
#' | Function | Component | Formula |
#' | :------- | :-------- | :------ |
#' | `sorensen()` | Sorensen dissimilarity (total) | \eqn{\frac{b + c}{2a + b + c}} |
#' | `sorensen_turnover()` | Simpson dissimilarity, \eqn{\beta_{sim}} | \eqn{\frac{\min(b, c)}{a + \min(b, c)}} |
#' | `sorensen_nestedness()` | nestedness-resultant, \eqn{\beta_{sne}} | \eqn{\beta_{sor} - \beta_{sim}} |
#' | `jaccard()` | Jaccard dissimilarity (total) | \eqn{\frac{b + c}{a + b + c}} |
#' | `jaccard_turnover()` | \eqn{\beta_{jtu}} | \eqn{\frac{2\min(b, c)}{a + 2\min(b, c)}} |
#' | `jaccard_nestedness()` | \eqn{\beta_{jne}} | \eqn{\beta_{jac} - \beta_{jtu}} |
#'
#' Turnover is 0 when one sample's features are a subset of the other's,
#' and 1 when the samples share none; nestedness is 0 when they have the
#' same number of features. Turnover and nestedness add up to the total,
#' [sorensen()] or [jaccard()].
#'
#' When a sample has no features, turnover and nestedness are undefined
#' (`NaN`), as in abdiv and betapart, although the total is 1 (or `NaN`
#' when both samples are empty).
#'
#' Neither component is a metric: a sample nested in two others has
#' turnover 0 to both, while they can have turnover 1 to each other.
#'
#' The partition is debated: Podani & Schmera (2011) and Carvalho et al.
#' (2012) split the same totals differently, into replacement and richness
#' difference (adespatial's `beta.div.comp()`). See Baselga & Leprieur
#' (2015) for a comparison.
#'
#' **Same names in other packages:** abdiv's functions are named
#' `sorenson_turnover()` and `sorenson_nestedness()` (with an "o"), and
#' `jaccard_turnover()` and `jaccard_nestedness()`, taking two vectors; they
#' give the same values. betapart's `beta.pair()` returns all three
#' matrices of a family, as [beta_partition()] does.
#'
#' @return A distance matrix: a `dist` object of class
#'   `c("ecodive_dist", "dist")`, as [beta_div()] returns.
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components
#' of beta diversity. *Global Ecology and Biogeography*, 19, 134-143.
#' \doi{10.1111/j.1466-8238.2009.00490.x}
#'
#' Baselga, A. (2012). The relationship between species replacement,
#' dissimilarity derived from nestedness, and nestedness. *Global Ecology
#' and Biogeography*, 21, 1223-1232. \doi{10.1111/j.1466-8238.2011.00756.x}
#'
#' Lennon, J. J., Koleff, P., Greenwood, J. J. D., & Gaston, K. J. (2001).
#' The geographical structure of British bird distributions: diversity,
#' spatial turnover and scale. *Journal of Animal Ecology*, 70, 966-979.
#' \doi{10.1046/j.0021-8790.2001.00563.x}
#'
#' Podani, J., & Schmera, D. (2011). A new conceptual and methodological
#' framework for exploring and explaining pattern in presence-absence data.
#' *Oikos*, 120, 1625-1638. \doi{10.1111/j.1600-0706.2011.19451.x}
#'
#' Carvalho, J. C., Cardoso, P., & Gomes, P. (2012). Determining the
#' relative roles of species replacement and species richness differences
#' in generating beta-diversity patterns. *Global Ecology and
#' Biogeography*, 21, 760-771. \doi{10.1111/j.1466-8238.2011.00694.x}
#'
#' Baselga, A., & Leprieur, F. (2015). Comparing methods to separate
#' components of beta diversity. *Methods in Ecology and Evolution*, 6,
#' 1069-1079. \doi{10.1111/2041-210X.12388}
#'
#' @export
#' @examples
#'     sorensen_turnover(ex_counts)
#'     sorensen_nestedness(ex_counts)
#'
#'     # Turnover + nestedness = total.
#'     all.equal(
#'       as.vector(jaccard_turnover(ex_counts) + jaccard_nestedness(ex_counts)),
#'       as.vector(jaccard(ex_counts)) )
sorensen_turnover <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {
  norm <- 'none'
  validate_args()
  part_metric(BDIV_SORENSEN_TURNOVER, 'sorensen_turnover')
}

#' @rdname sorensen_turnover
#' @export
sorensen_nestedness <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {
  norm <- 'none'
  validate_args()
  part_metric(BDIV_SORENSEN_NESTEDNESS, 'sorensen_nestedness')
}

#' @rdname sorensen_turnover
#' @export
jaccard_turnover <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {
  norm <- 'none'
  validate_args()
  part_metric(BDIV_JACCARD_TURNOVER, 'jaccard_turnover')
}

#' @rdname sorensen_turnover
#' @export
jaccard_nestedness <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {
  norm <- 'none'
  validate_args()
  part_metric(BDIV_JACCARD_NESTEDNESS, 'jaccard_nestedness')
}




#' Balanced variation and abundance gradients of abundance dissimilarity
#'
#' Baselga's partition of the Bray-Curtis and Ruzicka (Soergel)
#' dissimilarities between two samples into **balanced variation** in
#' abundance (some features gain what others lose) and **abundance
#' gradients** (every feature is less abundant in one sample).
#'
#' @inherit documentation
#'
#' @param counts   A numeric matrix of count data (samples \eqn{\times}
#'        features). Values must not be negative.
#'
#' @param norm   Normalize the incoming counts first: `'none'` (the
#'        default, as betapart and abdiv), `'percent'` (relative
#'        abundances), `'chord'`, `'hellinger'` or `'binary'`. `'clr'`,
#'        `'rclr'` and `'log'` give negative values and are refused.
#'
#' @family Partitioned metrics
#' @seealso [beta_partition()] for the three matrices of a family at once,
#'   [beta_partition_multi()] for whole sets of samples, [beta_div()].
#'
#' @details
#' With abundances \eqn{x_k} and \eqn{y_k} of each feature in the two
#' samples, \eqn{A = \sum_k \min(x_k, y_k)} (the abundance they share),
#' \eqn{B = \sum_k \max(x_k - y_k, 0)} and \eqn{C = \sum_k \max(y_k - x_k,
#' 0)} (what only the first and only the second has); Baselga (2013, 2017):
#'
#' | Function | Component | Formula |
#' | :------- | :-------- | :------ |
#' | `bray()` | Bray-Curtis dissimilarity (total) | \eqn{\frac{B + C}{2A + B + C}} |
#' | `bray_balanced()` | balanced variation, \eqn{d_{BC-bal}} | \eqn{\frac{\min(B, C)}{A + \min(B, C)}} |
#' | `bray_gradient()` | abundance gradient, \eqn{d_{BC-gra}} | \eqn{d_{BC} - d_{BC-bal}} |
#' | `soergel()` | Ruzicka dissimilarity (total) | \eqn{\frac{B + C}{A + B + C}} |
#' | `soergel_balanced()` | balanced variation, \eqn{d_{RUZ-bal}} | \eqn{\frac{2\min(B, C)}{A + 2\min(B, C)}} |
#' | `soergel_gradient()` | abundance gradient, \eqn{d_{RUZ-gra}} | \eqn{d_{RUZ} - d_{RUZ-bal}} |
#'
#' [soergel()] is the Ruzicka dissimilarity (vegan's quantitative Jaccard),
#' so its components are named after it; `"ruzicka_balanced"` and
#' `"ruzicka_gradient"` select them in [beta_div()]. The two components add
#' up to the total, [bray()] or [soergel()]. On presence/absence data
#' (`norm = 'binary'`) they are the turnover and nestedness of
#' [sorensen_turnover()].
#'
#' When a sample is empty (all zeros), both components are undefined
#' (`NaN`), as in abdiv and betapart. Neither is a metric.
#'
#' **Same names in other packages:** abdiv's functions are named
#' `bray_curtis_balanced()`, `bray_curtis_gradient()`, `ruzicka_balanced()`
#' and `ruzicka_gradient()`, taking two vectors; they give the same values.
#' betapart's `beta.pair.abund()` returns all three matrices of a family,
#' as [beta_partition()] does.
#'
#' @return A distance matrix: a `dist` object of class
#'   `c("ecodive_dist", "dist")`, as [beta_div()] returns.
#'
#' @references
#' Baselga, A. (2013). Separating the two components of abundance-based
#' dissimilarity: balanced changes in abundance vs. abundance gradients.
#' *Methods in Ecology and Evolution*, 4, 552-557.
#' \doi{10.1111/2041-210X.12029}
#'
#' Baselga, A. (2017). Partitioning abundance-based multiple-site
#' dissimilarity into components: balanced variation in abundance and
#' abundance gradients. *Methods in Ecology and Evolution*, 8, 799-808.
#' \doi{10.1111/2041-210X.12693}
#'
#' @export
#' @examples
#'     bray_balanced(ex_counts)
#'     bray_gradient(ex_counts)
#'
#'     # Of relative abundances.
#'     soergel_balanced(ex_counts, norm = 'percent')
bray_balanced <- function (counts, margin = 1L, norm = 'none', pairs = NULL, clust = NULL, cpus = n_cpus()) {
  part_norm(norm, FALSE)
  validate_args()
  part_metric(BDIV_BRAY_BALANCED, 'bray_balanced')
}

#' @rdname bray_balanced
#' @export
bray_gradient <- function (counts, margin = 1L, norm = 'none', pairs = NULL, clust = NULL, cpus = n_cpus()) {
  part_norm(norm, FALSE)
  validate_args()
  part_metric(BDIV_BRAY_GRADIENT, 'bray_gradient')
}

#' @rdname bray_balanced
#' @export
soergel_balanced <- function (counts, margin = 1L, norm = 'none', pairs = NULL, clust = NULL, cpus = n_cpus()) {
  part_norm(norm, FALSE)
  validate_args()
  part_metric(BDIV_SOERGEL_BALANCED, 'soergel_balanced')
}

#' @rdname bray_balanced
#' @export
soergel_gradient <- function (counts, margin = 1L, norm = 'none', pairs = NULL, clust = NULL, cpus = n_cpus()) {
  part_norm(norm, FALSE)
  validate_args()
  part_metric(BDIV_SOERGEL_GRADIENT, 'soergel_gradient')
}




#' Partition beta diversity into turnover and nestedness
#'
#' Baselga's partition of the dissimilarity between every pair of samples
#' into turnover and nestedness (presence/absence), or into balanced
#' variation and abundance gradients (abundances): all three distance
#' matrices of one family, computed in one pass, as betapart's
#' `beta.pair()` and `beta.pair.abund()` return them.
#'
#' @inherit documentation
#'
#' @param counts   A numeric matrix of count data (samples \eqn{\times}
#'        features). Values must not be negative.
#'
#' @param index   The family of dissimilarities, as betapart's
#'        `index.family`: `"sorensen"` (the default) or `"jaccard"`, of
#'        presence/absence, or `"bray"` (Bray-Curtis) or `"ruzicka"`
#'        (`"soergel"` is accepted too), of abundances.
#'
#' @param norm   Only for `index = "bray"` or `"ruzicka"`. Normalize the
#'        incoming counts first: `'none'` (the default, as betapart),
#'        `'percent'` (relative abundances), `'chord'`, `'hellinger'` or
#'        `'binary'`. `'clr'`, `'rclr'` and `'log'` give negative values and
#'        are refused.
#'
#' @param clust  Whether and how to cluster the samples, as for any metric
#'        (see [beta_div()]): the samples are ordered by the clustering of
#'        the **total** dissimilarity, and all three matrices are put in
#'        that order, so that they stay aligned. `NULL` (the default) keeps
#'        the input order.
#'
#' @details
#' The three matrices are those of the family's three metrics:
#'
#' | `index` | `turnover` / `balanced` | `nestedness` / `gradient` | `total` |
#' | :------ | :---------------------- | :------------------------ | :------ |
#' | `"sorensen"` | [sorensen_turnover()] | [sorensen_nestedness()] | [sorensen()] |
#' | `"jaccard"`  | [jaccard_turnover()]  | [jaccard_nestedness()]  | [jaccard()]  |
#' | `"bray"`     | [bray_balanced()]     | [bray_gradient()]       | [bray()]     |
#' | `"ruzicka"`  | [soergel_balanced()]  | [soergel_gradient()]    | [soergel()]  |
#'
#' and are identical to theirs. Their formulas, and the debate about the
#' partition, are on those functions' pages.
#'
#' **Differences from betapart:** the abundance indices can use
#' relative abundances (`norm`); the matrices are named `turnover`,
#' `nestedness` and `total` (or `balanced`, `gradient` and `total`) rather
#' than `beta.sim`, `beta.sne`, `beta.sor`, ...; and `pairs` computes only
#' some of the distances.
#'
#' @return A list of class `ecodive_beta_partition`, with three distance
#'   matrices (`ecodive_dist` objects): `turnover`, `nestedness` and `total`
#'   for the presence/absence indices, or `balanced`, `gradient` and
#'   `total` for the abundance indices. Attribute `index` is the family.
#'   `plot()` draws each pair's two components against each other (see
#'   [plot.ecodive_beta_partition()]).
#'
#' @family Partitioned metrics
#' @seealso [beta_partition_multi()] for whole sets of samples.
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components
#' of beta diversity. *Global Ecology and Biogeography*, 19, 134-143.
#' \doi{10.1111/j.1466-8238.2009.00490.x}
#'
#' Baselga, A. (2013). Separating the two components of abundance-based
#' dissimilarity: balanced changes in abundance vs. abundance gradients.
#' *Methods in Ecology and Evolution*, 4, 552-557.
#' \doi{10.1111/2041-210X.12029}
#'
#' Baselga, A., & Orme, C. D. L. (2012). betapart: an R package for the
#' study of beta diversity. *Methods in Ecology and Evolution*, 3, 808-812.
#' \doi{10.1111/j.2041-210X.2012.00224.x}
#'
#' @export
#' @examples
#'     bp <- beta_partition(ex_counts)
#'     bp
#'     bp$turnover
#'
#'     plot(bp)
#'
#'     beta_partition(ex_counts, 'bray', norm = 'percent')
beta_partition <- function (
    counts, index = 'sorensen', margin = 1L, norm = 'none',
    pairs = NULL, clust = NULL, cpus = n_cpus() ) {

  fam <- part_index(index)
  part_norm(norm, fam$incidence)

  validate_counts()
  part_nonneg(counts)
  validate_margin()
  validate_norm()
  validate_pairs()
  validate_clust()
  validate_cpus()

  res <- .Call(C_beta_partition, fam$code, counts, margin, norm, pairs, cpus)

  ids <- c(fam$metric1, fam$metric2, fam$metric3)
  res <- lapply(1:3, function (i) {
    d <- res[[i]]
    attr(d, 'method') <- ids[[i]]
    set_display(d, title = METRICS$name[match(ids[[i]], METRICS$id)])
  })

  # Cluster by the total, and put the other two in the same order.
  if (!is.null(clust)) {
    total <- dist_clust(res[[3]], clust, cpus)
    ord   <- attr(total, 'order', exact = TRUE)
    res[1:2] <- lapply(res[1:2], function (d) {
      if (!is.null(ord)) {
        d <- dist_reorder(d, ord)
        attr(d, 'order') <- ord
      }
      set_display(d, clust = clust)
    })
    res[[3]] <- total
  }

  names(res) <- c(fam$part1, fam$part2, 'total')
  structure(res, class = 'ecodive_beta_partition', index = fam$index)
}




#' Multiple-site partition of beta diversity
#'
#' Baselga's multiple-site dissimilarity of a whole set of samples, and its
#' partition into turnover and nestedness (presence/absence) or balanced
#' variation and abundance gradients (abundances), as betapart's
#' `beta.multi()` and `beta.multi.abund()` compute them; and, with
#' `n_sites`, the same for random subsets of a fixed number of samples, as
#' `beta.sample()` and `beta.sample.abund()`.
#'
#' @inherit beta_partition params
#' @inherit documentation params
#'
#' @param n_sites   `NULL` (the default) for the values of all the samples
#'        only; or a number of samples (at least 2, at most the number of
#'        samples) to draw, without replacement, `n_draws` times
#'        (betapart's `sites`).
#'
#' @param n_draws   How many random subsets of `n_sites` samples to draw
#'        (betapart's `samples`). Default: `100`.
#'
#' @param seed   Random seed for the draws: any integer. Results are the
#'        same for a given seed, whatever `cpus`.
#'
#' @details
#' With \eqn{N} samples, \eqn{S_i} the number of features of sample
#' \eqn{i}, \eqn{S_T} the number of features of all the samples together,
#' and \eqn{b_{ij}} the number of features in sample \eqn{i} but not in
#' sample \eqn{j}, Baselga (2010, 2012) defines, with sums over the pairs
#' \eqn{i < j}:
#'
#' * \eqn{\beta_{SOR} = \frac{\sum \min(b_{ij}, b_{ji}) + \sum \max(b_{ij}, b_{ji})}{2(\sum_i S_i - S_T) + \sum \min(b_{ij}, b_{ji}) + \sum \max(b_{ij}, b_{ji})}}
#' * \eqn{\beta_{SIM} = \frac{\sum \min(b_{ij}, b_{ji})}{\sum_i S_i - S_T + \sum \min(b_{ij}, b_{ji})}}, and \eqn{\beta_{SNE} = \beta_{SOR} - \beta_{SIM}}
#' * \eqn{\beta_{JAC} = \frac{\sum \min(b_{ij}, b_{ji}) + \sum \max(b_{ij}, b_{ji})}{\sum_i S_i - S_T + \sum \min(b_{ij}, b_{ji}) + \sum \max(b_{ij}, b_{ji})}}
#' * \eqn{\beta_{JTU} = \frac{2 \sum \min(b_{ij}, b_{ji})}{\sum_i S_i - S_T + 2 \sum \min(b_{ij}, b_{ji})}}, and \eqn{\beta_{JNE} = \beta_{JAC} - \beta_{JTU}}
#'
#' For abundances (Baselga 2017), \eqn{b_{ij}} is
#' \eqn{B_{ij} = \sum_k \max(x_{ik} - x_{jk}, 0)}, and
#' \eqn{\sum_i S_i - S_T} is \eqn{A = \sum_k (\sum_i x_{ik} - \max_i
#' x_{ik})}; the Bray-Curtis values then have the shapes of the Sorensen
#' ones, and the Ruzicka values those of the Jaccard ones. With two
#' samples, each value is the pairwise one ([beta_partition()]).
#'
#' **Multiple-site values depend on the number of samples:** they grow as
#' samples are added, so sets of different sizes can't be compared
#' directly. `n_sites` draws subsets of the same size from each set, as
#' Baselga (2010) recommends, and the means (or the distributions) of the
#' draws can be compared instead.
#'
#' **Differences from betapart:** draws use ecodive's own random number
#' generator (PCG; `seed`), not R's, so they differ from `beta.sample()`'s
#' while having the same distribution; the samples of each draw are taken
#' in their input order; `sd` is the standard deviation of the draws, as
#' `sd.values`; and the abundance indices can use relative abundances
#' (`norm`).
#'
#' @return A list of class `ecodive_beta_partition_multi`:
#'   * `value`: the three values of all the samples, named `turnover`,
#'     `nestedness` and `total` (or `balanced`, `gradient` and `total` for
#'     the abundance indices).
#'   * With `n_sites`: `draws`, a matrix of the three values (columns) of
#'     each draw (rows); and `mean` and `sd`, their means and standard
#'     deviations.
#'
#'   Attributes: `index`, the family; `n_samples`, the number of samples
#'   behind `value`; and with `n_sites`, `n_sites` and `n_draws`.
#'
#' @family Partitioned metrics
#' @seealso [beta_partition()] for every pair of samples.
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components
#' of beta diversity. *Global Ecology and Biogeography*, 19, 134-143.
#' \doi{10.1111/j.1466-8238.2009.00490.x}
#'
#' Baselga, A. (2012). The relationship between species replacement,
#' dissimilarity derived from nestedness, and nestedness. *Global Ecology
#' and Biogeography*, 21, 1223-1232. \doi{10.1111/j.1466-8238.2011.00756.x}
#'
#' Baselga, A. (2017). Partitioning abundance-based multiple-site
#' dissimilarity into components: balanced variation in abundance and
#' abundance gradients. *Methods in Ecology and Evolution*, 8, 799-808.
#' \doi{10.1111/2041-210X.12693}
#'
#' @export
#' @examples
#'     beta_partition_multi(ex_counts)
#'
#'     # 50 draws of 4 samples.
#'     bm <- beta_partition_multi(ex_counts, 'jaccard', n_sites = 4, n_draws = 50)
#'     bm
#'     plot(bm)
beta_partition_multi <- function (
    counts, index = 'sorensen', margin = 1L, norm = 'none',
    n_sites = NULL, n_draws = 100L, seed = 0L, cpus = n_cpus() ) {

  fam <- part_index(index)
  part_norm(norm, fam$incidence)

  validate_counts()
  part_nonneg(counts)
  validate_margin()
  validate_norm()
  validate_seed()
  validate_cpus()

  n <- if (margin == 1L) nrow(counts) else ncol(counts)
  if (n < 2)
    stop('`counts` must have at least two samples.', call. = FALSE)

  if (!is.null(n_sites)) {
    if (!is.numeric(n_sites) || length(n_sites) != 1 || is.na(n_sites) ||
        n_sites %% 1 != 0 || n_sites < 2 || n_sites > n)
      stop(sprintf(
        '`n_sites` must be NULL or a whole number from 2 to %d (the number of samples).', n ),
        call. = FALSE )
    if (!is.numeric(n_draws) || length(n_draws) != 1 || is.na(n_draws) ||
        n_draws %% 1 != 0 || n_draws < 1 || n_draws > .Machine$integer.max)
      stop('`n_draws` must be a whole number of at least 1.', call. = FALSE)
    n_sites <- as.integer(n_sites)
    n_draws <- as.integer(n_draws)
  }

  parts <- c(fam$part1, fam$part2, 'total')

  value <- .Call(C_beta_multi, fam$code, counts, margin, norm, 0L, 1L, seed, cpus)
  value <- stats::setNames(as.vector(value), parts)
  res   <- list(value = value)

  if (!is.null(n_sites)) {
    draws <- .Call(C_beta_multi, fam$code, counts, margin, norm, n_sites, n_draws, seed, cpus)
    colnames(draws) <- parts
    res$draws <- draws
    res$mean  <- colMeans(draws)
    res$sd    <- apply(draws, 2L, stats::sd)
  }

  structure(
    res, class = 'ecodive_beta_partition_multi', index = fam$index,
    n_samples = n,
    n_sites   = if (!is.null(n_sites)) n_sites,
    n_draws   = if (!is.null(n_sites)) n_draws )
}




#' @export
print.ecodive_beta_partition <- function (x, digits = 3L, ...) {
  fam  <- part_index(attr(x, 'index', exact = TRUE))
  size <- attr(x[[3]], 'Size', exact = TRUE)
  cat(sprintf(
    'Partition of %s: %d samples, %d pairs\n\n',
    fam$name, size, length(x[[3]]) ))
  labs <- c(fam$label1, fam$label2, 'Total')
  for (i in 1:3) {
    v   <- unclass(x[[i]])
    ok  <- !is.na(v)
    cat(sprintf(
      '  %-20s %-20s mean %s%s\n',
      paste0('$', names(x)[[i]]), labs[[i]],
      if (any(ok)) formatC(mean(v[ok]), digits = digits, format = 'g') else 'NA',
      if (any(!ok)) sprintf(' (%d NA or NaN)', sum(!ok)) else '' ))
  }
  invisible(x)
}


#' @export
print.ecodive_beta_partition_multi <- function (x, digits = 3L, ...) {
  fam <- part_index(attr(x, 'index', exact = TRUE))
  fmt <- function (v) formatC(v, digits = digits, format = 'g')
  cat(sprintf(
    'Multiple-site partition of %s: %d samples\n\n',
    fam$name, attr(x, 'n_samples', exact = TRUE) ))
  labs <- c(fam$label1, fam$label2, 'Total')
  if (is.null(x$draws)) {
    for (i in 1:3)
      cat(sprintf('  %-20s %s\n', labs[[i]], fmt(x$value[[i]])))
  }
  else {
    cat(sprintf(
      '  %-20s %-12s %s\n', '',  'All samples',
      sprintf('%d draws of %d: mean (sd)',
              attr(x, 'n_draws', exact = TRUE), attr(x, 'n_sites', exact = TRUE)) ))
    for (i in 1:3)
      cat(sprintf(
        '  %-20s %-12s %s (%s)\n', labs[[i]], fmt(x$value[[i]]),
        fmt(x$mean[[i]]), fmt(x$sd[[i]]) ))
  }
  invisible(x)
}




#' Plot a partition of beta diversity
#'
#' For [beta_partition()]: each pair of samples as a point, its turnover
#' (or balanced variation) against its nestedness (or abundance gradient),
#' with dashed lines where the total dissimilarity is 0.25, 0.5, 0.75 and
#' 1. For [beta_partition_multi()]: with draws, the density of each
#' component's values over the draws, as betapart's vignette plots them;
#' without, one stacked bar of the two components of all the samples.
#'
#' @param x  The result of [beta_partition()] or [beta_partition_multi()].
#'
#' @param groups  Optional, for [beta_partition()]: a vector or factor
#'        with each sample's group, in the input order of the samples.
#'        Pairs within a group are coloured by group; pairs between groups
#'        are grey.
#'
#' @param col  Colours: of the groups, or of the three components (for
#'        [beta_partition_multi()]). Defaults to [grDevices::hcl.colors()]'s
#'        `"Dark 3"` palette.
#'
#' @param legend  Where to put the legend: a position keyword for
#'        [graphics::legend()] (`"topright"`, ...), or `NULL` for none.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `xlab`, `pch` or `cex`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     bp <- beta_partition(ex_counts)
#'     plot(bp, groups = c('a', 'a', 'b', 'b'))
#'
#'     plot(beta_partition_multi(ex_counts, n_sites = 3, n_draws = 50))
plot.ecodive_beta_partition <- function (x, groups = NULL, col = NULL, legend = 'topright', ...) {

  fam    <- part_index(attr(x, 'index', exact = TRUE))
  legend <- plot_legend_pos(legend)
  tu     <- as.vector(x[[1]])
  ne     <- as.vector(x[[2]])
  n      <- attr(x[[3]], 'Size', exact = TRUE)

  # Each pair's group, NA between groups.
  pair_grp <- NULL
  if (!is.null(groups)) {
    if (length(groups) != n)
      stop(sprintf('`groups` must have one value per sample (%d).', n), call. = FALSE)
    ord <- attr(x[[3]], 'order', exact = TRUE)
    if (!is.null(ord)) groups <- groups[ord]
    groups <- as.factor(groups)
    ij     <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
    gi     <- groups[ij[, 2]]
    gj     <- groups[ij[, 1]]
    pair_grp <- factor(ifelse(as.character(gi) == as.character(gj), as.character(gi), NA),
                       levels = levels(groups))
  }

  args <- list(...)
  pcol <- PLOT_MUTED
  if (!is.null(pair_grp)) {
    k <- nlevels(pair_grp)
    if (is.null(col)) col <- plot_palette(k)
    col  <- rep_len(col, k)
    pcol <- ifelse(is.na(pair_grp), PLOT_MUTED, col[as.integer(pair_grp)])
  }
  else if (!is.null(col)) pcol <- col[[1]]

  ok  <- is.finite(tu) & is.finite(ne)
  lim <- c(0, max(1, tu[ok], ne[ok]))

  plot_do(
    graphics::plot.default, c(list(tu[ok], ne[ok]), plot_drop(args, 'col')),
    list(
      xlim = lim, ylim = lim, asp = 1,
      xlab = fam$label1, ylab = fam$label2, pch = 19,
      col  = if (length(pcol) > 1) pcol[ok] else pcol,
      main = paste('Partition of', fam$name),
      sub  = sprintf('%d pairs of samples%s', sum(ok),
                     if (any(!ok)) sprintf(' (%d undefined left out)', sum(!ok)) else '') ))

  for (tot in c(0.25, 0.5, 0.75, 1))
    graphics::segments(0, tot, tot, 0, lty = 2, col = 'grey70')

  if (!is.null(pair_grp))
    plot_legend(
      legend, c(levels(pair_grp), 'Between groups'),
      c(col, PLOT_MUTED), pch = 19 )

  invisible(x)
}


#' @rdname plot.ecodive_beta_partition
#' @export
plot.ecodive_beta_partition_multi <- function (x, col = NULL, legend = 'topright', ...) {

  fam    <- part_index(attr(x, 'index', exact = TRUE))
  legend <- plot_legend_pos(legend)
  labs   <- c(fam$label1, fam$label2, 'Total')
  if (is.null(col)) col <- plot_palette(3L)
  col  <- rep_len(col, 3L)
  args <- list(...)

  if (is.null(x$draws)) {
    v <- matrix(x$value[1:2], ncol = 1L)
    plot_do(
      graphics::barplot, c(list(v), plot_drop(args, 'col')),
      list(
        col = col[1:2], ylim = c(0, max(1, x$value[[3]], na.rm = TRUE)),
        ylab = 'Dissimilarity',
        main = paste('Multiple-site partition of', fam$name),
        sub  = sprintf('%d samples', attr(x, 'n_samples', exact = TRUE)) ))
    plot_legend(legend, labs[2:1], col[2:1], pch = 15)
    return (invisible(x))
  }

  # Each component's density over the draws; a single value (no spread)
  # as a vertical line.
  dens <- lapply(1:3, function (i) {
    v <- x$draws[, i]
    v <- v[is.finite(v)]
    if (length(v) < 2 || stats::sd(v) == 0) return (v[1])
    stats::density(v)
  })
  xs <- unlist(lapply(dens, function (d) if (is.list(d)) d$x else d))
  ys <- unlist(lapply(dens, function (d) if (is.list(d)) d$y else 0))

  plot_do(
    graphics::plot.default, c(list(NA, NA), plot_drop(args, 'col')),
    list(
      type = 'n', xlim = range(c(0, xs, 1), na.rm = TRUE),
      ylim = c(0, max(c(1, ys), na.rm = TRUE)),
      xlab = 'Dissimilarity', ylab = 'Density',
      main = paste('Multiple-site partition of', fam$name),
      sub  = sprintf('%d draws of %d samples',
                     attr(x, 'n_draws', exact = TRUE), attr(x, 'n_sites', exact = TRUE)) ))

  for (i in 1:3) {
    d <- dens[[i]]
    if (is.list(d)) graphics::lines(d, col = col[[i]], lwd = 2)
    else if (length(d) && !is.na(d)) graphics::abline(v = d, col = col[[i]], lwd = 2)
  }
  plot_legend(legend, labs, col, pch = NA, lty = 1)

  invisible(x)
}
