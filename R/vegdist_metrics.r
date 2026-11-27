# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

# The metrics of vegan's vegdist() that ecodive lacked (plans/
# vegdist_metrics.md). Written from vegdist()'s documentation, the cited
# papers, and vegan's outputs; vegan's source (GPL) was not read.


# Refuses negative values in `counts` (already validated).
vd_nonneg <- function (counts) {
  if (isTRUE(any(stored_values(counts) < 0)))
    stop('`counts` contains negative values.', call. = FALSE)
  invisible()
}


#' Kulczynski dissimilarity
#'
#' The quantitative Kulczynski dissimilarity: one minus the mean of the
#' shares of each sample's total that the two samples have in common.
#'
#' @inherit documentation
#' @inherit pseudocount_section
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Kulczynski dissimilarity is defined as:
#' \deqn{1 - \frac{1}{2}\left(\frac{\sum_{i=1}^{n}\min(X_i, Y_i)}{\sum_{i=1}^{n}X_i} + \frac{\sum_{i=1}^{n}\min(X_i, Y_i)}{\sum_{i=1}^{n}Y_i}\right)}
#'
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#'
#' It ranges from 0 (identical samples) to 1 (no shared features). With
#' `norm = 'binary'` it is \eqn{1 - (J/A + J/B)/2}, with \eqn{A} and
#' \eqn{B} the samples' numbers of features and \eqn{J} the number in
#' common, as vegan's `vegdist(x, "kulczynski", binary = TRUE)`.
#'
#' It equals vegan's `vegdist(x, "kulczynski")` and abdiv's
#' `weighted_kulczynski_second()`.
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,]
#' y <- ex_counts[2,]
#' 1 - (sum(pmin(x, y)) / sum(x) + sum(pmin(x, y)) / sum(y)) / 2
#' ```
#'
#' @references
#' Kulczynski, S. (1927). Die Pflanzenassoziationen der Pieninen. *Bulletin International de l'Academie Polonaise des Sciences et des Lettres, Classe des Sciences Mathematiques et Naturelles, Serie B*, Suppl. II, 57-203.
#'
#' Faith, D. P., Minchin, P. R., & Belbin, L. (1987). Compositional dissimilarity as a robust measure of ecological distance. *Vegetatio*, 69, 57-68. \doi{10.1007/BF00038687}
#'
#' @export
#' @examples
#'     kulczynski(ex_counts)
kulczynski <- function (counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_KULCZYNSKI, counts, margin, norm, pairs, cpus, pseudocount, NULL), 'kulczynski')
}


#' Binomial deviance dissimilarity
#'
#' The deviance of a binomial model of each feature's split between two
#' samples (Anderson & Millar 2004).
#'
#' @inherit documentation
#' @inherit pseudocount_section
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The binomial deviance dissimilarity is defined as:
#' \deqn{\sum_{i=1}^{n}\frac{1}{N_i}\left[X_i\ln\frac{X_i}{N_i} + Y_i\ln\frac{Y_i}{N_i} - N_i\ln\frac{1}{2}\right]}
#'
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{N_i = X_i + Y_i}, over the features with \eqn{N_i > 0}, and
#'   \eqn{0 \ln 0 = 0}.
#' * \eqn{n} : The number of features.
#'
#' It is 0 for identical samples (`NaN` for two empty samples) and has no
#' upper limit: each feature found in only one sample adds \eqn{\ln 2}. With `norm = 'binary'` it is
#' \eqn{\ln(2) (A + B - 2J)}, as vegan's `vegdist(x, "binomial", binary =
#' TRUE)`. The values depend on the scale of the counts: proportions give
#' other values than counts.
#'
#' It equals vegan's `vegdist(x, "binomial")` and abdiv's
#' `binomial_deviance()`. Named `binomial_deviance()`, not `binomial()`,
#' to leave [stats::binomial()] unmasked; `beta_div(x, 'binomial')` finds
#' it.
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,]
#' y <- ex_counts[2,]
#' k <- x + y > 0; x <- x[k]; y <- y[k]; n <- x + y
#' xl <- ifelse(x > 0, x * log(x / n), 0)
#' yl <- ifelse(y > 0, y * log(y / n), 0)
#' sum((xl + yl - n * log(1/2)) / n)
#' ```
#'
#' @references
#' Anderson, M. J., & Millar, R. B. (2004). Spatial variation and effects of habitat on temperate reef fish assemblages in northeastern New Zealand. *Journal of Experimental Marine Biology and Ecology*, 305, 191-221. \doi{10.1016/j.jembe.2003.12.011}
#'
#' @export
#' @examples
#'     binomial_deviance(ex_counts)
binomial_deviance <- function (counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_BINOMIAL, counts, margin, norm, pairs, cpus, pseudocount, NULL), 'binomial_deviance')
}


#' Alternative Gower distance
#'
#' The mean absolute difference over the features found in either sample
#' (Anderson et al. 2006), vegan's `"altGower"`.
#'
#' @inherit documentation
#' @inherit pseudocount_section
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `gower()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The alternative Gower distance is defined as:
#' \deqn{\frac{1}{N_Z}\sum_{i=1}^{n}|X_i - Y_i|}
#'
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{N_Z} : The number of features that are not zero in both samples.
#' * \eqn{n} : The number of features.
#'
#' Unlike [gower()], it doesn't scale each feature by its range, and it
#' leaves out the features absent from both samples. With `norm = 'binary'`
#' it is \eqn{(A + B - 2J) / (A + B - J)}, the Jaccard distance, as vegan's
#' `vegdist(x, "altGower", binary = TRUE)`. With `norm = 'clr'` no value is
#' zero, so \eqn{N_Z = n}.
#'
#' Anderson et al. (2006) apply it to log-transformed counts, as vegan's
#' example does with `decostand(x, "log")`: \eqn{\log_2(x) + 1} for
#' \eqn{x > 0}, and 0 for zeros. That is not ecodive's `norm = 'log'`
#' (the log of proportions, see [normalize()]); transform the counts first
#' to follow them.
#'
#' It equals vegan's `vegdist(x, "altGower")` and abdiv's
#' `modified_mean_character_difference()`.
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,]
#' y <- ex_counts[2,]
#' sum(abs(x - y)) / sum(x != 0 | y != 0)
#' ```
#'
#' @references
#' Anderson, M. J., Ellingsen, K. E., & McArdle, B. H. (2006). Multivariate dispersion as a measure of beta diversity. *Ecology Letters*, 9, 683-693. \doi{10.1111/j.1461-0248.2006.00926.x}
#'
#' @export
#' @examples
#'     alt_gower(ex_counts)
alt_gower <- function (counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_ALT_GOWER, counts, margin, norm, pairs, cpus, pseudocount, NULL), 'alt_gower')
}


#' Cao dissimilarity (CYd)
#'
#' The Cao (CY) index of dissimilarity (Cao et al. 1997), for counts.
#'
#' @inherit documentation
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Cao dissimilarity is defined as:
#' \deqn{\frac{1}{S}\sum_{i=1}^{S}\left[\ln\frac{X_i + Y_i}{2} - \frac{X_i\ln Y_i + Y_i\ln X_i}{X_i + Y_i}\right]}
#'
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature, each
#'   raised to at least 0.1 (so zeros count as 0.1).
#' * \eqn{S} : The number of features found in either sample (the sum is
#'   over those).
#'
#' It is 0 for identical samples, has no upper limit, and uses natural
#' logarithms. It is intended for counts: the 0.1 that stands in for zeros
#' assumes that the smallest non-zero value is 1. ecodive doesn't check
#' for whole numbers (vegan warns).
#'
#' It equals vegan's `vegdist(x, "cao")`. abdiv's `cy_dissimilarity()` is a
#' different scaling (base-10 logarithms).
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,]
#' y <- ex_counts[2,]
#' k <- x > 0 | y > 0
#' x <- pmax(x[k], 0.1); y <- pmax(y[k], 0.1); n <- x + y
#' mean(log(n / 2) - (x * log(y) + y * log(x)) / n)
#' ```
#'
#' @references
#' Cao, Y., Williams, W. P., & Bark, A. W. (1997). Similarity measure bias in river benthic Aufwuchs community analysis. *Water Environment Research*, 69, 95-106. \doi{10.2175/106143097X125227}
#'
#' @export
#' @examples
#'     cao(ex_counts)
cao <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  norm <- 'none'
  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_CAO, counts, margin, norm, pairs, cpus, 0, NULL), 'cao')
}


#' Chao-Jaccard dissimilarity
#'
#' Chao et al.'s (2005) abundance-based Jaccard dissimilarity, which
#' estimates the shares of individuals that belong to shared species,
#' including shared species not seen.
#'
#' @inherit documentation
#' @inherit bdiv_assert_integer
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Chao-Jaccard dissimilarity is defined as:
#' \deqn{1 - \frac{UV}{U + V - UV}}
#'
#' with
#' \deqn{U = \min\left(1,\ \frac{C_X}{N_X} + \frac{N_Y - 1}{N_Y}\cdot\frac{f_{1}^{Y}}{2f_{2}^{Y}}\cdot\frac{S_X}{N_X}\right)}
#' and \eqn{V} the same with \eqn{X} and \eqn{Y} swapped.
#'
#' Where:
#' * \eqn{N_X}, \eqn{N_Y} : The total counts of each sample.
#' * \eqn{C_X} : The counts in \eqn{X} of the features shared with \eqn{Y}.
#' * \eqn{f_{1}^{Y}}, \eqn{f_{2}^{Y}} : The number of shared features with 1
#'   and with 2 counts in \eqn{Y}. When \eqn{f_{2}^{Y} = 0}, 1 is used.
#' * \eqn{S_X} : The counts in \eqn{X} of the shared features with 1 count
#'   in \eqn{Y}.
#'
#' It ranges from 0 to 1: 1 when the samples share no features (and `NaN`
#' when both are empty). It is not a metric.
#'
#' It equals vegan's `vegdist(x, "chao")`. It is the Jaccard-type index;
#' vegan's `chaodist()` offers the Sorensen type and others. Named
#' `chao_jaccard()`, not `chao()`, to keep it apart from the alpha metric
#' [chao1()]; `beta_div(x, 'chao')` finds it.
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,]
#' y <- ex_counts[2,]
#' u <- function (x, y) {
#'   sh <- x > 0 & y > 0
#'   f1 <- sum(sh & y == 1); f2 <- max(1, sum(sh & y == 2))
#'   min(1, sum(x[sh]) / sum(x) +
#'     (sum(y) - 1) / sum(y) * f1 / (2 * f2) * sum(x[sh & y == 1]) / sum(x))
#' }
#' U <- u(x, y); V <- u(y, x)
#' 1 - U * V / (U + V - U * V)
#' ```
#'
#' @references
#' Chao, A., Chazdon, R. L., Colwell, R. K., & Shen, T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. *Ecology Letters*, 8, 148-159. \doi{10.1111/j.1461-0248.2004.00707.x}
#'
#' @export
#' @examples
#'     chao_jaccard(ex_counts)
chao_jaccard <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  norm <- 'none'
  validate_args()
  vd_nonneg(counts)
  assert_integer_counts()
  dist_result(.Call(C_beta_div, BDIV_CHAO_JACCARD, counts, margin, norm, pairs, cpus, 0, NULL), 'chao_jaccard')
}


#' Mahalanobis distance
#'
#' The Euclidean distance after the features are centred, scaled to unit
#' variance, and decorrelated: distances that account for the features'
#' covariance.
#'
#' @inherit documentation
#' @inherit pseudocount_section
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Mahalanobis distance is defined as:
#' \deqn{\sqrt{(X - Y)^{T}\Sigma^{-1}(X - Y)}}
#'
#' Where:
#' * \eqn{X}, \eqn{Y} : The two samples' vectors of counts (after `norm`).
#' * \eqn{\Sigma} : The features' sample covariance matrix (with
#'   \eqn{n - 1}), over all the samples of `counts`.
#'
#' It is a metric. It depends on the whole table: removing samples changes
#' the covariance and so every distance. ecodive computes it without
#' inverting \eqn{\Sigma}: from a QR decomposition of the centred table
#' \eqn{X_c = QR}, the whitened table \eqn{\sqrt{n - 1}\,X_c R^{-1}}, whose
#' Euclidean distances are the Mahalanobis distances.
#'
#' It needs more samples than features, and no feature that is constant or
#' a linear combination of others (a diagonal element of \eqn{R} below
#' \eqn{10^{-9}} times the largest column norm counts as zero). Otherwise
#' it stops with an error: reduce the features first, e.g. to the leading
#' axes of [pca()]. `norm = 'clr'`, `'percent'` and `'chisq'` make the
#' features collinear (each sample's values have a fixed sum), so they
#' always stop it.
#'
#' It equals vegan's `vegdist(x, "mahalanobis")`. Differences from vegan:
#' where \eqn{\Sigma} is singular, vegan returns nearly constant distances
#' with a warning, and ecodive an error. Named `mahalanobis_distance()`,
#' not `mahalanobis()`, to leave [stats::mahalanobis()] unmasked;
#' `beta_div(x, 'mahalanobis')` finds it. Unlike the other metrics, it
#' accepts negative values.
#'
#' Base R Equivalent:
#' ```r
#' x <- as.matrix(iris[, 1:4])
#' dist(scale(x, scale = FALSE) %*% solve(chol(cov(x))))
#' ```
#'
#' @references
#' Mahalanobis, P. C. (1936). On the generalized distance in statistics. *Proceedings of the National Institute of Sciences of India*, 2, 49-55.
#'
#' Mardia, K. V., Kent, J. T., & Bibby, J. M. (1979). *Multivariate Analysis*. Academic Press.
#'
#' @export
#' @examples
#'     # It needs more samples than features: 150 flowers, 4 measurements.
#'     dm <- mahalanobis_distance(iris[, 1:4])
#'     dm
mahalanobis_distance <- function (counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  validate_args()
  dist_result(.Call(C_beta_div, BDIV_MAHALANOBIS, counts, margin, norm, pairs, cpus, pseudocount, NULL), 'mahalanobis_distance')
}


#' Chi-square distance
#'
#' The Euclidean distance after the chi-square transformation: the distance
#' that correspondence analysis ([ca()]) preserves.
#'
#' @inherit documentation
#'
#' @family Abundance metrics
#' @seealso `beta_div()`, `normalize()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The chi-square distance is defined as:
#' \deqn{\sqrt{T\sum_{i=1}^{n}\frac{1}{c_i}\left(\frac{X_i}{X_T} - \frac{Y_i}{Y_T}\right)^2}}
#'
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{X_T}, \eqn{Y_T} : Total counts in each sample.
#' * \eqn{c_i} : The total count of the \eqn{i}-th feature over all the
#'   samples of `counts` (features with \eqn{c_i = 0} are left out).
#' * \eqn{T} : The total count of the table.
#' * \eqn{n} : The number of features.
#'
#' It is a metric. It depends on the whole table, through \eqn{c_i} and
#' \eqn{T}. It is `euclidean(counts, norm = 'chisq')` (see [normalize()]).
#' It is not [squared_chisq()] or [psym_chisq()], which compare each pair
#' of samples on its own.
#'
#' It equals vegan's `vegdist(x, "chisq")`, which is Euclidean distance
#' after `decostand(x, "chi.square")` (Legendre & Gallagher 2001).
#' Differences from vegan: an empty sample is all zeros after the
#' transformation, so its distances are finite (vegan gives `NA`).
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts
#' x <- x[, colSums(x) > 0]
#' y <- sqrt(sum(x)) * sweep(x / rowSums(x), 2, sqrt(colSums(x)), '/')
#' dist(y)
#' ```
#'
#' @references
#' Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful transformations for ordination of species data. *Oecologia*, 129, 271-280. \doi{10.1007/s004420100716}
#'
#' @export
#' @examples
#'     chisq_distance(ex_counts)
chisq_distance <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  norm <- 'chisq'
  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, norm, pairs, cpus, 0, NULL), 'chisq_distance')
}


#' Mountford dissimilarity
#'
#' Mountford's (1962) index of similarity, from Fisher's log-series, as a
#' dissimilarity.
#'
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#'
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Mountford dissimilarity is \eqn{\ln 2 - M}, where \eqn{M} is the
#' positive root of
#' \deqn{e^{AM} + e^{BM} = 1 + e^{(A + B - J)M}}
#'
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#'
#' It ranges from 0 to \eqn{\ln 2}: \eqn{\ln 2} when the samples share no
#' features (\eqn{J = 0}), and 0 when one sample's features are all in the
#' other (\eqn{J = A} or \eqn{J = B}), so not only for identical samples.
#' It is `NaN` when either sample is empty. It is not a metric.
#'
#' \eqn{M} is solved by Newton's method, safeguarded by bisection, to full
#' double precision, in a form that can't overflow. It equals vegan's
#' `vegdist(x, "mountford")` to about \eqn{10^{-10}} while \eqn{A} and
#' \eqn{B} are at most about 50. Differences from vegan: for richer
#' samples, vegan's iterations can stop far from the root (by up to 0.65
#' at a few hundred features), and ecodive's value is the root.
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,] > 0
#' y <- ex_counts[2,] > 0
#' A <- sum(x); B <- sum(y); J <- sum(x & y)
#' g <- function (M) expm1(-(B - J) * M) + expm1(-(A - J) * M) - expm1(-(A + B - J) * M)
#' log(2) - uniroot(g, c(1e-9, 10), tol = 1e-12)$root  # for 0 < J < min(A, B)
#' ```
#'
#' @references
#' Mountford, M. D. (1962). An index of similarity and its application to classification problems. In P. W. Murphy (Ed.), *Progress in Soil Zoology* (pp. 43-50). Butterworths.
#'
#' @export
#' @examples
#'     mountford(ex_counts)
mountford <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  norm <- 'none'
  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_MOUNTFORD, counts, margin, norm, pairs, cpus, 0, NULL), 'mountford')
}


#' Raup-Crick dissimilarity
#'
#' The probability, under random placement of features, of the two samples
#' sharing at least as many features as they do (Raup & Crick 1979):
#' the analytic, presence/absence form of vegan's `vegdist()`.
#'
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#'
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#'
#' @details
#' The Raup-Crick dissimilarity is:
#' \deqn{P(X \ge J), \quad X \sim \mathrm{Hypergeometric}(S, A, B)}
#'
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#' * \eqn{S} : The number of features (columns) of `counts`, **including
#'   features absent from every sample**: adding an empty feature changes
#'   every value.
#' * \eqn{X} : The number of features two random samples of \eqn{A} and
#'   \eqn{B} of the \eqn{S} features have in common, with all features
#'   equally likely.
#'
#' It ranges from 0 to 1, and is not a metric. Samples that share no
#' features have a dissimilarity of exactly 1, while identical samples can
#' be slightly above 0. The smaller tail of the hypergeometric distribution
#' is summed directly, so small values keep their precision.
#'
#' It equals vegan's `vegdist(x, "raup")`, which computes \eqn{1 - P(X <
#' J)} and so returns 0 where the true value is below about
#' \eqn{10^{-16}}.
#'
#' This is not vegan's `raupcrick()` (Chase et al. 2011), which simulates
#' null communities with features drawn in proportion to their
#' frequencies, nor Stegen et al.'s abundance-based Raup-Crick
#' (\eqn{RC_{bray}}).
#'
#' Base R Equivalent:
#' ```r
#' x <- ex_counts[1,] > 0
#' y <- ex_counts[2,] > 0
#' S <- length(x); A <- sum(x); B <- sum(y); J <- sum(x & y)
#' phyper(J - 1, A, S - A, B, lower.tail = FALSE)
#' ```
#'
#' @references
#' Raup, D. M., & Crick, R. E. (1979). Measurement of faunal similarity in paleontology. *Journal of Paleontology*, 53, 1213-1227.
#'
#' Chase, J. M., Kraft, N. J. B., Smith, K. G., Vellend, M., & Inouye, B. D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. *Ecosphere*, 2, art24. \doi{10.1890/ES10-00117.1}
#'
#' @export
#' @examples
#'     raup_crick(ex_counts)
raup_crick <- function (counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus()) {

  norm <- 'none'
  validate_args()
  vd_nonneg(counts)
  dist_result(.Call(C_beta_div, BDIV_RAUP_CRICK, counts, margin, norm, pairs, cpus, 0, NULL), 'raup_crick')
}
