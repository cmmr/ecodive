# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Integer IDs for C code.
BDIV_BHATTACHARYYA <-  1L
BDIV_BRAY          <-  2L
BDIV_CANBERRA      <-  3L
BDIV_CHEBYSHEV     <-  4L
BDIV_CLARK         <-  6L
BDIV_DIVERGENCE    <-  7L
BDIV_EUCLIDEAN     <-  8L
BDIV_GOWER         <-  9L
BDIV_HAMMING       <- 10L
BDIV_HORN          <- 11L
BDIV_JACCARD       <- 12L
BDIV_JSD           <- 13L
BDIV_LORENTZIAN    <- 14L
BDIV_MANHATTAN     <- 15L
BDIV_MINKOWSKI     <- 16L
BDIV_MORISITA      <- 17L
BDIV_MOTYKA        <- 18L
BDIV_OCHIAI        <- 19L
BDIV_SOERGEL       <- 20L
BDIV_SORENSEN      <- 21L
BDIV_SQUARED_CHISQ <- 22L
BDIV_SQUARED_CHORD <- 23L
BDIV_WAVE_HEDGES   <- 24L

U_UNIFRAC <- 1L
W_UNIFRAC <- 2L
N_UNIFRAC <- 3L
G_UNIFRAC <- 4L
V_UNIFRAC <- 5L



#' Beta Diversity Wrapper Function
#' 
#' @inherit documentation
#' 
#' @name beta_div
#'        
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features). 
#'        Typically contains absolute abundances (integer counts), though 
#'        proportions are also accepted by some diversity metrics.
#' 
#' @param metric   The name of a beta diversity metric. One of `c('aitchison',
#'   'bhattacharyya', 'bray', 'canberra', 'chebyshev', 'chord', 'clark',
#'   'divergence', 'euclidean', 'generalized_unifrac', 'gower', 'hamming',
#'   'hellinger', 'horn', 'jaccard', 'jensen', 'jsd', 'lorentzian', 'manhattan',
#'   'matusita', 'minkowski', 'morisita', 'motyka', 'normalized_unifrac',
#'   'ochiai', 'psym_chisq', 'soergel', 'sorensen', 'squared_chisq',
#'   'squared_chord', 'squared_euclidean', 'topsoe', 'unweighted_unifrac',
#'   'variance_adjusted_unifrac', 'wave_hedges', 'weighted_unifrac')`. Flexible
#'   matching is supported (see below). Programmatic access via
#'   `list_metrics('beta')`.
#' 
#' @return A numeric vector.
#' 
#' 
#' @details
#' 
#' **List of Beta Diversity Metrics**
#' 
#' | Option / Function Name      | Metric Name                                      |
#' | :-------------------------- | :----------------------------------------------- |
#' | `aitchison`                 | Aitchison distance                               |
#' | `bhattacharyya`             | Bhattacharyya distance                           |
#' | `bray`                      | Bray-Curtis dissimilarity                        |
#' | `canberra`                  | Canberra distance                                |
#' | `chebyshev`                 | Chebyshev distance                               |
#' | `chord`                     | Chord distance                                   |
#' | `clark`                     | Clark's divergence distance                      |
#' | `divergence`                | Divergence                                       |
#' | `euclidean`                 | Euclidean distance                               |
#' | `generalized_unifrac`       | Generalized UniFrac (GUniFrac)                   |
#' | `gower`                     | Gower distance                                   |
#' | `hamming`                   | Hamming distance                                 |
#' | `hellinger`                 | Hellinger distance                               |
#' | `horn`                      | Horn-Morisita dissimilarity                      |
#' | `jaccard`                   | Jaccard distance                                 |
#' | `jensen`                    | Jensen-Shannon distance                          |
#' | `jsd`                       | Jesen-Shannon divergence (JSD)                   |
#' | `lorentzian`                | Lorentzian distance                              |
#' | `manhattan`                 | Manhattan distance                               |
#' | `matusita`                  | Matusita distance                                |
#' | `minkowski`                 | Minkowski distance                               |
#' | `morisita`                  | Morisita dissimilarity                           |
#' | `motyka`                    | Motyka dissimilarity                             |
#' | `normalized_unifrac`        | Normalized Weighted UniFrac                      |
#' | `ochiai`                    | Otsuka-Ochiai dissimilarity                      |
#' | `psym_chisq`                | Probabilistic Symmetric Chi-Squared distance     |
#' | `soergel`                   | Soergel distance                                 |
#' | `sorensen`                  | Dice-Sorensen dissimilarity                      |
#' | `squared_chisq`             | Squared Chi-Squared distance                     |
#' | `squared_chord`             | Squared Chord distance                           |
#' | `squared_euclidean`         | Squared Euclidean distance                       |
#' | `topsoe`                    | Topsoe distance                                  |
#' | `unweighted_unifrac`        | Unweighted UniFrac                               |
#' | `variance_adjusted_unifrac` | Variance-Adjusted Weighted UniFrac (VAW-UniFrac) |
#' | `wave_hedges`               | Wave Hedges distance                             |
#' | `weighted_unifrac`          | Weighted UniFrac                                 |
#' 
#' 
#' 
#' **Flexible name matching**
#' 
#' Case insensitive and partial matching. Any runs of non-alpha characters are
#' converted to underscores. E.g. `metric = 'Weighted UniFrac` selects
#' `weighted_unifrac`.
#' 
#' UniFrac names can be shortened to the first letter plus "unifrac". E.g. 
#' `uunifrac`, `w_unifrac`, or `V UniFrac`. These also support partial matching.
#' 
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the metrics list.
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Bray-Curtis distances
#'     beta_div(ex_counts, 'bray')
#'     
#'     # Generalized UniFrac distances
#'     beta_div(ex_counts, 'GUniFrac', tree = ex_tree)
#'     
beta_div <- function (
    counts, 
    metric, 
    norm        = 'none', 
    power       = 1.5, 
    pseudocount = NULL, 
    alpha       = 0.5, 
    tree        = NULL, 
    pairs       = NULL, 
    margin      = 1L, 
    cpus        = n_cpus() ) {
  
  metric <- match_metric(metric, div = 'beta')
  args   <- mget(metric$params, environment())
  
  do.call(metric$func, args)
}


#' Aitchison distance
#' 
#' Calculates the Euclidean distance between centered log-ratio (CLR) transformed abundances.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Aitchison distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n} [(\ln{X_i} - X_L) - (\ln{Y_i} - Y_L)]^2}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts for the \eqn{i}-th feature.
#' * \eqn{X_L}, \eqn{Y_L} : Mean log of abundances. \eqn{X_L = \frac{1}{n}\sum_{i=1}^{n} \ln{X_i}}.
#' * \eqn{n} : The number of features.
#' 
#' **Parameter: pseudocount**
#' 
#' Because the formula uses logarithms, zeros in the data must be handled. 
#' The `pseudocount` argument adds a small value to all counts prior to calculation.
#' 
#' Base R Equivalent: 
#' ```r
#' x <- log((x + pseudocount) / exp(mean(log(x + pseudocount))))
#' y <- log((y + pseudocount) / exp(mean(log(y + pseudocount))))
#' sqrt(sum((x-y)^2)) # Euclidean distance
#' ```
#' 
#' @references
#' Aitchison, J. (1986). The statistical analysis of compositional data. Chapman and Hall. \doi{10.1007/978-94-009-4109-3}
#' 
#' Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society: Series B (Methodological)*, 44(2), 139-160. \doi{10.1111/j.2517-6161.1982.tb01195.x}
#'
#' @export
#' @examples
#'     aitchison(ex_counts, pseudocount = 1)
aitchison <- function (counts, pseudocount = NULL, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- structure('clr', pseudocount = pseudocount)
  validate_args()
  
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, norm, pairs, cpus, NULL)
}


#' Bhattacharyya distance
#' 
#' Measures the similarity of two probability distributions.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Bhattacharyya distance is defined as:
#' \deqn{-\ln{\sum_{i=1}^{n}\sqrt{P_{i}Q_{i}}}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' -log(sum(sqrt(p * q)))
#' ```
#' 
#' @references
#' Bhattacharyya, A. (1943). On a measure of divergence between two statistical populations defined by their probability distributions. *Bulletin of the Calcutta Mathematical Society*, 35, 99-109.
#' 
#' @export
#' @examples
#'     bhattacharyya(ex_counts)
bhattacharyya <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_beta_div, BDIV_BHATTACHARYYA, counts, margin, norm, pairs, cpus, NULL)
}


#' Bray-Curtis dissimilarity
#' 
#' A standard ecological metric quantifying the dissimilarity between communities.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Bray-Curtis dissimilarity is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} |X_i - Y_i|}{\sum_{i=1}^{n} (X_i + Y_i)}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(abs(x-y)) / sum(x+y)
#' ```
#' 
#' @references
#' Bray, J. R., & Curtis, J. T. (1957). An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*, 27(4), 325-349. \doi{10.2307/1942268}
#' 
#' @export
#' @examples
#'     bray(ex_counts)
bray <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_BRAY, counts, margin, norm, pairs, cpus, NULL)
}


#' Canberra distance
#' 
#' A weighted version of the Manhattan distance, sensitive to differences when both values are small.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Canberra distance is defined as:
#' \deqn{\sum_{i=1}^{n} \frac{|X_i - Y_i|}{X_i + Y_i}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(abs(x-y) / (x+y))
#' ```
#' 
#' @references
#' Lance, G. N., & Williams, W. T. (1966). Computer programs for hierarchical polythetic classification ("similarity analyses"). *The Computer Journal*, 9(1), 60-64. \doi{10.1093/comjnl/9.1.60}
#' 
#' @export
#' @examples
#'     canberra(ex_counts)
canberra <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CANBERRA, counts, margin, norm, pairs, cpus, NULL)
}


#' Chebyshev distance
#' 
#' The maximum difference between any single feature across two samples.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Chebyshev distance is defined as:
#' \deqn{\max(|X_i - Y_i|)}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' 
#' Base R Equivalent: 
#' ```r
#' max(abs(x-y))
#' ```
#' 
#' @references
#' Cantrell, C. D. (2000). Modern mathematical methods for physicists and engineers. Cambridge University Press.
#' 
#' @export
#' @examples
#'     chebyshev(ex_counts)
chebyshev <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CHEBYSHEV, counts, margin, norm, pairs, cpus, NULL)
}


#' Chord distance
#' 
#' Euclidean distance between normalized vectors.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Chord distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n} \left(\frac{X_i}{\sqrt{\sum_{j=1}^{n} X_j^2}} - \frac{Y_i}{\sqrt{\sum_{j=1}^{n} Y_j^2}}\right)^2}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sqrt(sum(((x / sqrt(sum(x ^ 2))) - (y / sqrt(sum(y ^ 2))))^2))
#' ```
#' 
#' @references
#' Orlóci, L. (1967). An agglomerative method for classification of plant communities. *Journal of Ecology*, 55(1), 193-206. \doi{10.2307/2257725}
#' 
#' @export
#' @examples
#'     chord(ex_counts)
chord <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'chord'
  validate_args()
  
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, norm, pairs, cpus, NULL)
}


#' Clark's divergence distance
#' 
#' Also known as the coefficient of divergence.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' Clark's divergence distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n}\left(\frac{X_i - Y_i}{X_i + Y_i}\right)^{2}}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sqrt(sum((abs(x - y) / (x + y)) ^ 2))
#' ```
#' 
#' @references
#' Clark, P. J. (1952). An extension of the coefficient of divergence for use with multiple characters. *Copeia*, 1952(2), 61-64. \doi{10.2307/1438598}
#' 
#' @export
#' @examples
#'     clark(ex_counts)
clark <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CLARK, counts, margin, norm, pairs, cpus, NULL)
}


#' Divergence
#' 
#' A probabilistic divergence metric.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' Divergence is defined as:
#' \deqn{2\sum_{i=1}^{n} \frac{(P_i - Q_i)^2}{(P_i + Q_i)^2}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' 2 * sum((p - q)^2 / (p + q)^2)
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     divergence(ex_counts)
divergence <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_beta_div, BDIV_DIVERGENCE, counts, margin, norm, pairs, cpus, NULL)
}


#' Euclidean distance
#' 
#' The straight-line distance between two points in multidimensional space.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Euclidean distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n} (X_i - Y_i)^2}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sqrt(sum((x-y)^2))
#' ```
#' 
#' @references
#' Legendre, P., & Legendre, L. (2012). Numerical ecology. Elsevier.
#' 
#' @export
#' @examples
#'     euclidean(ex_counts)
euclidean <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, norm, pairs, cpus, NULL)
}


#' Gower distance
#' 
#' A distance metric that normalizes differences by the range of the feature.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Gower distance is defined as:
#' \deqn{\frac{1}{n}\sum_{i=1}^{n}\frac{|X_i - Y_i|}{R_i}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{R_i} : The range of the \eqn{i}-th feature across all samples (max - min).
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' r <- abs(x - y)
#' n <- length(x) # <-- not `sum(x|y)`
#' sum(abs(x-y) / r) / n
#' ```
#' 
#' @references
#' Gower, J. C. (1971). A general coefficient of similarity and some of its properties. *Biometrics*, 27(4), 857-871. \doi{10.2307/2528823}
#' 
#' @export
#' @examples
#'     gower(ex_counts)
gower <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  # range_vec <- apply(counts, 2L, function (x) diff(range(x)))
  
  .Call(C_beta_div, BDIV_GOWER, counts, margin, norm, pairs, cpus, NULL)
}


#' Hellinger distance
#' 
#' A distance metric related to the Bhattacharyya distance, often used for community data with many zeros.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Hellinger distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n}(\sqrt{P_i} - \sqrt{Q_i})^{2}}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sqrt(sum((sqrt(p) - sqrt(q)) ^ 2))
#' ```
#' 
#' @references
#' Rao, C. R. (1995). A review of canonical coordinates and an alternative to correspondence analysis using Hellinger distance. *Qüestiió*, 19, 23-63.
#' 
#' Hellinger, E. (1909). Neue Begründung der Theorie quadratischer Formen von unendlichvielen Veränderlichen. *Journal für die reine und angewandte Mathematik*, 136, 210–271. \doi{10.1515/crll.1909.136.210}
#'
#' @export
#' @examples
#'     hellinger(ex_counts)
hellinger <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  sqc <- .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, margin, norm, pairs, cpus, NULL)
  
  sqrt(sqc)
}


#' Horn-Morisita dissimilarity
#' 
#' A similarity index based on Simpson's diversity index, suitable for abundance data.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Horn-Morisita dissimilarity is defined as:
#' \deqn{1 - \frac{2\sum_{i=1}^{n}P_{i}Q_{i}}{\sum_{i=1}^{n}P_i^2 + \sum_{i=1}^{n}Q_i^2}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
#' 1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))
#' ```
#' 
#' @references
#' Horn, H. S. (1966). Measurement of "overlap" in comparative ecological studies. *The American Naturalist*, 100(914), 419-424. \doi{10.1086/282436}
#' 
#' @export
#' @examples
#'     horn(ex_counts)
horn <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_HORN, counts, margin, norm, pairs, cpus, NULL)
}


#' Jensen-Shannon distance
#' 
#' The square root of the Jensen-Shannon divergence.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Jensen-Shannon distance is defined as:
#' \deqn{\sqrt{\frac{1}{2}\left[\sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)\right]}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sqrt(sum(p * log(2 * p / (p+q)), q * log(2 * q / (p+q))) / 2)
#' ```
#' 
#' @references
#' Endres, D. M., & Schindelin, J. E. (2003). A new metric for probability distributions. *IEEE Transactions on Information Theory*, 49(7), 1858-1860. \doi{10.1109/TIT.2003.813506}
#' 
#' @export
#' @examples
#'     jensen(ex_counts)
jensen <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  jsd <- .Call(C_beta_div, BDIV_JSD, counts, margin, norm, pairs, cpus, NULL)
  
  sqrt(jsd)
}


#' Jensen-Shannon divergence (JSD)
#' 
#' A symmetrized and smoothed version of the Kullback-Leibler divergence.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Jensen-Shannon divergence (JSD) is defined as:
#' \deqn{\frac{1}{2}\left[\sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)\right]}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sum(p * log(2 * p / (p+q)), q * log(2 * q / (p+q))) / 2
#' ```
#' 
#' @references
#' Lin, J. (1991). Divergence measures based on the Shannon entropy. *IEEE Transactions on Information Theory*, 37(1), 145-151. \doi{10.1109/18.61115}
#' 
#' @export
#' @examples
#'     jsd(ex_counts)
jsd <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_beta_div, BDIV_JSD, counts, margin, norm, pairs, cpus, NULL)
}


#' Lorentzian distance
#' 
#' A log-based distance metric that is robust to outliers.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Lorentzian distance is defined as:
#' \deqn{\sum_{i=1}^{n}\ln{(1 + |X_i - Y_i|)}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(log(1 + abs(x - y)))
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     lorentzian(ex_counts)
lorentzian <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_LORENTZIAN, counts, margin, norm, pairs, cpus, NULL)
}


#' Manhattan distance
#' 
#' The sum of absolute differences, also known as the taxicab geometry.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Manhattan distance is defined as:
#' \deqn{\sum_{i=1}^{n} |X_i - Y_i|}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(abs(x-y))
#' ```
#' 
#' @references
#' Krause, E. F. (1987). Taxicab geometry: An adventure in non-Euclidean geometry. Dover Publications.
#' 
#' @export
#' @examples
#'     manhattan(ex_counts)
manhattan <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MANHATTAN, counts, margin, norm, pairs, cpus, NULL)
}


#' Matusita distance
#' 
#' A distance measure closely related to the Hellinger distance.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Matusita distance is defined as:
#' \deqn{\sqrt{\sum_{i=1}^{n}\left(\sqrt{P_i} - \sqrt{Q_i}\right)^2}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sqrt(sum((sqrt(p) - sqrt(q)) ^ 2))
#' ```
#' 
#' @references
#' Matusita, K. (1955). Decision rules, based on the distance, for problems of fit, two samples, and estimation. *The Annals of Mathematical Statistics*, 26(4), 631-640. \doi{10.1214/aoms/1177728422}
#' 
#' @export
#' @examples
#'     matusita(ex_counts)
matusita <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  sqc <- .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, margin, norm, pairs, cpus, NULL)
  
  sqrt(sqc)
}


#' Minkowski distance
#' 
#' A generalized metric that includes Euclidean and Manhattan distance as special cases.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Minkowski distance is defined as:
#' \deqn{\sqrt[p]{\sum_{i=1}^{n} (X_i - Y_i)^p}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' * \eqn{p} : The geometry of the space (power parameter).
#' 
#' **Parameter: power**
#' 
#' The `power` parameter (default 1.5) determines the value of \eqn{p} in the equation.
#' 
#' **Special Cases**
#' 
#' * **Manhattan distance**: When \eqn{p = 1}, the formula reduces to the sum of absolute differences.
#' * **Euclidean distance**: When \eqn{p = 2}, the formula reduces to the standard straight-line distance.
#' * **Chebyshev distance**: When \eqn{p \to \infty}, the formula reduces to the maximum absolute difference.
#'
#' Base R Equivalent: 
#' ```r
#' p <- 1.5
#' sum(abs(x - y)^p) ^ (1/p)
#' ```
#' 
#' @references
#' Deza, M. M., & Deza, E. (2009). Encyclopedia of distances. Springer.
#' 
#' Minkowski, H. (1896). *Geometrie der Zahlen*. Teubner.
#'
#' @export
#' @examples
#'     minkowski(ex_counts, power = 2) # Equivalent to Euclidean
minkowski <- function (counts, norm = 'none', power = 1.5, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MINKOWSKI, counts, margin, norm, pairs, cpus, power)
}


#' Morisita dissimilarity
#' 
#' A measure of overlap between samples that is independent of sample size. Requires integer counts.
#' 
#' @inherit documentation
#' @inherit bdiv_assert_integer
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Morisita dissimilarity is defined as:
#' \deqn{1 - \frac{2\sum_{i=1}^{n}X_{i}Y_{i}}{\left(\frac{\sum_{i=1}^{n}X_i(X_i - 1)}{X_T(X_T - 1)} + \frac{\sum_{i=1}^{n}Y_i(Y_i - 1)}{Y_T(Y_T - 1)}\right)X_{T}Y_{T}}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts of the \eqn{i}-th feature.
#' * \eqn{X_T}, \eqn{Y_T} : Total counts in each sample. \eqn{X_T = \sum_{i=1}^{n} X_i}.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' simp_x <- sum(x * (x - 1)) / (sum(x) * (sum(x) - 1))
#' simp_y <- sum(y * (y - 1)) / (sum(y) * (sum(y) - 1))
#' 1 - ((2 * sum(x * y)) / ((simp_x + simp_y) * sum(x) * sum(y)))
#' ```
#' 
#' @references
#' Morisita, M. (1959). Measuring of interspecific association and similarity between communities. *Memoirs of the Faculty of Science, Kyushu University, Series E (Biology)*, 3, 65-80.
#' 
#' @export
#' @examples
#'     morisita(ex_counts)
morisita <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  assert_integer_counts()
  
  .Call(C_beta_div, BDIV_MORISITA, counts, margin, norm, pairs, cpus, NULL)
}


#' Motyka dissimilarity
#' 
#' Also known as the Bray-Curtis dissimilarity when applied to abundance data, but formulated slightly differently.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Motyka dissimilarity is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} \max(X_i, Y_i)}{\sum_{i=1}^{n} (X_i + Y_i)}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(pmax(x, y)) / sum(x, y)
#' ```
#' 
#' @references
#' Motyka, J. (1947). O celach i metodach badan geobotanicznych. *Annales Universitatis Mariae Curie-Sklodowska, Sectio C*, 3, 1-168.
#' 
#' @export
#' @examples
#'     motyka(ex_counts)
motyka <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MOTYKA, counts, margin, norm, pairs, cpus, NULL)
}


#' Probabilistic Symmetric Chi-Squared distance
#' 
#' A chi-squared based distance metric for comparing probability distributions.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Probabilistic Symmetric \eqn{\chi^2} distance is defined as:
#' \deqn{2\sum_{i=1}^{n}\frac{(P_i - Q_i)^2}{P_i + Q_i}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' 2 * sum((p - q)^2 / (p + q))
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     psym_chisq(ex_counts)
psym_chisq <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  scs <- .Call(C_beta_div, BDIV_SQUARED_CHISQ, counts, margin, norm, pairs, cpus, NULL)
  
  2 * scs
}


#' Soergel distance
#' 
#' A distance metric related to the Bray-Curtis and Jaccard indices.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Soergel distance is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} |X_i - Y_i|}{\sum_{i=1}^{n} \max(X_i, Y_i)}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(abs(x - y)) / sum(pmax(x, y))
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     soergel(ex_counts)
soergel <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_SOERGEL, counts, margin, norm, pairs, cpus, NULL)
}


#' Squared Chi-Squared distance
#' 
#' The squared version of the Chi-Squared distance.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Squared \eqn{\chi^2} distance is defined as:
#' \deqn{\sum_{i=1}^{n}\frac{(P_i - Q_i)^2}{P_i + Q_i}}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sum((p - q)^2 / (p + q))
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     squared_chisq(ex_counts)
squared_chisq <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_beta_div, BDIV_SQUARED_CHISQ, counts, margin, norm, pairs, cpus, NULL)
}


#' Squared Chord distance
#' 
#' The squared version of the Chord distance.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Squared Chord distance is defined as:
#' \deqn{\sum_{i=1}^{n}\left(\sqrt{P_i} - \sqrt{Q_i}\right)^2}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sum((sqrt(x) - sqrt(y)) ^ 2)
#' ```
#' 
#' @references
#' Legendre, P., & Legendre, L. (2012). Numerical ecology. Elsevier.
#' 
#' @export
#' @examples
#'     squared_chord(ex_counts)
squared_chord <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, margin, norm, pairs, cpus, NULL)
}


#' Squared Euclidean distance
#' 
#' The squared Euclidean distance between two vectors.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Squared Euclidean distance is defined as:
#' \deqn{\sum_{i=1}^{n} (X_i - Y_i)^2}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum((x-y)^2)
#' ```
#' 
#' @references
#' Legendre, P., & Legendre, L. (2012). Numerical ecology. Elsevier.
#' 
#' @export
#' @examples
#'     squared_euclidean(ex_counts)
squared_euclidean <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  euc <- .Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, norm, pairs, cpus, NULL)
  
  euc ^ 2
}


#' Topsoe distance
#' 
#' A symmetric divergence measure based on the Jensen-Shannon divergence.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Topsoe distance is defined as:
#' \deqn{\sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)}
#' 
#' Where:
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' p <- x / sum(x)
#' q <- y / sum(y)
#' sum(p * log(2 * p / (p+q)), q * log(2 * y / (p+q)))
#' ```
#' 
#' @references
#' Topsoe, F. (2000). Some inequalities for information divergence and related measures of discrimination. *IEEE Transactions on Information Theory*, 46(4), 1602-1609. \doi{10.1109/18.850703}
#' 
#' @export
#' @examples
#'     topsoe(ex_counts)
topsoe <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  jsd <- .Call(C_beta_div, BDIV_JSD, counts, margin, norm, pairs, cpus, NULL)
  
  2 * jsd
}


#' Wave Hedges distance
#' 
#' A distance metric derived from the Hedges' distance.
#' 
#' @inherit documentation
#' 
#' @family Abundance metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Wave Hedges distance is defined as:
#' \deqn{\sum_{i=1}^{n}\frac{|X_i - Y_i|}{\max(X_i, Y_i)}}
#' 
#' Where:
#' * \eqn{X_i}, \eqn{Y_i} : Absolute abundances of the \eqn{i}-th feature.
#' * \eqn{n} : The number of features.
#' 
#' Base R Equivalent: 
#' ```r
#' sum(abs(x - y) / pmax(x, y))
#' ```
#' 
#' @references
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures between probability density functions. *International Journal of Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' @export
#' @examples
#'     wave_hedges(ex_counts)
wave_hedges <- function (counts, norm = 'none', margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_WAVE_HEDGES, counts, margin, norm, pairs, cpus, NULL)
}




# Presence/Absence -----------


#' Hamming distance
#' 
#' Measures the minimum number of substitutions required to change one string into the other.
#' 
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#' 
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Hamming distance is defined as:
#' \deqn{(A + B) - 2J}
#' 
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#' 
#' Base R Equivalent: 
#' ```r
#' sum(xor(x, y))
#' ```
#' 
#' @references
#' Hamming, R. W. (1950). Error detecting and error correcting codes. *Bell System Technical Journal*, 29(2), 147-160. \doi{10.1002/j.1538-7305.1950.tb00463.x}
#' 
#' @export
#' @examples
#'     hamming(ex_counts)
hamming <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_beta_div, BDIV_HAMMING, counts, margin, norm, pairs, cpus, NULL)
}


#' Jaccard distance
#' 
#' Measures dissimilarity between sample sets.
#' 
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#' 
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Jaccard distance is defined as:
#' \deqn{1 - \frac{J}{(A + B - J)}}
#' 
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#' 
#' Base R Equivalent: 
#' ```r
#' 1 - sum(x & y) / sum(x | y)
#' ```
#' 
#' @references
#' Jaccard, P. (1912). The distribution of the flora in the alpine zone. *New Phytologist*, 11(2), 37-50. \doi{10.1111/j.1469-8137.1912.tb05611.x}
#' 
#' Jaccard, P. (1908). Nouvelles recherches sur la distribution florale. *Bulletin de la Societe Vaudoise des Sciences Naturelles*, 44(163), 223-270. \doi{10.5169/seals-268384}
#'
#' @export
#' @examples
#'     jaccard(ex_counts)
jaccard <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_beta_div, BDIV_JACCARD, counts, margin, norm, pairs, cpus, NULL)
}


#' Otsuka-Ochiai dissimilarity
#' 
#' Also known as the cosine similarity for binary data.
#' 
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#' 
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Otsuka-Ochiai dissimilarity is defined as:
#' \deqn{1 - \frac{J}{\sqrt{AB}}}
#' 
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#' 
#' Base R Equivalent: 
#' ```r
#' 1 - sum(x & y) / sqrt(sum(x>0) * sum(y>0)) 
#' ```
#' 
#' @references
#' Ochiai, A. (1957). Zoogeographic studies on the soleoid fishes found in Japan and its neighbouring regions. *Bulletin of the Japanese Society of Scientific Fisheries*, 22, 526-530.
#' 
#' @export
#' @examples
#'     ochiai(ex_counts)
ochiai <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_beta_div, BDIV_OCHIAI, counts, margin, norm, pairs, cpus, NULL)
}


#' Dice-Sorensen dissimilarity
#' 
#' A statistic used for comparing the similarity of two samples.
#' 
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#' 
#' @family Presence/Absence metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Dice-Sorensen dissimilarity is defined as:
#' \deqn{\frac{2J}{(A + B)}}
#' 
#' Where:
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common (intersection).
#' 
#' Base R Equivalent: 
#' ```r
#' 2 * sum(x & y) / sum(x>0, y>0)
#' ```
#' 
#' @references
#' Sørensen, T. (1948). A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. *Kongelige Danske Videnskabernes Selskab, Biologiske Skrifter*, 5, 1-34.
#' 
#' Dice, L. R. (1945). Measures of the amount of ecologic association between species. *Ecology*, 26(3), 297–302. \doi{10.2307/1932409}
#'
#' @export
#' @examples
#'     sorensen(ex_counts)
sorensen <- function (counts, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_beta_div, BDIV_SORENSEN, counts, margin, norm, pairs, cpus, NULL)
}




# UniFrac Family -----------

#' Unweighted UniFrac
#' 
#' A phylogenetic distance metric that accounts for the presence/absence of lineages.
#' 
#' @inherit documentation
#' @inherit bdiv_binary_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Unweighted UniFrac distance is defined as:
#' \deqn{\frac{1}{n}\sum_{i=1}^{n} L_i|A_i - B_i|}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{A_i}, \eqn{B_i} : Binary values (0 or 1) indicating if descendants of branch \eqn{i} are present in sample A or B.
#' 
#' @references
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology*, 71(12), 8228-8235. \doi{10.1128/AEM.71.12.8228-8235.2005}
#' 
#' @export
#' @examples
#'     unweighted_unifrac(ex_counts, tree = ex_tree)
unweighted_unifrac <- function (counts, tree = NULL, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  .Call(C_unifrac, U_UNIFRAC, counts, tree, margin, pairs, cpus, NULL)
}


#' Weighted UniFrac
#' 
#' A phylogenetic distance metric that accounts for the relative abundance of lineages.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Weighted UniFrac distance is defined as:
#' \deqn{\sum_{i=1}^{n} L_i|P_i - Q_i|}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{P_i}, \eqn{Q_i} : The proportion of the community descending from branch \eqn{i} in sample P and Q.
#' 
#' @references
#' Lozupone, C. A., Hamady, M., Kelley, S. T., & Knight, R. (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. *Applied and Environmental Microbiology*, 73(5), 1576-1585. \doi{10.1128/AEM.01996-06}
#' 
#' @export
#' @examples
#'     weighted_unifrac(ex_counts, tree = ex_tree)
weighted_unifrac <- function (counts, tree = NULL, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  .Call(C_unifrac, W_UNIFRAC, counts, tree, margin, pairs, cpus, NULL)
}


#' Normalized Weighted UniFrac
#' 
#' Weighted UniFrac normalized by the tree length to allow comparison between trees.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Normalized Weighted UniFrac distance is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} L_i|P_i - Q_i|}{\sum_{i=1}^{n} L_i(P_i + Q_i)}}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{P_i}, \eqn{Q_i} : The proportion of the community descending from branch \eqn{i} in sample P and Q.
#' 
#' @references
#' Lozupone, C. A., Hamady, M., Kelley, S. T., & Knight, R. (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. *Applied and Environmental Microbiology*, 73(5), 1576-1585. \doi{10.1128/AEM.01996-06}
#' 
#' @export
#' @examples
#'     normalized_unifrac(ex_counts, tree = ex_tree)
normalized_unifrac <- function (counts, tree = NULL, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  .Call(C_unifrac, N_UNIFRAC, counts, tree, margin, pairs, cpus, NULL)
} 


#' Generalized UniFrac (GUniFrac)
#' 
#' A unified UniFrac distance that balances the weight of abundant and rare lineages.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Generalized UniFrac distance is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}\left|\frac{P_i - Q_i}{P_i + Q_i}\right|}{\sum_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}}}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{P_i}, \eqn{Q_i} : The proportion of the community descending from branch \eqn{i} in sample P and Q.
#' * \eqn{\alpha} : A scalable weighting factor.
#' 
#' **Parameter: alpha**
#' 
#' The `alpha` parameter controls the weight given to abundant lineages. \eqn{\alpha = 1} corresponds to Weighted UniFrac, while \eqn{\alpha = 0} corresponds to Unweighted UniFrac.
#' 
#' @references
#' Chen, J., Bittinger, K., Charlson, E. S., Hoffmann, C., Lewis, J., Wu, G. D., ... & Li, H. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. *Bioinformatics*, 28(16), 2106-2113. \doi{10.1093/bioinformatics/bts342}
#' 
#' @export
#' @examples
#'     generalized_unifrac(ex_counts, tree = ex_tree, alpha = 0.5)
generalized_unifrac <- function (counts, tree = NULL, alpha = 0.5, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  .Call(C_unifrac, G_UNIFRAC, counts, tree, margin, pairs, cpus, alpha)
}


#' Variance-Adjusted Weighted UniFrac
#' 
#' A weighted UniFrac that adjusts for the expected variance of the metric.
#' 
#' @inherit documentation
#' @inherit bdiv_percent_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `beta_div()`, `vignette('bdiv')`, `vignette('bdiv_guide')`
#' 
#' @details
#' The Variance-Adjusted Weighted UniFrac distance is defined as:
#' \deqn{\frac{\sum_{i=1}^{n} L_i\frac{|P_i - Q_i|}{\sqrt{(P_i + Q_i)(2 - P_i - Q_i)}} }{\sum_{i=1}^{n} L_i\frac{P_i + Q_i}{\sqrt{(P_i + Q_i)(2 - P_i - Q_i)}} }}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{P_i}, \eqn{Q_i} : The proportion of the community descending from branch \eqn{i} in sample P and Q.
#' 
#' @references
#' Chang, Q., Luan, Y., & Sun, F. (2011). Variance adjusted weighted UniFrac: a powerful beta diversity measure for comparing communities based on phylogeny. *BMC Bioinformatics*, 12, 118. \doi{10.1186/1471-2105-12-118}
#' 
#' @export
#' @examples
#'     variance_adjusted_unifrac(ex_counts, tree = ex_tree)
variance_adjusted_unifrac <- function (counts, tree = NULL, margin = 1L, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  .Call(C_unifrac, V_UNIFRAC, counts, tree, margin, pairs, cpus, NULL)
}
