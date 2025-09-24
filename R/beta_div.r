# Copyright (c) 2025 ecodive authors
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
#' @name beta_div
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
#' @param ...  Additional options to pass through to the called function. I.e.
#'   `tree`, `pairs`, `alpha`, or `cpus`. Named options will only be passed 
#'   though if the function expects them.
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
    rescale     = TRUE, 
    power       = 1.5, 
    pseudocount = NULL, 
    alpha       = 0.5, 
    tree        = NULL, 
    pairs       = NULL, 
    cpus        = n_cpus() ) {
  
  metric <- match_metric(metric, div = 'beta')
  args   <- mget(metric$params, environment())
  
  do.call(metric$func, args)
}





#' Beta Diversity Metrics
#' 
#' @inherit documentation
#' @name bdiv_functions
#' @family bdiv_functions
#' 
#' @return A `dist` object.
#' 
#' 
#' @section Formulas:
#' 
#' Given:
#' 
#' * \eqn{n} : The number of features.
#' * \eqn{X_i}, \eqn{Y_i} : Absolute counts for the \eqn{i}-th feature in samples \eqn{X} and \eqn{Y}.
#' * \eqn{X_T}, \eqn{Y_T} : Total counts in each sample. \eqn{X_T = \sum_{i=1}^{n} X_i}
#' * \eqn{P_i}, \eqn{Q_i} : Proportional abundances of \eqn{X_i} and \eqn{Y_i}. \eqn{P_i = X_i / X_T}
#' * \eqn{X_L}, \eqn{Y_L} : Mean log of abundances. \eqn{X_L = \frac{1}{n}\sum_{i=1}^{n} \ln{X_i}}
#' * \eqn{R_i} : The range of the \eqn{i}-th feature across all samples (max - min).
#' 
#' |              |                                    |
#' | :----------- | :--------------------------------- |
#' | **Aitchison distance**                            <br> `aitchison()`         | \eqn{\sqrt{\sum_{i=1}^{n} [(\ln{X_i} - X_L) - (\ln{Y_i} - Y_L)]^2}} |
#' | **Bhattacharyya distance**                        <br> `bhattacharyya()`     | \eqn{-\ln{\sum_{i=1}^{n}\sqrt{P_{i}Q_{i}}}} |
#' | **Bray-Curtis dissimilarity**                     <br> `bray()`              | \eqn{\displaystyle \frac{\sum_{i=1}^{n} |P_i - Q_i|}{\sum_{i=1}^{n} (P_i + Q_i)}} |
#' | **Canberra distance**                             <br> `canberra()`          | \eqn{\displaystyle \sum_{i=1}^{n} \frac{|P_i - Q_i|}{P_i + Q_i}} |
#' | **Chebyshev distance**                            <br> `chebyshev()`         | \eqn{\max(|P_i - Q_i|)} |
#' | **Chord distance**                                <br> `chord()`             | \eqn{\displaystyle \sqrt{\sum_{i=1}^{n} \left(\frac{X_i}{\sqrt{\sum_{j=1}^{n} X_j^2}} - \frac{Y_i}{\sqrt{\sum_{j=1}^{n} Y_j^2}}\right)^2}} |
#' | **Clark's divergence distance**                   <br> `clark()`             | \eqn{\displaystyle \sqrt{\sum_{i=1}^{n}\left(\frac{P_i - Q_i}{P_i + Q_i}\right)^{2}}} |
#' | **Divergence**                                    <br> `divergence()`        | \eqn{\displaystyle 2\sum_{i=1}^{n} \frac{(P_i - Q_i)^2}{(P_i + Q_i)^2}} |
#' | **Euclidean distance**                            <br> `euclidean()`         | \eqn{\sqrt{\sum_{i=1}^{n} (P_i - Q_i)^2}} |
#' | **Gower distance**                                <br> `gower()`             | \eqn{\displaystyle \frac{1}{n}\sum_{i=1}^{n}\frac{|P_i - Q_i|}{R_i}} |
#' | **Hellinger distance**                            <br> `hellinger()`         | \eqn{\sqrt{\sum_{i=1}^{n}(\sqrt{P_i} - \sqrt{Q_i})^{2}}} |
#' | **Horn-Morisita dissimilarity**                   <br> `horn()`              | \eqn{\displaystyle 1 - \frac{2\sum_{i=1}^{n}P_{i}Q_{i}}{\sum_{i=1}^{n}P_i^2 + \sum_{i=1}^{n}Q_i^2}} |
#' | **Jensen-Shannon distance**                       <br> `jensen()`            | \eqn{\displaystyle \sqrt{\frac{1}{2}\left[\sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)\right]}} |
#' | **Jensen-Shannon divergence (JSD)**               <br> `jsd()`               | \eqn{\displaystyle \frac{1}{2}\left[\sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)\right]} |
#' | **Lorentzian distance**                           <br> `lorentzian()`        | \eqn{\sum_{i=1}^{n}\ln{(1 + |P_i - Q_i|)}} |
#' | **Manhattan distance**                            <br> `manhattan()`         | \eqn{\sum_{i=1}^{n} |P_i - Q_i|} |
#' | **Matusita distance**                             <br> `matusita()`          | \eqn{\sqrt{\sum_{i=1}^{n}\left(\sqrt{P_i} - \sqrt{Q_i}\right)^2}} |
#' | **Minkowski distance**                            <br> `minkowski()`         | \eqn{\sqrt[p]{\sum_{i=1}^{n} (P_i - Q_i)^p}} <br> Where \eqn{p} is the geometry of the space. |
#' | **Morisita dissimilarity** <br> * Integers Only   <br> `morisita()`          | \eqn{\displaystyle 1 - \frac{2\sum_{i=1}^{n}X_{i}Y_{i}}{\displaystyle \left(\frac{\sum_{i=1}^{n}X_i(X_i - 1)}{X_T(X_T - 1)} + \frac{\sum_{i=1}^{n}Y_i(Y_i - 1)}{Y_T(Y_T - 1)}\right)X_{T}Y_{T}}} |
#' | **Motyka dissimilarity**                          <br> `motyka()`            | \eqn{\displaystyle \frac{\sum_{i=1}^{n} \max(P_i, Q_i)}{\sum_{i=1}^{n} (P_i + Q_i)}} |
#' | **Probabilistic Symmetric \eqn{\chi^2} distance** <br> `psym_chisq()`        | \eqn{\displaystyle 2\sum_{i=1}^{n}\frac{(P_i - Q_i)^2}{P_i + Q_i}} |
#' | **Soergel distance**                              <br> `soergel()`           | \eqn{\displaystyle \frac{\sum_{i=1}^{n} |P_i - Q_i|}{\sum_{i=1}^{n} \max(P_i, Q_i)}} |
#' | **Squared \eqn{\chi^2} distance**                 <br> `squared_chisq()`     | \eqn{\displaystyle \sum_{i=1}^{n}\frac{(P_i - Q_i)^2}{P_i + Q_i}} |
#' | **Squared Chord distance**                        <br> `squared_chord()`     | \eqn{\sum_{i=1}^{n}\left(\sqrt{P_i} - \sqrt{Q_i}\right)^2} |
#' | **Squared Euclidean distance**                    <br> `squared_euclidean()` | \eqn{\sum_{i=1}^{n} (P_i - Q_i)^2} |
#' | **Topsoe distance**                               <br> `topsoe()`            | \eqn{\displaystyle \sum_{i=1}^{n}P_i\ln\left(\frac{2P_i}{P_i + Q_i}\right) + \sum_{i=1}^{n}Q_i\ln\left(\frac{2Q_i}{P_i + Q_i}\right)} |
#' | **Wave Hedges distance**                          <br> `wave_hedges()`       | \eqn{\displaystyle \frac{\sum_{i=1}^{n} |P_i - Q_i|}{\sum_{i=1}^{n} \max(P_i, Q_i)}} |
#' 
#' 
#' ## Presence / Absence
#' 
#' Given:
#' 
#' * \eqn{A}, \eqn{B} : Number of features in each sample.
#' * \eqn{J} : Number of features in common.
#' 
#' |                   |                             |
#' | :---------------- | :-------------------------- |
#' | **Dice-Sorensen dissimilarity** <br> `sorensen()` | \eqn{\displaystyle \frac{2J}{(A + B)}}         |
#' | **Hamming distance**            <br> `hamming()`  | \eqn{\displaystyle (A + B) - 2J}               |
#' | **Jaccard distance**            <br> `jaccard()`  | \eqn{\displaystyle 1 - \frac{J}{(A + B - J)]}} |
#' | **Otsuka-Ochiai dissimilarity** <br> `ochiai()`   | \eqn{\displaystyle 1 - \frac{J}{\sqrt{AB}}}   |
#' 
#' 
#' ## Phylogenetic
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a pair of samples' binary
#' (\eqn{A} and \eqn{B}) or proportional abundances (\eqn{P} and \eqn{Q}) on
#' each of those branches.
#' 
#' |              |              |
#' | :----------- | :----------- |
#' | **Unweighted UniFrac**                  <br> `unweighted_unifrac()`        | \eqn{\displaystyle \frac{1}{n}\sum_{i=1}^{n} L_i|A_i - B_i|} |
#' | **Weighted UniFrac**                    <br> `weighted_unifrac()`          | \eqn{\displaystyle \sum_{i=1}^{n} L_i|P_i - Q_i|} |
#' | **Normalized Weighted UniFrac**         <br> `normalized_unifrac()`        | \eqn{\displaystyle \frac{\sum_{i=1}^{n} L_i|P_i - Q_i|}{\sum_{i=1}^{n} L_i(P_i + Q_i)}} |
#' | **Generalized UniFrac (GUniFrac)**      <br> `generalized_unifrac()`       | \eqn{\displaystyle \frac{\sum_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}\left|\displaystyle \frac{P_i - Q_i}{P_i + Q_i}\right|}{\sum_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}}} <br> Where \eqn{\alpha} is a scalable weighting factor. |
#' | **Variance-Adjusted Weighted UniFrac**  <br> `variance_adjusted_unifrac()` | \eqn{\displaystyle \frac{\displaystyle \sum_{i=1}^{n} L_i\displaystyle \frac{|P_i - Q_i|}{\sqrt{(P_i + Q_i)(2 - P_i - Q_i)}} }{\displaystyle \sum_{i=1}^{n} L_i\displaystyle \frac{P_i + Q_i}{\sqrt{(P_i + Q_i)(2 - P_i - Q_i)}} }} |
#' 
#' See `vignette('unifrac')` for detailed example UniFrac calculations.
#' 
#' 
#' @references
#' 
#' Levy, A., Shalom, B. R., & Chalamish, M. (2024). A guide to similarity
#' measures. *arXiv*.
#' 
#' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures
#' between probability density functions. *International Journal of Mathematical
#' Models and Methods in Applied Sciences*, 1(4), 300–307.
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     t(ex_counts)
#'     
#'     bray(ex_counts)
#'     
#'     jaccard(ex_counts)
#'     
#'     generalized_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for Saliva vs all.
#'     bray(ex_counts, pairs = 1:3)
#'     
NULL

#  Aitchison
#   x <- log((x + pseudocount) / exp(mean(log(x + pseudocount))))
#   y <- log((y + pseudocount) / exp(mean(log(y + pseudocount))))
#   sqrt(sum((x-y)^2)) # Euclidean distance
#' @export
#' @rdname bdiv_functions
aitchison <- function (counts, pseudocount = NULL, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  counts <- transform_clr(counts, pseudocount)
  
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, pairs, cpus, NULL)
}



#  Bhattacharyya
#  -log(sum(sqrt(x * y)))
#' @export
#' @rdname bdiv_functions
bhattacharyya <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_BHATTACHARYYA, counts, pairs, cpus, NULL)
}



# Bray-Curtis  
# sum(abs(x-y)) / sum(x+y)
#' @export
#' @rdname bdiv_functions
bray <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_BRAY, counts, pairs, cpus, NULL)
}


#  Canberra
#  sum(abs(x-y) / (x+y))
#' @export
#' @rdname bdiv_functions
canberra <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CANBERRA, counts, pairs, cpus, NULL)
}


#  Chebyshev
#  max(abs(x-y))
#' @export
#' @rdname bdiv_functions
chebyshev <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CHEBYSHEV, counts, pairs, cpus, NULL)
}


#  Chord
#  sqrt(sum(((x / sqrt(sum(x ^ 2))) - (y / sqrt(sum(y ^ 2))))^2))
#' @export
#' @rdname bdiv_functions
chord <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  counts <- transform_chord(counts)
  
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, pairs, cpus, NULL)
}


#  Clark
#  sqrt(sum((abs(x - y) / (x + y)) ^ 2))  
#' @export
#' @rdname bdiv_functions
clark <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_CLARK, counts, pairs, cpus, NULL)
}


#  Divergence
#  2 * sum((x - y)^2 / (x + y)^2)
#' @export
#' @rdname bdiv_functions
divergence <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_DIVERGENCE, counts, pairs, cpus, NULL)
}


#  Euclidean
#  sqrt(sum((x-y)^2))
#' @export
#' @rdname bdiv_functions
euclidean <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_EUCLIDEAN, counts, pairs, cpus, NULL)
}


#  Gower
#  r <- abs(x - y)
#  n <- length(x) # <-- not `sum(x|y)`
#  sum(abs(x-y) / r) / n
#' @export
#' @rdname bdiv_functions
gower <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  range_vec <- apply(counts, 2L, function (x) diff(range(x)))
  
  .Call(C_beta_div, BDIV_GOWER, counts, pairs, cpus, range_vec)
}


#  Hellinger
#  sqrt(sum((sqrt(x) - sqrt(y)) ^ 2))
#' @export
#' @rdname bdiv_functions
hellinger <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  sqc <- .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, pairs, cpus, NULL)
  
  sqrt(sqc)
}


#  Horn-Morisita
#  z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
#  1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))
#' @export
#' @rdname bdiv_functions
horn <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_HORN, counts, pairs, cpus, NULL)
}


#  Jensen-Shannon distance
#  sqrt(sum(x * log(2 * x / (x+y)), y * log(2 * y / (x+y))) / 2)  
#' @export
#' @rdname bdiv_functions
jensen <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  jsd <- .Call(C_beta_div, BDIV_JSD, counts, pairs, cpus, NULL)
  
  sqrt(jsd)
}


#  Jensen-Shannon Divergence (JSD)
#  sum(x * log(2 * x / (x+y)), y * log(2 * y / (x+y))) / 2  
#' @export
#' @rdname bdiv_functions
jsd <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_JSD, counts, pairs, cpus, NULL)
}


#  Lorentzian
#  sum(log(1 + abs(x - y)))
#' @export
#' @rdname bdiv_functions
lorentzian <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_LORENTZIAN, counts, pairs, cpus, NULL)
}


#  Manhattan
#  sum(abs(x-y))
#' @export
#' @rdname bdiv_functions
manhattan <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MANHATTAN, counts, pairs, cpus, NULL)
}


#  Matusita
#  sqrt(sum((sqrt(x) - sqrt(y)) ^ 2))
#' @export
#' @rdname bdiv_functions
matusita <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  sqc <- .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, pairs, cpus, NULL)
  
  sqrt(sqc)
}


#  Minkowski
#  p <- 1.5
#  sum(abs(x - y)^p) ^ (1/p)
#' @export
#' @rdname bdiv_functions
minkowski <- function (counts, rescale = TRUE, power = 1.5, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MINKOWSKI, counts, pairs, cpus, power)
}


#  Morisita
#  simp_x <- sum(x * (x - 1)) / (sum(x) * (sum(x) - 1))
#  simp_y <- sum(y * (y - 1)) / (sum(y) * (sum(y) - 1))
#  1 - ((2 * sum(x * y)) / ((simp_x + simp_y) * sum(x) * sum(y))) 
#' @export
#' @rdname bdiv_functions
morisita <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_beta_div, BDIV_MORISITA, counts, pairs, cpus, NULL)
}


#  Motyka
#  sum(pmax(x, y)) / sum(x, y)
#' @export
#' @rdname bdiv_functions
motyka <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_MOTYKA, counts, pairs, cpus, NULL)
}


#  Probabilistic Symmetric Chi-Squared  
#  2 * sum((x - y)^2 / (x + y))
#' @export
#' @rdname bdiv_functions
psym_chisq <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  scs <- .Call(C_beta_div, BDIV_SQUARED_CHISQ, counts, pairs, cpus, NULL)
  
  2 * scs
}


#  Soergel
#  sum(abs(x - y)) / sum(pmax(x, y))
#' @export
#' @rdname bdiv_functions
soergel <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_SOERGEL, counts, pairs, cpus, NULL)
}


#  Squared Chi-Squared
#  sum((x - y)^2 / (x + y))
#' @export
#' @rdname bdiv_functions
squared_chisq <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_SQUARED_CHISQ, counts, pairs, cpus, NULL)
}


#  Squared Chord
#  sum((sqrt(x) - sqrt(y)) ^ 2)
#' @export
#' @rdname bdiv_functions
squared_chord <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_SQUARED_CHORD, counts, pairs, cpus, NULL)
}


#  Squared Euclidean
#  sum((x-y)^2)
#' @export
#' @rdname bdiv_functions
squared_euclidean <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  euc <- .Call(C_beta_div, BDIV_EUCLIDEAN, counts, pairs, cpus, NULL)
  
  euc ^ 2
}


#  Topsoe
#  sum(x * log(2 * x / (x+y)), y * log(2 * y / (x+y)))
#' @export
#' @rdname bdiv_functions
topsoe <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  
  jsd <- .Call(C_beta_div, BDIV_JSD, counts, pairs, cpus, NULL)
  
  2 * jsd
}


#  Wave Hedges
#  sum(abs(x - y) / pmax(x, y))
#' @export
#' @rdname bdiv_functions
wave_hedges <- function (counts, rescale = TRUE, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_WAVE_HEDGES, counts, pairs, cpus, NULL)
}




# Presence/Absence -----------


#  Hamming
#  sum(xor(x, y))
#' @export
#' @rdname bdiv_functions
hamming <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_HAMMING, counts, pairs, cpus, NULL)
}


#  Jaccard
#  1 - sum(x & y) / sum(x | y)
#' @export
#' @rdname bdiv_functions
jaccard <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_JACCARD, counts, pairs, cpus, NULL)
}


#  Otsuka-Ochiai
#  1 - sum(x & y) / sqrt(sum(x>0) * sum(y>0)) 
#' @export
#' @rdname bdiv_functions
ochiai <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_OCHIAI, counts, pairs, cpus, NULL)
}


#  Dice-Sorensen
#  2 * sum(x & y) / sum(x>0, y>0)
#' @export
#' @rdname bdiv_functions
sorensen <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_beta_div, BDIV_SORENSEN, counts, pairs, cpus, NULL)
}




# UniFrac Family -----------

#' @export
#' @rdname bdiv_functions
unweighted_unifrac <- function (counts, tree = NULL, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_unifrac, U_UNIFRAC, counts, tree, pairs, cpus, NULL)
}


#' @export
#' @rdname bdiv_functions
weighted_unifrac <- function (counts, tree = NULL, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_unifrac, W_UNIFRAC, counts, tree, pairs, cpus, NULL)
}


#' @export
#' @rdname bdiv_functions
normalized_unifrac <- function (counts, tree = NULL, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_unifrac, N_UNIFRAC, counts, tree, pairs, cpus, NULL)
} 


#' @export
#' @rdname bdiv_functions
generalized_unifrac <- function (counts, tree = NULL, alpha = 0.5, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_unifrac, G_UNIFRAC, counts, tree, pairs, cpus, alpha)
}


#' @export
#' @rdname bdiv_functions
variance_adjusted_unifrac <- function (counts, tree = NULL, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_unifrac, V_UNIFRAC, counts, tree, pairs, cpus, NULL)
}
