# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Formulas as given by
# help('vegdist', 'vegan')

# References as given by
# https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282



#' Beta diversity wrapper
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @param metric   The name of a beta diversity metric. Current options are
#'   `c('bray_curtis', 'sorenson', 'canberra', 'euclidean', 'gower', 'jaccard',
#'   'kulczynski', 'manhattan', 'unweighted_unifrac', 'weighted_unifrac',
#'   'weighted_normalized_unifrac', 'generalized_unifrac',
#'   'variance_adjusted_unifrac')`. Supports partial name matching. Options are
#'   also available via `names(metrics$beta)`.
#'   
#' @param ...  Additional options to pass through to the called function. I.e.
#'   `tree`, `weighted`, `pairs`, `alpha`, or `cpus`.
#' 
#' @return A numeric vector.
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Shannon diversity values
#'     alpha_div(ex_counts, 'Shannon')
#'     
#'     # Chao1 diversity values
#'     alpha_div(ex_counts, 'c')
#'     
#'     # Faith PD values
#'     alpha_div(ex_counts, 'faith', tree = ex_tree)
#'     
#'     
beta_div <- function (counts, metric, ...) {
  metric <- match.arg(tolower(metric), names(metrics$beta))
  do.call(metrics$beta[[metric]], list(counts = counts, ...))
}



#' Bray-Curtis
#' 
#' Bray-Curtis beta diversity metric. Sorenson is an alias for unweighted Bray-Curtis.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{D = \displaystyle \frac{\sum_{i = 1}^{n} |x_i - y_i|}{\sum_{i = 1}^{n} (x_i + y_i)}}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)
#'   y <- c(0, 8, 0, 0, 5)
#'   sum(abs(x-y)) / sum(x+y)  
#'   #>  0.6428571
#' ```
#' 
#' @references
#' 
#' Sorenson T 1948.
#' A method of establishing groups of equal amplitude in plant sociology based on similarity of species content.
#' Kongelige Danske Videnskabernes Selskab, 5.
#' 
#' Bray JR and Curtis JT 1957.
#' An ordination of the upland forest communities of southern Wisconsin.
#' Ecological Monographs, 27(4).
#' \doi{10.2307/1942268}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Bray-Curtis weighted distance matrix
#'     bray_curtis(ex_counts)
#'     
#'     # Bray-Curtis unweighted distance matrix
#'     bray_curtis(ex_counts, weighted = FALSE)
#'     
#'     # Sorenson is the same as unweighted Bray-Curtis
#'     sorenson(ex_counts)
#'     
#'     # Only calculate distances for A vs all.
#'     bray_curtis(ex_counts, pairs = 1:3)
#'     
bray_curtis <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus()) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 1L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' @rdname bray_curtis
#' @export
sorenson <- function (counts, pairs = NULL, cpus = n_cpus()) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  weighted <- FALSE
  
  .Call(
    C_beta_div, 1L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Canberra
#' 
#' Canberra beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' OTUs must be removed if they are absent from both samples.
#' 
#' \deqn{D = \displaystyle \frac{1}{n}\sum_{i = 1}^{n} \frac{|x_i - y_i|}{x_i + y_i}}
#' 
#' ```
#'   x <- c(4, 0, 3, 0, 6)[-4]
#'   y <- c(0, 8, 0, 0, 5)[-4]
#'   sum(abs(x-y) / (x+y)) / length(x)  
#'   #>  0.7727273
#' ```
#' 
#' @references
#' 
#' Lance GN and Williams WT 1967.
#' A general theory of classificatory sorting strategies II. Clustering systems.
#' The computer journal, 10(3).
#' \doi{10.1093/comjnl/10.3.271}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Gower weighted distance matrix
#'     canberra(ex_counts)
#'     
#'     # Gower unweighted distance matrix
#'     canberra(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     canberra(ex_counts, pairs = 1:3)
#'     
canberra <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 2L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Euclidean
#' 
#' Euclidean beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{D = \displaystyle \sqrt{\sum_{i = 1}^{n} (x_i - y_i)^{2}}}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)  
#'   y <- c(0, 8, 0, 0, 5)  
#'   sqrt(sum((x-y)^2))
#'   #>  9.69536
#' ```
#' 
#' @references
#' 
#' Gower JC, Legendre P 1986.
#' Metric and Euclidean Properties of Dissimilarity Coefficients.
#' Journal of Classification. 3.
#'  \doi{10.1007/BF01896809}
#' 
#' Legendre P, Caceres M 2013.
#' Beta diversity as the variance of community data: dissimilarity coefficients and partitioning.
#' Ecology Letters. 16(8).
#' \doi{10.1111/ele.12141}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Euclidean weighted distance matrix
#'     euclidean(ex_counts)
#'     
#'     # Euclidean unweighted distance matrix
#'     euclidean(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     euclidean(ex_counts, pairs = 1:3)
#'     
euclidean <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 3L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Gower
#' 
#' Gower beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Each row (OTU) of `counts` is rescaled to the range 0-1. In cases where a 
#' row is all the same value, those values are replaced with `0`.
#' 
#'     counts                 scaled recounts
#'          A B C  D                 A   B   C D  
#'     OTU1 0 0 0  0    ->    OTU1 0.0 0.0 0.0 0
#'     OTU2 0 8 9 10    ->    OTU2 0.0 0.8 0.9 1
#'     OTU3 5 5 5  5    ->    OTU3 0.0 0.0 0.0 0
#'     OTU4 2 0 0  0    ->    OTU4 1.0 0.0 0.0 0
#'     OTU5 4 6 4  1    ->    OTU5 0.6 1.0 0.6 0
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from the scaled 
#' counts. `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{D = \displaystyle \frac{1}{n}\sum_{i = 1}^{n} |x_i - y_i|}
#' 
#' ```
#'   x <- c(0, 0, 0, 1, 0.6)
#'   y <- c(0, 0.8, 0, 0, 1)
#'   sum(abs(x-y)) / length(x)  
#'   #>  0.44
#' ```
#' 
#' @references
#' 
#' Gower JC 1971.
#' A general coefficient of similarity and some of its properties.
#' Biometrics. 27(4).
#'  \doi{10.2307/2528823}
#' 
#' Gower JC, Legendre P 1986.
#' Metric and Euclidean Properties of Dissimilarity Coefficients.
#' Journal of Classification. 3.
#'  \doi{10.1007/BF01896809}
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Gower weighted distance matrix
#'     gower(ex_counts)
#'     
#'     # Gower unweighted distance matrix
#'     gower(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     gower(ex_counts, pairs = 1:3)
#'     
gower <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_gower, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Jaccard
#' 
#' Jaccard beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{b = \displaystyle \frac{\sum_{i = 1}^{n} |x_i - y_i|}{\sum_{i = 1}^{n} x_i + y_i}}
#' \deqn{D = \displaystyle \frac{2b}{1 + b}}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)  
#'   y <- c(0, 8, 0, 0, 5)  
#'   bray <- sum(abs(x-y)) / sum(x+y)  
#'   2 * bray / (1 + bray)
#'   #>  0.7826087
#' ```
#' 
#' @references
#' 
#' Jaccard P 1908.
#' Nouvellesrecherches sur la distribution florale.
#' Bulletin de la Societe Vaudoise des Sciences Naturelles, 44(163).
#' \doi{10.5169/seals-268384}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Jaccard weighted distance matrix
#'     jaccard(ex_counts)
#'     
#'     # Jaccard unweighted distance matrix
#'     jaccard(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     jaccard(ex_counts, pairs = 1:3)
#'     
jaccard <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 4L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Kulczynski
#' 
#' Kulczynski beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{t = \displaystyle \sum_{i = 1}^{n} min(x_i,y_i)}
#' \deqn{D = \displaystyle 1 - 0.5(\frac{t}{\sum_{i = 1}^{n} x_i} + \frac{t}{\sum_{i = 1}^{n} y_i})}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)
#'   y <- c(0, 8, 0, 0, 5)
#'   t <- sum(pmin(x,y))
#'   1 - (t/sum(x) + t/sum(y)) / 2  
#'   #>  0.6410256
#' ```
#' 
#' @references
#' 
#' Kulcynski S 1927.
#' Die Pflanzenassoziationen der Pieninen.
#' Bulletin International de l'Académie Polonaise des Sciences et des Lettres, Classe des Sciences Mathématiques et Naturelles, Série B: Sciences Naturelles.
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Kulczynski weighted distance matrix
#'     kulczynski(ex_counts)
#'     
#'     # Kulczynski unweighted distance matrix
#'     kulczynski(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     kulczynski(ex_counts, pairs = 1:3)
#'     
kulczynski <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 5L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Manhattan
#' 
#' Manhattan beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' In the formulas below, `x` and `y` are two columns (samples) from `counts`. 
#' `n` is the number of rows (OTUs) in `counts`.
#' 
#' \deqn{D = \displaystyle \sum_{i = 1}^{n} |x_i - y_i|}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)  
#'   y <- c(0, 8, 0, 0, 5)  
#'   sum(abs(x-y))
#'   #>  18
#' ```
#' 
#' @references
#' 
#' Paul EB 2006.
#' Manhattan distance.
#' Dictionary of Algorithms and Data Structures.
#' <https://xlinux.nist.gov/dads/HTML/manhattanDistance.html>
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Manhattan weighted distance matrix
#'     manhattan(ex_counts)
#'     
#'     # Manhattan unweighted distance matrix
#'     manhattan(ex_counts, weighted = FALSE)
#'     
#'     # Only calculate distances for A vs all.
#'     manhattan(ex_counts, pairs = 1:3)
#'     
manhattan <- function (
    counts, 
    weighted = TRUE, 
    pairs    = NULL, 
    cpus     = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_beta_div, 6L, 
    counts, weighted, pairs, cpus, result_dist )
}


#' Unweighted UniFrac
#' 
#' Unweighted UniFrac beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a pair of samples' 
#' abundances (\eqn{A} and \eqn{B}) on each of those branches:
#' 
#' \deqn{D = \displaystyle \frac{\sum_{i = 1}^{n} L_i(|A_i - B_i|)}{\sum_{i = 1}^{n} L_i(max(A_i,B_i))}}
#' 
#' Abundances in \eqn{A} and \eqn{B} are coded as `1` or `0` to indicate their 
#' presence or absence, respectively, on each branch.
#' 
#' See <https://cmmr.github.io/ecodive/articles/unifrac.html> for details and 
#' a worked example.
#' 
#' @references
#' 
#' Lozupone C, Knight R 2005.
#' UniFrac: A new phylogenetic method for comparing microbial communities.
#' Applied and Environmental Microbiology, 71(12).
#' \doi{10.1128/AEM.71.12.8228-8235.2005}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Unweighted UniFrac distance matrix
#'     unweighted_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for A vs all.
#'     unweighted_unifrac(ex_counts, tree = ex_tree, pairs = 1:3)
#'     
unweighted_unifrac <- function (
    counts, 
    tree   = NULL, 
    pairs  = NULL, 
    cpus   = n_cpus() ) {
  
  validate_args()
  alpha <- NA_real_
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_unifrac, 1L, 
    counts, tree, alpha, pairs, cpus, result_dist )
}


#' Weighted UniFrac
#' 
#' Weighted UniFrac beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a pair of samples' 
#' abundances (\eqn{A} and \eqn{B}) on each of those branches:
#' 
#' \deqn{D = \sum_{i = 1}^{n} L_i|\frac{A_i}{A_T} - \frac{B_i}{B_T}|}
#' 
#' See `vignette('unifrac')` for details and a worked example.
#' 
#' @references
#' 
#' Lozupone CA, Hamady M, Kelley ST, Knight R 2007.
#' Quantitative and Qualitative \eqn{\beta} Diversity Measures Lead to Different Insights into Factors That Structure Microbial Communities. 
#' Applied and Environmental Microbiology, 73(5).
#' \doi{10.1128/AEM.01996-06}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Weighted UniFrac distance matrix
#'     weighted_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for A vs all.
#'     weighted_unifrac(ex_counts, tree = ex_tree, pairs = 1:3)
#'     
weighted_unifrac <- function (
    counts, 
    tree   = NULL, 
    pairs  = NULL, 
    cpus   = n_cpus() ) {
  
  validate_args()
  alpha <- NA_real_
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_unifrac, 2L, 
    counts, tree, alpha, pairs, cpus, result_dist )
}


#' Normalized UniFrac
#' 
#' Normalized UniFrac beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a pair of samples' 
#' abundances (\eqn{A} and \eqn{B}) on each of those branches:
#' 
#' \deqn{D = \displaystyle \frac{\sum_{i = 1}^{n} L_i|\frac{A_i}{A_T} - \frac{B_i}{B_T}|}{\sum_{i = 1}^{n} L_i(\frac{A_i}{A_T} + \frac{B_i}{B_T})}}
#' 
#' See `vignette('unifrac')` for details and a worked example.
#' 
#' @references
#' 
#' Lozupone CA, Hamady M, Kelley ST, Knight R 2007.
#' Quantitative and Qualitative \eqn{\beta} Diversity Measures Lead to Different Insights into Factors That Structure Microbial Communities. 
#' Applied and Environmental Microbiology, 73(5).
#' \doi{10.1128/AEM.01996-06}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # UniFrac weighted distance matrix
#'     weighted_normalized_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for A vs all.
#'     weighted_normalized_unifrac(ex_counts, tree = ex_tree, pairs = 1:3)
#'     
weighted_normalized_unifrac <- function (
    counts, 
    tree   = NULL, 
    pairs  = NULL, 
    cpus   = n_cpus() ) {
  
  validate_args()
  alpha <- NA_real_
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_unifrac, 3L, 
    counts, tree, alpha, pairs, cpus, result_dist )
} 


#' Generalized UniFrac
#' 
#' Generalized UniFrac beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L}, a pair of samples' 
#' abundances (\eqn{A} and \eqn{B}) on each of those branches, and 
#' abundance weighting \eqn{0 \le \alpha \le 1}:
#' 
#' \deqn{D = \displaystyle \frac{\sum_{i = 1}^{n} L_i(\frac{A_i}{A_T} + \frac{B_i}{B_T})^{\alpha}|\displaystyle \frac{\frac{A_i}{A_T} - \frac{B_i}{B_T}}{\frac{A_i}{A_T} + \frac{B_i}{B_T}} |}{\sum_{i = 1}^{n} L_i(\frac{A_i}{A_T} + \frac{B_i}{B_T})^{\alpha}}}
#' 
#' See `vignette('unifrac')` for details and a worked example.
#' 
#' @references
#' 
#' Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD, Collman RG, Bushman FD, Li H 2012.
#' Associating microbiome composition with environmental covariates using generalized UniFrac distances. 
#' Bioinformatics, 28(16).
#' \doi{10.1093/bioinformatics/bts342}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Generalized UniFrac distance matrix
#'     generalized_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for A vs all.
#'     generalized_unifrac(ex_counts, tree = ex_tree, pairs = 1:3)
#'     
generalized_unifrac <- function (
    counts, 
    tree   = NULL, 
    alpha  = 0.5, 
    pairs  = NULL, 
    cpus   = n_cpus() ) {
  
  validate_args()
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_unifrac, 4L, 
    counts, tree, alpha, pairs, cpus, result_dist )
}


#' Variance Adjusted UniFrac
#' 
#' Variance Adjusted UniFrac beta diversity metric.
#' 
#' @inherit documentation
#' @family beta_diversity
#' 
#' @return A `dist` object.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a pair of samples' 
#' abundances (\eqn{A} and \eqn{B}) on each of those branches:
#' 
#' \deqn{D = \displaystyle \frac{\sum_{i = 1}^{n} L_i\displaystyle \frac{|\frac{A_i}{A_T} - \frac{B_i}{B_T}|}{\sqrt{(A_i + B_i)(A_T + B_T - A_i - B_i)}} }{\sum_{i = 1}^{n} L_i\displaystyle \frac{\frac{A_i}{A_T} + \frac{B_i}{B_T}}{\sqrt{(A_i + B_i)(A_T + B_T - A_i - B_i)}} }}
#' 
#' See `vignette('unifrac')` for details and a worked example.
#' 
#' @references
#' 
#' Chang Q, Luan Y, Sun F 2011.
#' Variance adjusted weighted UniFrac: a powerful beta diversity measure for comparing communities based on phylogeny. 
#' BMC Bioinformatics, 12.
#' \doi{10.1186/1471-2105-12-118}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Variance Adjusted UniFrac distance matrix
#'     variance_adjusted_unifrac(ex_counts, tree = ex_tree)
#'     
#'     # Only calculate distances for A vs all.
#'     variance_adjusted_unifrac(ex_counts, tree = ex_tree, pairs = 1:3)
#'     
variance_adjusted_unifrac <- function (
    counts, 
    tree   = NULL, 
    pairs  = NULL, 
    cpus   = n_cpus() ) {
  
  validate_args()
  alpha <- NA_real_
  result_dist <- init_result_dist(counts)
  
  .Call(
    C_unifrac, 5L, 
    counts, tree, alpha, pairs, cpus, result_dist )
}



init_result_dist <- function (counts) {
  n <- ncol(counts)
  structure(
    rep(NA_real_, n * (n - 1) / 2),
    class  = 'dist',
    Labels = colnames(counts),
    Size   = n,
    Diag   = FALSE,
    Upper  = FALSE )
}

