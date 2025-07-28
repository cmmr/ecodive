# References as given by
# https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282


#' Chao1
#' 
#' Chao1 alpha diversity metric.\cr\cr
#' A non-parametric estimator of the number of unobserved species in a sample.
#' The Chao1 index estimates total species richness based on the number of 
#' species that occur only once (singletons) and twice (doubletons) in the 
#' sample.
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @return A numeric vector.
#' 
#' @section Calculation:
#' 
#' Prerequisite: all counts are whole numbers.
#' 
#' In the formulas below, `x` is a single column (sample) from `counts`. 
#' \eqn{n} is the total number of non-zero OTUs, \eqn{a} is the number of 
#' singletons, and \eqn{b} is the number of doubletons.
#' 
#' \tabular{lll}{
#' \strong{Formula} \tab  \tab \strong{Example} \cr
#' \deqn{D = \displaystyle n + \frac{a^{2}}{2b}}
#' \tab
#' \tab
#' ```
#'   x <- c(1, 0, 3, 2, 6)
#'   sum(x>0) + (sum(x==1) ^ 2) / (2 * sum(x==2))  
#'   #>  4.5
#' ```
#' }
#' 
#' Note that when \eqn{x} does not have any singletons or doubletons 
#' (\eqn{a = 0, b = 0}), the result will be `NaN`. When \eqn{x} has singletons
#' but no doubletons (\eqn{a > 0, b = 0}), the result will be `Inf`.
#' 
#' @references
#' 
#' Chao A 1984.
#' Non-parametric estimation of the number of classes in a population.
#' Scandinavian Journal of Statistics, 11:265-270.
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Chao1 diversity values
#'     chao1(ex_counts)
#'     
#'     # Low diversity
#'     chao1(c(100, 1, 1, 1, 1)) # Inf
#'     
#'     # High diversity
#'     chao1(c(20, 20, 20, 20, 20)) # NaN
#'     
#'     # Low richness
#'     chao1(1:3) # 3.5
#'     
#'     # High richness
#'     chao1(1:100) # 100.5
#'     
chao1 <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  result_vec <- init_result_vec(counts)
  
  .Call(
    C_alpha_div, 1L, 
    counts, cpus, result_vec )
}


#' Faith's PD
#' 
#' Faith's phylogenetic diversity metric.\cr\cr
#' A higher value indicates a greater amount of evolutionary history 
#' represented within the community, suggesting higher biodiversity in terms of 
#' evolutionary relationships.
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @return A numeric vector.
#' 
#' @section Calculation:
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a sample's
#' abundances on each of those branches coded as 1 for present or 0 for absent:
#' 
#' \deqn{\sum_{i = 1}^{n} P_i \times L_i}
#' 
#' @references
#' 
#' Faith DP 1992.
#' Conservation evaluation and phylogenetic diversity.
#' Biological Conservation, 61:1-10.
#' <https://doi.org/10.1016/0006-3207(92)91201-3>
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Faith diversity values
#'     faith(ex_counts, tree = ex_tree)
#'     
faith <- function (counts, tree, cpus = n_cpus()) {
  
  validate_args()
  result_vec <- init_result_vec(counts)
  
  .Call(
    C_faith, 
    counts, tree, cpus, result_vec )
}


#' Inverse Simpson
#' 
#' Inverse Simpson alpha diversity metric.
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @return A numeric vector.
#' 
#' @section Calculation:
#' 
#' Pre-transformation: drop all OTUs with zero abundance.
#' 
#' In the formulas below, \eqn{x} is a single column (sample) from `counts`.
#' \eqn{p} are the relative abundances.
#' 
#' \tabular{lll}{
#' \strong{Formula} \tab  \tab \strong{Example} \cr
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle 1 / \sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' \tab
#' \tab
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   1 / sum(p * log(p))
#'   #>  -0.7636352
#' ```
#' }
#' 
#' @references
#' 
#' Simpson EH 1949.
#' Measurement of diversity.
#' Nature, 163.
#' <https://doi.org/10.1038/163688a0>
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Inverse Simpson diversity values
#'     inv_simpson(ex_counts)
#'     
#'     # Low diversity
#'     inv_simpson(c(100, 1, 1, 1, 1)) # 1.08
#'     
#'     # High diversity
#'     inv_simpson(c(20, 20, 20, 20, 20)) # 5
#'     
#'     # Low richness
#'     inv_simpson(1:3) # 2.57
#'     
#'     # High richness
#'     inv_simpson(1:100) # 75.37
#'     
inv_simpson <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  result_vec <- init_result_vec(counts)
  
  .Call(
    C_alpha_div, 2L, 
    counts, cpus, result_vec )
}


#' Shannon
#' 
#' Shannon alpha diversity metric.\cr\cr
#' The index considers both the number of different OTUs (richness) and how 
#' evenly the observations are distributed among those OTUs (evenness).
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @return A numeric vector.
#' 
#' @section Calculation:
#' 
#' Pre-transformation: drop all OTUs with zero abundance.
#' 
#' In the formulas below, \eqn{x} is a single column (sample) from `counts`.
#' \eqn{p_i} is the proportion of the \eqn{i}-th OTU in the total community.
#' 
#' \tabular{lll}{
#' \strong{Formula} \tab  \tab \strong{Example} \cr
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle -\sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' \tab
#' \tab
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   -sum(p * log(p))
#'   #>  1.309526
#' ```
#' }
#' 
#' @references
#' 
#' Shannon CE, Weaver W 1949.
#' The Mathematical Theory of Communication.
#' University of Illinois Press.
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Shannon diversity values
#'     shannon(ex_counts)
#'     
#'     # Low diversity
#'     shannon(c(100, 1, 1, 1, 1)) # 0.22
#'     
#'     # High diversity
#'     shannon(c(20, 20, 20, 20, 20)) # 1.61
#'     
#'     # Low richness
#'     shannon(1:3) # 1.01
#'     
#'     # High richness
#'     shannon(1:100) # 4.42
#'     
shannon <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  result_vec <- init_result_vec(counts)
  
  .Call(
    C_alpha_div, 3L, 
    counts, cpus, result_vec )
}


#' Simpson
#' 
#' Simpson alpha diversity metric.\cr\cr
#' Gauges the uniformity of species within a community. A Simpson index of `1` 
#' indicates that one or a few high abundance OTUs dominate the community, 
#' which is indicative of low diversity.
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @return A numeric vector.
#' 
#' @section Calculation:
#' 
#' Pre-transformation: drop all OTUs with zero abundance.
#' 
#' In the formulas below, \eqn{x} is a single column (sample) from `counts`.
#' \eqn{p} are the relative abundances.
#' 
#' \tabular{lll}{
#' \strong{Formula} \tab  \tab \strong{Example} \cr
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle 1 - \sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' \tab
#' \tab
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   1 - sum(p * log(p))
#'   #>  2.309526
#' ```
#' }
#' 
#' @references
#' 
#' Simpson EH 1949.
#' Measurement of diversity.
#' Nature, 163.
#' <https://doi.org/10.1038/163688a0>
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Simpson diversity values
#'     simpson(ex_counts)
#'     
#'     # Low diversity
#'     simpson(c(100, 1, 1, 1, 1)) # 0.075
#'     
#'     # High diversity
#'     simpson(c(20, 20, 20, 20, 20)) # 0.8
#'     
#'     # Low richness
#'     simpson(1:3) # 0.61
#'     
#'     # High richness
#'     simpson(1:100) # 0.99
#'     
simpson <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  result_vec <- init_result_vec(counts)
  
  .Call(
    C_alpha_div, 4L, 
    counts, cpus, result_vec )
}






init_result_vec <- function (counts) {
  result_vec        <- rep(NA_real_, ncol(counts))
  names(result_vec) <- colnames(counts)
  return (result_vec)
}
