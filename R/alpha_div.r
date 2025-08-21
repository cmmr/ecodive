# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# References as given by
# https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282



#' Alpha diversity wrapper
#' 
#' @inherit documentation
#' @family alpha_diversity
#' 
#' @param metric   The name of an alpha diversity metric. Current options are
#'   `c('observed', 'chao1', 'shannon', 'simpson', 'inv_simpson', 'faith')`.
#'   Supports case-insensitive and partial name matching. Options are also
#'   available via `names(metrics$alpha)`.
#'   
#' @param ...  Additional options to pass through to the called function. I.e.
#'   `cpus` or `tree`.
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
alpha_div <- function (counts, metric, ...) {
  metric <- as_alpha_metric(metric)
  do.call(metrics$alpha[[metric]], list(counts = counts, ...))
}



#' Chao1
#' 
#' Chao1 alpha diversity metric.
#' 
#' The Chao1 index is a non-parametric estimator that seeks to predict the true
#' species richness of a community, including species that were likely missed
#' due to undersampling. It works by using the counts of the rarest observed
#' taxa — specifically "singletons" (taxa seen only once) and "doubletons" (taxa
#' seen twice) — to estimate how many species were not detected at all. The
#' logic is that if you find many species represented by only one or two
#' individuals, it is highly probable that many other rare species were missed
#' entirely.
#' 
#' **Important Caveat:** The Chao1 estimator is mathematically dependent on
#' singleton counts. However, modern bioinformatic pipelines that generate ASVs
#' (like DADA2) are designed to remove singletons, as they are often
#' indistinguishable from sequencing errors. Using Chao1 on data that has had
#' singletons removed will produce a scientifically meaningless result that is
#' often just the same as the observed richness. Therefore, this metric is
#' considered methodologically unsound for most modern ASV-based workflows and
#' should be used with extreme caution.
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
#' \deqn{D = \displaystyle n + \frac{a^{2}}{2b}}
#' 
#' ```
#'   x <- c(1, 0, 3, 2, 6)
#'   sum(x>0) + (sum(x==1) ^ 2) / (2 * sum(x==2))  
#'   #>  4.5
#' ```
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
#' Faith's phylogenetic diversity metric.
#' 
#' Faith's Phylogenetic Diversity (PD) is a unique richness metric that
#' incorporates the evolutionary relationships between the members of a
#' community. Instead of simply counting the number of unique taxa, Faith's PD
#' is calculated as the sum of the lengths of all the branches on a phylogenetic
#' tree that connect all the species present in a sample. The rationale is that
#' a community composed of distantly related organisms (e.g., from different
#' phyla) is more "diverse" in an evolutionary sense than a community with the
#' same number of species that are all closely related (e.g., from the same
#' genus). A higher PD value indicates greater phylogenetic diversity. Because
#' evolutionary relatedness often correlates with function, Faith's PD can also
#' serve as a valuable proxy for the unmeasured functional diversity of a
#' community. This metric requires a phylogenetic tree as an input for its
#' calculation.
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
#' \doi{10.1016/0006-3207(92)91201-3}
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Faith diversity values
#'     faith(ex_counts, tree = ex_tree)
#'     
faith <- function (counts, tree = NULL, cpus = n_cpus()) {
  
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
#' The Inverse Simpson index is another way to present the information from the
#' Simpson index, but in a more intuitive format. While the standard Simpson
#' index calculates the probability of two randomly selected individuals
#' belonging to the same species (where a lower value means more diversity), the
#' Inverse Simpson index is the reciprocal of that value (1/D). This simple
#' transformation makes the index easier to interpret: the value increases as
#' diversity increases.
#' 
#' The primary advantage of the Inverse Simpson index is that its value can be
#' understood as the "effective number of species". This means a community with
#' an Inverse Simpson index of 10 has a diversity that is equivalent to a
#' community composed of 10 equally abundant species. This provides a more
#' direct and biologically meaningful interpretation compared to more abstract
#' indices. Like the standard Simpson index, it is more heavily weighted by the
#' most abundant species and is less sensitive to sampling depth than richness
#' metrics like Observed Features or Chao1.
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
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle 1 / \sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   1 / sum(p * log(p))
#'   #>  -0.7636352
#' ```
#' 
#' @references
#' 
#' Simpson EH 1949.
#' Measurement of diversity.
#' Nature, 163.
#' \doi{10.1038/163688a0}
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


#' Observed Features
#' 
#' Observed Features alpha diversity metric.
#' 
#' This is the most straightforward and intuitive measure of diversity. It is a
#' simple count of the number of unique microbial taxa (such as Amplicon
#' Sequence Variants, or ASVs) detected in a sample. A higher value indicates
#' greater richness. While easy to understand, this metric is highly sensitive
#' to the number of sequences per sample (sequencing depth). A sample with more
#' sequences is more likely to detect rare taxa by chance, leading to an
#' inflated richness value. Therefore, it is not appropriate to directly compare
#' the Observed Features of samples with different sequencing depths without
#' first normalizing the data, typically through a process called rarefaction
#' (subsampling all samples to an equal depth).
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
#' 
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle \sum_{i = 1}^{n} 1}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]
#'   length(x)
#'   #>  4
#' ```
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Observed features
#'     observed(ex_counts)
#'     
observed <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  counts[] <- as.logical(counts)
  colSums(counts)
}


#' Shannon
#' 
#' Shannon alpha diversity metric.
#' 
#' The Shannon index is a widely used metric that quantifies diversity by
#' considering both the number of species (richness) and their abundance
#' distribution (evenness). Borrowed from information theory, it measures the
#' "uncertainty" or entropy in predicting the identity of a microbe drawn
#' randomly from the sample. A community with many different species that are
#' present in similar proportions will have high uncertainty and thus a high
#' Shannon index value. Compared to the Simpson index, the Shannon index gives
#' more equitable weight to both rare and abundant species, making it more
#' sensitive to changes in richness. Higher values indicate greater community
#' diversity.
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
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle -\sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   -sum(p * log(p))
#'   #>  1.309526
#' ```
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
#' Simpson alpha diversity metric.
#' 
#' The Simpson index is a popular metric that incorporates both richness and
#' evenness to describe community diversity. The most common version, the
#' Gini-Simpson index (implemented here), measures the probability that two
#' individuals selected randomly from the community will belong to  different
#' species. The value ranges from 0 to 1, where higher values indicate greater
#' diversity. Because the calculation involves squaring the proportional
#' abundances of each species, the index is heavily weighted by the most
#' abundant (dominant) taxa and is less sensitive to the presence of rare
#' species. A low Simpson index suggests that the community is dominated by one
#' or a few species, making it a strong measure of community dominance.
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
#' \deqn{p_{i} = \displaystyle \frac{x_i}{\sum x}}
#' \deqn{D = \displaystyle 1 - \sum_{i = 1}^{n} p_{i}\times\ln(p_{i})}
#' 
#' ```
#'   x <- c(4, 0, 3, 2, 6)[-2]  
#'   p <- x / sum(x)
#'   1 - sum(p * log(p))
#'   #>  2.309526
#' ```
#' 
#' @references
#' 
#' Simpson EH 1949.
#' Measurement of diversity.
#' Nature, 163.
#' \doi{10.1038/163688a0}
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
