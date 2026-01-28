# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Integer IDs for C code.
ADIV_ACE         <-  1L
ADIV_BERGER      <-  2L
ADIV_BRILLOUIN   <-  3L
ADIV_CHAO1       <-  4L
ADIV_FAITH       <-  5L
ADIV_FISHER      <-  6L
ADIV_INV_SIMPSON <-  7L
ADIV_MARGALEF    <-  8L
ADIV_MCINTOSH    <-  9L
ADIV_MENHINICK   <- 10L
ADIV_OBSERVED    <- 11L
ADIV_SHANNON     <- 12L
ADIV_SIMPSON     <- 13L
ADIV_SQUARES     <- 14L


#' Alpha Diversity Wrapper Function
#' 
#' @inherit documentation
#' @name alpha_div
#' 
#' @param metric   The name of an alpha diversity metric. One of `c('ace',
#'   'berger', 'brillouin', 'chao1', 'faith', 'fisher', 'inv_simpson',
#'   'margalef', 'mcintosh', 'menhinick', 'observed', 'shannon', 'simpson',
#'   'squares')`. Case-insensitive and partial name matching is supported.
#'   Programmatic access via `list_metrics('alpha')`.
#' 
#' @return A numeric vector.
#' 
#' @details
#' 
#' ## Integer Count Requirements
#' 
#' A frequent and critical error in alpha diversity analysis is providing the
#' wrong type of data to a metric's formula. Some indices are mathematically
#' defined based on counts of individuals and require raw, integer abundance
#' data. Others are based on proportional abundances and can accept either
#' integer counts (which are then converted to proportions) or pre-normalized
#' proportional data. Using proportional data with a metric that requires
#' integer counts will return an error message.
#' 
#' ### Requires Integer Counts Only
#' 
#' * Chao1  
#' * ACE
#' * Squares Richness Estimator
#' * Margalef's Index
#' * Menhinick's Index
#' * Fisher's Alpha
#' * Brillouin Index
#' 
#' ### Can Use Proportional Data
#' 
#' * Observed Features
#' * Shannon Index
#' * Gini-Simpson Index
#' * Inverse Simpson Index
#' * Berger-Parker Index
#' * McIntosh Index
#' * Faith's PD
#' 
#' 
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
alpha_div <- function (
    counts, 
    metric, 
    norm   = 'percent', 
    cutoff = 10L, 
    digits = 3L, 
    tree   = NULL, 
    margin = 1L,
    cpus   = n_cpus() ) {
  
  metric <- match_metric(metric, div = 'alpha')
  args   <- mget(metric$params, environment())
  
  do.call(metric$func, args)
}


#' Abundance-based Coverage Estimator (ACE)
#' 
#' A non-parametric estimator of species richness that separates features into abundant and rare groups.
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' The ACE metric separates features into "abundant" and "rare" groups based on a cutoff (usually 10 counts). It assumes that the presence of abundant species is certain, while the true number of rare species must be estimated.
#' 
#' **Equations:**
#' 
#' \deqn{C_{ace} = 1 - \frac{F_1}{X_{rare}}}
#' 
#' \deqn{\gamma_{ace}^2 = \max\left[\frac{F_{rare} \sum_{i=1}^{r}i(i-1)F_i}{C_{ace}X_{rare}(X_{rare} - 1)} - 1, 0\right]}
#' 
#' \deqn{D_{ace} = F_{abund} + \frac{F_{rare}}{C_{ace}} + \frac{F_1}{C_{ace}}\gamma_{ace}^2}
#' 
#' Where:
#' * \eqn{r} : Rare cutoff (default 10). Features with \eqn{\le r} counts are considered rare.
#' * \eqn{F_i} : Number of features with exactly \eqn{i} counts.
#' * \eqn{F_1} : Number of features where \eqn{X_i = 1} (singletons).
#' * \eqn{F_{rare}} : Number of rare features where \eqn{X_i \le r}.
#' * \eqn{F_{abund}} : Number of abundant features where \eqn{X_i > r}.
#' * \eqn{X_{rare}} : Total counts belonging to rare features.
#' * \eqn{C_{ace}} : The sample abundance coverage estimator.
#' * \eqn{\gamma_{ace}^2} : The estimated coefficient of variation.
#' 
#' **Parameter: cutoff**
#' The integer threshold distinguishing rare from abundant species. Standard practice is to use 10.
#' 
#' @references
#' Chao, A., & Lee, S. M. (1992). Estimating the number of classes via sample coverage. *Journal of the American Statistical Association*, 87(417), 210-217. \doi{10.1080/01621459.1992.10475194}
#' 
#' @export
#' @examples
#'     ace(ex_counts)
ace <- function (counts, cutoff = 10L, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_ACE, counts, margin, norm, cpus, cutoff)
}


#' Berger-Parker Index
#' 
#' A measure of the numerical importance of the most abundant species.
#' 
#' @inherit documentation
#' @family Dominance metrics
#' @seealso alpha_div
#' 
#' @details
#' The Berger-Parker index is defined as the proportional abundance of the most dominant feature:
#' \deqn{\max(P_i)}
#' 
#' Where:
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' max(x / sum(x))
#' ```
#' 
#' @references
#' Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic foraminifera in deep-sea sediments. *Science*, 168(3937), 1345-1347. \doi{10.1126/science.168.3937.1345}
#' 
#' @export
#' @examples
#'     berger(ex_counts)
berger <- function (counts, norm = 'percent', margin = 1L, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_BERGER, counts, margin, norm, cpus, NULL)
}


#' Brillouin Index
#' 
#' A diversity index derived from information theory, appropriate for fully censused communities.
#' 
#' @inherit documentation
#' @family Diversity metrics
#' @seealso alpha_div
#' 
#' @details
#' The Brillouin index is defined as:
#' \deqn{\frac{\ln{[(\sum_{i = 1}^{n} X_i)!]} - \sum_{i = 1}^{n} \ln{(X_i!)}}{\sum_{i = 1}^{n} X_i}}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' # note: lgamma(x + 1) == log(x!)
#' (lgamma(sum(x) + 1) - sum(lgamma(x + 1))) / sum(x)
#' ```
#' 
#' @references
#' Brillouin, L. (1956). Science and information theory. Academic Press.
#' 
#' @export
#' @examples
#'     brillouin(ex_counts)
brillouin <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_BRILLOUIN, counts, margin, norm, cpus, NULL)
}


#' Chao1 Richness Estimator
#' 
#' A non-parametric estimator of the lower bound of species richness.
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' The Chao1 estimator uses the ratio of singletons to doubletons to estimate the number of missing species:
#' \deqn{n + \frac{(F_1)^2}{2 F_2}}
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{F_1} : Number of features observed once (singletons).
#' * \eqn{F_2} : Number of features observed twice (doubletons).
#' 
#' **Base R Equivalent:**
#' ```r
#' sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))
#' ```
#' 
#' @references
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scandinavian Journal of Statistics*, 11, 265-270.
#' 
#' @export
#' @examples
#'     chao1(ex_counts)
chao1 <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_CHAO1, counts, margin, norm, cpus, NULL)
}


#' Faith's Phylogenetic Diversity (PD)
#' 
#' Calculates the sum of the branch lengths for all species present in a sample.
#' 
#' @inherit documentation
#' @family Phylogenetic metrics
#' @seealso alpha_div
#' 
#' @details
#' Faith's PD is defined as:
#' \deqn{\sum_{i = 1}^{n} L_i A_i}
#' 
#' Where:
#' * \eqn{n} : The number of branches in the phylogenetic tree.
#' * \eqn{L_i} : The length of the \eqn{i}-th branch.
#' * \eqn{A_i} : A binary value (1 if any descendants of branch \eqn{i} are present in the sample, 0 otherwise).
#' 
#' @references
#' Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61(1), 1-10. \doi{10.1016/0006-3207(92)91201-3}
#' 
#' @export
#' @examples
#'     faith(ex_counts, tree = ex_tree)
faith <- function (counts, tree = NULL, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_alpha_div, ADIV_FAITH, counts, margin, norm, cpus, tree)
}


#' Fisher's Alpha
#' 
#' A parametric diversity index assuming species abundances follow a log-series distribution.
#' 
#' @inherit documentation
#' @family Diversity metrics
#' @seealso alpha_div
#' 
#' @details
#' Fisher's Alpha (\eqn{\alpha}) is the parameter in the equation:
#' \deqn{\frac{n}{\alpha} = \ln{\left(1 + \frac{X_T}{\alpha}\right)}}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{X_T} : Total of all counts (sequencing depth).
#' 
#' The value of \eqn{\alpha} is solved for iteratively.
#' 
#' **Parameter: digits**
#' The precision (number of decimal places) to use when solving the equation.
#' 
#' @references
#' Fisher, R. A., Corbet, A. S., & Williams, C. B. (1943). The relation between the number of species and the number of individuals in a random sample of an animal population. *Journal of Animal Ecology*, 12, 42-58. \doi{10.2307/1411}
#' 
#' @export
#' @examples
#'     fisher(ex_counts)
fisher <- function (counts, digits = 3L, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_FISHER, counts, margin, norm, cpus, digits)
}


#' Inverse Simpson Index
#' 
#' A transformation of the Simpson index that represents the "effective number of species".
#' 
#' @inherit documentation
#' @family Diversity metrics
#' @seealso alpha_div
#' 
#' @details
#' The Inverse Simpson index is defined as:
#' \deqn{1 / \sum_{i = 1}^{n} P_i^2}
#' 
#' Where:
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' p <- x / sum(x)
#' 1 / sum(p ** 2)
#' ```
#' 
#' @references
#' Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688. \doi{10.1038/163688a0}
#' 
#' @export
#' @examples
#'     inv_simpson(ex_counts)
inv_simpson <- function (counts, norm = 'percent', margin = 1L, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_INV_SIMPSON, counts, margin, norm, cpus, NULL)
}


#' Margalef's Richness Index
#' 
#' A richness metric that normalizes the number of species by the log of the total sample size.
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' Margalef's index is defined as:
#' \deqn{\frac{n - 1}{\ln{X_T}}}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{X_T} : Total of all counts (sequencing depth).
#' 
#' **Base R Equivalent:**
#' ```r
#' (sum(x > 0) - 1) / log(sum(x))
#' ```
#' 
#' @references
#' Margalef, R. (1958). Information theory in ecology. *General Systems*, 3, 36-71.
#' 
#' @export
#' @examples
#'     margalef(ex_counts)
margalef <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MARGALEF, counts, margin, norm, cpus, NULL)
}


#' McIntosh Index
#' 
#' A dominance index based on the Euclidean distance from the origin.
#' 
#' @inherit documentation
#' @family Dominance metrics
#' @seealso alpha_div
#' 
#' @details
#' The McIntosh index is defined as:
#' \deqn{\frac{X_T - \sqrt{\sum_{i = 1}^{n} (X_i)^2}}{X_T - \sqrt{X_T}}}
#' 
#' Where:
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' * \eqn{X_T} : Total of all counts.
#' 
#' **Base R Equivalent:**
#' ```r
#' (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
#' ```
#' 
#' @references
#' McIntosh, R. P. (1967). An index of diversity and the relation of certain concepts to diversity. *Ecology*, 48(3), 392-404. \doi{10.2307/1932674}
#' 
#' @export
#' @examples
#'     mcintosh(ex_counts)
mcintosh <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MCINTOSH, counts, margin, norm, cpus, NULL)
}


#' Menhinick's Richness Index
#' 
#' A richness metric that normalizes the number of species by the square root of the total sample size.
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' Menhinick's index is defined as:
#' \deqn{\frac{n}{\sqrt{X_T}}}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{X_T} : Total of all counts.
#' 
#' **Base R Equivalent:**
#' ```r
#' sum(x > 0) / sqrt(sum(x))
#' ```
#' 
#' @references
#' Menhinick, E. F. (1964). A comparison of some species-individuals diversity indices applied to samples of field insects. *Ecology*, 45(4), 859-861. \doi{10.2307/1934933}
#' 
#' @export
#' @examples
#'     menhinick(ex_counts)
menhinick <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MENHINICK, counts, margin, norm, cpus, NULL)
}


#' Observed Features
#' 
#' The count of unique features (richness) in a sample.
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' Observed features is defined simply as the number of features with non-zero abundance:
#' \deqn{n}
#' 
#' **Base R Equivalent:**
#' ```r
#' sum(x > 0)
#' ```
#' 
#' @export
#' @examples
#'     observed(ex_counts)
observed <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  
  .Call(C_alpha_div, ADIV_OBSERVED, counts, margin, norm, cpus, NULL)
}


#' Shannon Diversity Index
#' 
#' A commonly used diversity index accounting for both abundance and evenness.
#' 
#' @inherit documentation
#' @family Diversity metrics
#' @seealso alpha_div
#' 
#' @details
#' The Shannon index (entropy) is defined as:
#' \deqn{-\sum_{i = 1}^{n} P_i \times \ln(P_i)}
#' 
#' Where:
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' p <- x / sum(x)
#' -sum(p * log(p))
#' ```
#' 
#' @references
#' Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27, 379-423.
#' 
#' @export
#' @examples
#'     shannon(ex_counts)
shannon <- function (counts, norm = 'percent', margin = 1L, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_SHANNON, counts, margin, norm, cpus, NULL)
}


#' Gini-Simpson Index
#' 
#' The probability that two entities taken at random from the dataset represent different types.
#' 
#' @inherit documentation
#' @family Diversity metrics
#' @seealso alpha_div
#' 
#' @details
#' The Gini-Simpson index is defined as:
#' \deqn{1 - \sum_{i = 1}^{n} P_i^2}
#' 
#' Where:
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' p <- x / sum(x)
#' 1 - sum(p ** 2)
#' ```
#' 
#' @references
#' Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688. \doi{10.1038/163688a0}
#' 
#' @export
#' @examples
#'     simpson(ex_counts)
simpson <- function (counts, norm = 'percent', margin = 1L, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_SIMPSON, counts, margin, norm, cpus, NULL)
}


#' Squares Richness Estimator
#' 
#' A richness estimator based on the concept of "squares" (counts of species observed once or twice).
#' 
#' @inherit documentation
#' @family Richness metrics
#' @seealso alpha_div
#' 
#' @details
#' The Squares estimator is defined as:
#' \deqn{n + \frac{(F_1)^2 \sum_{i=1}^{n} (X_i)^2}{X_T^2 - nF_1}}
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{X_T} : Total of all counts.
#' * \eqn{F_1} : Number of features observed once (singletons).
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' N  <- sum(x)      # sampling depth
#' S  <- sum(x > 0)  # observed features
#' F1 <- sum(x == 1) # singletons
#' S + ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S))
#' ```
#' 
#' @references
#' Alroy, J. (2018). Limits to species richness estimates based on subsampling. *Paleobiology*, 44(2), 177-194. \doi{10.1017/pab.2017.38}
#' 
#' @export
#' @examples
#'     squares(ex_counts)
squares <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- NULL
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_SQUARES, counts, margin, norm, cpus, NULL)
}
