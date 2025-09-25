# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Integer IDs for C code.
ADIV_ACE         <-  1L
ADIV_BERGER      <-  2L
ADIV_BRILLOUIN   <-  3L
ADIV_CHAO1       <-  4L
ADIV_FISHER      <-  5L
ADIV_INV_SIMPSON <-  6L
ADIV_MARGALEF    <-  7L
ADIV_MCINTOSH    <-  8L
ADIV_MENHINICK   <-  9L
ADIV_OBSERVED    <- 10L
ADIV_SHANNON     <- 11L
ADIV_SIMPSON     <- 12L
ADIV_SQUARES     <- 13L


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
    norm = 'percent', 
    cutoff  = 10, 
    digits  = 3L, 
    tree    = NULL, 
    cpus    = n_cpus() ) {
  
  metric <- match_metric(metric, div = 'alpha')
  args   <- mget(metric$params, environment())
  
  do.call(metric$func, args)
}





#' Alpha Diversity Metrics
#' 
#' 
#' @inherit documentation
#' @name adiv_functions
#' @family adiv_functions
#' 
#' @return A numeric vector.
#' 
#' 
#' @section Formulas:
#' 
#' Prerequisite: all counts are whole numbers.
#' 
#' Given:
#' 
#' * \eqn{n} : The number of features (e.g. species, OTUs, ASVs, etc).
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' * \eqn{X_T} : Total of all counts (i.e. sequencing depth). \eqn{X_T = \sum_{i=1}^{n} X_i}
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature. \eqn{P_i = X_i / X_T}
#' * \eqn{F_1} : Number of features where \eqn{X_i = 1} (i.e. singletons).
#' * \eqn{F_2} : Number of features where \eqn{X_i = 2} (i.e. doubletons).
#' 
#' |              |                                    |
#' | :----------- | :--------------------------------- |
#' | **Abundance-based Coverage Estimator (ACE)** <br> `ace()`         | See below. |
#' | **Berger-Parker Index**                      <br> `berger()`      | \eqn{\max(P_i)} |
#' | **Brillouin Index**                          <br> `brillouin()`   | \eqn{\displaystyle \frac{\ln{[(\sum_{i = 1}^{n} X_i)!]} - \sum_{i = 1}^{n} \ln{(X_i!)}}{\sum_{i = 1}^{n} X_i}} |
#' | **Chao1**                                    <br> `chao1()`       | \eqn{\displaystyle n + \frac{(F_1)^2}{2 F_2}} |
#' | **Faith's Phylogenetic Diversity**           <br> `faith()`       | See below. |
#' | **Fisher's Alpha (\eqn{\alpha})**            <br> `fisher()`      | \eqn{\displaystyle \frac{n}{\alpha} = \ln{\left(1 + \frac{X_T}{\alpha}\right)}} <br> The value of \eqn{\alpha} must be solved for iteratively. |
#' | **Gini-Simpson Index**                       <br> `simpson()`     | \eqn{1 - \sum_{i = 1}^{n} P_i^2} |
#' | **Inverse Simpson Index**                    <br> `inv_simpson()` | \eqn{1 / \sum_{i = 1}^{n} P_i^2} |
#' | **Margalef's Richness Index**                <br> `margalef()`    | \eqn{\displaystyle \frac{n - 1}{\ln{X_T}}} |
#' | **McIntosh Index**                           <br> `mcintosh()`    | \eqn{\displaystyle \frac{X_T - \sqrt{\sum_{i = 1}^{n} (X_i)^2}}{X_T - \sqrt{X_T}}} |
#' | **Menhinick's Richness Index**               <br> `menhinick()`   | \eqn{\displaystyle \frac{n}{\sqrt{X_T}}} |
#' | **Observed Features**                        <br> `observed()`    | \eqn{n} |
#' | **Shannon Diversity Index**                  <br> `shannon()`     | \eqn{-\sum_{i = 1}^{n} P_i \times \ln(P_i)} |
#' | **Squares Richness Estimator**               <br> `squares()`     | \eqn{\displaystyle n + \frac{(F_1)^2 \sum_{i=1}^{n} (X_i)^2}{X_T^2 - nF_1}} |
#' 
#' 
#' ## Abundance-based Coverage Estimator (ACE)
#' 
#' Given:
#' * \eqn{n} : The number of features (e.g. species, OTUs, ASVs, etc).
#' * \eqn{r} : Rare cutoff. Features with \eqn{\le r} counts are considered rare.
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' * \eqn{F_i} : Number of features with exactly \eqn{i} counts.
#' * \eqn{F_1} : Number of features where \eqn{X_i = 1} (i.e. singletons).
#' * \eqn{F_{rare}} : Number of rare features where \eqn{X_i \le r}.
#' * \eqn{F_{abund}} : Number of abundant features where \eqn{X_i > r}.
#' * \eqn{X_{rare}} : Total counts belonging to rare features.
#' * \eqn{C_{ace}} : The sample abundance coverage estimator, defined below.
#' * \eqn{\gamma_{ace}^2} : The estimated coefficient of variation, defined below.
#' * \eqn{D_{ace}} : Estimated number of features in the sample.
#' 
#' \eqn{\displaystyle C_{ace} = 1 - \frac{F_1}{X_{rare}}}
#' 
#' \eqn{\displaystyle \gamma_{ace}^2 = \max\left[\frac{F_{rare} \sum_{i=1}^{r}i(i-1)F_i}{C_{ace}X_{rare}(X_{rare} - 1)} - 1, 0\right]}
#' 
#' \eqn{\displaystyle D_{ace} = F_{abund} + \frac{F_{rare}}{C_{ace}} + \frac{F_1}{C_{ace}}\gamma_{ace}^2 }
#' 
#' 
#' 
#' 
#' ## Faith's Phylogenetic Diversity (Faith's PD)
#' 
#' Given \eqn{n} branches with lengths \eqn{L} and a sample's abundances 
#' \eqn{A} on each of those branches coded as 1 for present or 0 for absent:
#' 
#' \eqn{\sum_{i = 1}^{n} L_i A_i}
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     t(ex_counts)
#'     
#'     ace(ex_counts)
#'     
#'     chao1(ex_counts)
#'     
#'     squares(ex_counts)
#'     
NULL


#  Abundance-based Coverage Estimator (ACE)
#' @export
#' @rdname adiv_functions
ace <- function (counts, cutoff = 10, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_ACE, counts, cpus, cutoff)
}


#  Berger-Parker
#  max(x / sum(x))
#' @export
#' @rdname adiv_functions
berger <- function (counts, norm = 'percent', cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_BERGER, counts, cpus, NULL)
}


#  Brillouin Index
#  note: lgamma(x + 1) == log(x!)
#  (lgamma(sum(x) + 1) - sum(lgamma(x + 1))) / sum(x)
#' @export
#' @rdname adiv_functions
brillouin <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_BRILLOUIN, counts, cpus, NULL)
}


#  Chao1
#  sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))
#' @export
#' @rdname adiv_functions
chao1 <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_CHAO1, counts, cpus, NULL)
}


#  Faith's Phylogenetic Diversity
#' @export
#' @rdname adiv_functions
faith <- function (counts, tree = NULL, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_faith, counts, tree, cpus)
}


#  Fisher's Alpha
#' @export
#' @rdname adiv_functions
fisher <- function (counts, digits = 3L, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_FISHER, counts, cpus, digits)
}


#  Inverse Simpson Index
#  p <- x / sum(x)
#  1 / sum(p ** 2)
#' @export
#' @rdname adiv_functions
inv_simpson <- function (counts, norm = 'percent', cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_INV_SIMPSON, counts, cpus, NULL)
}


#  Margalef's Richness Index
#  (sum(x > 0) - 1) / log(sum(x))
#' @export
#' @rdname adiv_functions
margalef <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MARGALEF, counts, cpus, NULL)
}


#  McIntosh Index
#  (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
#' @export
#' @rdname adiv_functions
mcintosh <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MCINTOSH, counts, cpus, NULL)
}


#  Menhinick's Richness Index
#  sum(x > 0) / sqrt(sum(x))
#' @export
#' @rdname adiv_functions
menhinick <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MENHINICK, counts, cpus, NULL)
}


#  Observed Features
#  sum(x>0)
#' @export
#' @rdname adiv_functions
observed <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_OBSERVED, counts, cpus, NULL)
}


#  Shannon Diversity Index
#  p <- x / sum(x)
#  -sum(p * log(p))
#' @export
#' @rdname adiv_functions
shannon <- function (counts, norm = 'percent', cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_SHANNON, counts, cpus, NULL)
}


#  Gini-Simpson Index
#  p <- x / sum(x)
#  1 - sum(p ** 2)
#' @export
#' @rdname adiv_functions
simpson <- function (counts, norm = 'percent', cpus = n_cpus()) {
  
  validate_args()
  .Call(C_alpha_div, ADIV_SIMPSON, counts, cpus, NULL)
}


#  Squares Richness Estimator
#  N  <- sum(x)      # sampling depth
#  S  <- sum(x > 0)  # observed features
#  F1 <- sum(x == 1) # singletons
#  S + ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S))
#' @export
#' @rdname adiv_functions
squares <- function (counts, cpus = n_cpus()) {
  
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_SQUARES, counts, cpus, NULL)
}
