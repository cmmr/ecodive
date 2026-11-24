# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The metrics dominance(), gini(), goods_coverage(), heip(), hill(),
# mcintosh_evenness(), pielou(), renyi(), simpson_evenness(), strong() and
# tsallis(), and chao1()'s bias correction, follow skbio.diversity.alpha
# (scikit-bio 0.7.4), whose formulas and edge cases src/alpha_div.c ports,
# under this licence:
#
# Copyright (c) 2013--, scikit-bio development team.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the names scikit-bio, skbio, or biocore nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
ADIV_DOMINANCE   <- 15L
ADIV_GINI        <- 16L
ADIV_GOODS       <- 17L
ADIV_HEIP        <- 18L
ADIV_HILL        <- 19L
ADIV_MCINTOSH_E  <- 20L
ADIV_PIELOU      <- 21L
ADIV_RENYI       <- 22L
ADIV_SIMPSON_E   <- 23L
ADIV_STRONG      <- 24L
ADIV_TSALLIS     <- 25L


#' Alpha Diversity Wrapper Function
#' 
#' @inherit documentation
#' @name alpha_div
#' 
#' @param metric   The name of an alpha diversity metric. One of `c('ace',
#'   'berger', 'brillouin', 'chao1', 'dominance', 'faith', 'fisher', 'gini',
#'   'goods_coverage', 'heip', 'hill', 'inv_simpson', 'margalef',
#'   'mcintosh', 'mcintosh_evenness', 'menhinick', 'observed', 'pielou',
#'   'renyi', 'shannon', 'simpson', 'simpson_evenness', 'squares', 'strong',
#'   'tsallis')`. Case-insensitive and partial name matching is supported,
#'   as are scikit-bio's names (such as `'pielou_e'`, `'simpson_d'`,
#'   `'observed_features'` and `'enspie'`). Programmatic access via
#'   `list_metrics('alpha')`.
#' 
#' @return A numeric vector, one value per sample, named by the samples.
#'   To compute several metrics, call `alpha_div()` once for each (see the
#'   examples).
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
#' * McIntosh Index
#' * McIntosh's Evenness
#' * Good's Coverage
#' 
#' ### Can Use Proportional Data
#' 
#' * Observed Features
#' * Shannon Index
#' * Gini-Simpson Index
#' * Inverse Simpson Index
#' * Hill Numbers, Rényi and Tsallis Entropies
#' * Pielou's, Heip's and Simpson's Evenness
#' * Berger-Parker Index
#' * Simpson's and Strong's Dominance Indices
#' * Gini Index
#' * Faith's PD
#' 
#' 
#' ## Differences from scikit-bio
#' 
#' The metrics follow `skbio.diversity.alpha` (scikit-bio 0.7.4), which
#' QIIME 2's `diversity alpha` runs, with these differences:
#' 
#' * Logarithms are natural. scikit-bio's `shannon()`, `pielou_e()` and
#'   `renyi()` take a `base`; divide ecodive's [shannon()] and [renyi()] by
#'   `log(base)` instead. Pielou's evenness doesn't depend on the base.
#' * No `finite` option for [simpson()], [dominance()] and [inv_simpson()]
#'   (scikit-bio's estimator without replacement).
#' * The [gini()] index is computed by rectangles only (scikit-bio's default
#'   `method = "rectangles"`, not `"trapezoids"`).
#' * A sample with no observations is `NA` for every metric. scikit-bio
#'   gives `nan` for most metrics, but `0` for `chao1()` and `gini_index()`.
#' * `chao1(bias_corrected = FALSE)` is the classic formula, falling back on
#'   the bias-corrected one when a sample has no singletons or no doubletons,
#'   exactly as scikit-bio does.
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
#'     # Hill number of order 1 (the exponential of Shannon's entropy)
#'     alpha_div(ex_counts, 'hill', q = 1)
#'     
#'     # Several metrics: one call each
#'     sapply(c('observed', 'shannon', 'pielou'), alpha_div, counts = ex_counts)
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
    q      = 2,
    bias_corrected = TRUE,
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
#' @inherit adiv_assert_integer
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_ACE, counts, margin, norm, cpus, cutoff)
}


#' Berger-Parker Index
#' 
#' A measure of the numerical importance of the most abundant species.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Dominance metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' x <- ex_counts[1,]
#' p <- x / sum(x)
#' max(p)
#' ```
#' 
#' @references
#' Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic foraminifera in deep-sea sediments. *Science*, 168(3937), 1345-1347. \doi{10.1126/science.168.3937.1345}
#' 
#' @export
#' @examples
#'     berger(ex_counts)
berger <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_BERGER, counts, margin, norm, cpus, NULL)
}


#' Brillouin Index
#' 
#' A diversity index derived from information theory, appropriate for fully censused communities.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' x <- ex_counts[1,]
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
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_BRILLOUIN, counts, margin, norm, cpus, NULL)
}


#' Chao1 Richness Estimator
#' 
#' A non-parametric estimator of the lower bound of species richness.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The Chao1 estimator uses the ratio of singletons to doubletons to
#' estimate the number of missing species. By default it is the
#' bias-corrected form (Chao 1987), as in scikit-bio and QIIME 2:
#' \deqn{n + \frac{F_1 (F_1 - 1)}{2 (F_2 + 1)}}
#' 
#' With `bias_corrected = FALSE`, it is the classic form (Chao 1984):
#' \deqn{n + \frac{(F_1)^2}{2 F_2}}
#' which is infinite when there are no doubletons. So, as scikit-bio's
#' `chao1(bias_corrected = False)`, the classic form is used only when
#' \eqn{F_1 > 0} and \eqn{F_2 > 0}, and the bias-corrected form otherwise.
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{F_1} : Number of features observed once (singletons).
#' * \eqn{F_2} : Number of features observed twice (doubletons).
#' 
#' `vegan::estimateR()` and iNEXT multiply the bias-corrected term by
#' \eqn{(N - 1)/N}, where \eqn{N} is the number of observations, so their
#' estimates are slightly smaller than these for shallow samples with
#' singletons.
#' 
#' ecodive before version 2.4.0 used the classic form without the fallback,
#' which gave `Inf` without doubletons and `NaN` without singletons and
#' doubletons.
#' 
#' **Base R Equivalent:**
#' ```r
#' x  <- ex_counts[1,]
#' F1 <- sum(x == 1)
#' F2 <- sum(x == 2)
#' sum(x > 0) + F1 * (F1 - 1) / (2 * (F2 + 1))
#' ```
#' 
#' @references
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scandinavian Journal of Statistics*, 11, 265-270.
#' 
#' Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. *Biometrics*, 43(4), 783-791. \doi{10.2307/2531532}
#' 
#' @export
#' @examples
#'     chao1(ex_counts)
#'     
#'     # The classic form
#'     chao1(ex_counts, bias_corrected = FALSE)
chao1 <- function (counts, bias_corrected = TRUE, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_CHAO1, counts, margin, norm, cpus, bias_corrected)
}


#' Faith's Phylogenetic Diversity (PD)
#' 
#' Calculates the sum of the branch lengths for all species present in a sample.
#' 
#' @inherit documentation
#' @inherit adiv_binary_normalized
#' 
#' @family Phylogenetic metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
  
  norm <- 'none'
  validate_args()
  
  .Call(C_alpha_div, ADIV_FAITH, counts, margin, norm, cpus, tree)
}


#' Fisher's Alpha
#' 
#' A parametric diversity index assuming species abundances follow a log-series distribution.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' 
#' The precision (number of decimal places) to use when solving the equation.
#' 
#' @references
#' Fisher, R. A., Corbet, A. S., & Williams, C. B. (1943). The relation between the number of species and the number of individuals in a random sample of an animal population. *Journal of Animal Ecology*, 12, 42-58. \doi{10.2307/1411}
#' 
#' @export
#' @examples
#'     fisher(ex_counts)
fisher <- function (counts, digits = 3L, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_FISHER, counts, margin, norm, cpus, digits)
}


#' Inverse Simpson Index
#' 
#' A transformation of the Simpson index that represents the "effective number of species".
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The Inverse Simpson index is defined as:
#' \deqn{1 / \sum_{i = 1}^{n} P_i^2}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
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
inv_simpson <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_INV_SIMPSON, counts, margin, norm, cpus, NULL)
}


#' Margalef's Richness Index
#' 
#' A richness metric that normalizes the number of species by the log of the total sample size.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' x <- ex_counts[1,]
#' (sum(x > 0) - 1) / log(sum(x))
#' ```
#' 
#' @references
#' Margalef, R. (1958). Information theory in ecology. *General Systems*, 3, 36-71.
#' 
#' Gamito, S. (2010). Caution is needed when applying Margalef diversity index. *Ecological Indicators*, 10(2), 550-551. \doi{10.1016/j.ecolind.2009.07.006}
#' 
#' @export
#' @examples
#'     margalef(ex_counts)
margalef <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MARGALEF, counts, margin, norm, cpus, NULL)
}


#' McIntosh Index
#' 
#' A dominance index based on the Euclidean distance from the origin.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Dominance metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The McIntosh index is defined as:
#' \deqn{\frac{X_T - \sqrt{\sum_{i = 1}^{n} (X_i)^2}}{X_T - \sqrt{X_T}}}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' * \eqn{X_T} : Total of all counts.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
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
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MCINTOSH, counts, margin, norm, cpus, NULL)
}


#' Menhinick's Richness Index
#' 
#' A richness metric that normalizes the number of species by the square root of the total sample size.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' x <- ex_counts[1,]
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
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MENHINICK, counts, margin, norm, cpus, NULL)
}


#' Observed Features
#' 
#' The count of unique features (richness) in a sample.
#' 
#' @inherit documentation
#' @inherit adiv_binary_normalized
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' Observed features is defined simply as the number of features with non-zero abundance:
#' \deqn{n}
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' sum(x > 0)
#' ```
#' 
#' @export
#' @examples
#'     observed(ex_counts)
observed <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  
  .Call(C_alpha_div, ADIV_OBSERVED, counts, margin, norm, cpus, NULL)
}


#' Shannon Diversity Index
#' 
#' A commonly used diversity index accounting for both abundance and evenness.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The Shannon index (entropy) is defined as:
#' \deqn{-\sum_{i = 1}^{n} P_i \times \ln(P_i)}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' p <- x / sum(x)
#' -sum(p * log(p))
#' ```
#' 
#' @references
#' Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27, 379-423.
#' 
#' Shannon, C. E., & Weaver, W. (1949). *The Mathematical Theory of Communication*. University of Illinois Press.
#'
#' @export
#' @examples
#'     shannon(ex_counts)
shannon <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_SHANNON, counts, margin, norm, cpus, NULL)
}


#' Gini-Simpson Index
#' 
#' The probability that two entities taken at random from the dataset represent different types.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The Gini-Simpson index is defined as:
#' \deqn{1 - \sum_{i = 1}^{n} P_i^2}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
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
simpson <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_SIMPSON, counts, margin, norm, cpus, NULL)
}


#' Squares Richness Estimator
#' 
#' A richness estimator based on the concept of "squares" (counts of species observed once or twice).
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Richness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
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
#' x  <- ex_counts[1,]
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
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_SQUARES, counts, margin, norm, cpus, NULL)
}


#' Simpson's Dominance Index
#' 
#' The probability that two entities taken at random (with replacement) are
#' of the same type.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Dominance metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' Simpson's dominance index \eqn{D} is defined as:
#' \deqn{\sum_{i = 1}^{n} P_i^2}
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' It is `1 - simpson()` (the Gini-Simpson index) and `1 / inv_simpson()`.
#' scikit-bio calls it `dominance()` and `simpson_d()`; both names are
#' aliases in [alpha_div()].
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' p <- x / sum(x)
#' sum(p ** 2)
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `dominance()` is the same index.
#' 
#' @references
#' Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688. \doi{10.1038/163688a0}
#' 
#' @export
#' @examples
#'     dominance(ex_counts)
dominance <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_DOMINANCE, counts, margin, norm, cpus, NULL)
}


#' Gini Index
#' 
#' The inequality of the features' abundances: 0 when all features of the
#' table are equally abundant, approaching 1 when one feature has them all.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Dominance metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' The Gini index is computed from the Lorenz curve of all \eqn{K} features
#' of the table, by rectangles, as scikit-bio's `gini_index()`:
#' \deqn{\max\left(0,\ 1 - \frac{2}{K} \sum_{j = 1}^{K} L_j\right)}
#' 
#' Where:
#' * \eqn{K} : The number of features in the table, including those absent
#'   from the sample.
#' * \eqn{L_j} : The Lorenz curve: the share of the sample's total held by
#'   the \eqn{j} least abundant features, \eqn{\sum_{i \le j} X_{(i)} / X_T},
#'   with \eqn{X_{(1)} \le \ldots \le X_{(K)}} the sorted abundances.
#' 
#' **The Gini index depends on how many features the table has,** including
#' features the sample lacks: each absent feature makes the sample look less
#' even. Filtering features out of a table (or adding all-zero ones)
#' changes it, so compare samples only from the same table. A sample of one
#' table with \eqn{K} features, all equally abundant, has index 0; a sample
#' with all its counts in one feature has index \eqn{1 - 2/K} (by
#' rectangles, the Lorenz curve's area is overestimated by \eqn{1/(2K)}).
#' 
#' `alpha_div(counts, 'gini')` is this index; the Gini-Simpson index is
#' `'simpson'`.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' L <- cumsum(sort(x)) / sum(x)
#' max(0, 1 - 2 * sum(L) / length(x))
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `gini_index()` is the same index (with its default
#' `method = "rectangles"`).
#' 
#' @references
#' Gini, C. (1912). *Variabilità e mutabilità*. Cuppini, Bologna.
#' 
#' Wittebolle, L., Marzorati, M., Clement, L., Balloi, A., Daffonchio, D., Heylen, K., De Vos, P., Verstraete, W., & Boon, N. (2009). Initial community evenness favours functionality under selective stress. *Nature*, 458, 623-626. \doi{10.1038/nature07840}
#' 
#' @export
#' @examples
#'     gini(ex_counts)
gini <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  
  .Call(C_alpha_div, ADIV_GINI, counts, margin, norm, cpus, NULL)
}


#' Good's Coverage
#' 
#' The estimated fraction of a community's individuals that belong to the
#' features observed in the sample: a check that sequencing was deep
#' enough.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Coverage metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' Good's coverage is defined as:
#' \deqn{1 - \frac{F_1}{X_T}}
#' 
#' Where:
#' * \eqn{F_1} : Number of features observed once (singletons).
#' * \eqn{X_T} : Total of all counts (sequencing depth).
#' 
#' It is Good's (1953) estimate of the sample coverage: the probability
#' that the next individual drawn belongs to a feature already seen. A
#' sample of only singletons has coverage 0. Chao & Jost (2012) give a less
#' biased estimator, which uses doubletons too.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' 1 - sum(x == 1) / sum(x)
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `goods_coverage()` is the same estimate.
#' 
#' @references
#' Good, I. J. (1953). The population frequencies of species and the estimation of population parameters. *Biometrika*, 40, 237-264. \doi{10.1093/biomet/40.3-4.237}
#' 
#' Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. *Ecology*, 93, 2533-2547. \doi{10.1890/11-1952.1}
#' 
#' @export
#' @examples
#'     goods_coverage(ex_counts)
goods_coverage <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_GOODS, counts, margin, norm, cpus, NULL)
}


#' Heip's Evenness
#' 
#' An evenness index derived from the exponential of Shannon's entropy.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Evenness metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' Heip's evenness is defined as:
#' \deqn{\frac{e^H - 1}{n - 1}}
#' 
#' Where:
#' * \eqn{H} : Shannon's entropy, `shannon()`.
#' * \eqn{n} : The number of observed features.
#' 
#' A sample with one feature has evenness 1, as in scikit-bio.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' p <- x[x > 0] / sum(x)
#' (exp(-sum(p * log(p))) - 1) / (length(p) - 1)
#' ```
#' 
#' @references
#' Heip, C. (1974). A new index measuring evenness. *Journal of the Marine Biological Association of the UK*, 54, 555-557. \doi{10.1017/S0025315400022736}
#' 
#' @export
#' @examples
#'     heip(ex_counts)
heip <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_HEIP, counts, margin, norm, cpus, NULL)
}


#' Hill Numbers
#' 
#' The effective number of species of order \eqn{q}: how many equally
#' abundant features would give the same diversity.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `renyi()`, `tsallis()`, `vignette('adiv')`
#' 
#' @details
#' The Hill number of order \eqn{q} is defined as:
#' \deqn{{}^qD = \left(\sum_{i = 1}^{n} P_i^q\right)^{1/(1-q)}}
#' 
#' with its limits \eqn{{}^1D = e^H} (\eqn{H}: Shannon's entropy) and
#' \eqn{{}^\infty D = 1 / \max(P_i)}.
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' The order sets how much weight abundant features get. The Hill numbers
#' are a profile that unifies the common indices:
#' * \eqn{q = 0}: richness, `observed()`.
#' * \eqn{q = 1}: the exponential of Shannon's entropy, `exp(shannon())`,
#'   which weighs every feature by its frequency (Jost 2006).
#' * \eqn{q = 2}: the inverse Simpson index, `inv_simpson()`.
#' * \eqn{q = \infty}: the inverse Berger-Parker index, `1 / berger()`.
#' 
#' The default order, 2, is scikit-bio's, which makes `hill()` equal to
#' `inv_simpson()`; give `q` for any other order.
#' [renyi()] is `log(hill())`, and [tsallis()] is another function of it.
#' 
#' Observed Hill numbers grow with sequencing depth, most for small
#' \eqn{q}; compare samples at equal depth, for instance with
#' [alpha_div_rarefied()].
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' q <- 2
#' p <- x[x > 0] / sum(x)
#' sum(p ** q) ** (1 / (1 - q))
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `hill()` is the same (its `order` is `q`).
#' `vegan::renyi(x, scales, hill = TRUE)` gives a profile of Hill numbers
#' over several orders.
#' 
#' @references
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. *Ecology*, 54, 427-432. \doi{10.2307/1934352}
#' 
#' Jost, L. (2006). Entropy and diversity. *Oikos*, 113, 363-375. \doi{10.1111/j.2006.0030-1299.14714.x}
#' 
#' @export
#' @examples
#'     hill(ex_counts)
#'     
#'     # Richness, exp(Shannon), inverse Simpson
#'     sapply(c(0, 1, 2), function (q) hill(ex_counts, q = q))
hill <- function (counts, q = 2, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_HILL, counts, margin, norm, cpus, as.double(q))
}


#' McIntosh's Evenness
#' 
#' An evenness index derived from the McIntosh index.
#' 
#' @inherit documentation
#' @inherit adiv_assert_integer
#' 
#' @family Evenness metrics
#' @seealso `alpha_div()`, `mcintosh()`, `vignette('adiv')`
#' 
#' @details
#' McIntosh's evenness is defined as:
#' \deqn{\frac{\sqrt{\sum_{i = 1}^{n} (X_i)^2}}{\sqrt{(X_T - n + 1)^2 + n - 1}}}
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{X_i} : Integer count of the \eqn{i}-th feature.
#' * \eqn{X_T} : Total of all counts.
#' 
#' This follows scikit-bio's `mcintosh_e()` (Heip & Engels 1974).
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' n <- sum(x > 0)
#' sqrt(sum(x^2)) / sqrt((sum(x) - n + 1)^2 + n - 1)
#' ```
#' 
#' @references
#' Heip, C., & Engels, P. (1974). Comparing species diversity and evenness indices. *Journal of the Marine Biological Association of the UK*, 54, 559-563. \doi{10.1017/S0025315400022748}
#' 
#' @export
#' @examples
#'     mcintosh_evenness(ex_counts)
mcintosh_evenness <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  assert_integer_counts()
  
  .Call(C_alpha_div, ADIV_MCINTOSH_E, counts, margin, norm, cpus, NULL)
}


#' Pielou's Evenness
#' 
#' Shannon's entropy relative to its largest possible value for the number
#' of observed features.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Evenness metrics
#' @seealso `alpha_div()`, `shannon()`, `vignette('adiv')`
#' 
#' @details
#' Pielou's evenness \eqn{J} is defined as:
#' \deqn{\frac{H}{\ln n}}
#' 
#' Where:
#' * \eqn{H} : Shannon's entropy, `shannon()`.
#' * \eqn{n} : The number of observed features.
#' 
#' It is one of the alpha diversity outputs of QIIME 2's `core-metrics`
#' (`evenness_vector`). It doesn't depend on the base of the logarithm.
#' 
#' A sample with one feature has evenness 1, as in scikit-bio (where
#' `vegan::diversity(x) / log(vegan::specnumber(x))` gives `NaN`).
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' p <- x[x > 0] / sum(x)
#' -sum(p * log(p)) / log(length(p))
#' ```
#' 
#' @references
#' Pielou, E. C. (1966). The measurement of diversity in different types of biological collections. *Journal of Theoretical Biology*, 13, 131-144. \doi{10.1016/0022-5193(66)90013-0}
#' 
#' @export
#' @examples
#'     pielou(ex_counts)
pielou <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_PIELOU, counts, margin, norm, cpus, NULL)
}


#' Rényi Entropy
#' 
#' The generalized entropy of order \eqn{q}: the logarithm of the Hill
#' number of the same order.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `hill()`, `tsallis()`, `vignette('adiv')`
#' 
#' @details
#' Rényi's entropy of order \eqn{q} is defined as:
#' \deqn{{}^qH = \frac{1}{1 - q} \ln \sum_{i = 1}^{n} P_i^q}
#' 
#' with its limits \eqn{{}^1H = H} (Shannon's entropy) and
#' \eqn{{}^\infty H = -\ln \max(P_i)}; \eqn{{}^0H} is the log of the number
#' of observed features.
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' It is `log(hill(counts, q))`. A sample with one feature has entropy 0.
#' Logarithms are natural; for another base \eqn{b}, divide by `log(b)`.
#' The default order, 2, is scikit-bio's.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' q <- 2
#' p <- x[x > 0] / sum(x)
#' log(sum(p ** q)) / (1 - q)
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `renyi()` is the same (its `order` is `q`).
#' `vegan::renyi()` gives a profile over several orders (`scales`).
#' 
#' @references
#' Rényi, A. (1961). On measures of entropy and information. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 547-561.
#' 
#' @export
#' @examples
#'     renyi(ex_counts)
#'     renyi(ex_counts, q = 1) # Shannon's entropy
renyi <- function (counts, q = 2, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_RENYI, counts, margin, norm, cpus, as.double(q))
}


#' Simpson's Evenness
#' 
#' The inverse Simpson index relative to the number of observed features.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Evenness metrics
#' @seealso `alpha_div()`, `inv_simpson()`, `vignette('adiv')`
#' 
#' @details
#' Simpson's evenness is defined as:
#' \deqn{\frac{1}{n \sum_{i = 1}^{n} P_i^2}}
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' It is `inv_simpson() / observed()`.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' p <- x[x > 0] / sum(x)
#' 1 / (length(p) * sum(p ** 2))
#' ```
#' 
#' @references
#' Smith, B., & Wilson, J. B. (1996). A consumer's guide to evenness indices. *Oikos*, 76, 70-82. \doi{10.2307/3545749}
#' 
#' @export
#' @examples
#'     simpson_evenness(ex_counts)
simpson_evenness <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_SIMPSON_E, counts, margin, norm, cpus, NULL)
}


#' Strong's Dominance Index
#' 
#' The largest gap between the cumulative share of the most abundant
#' features and the share they would have if all were equally abundant.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Dominance metrics
#' @seealso `alpha_div()`, `vignette('adiv')`
#' 
#' @details
#' Strong's dominance index is defined as:
#' \deqn{\max_{i} \left( \frac{c_i}{X_T} - \frac{i}{n} \right)}
#' 
#' Where:
#' * \eqn{n} : The number of observed features.
#' * \eqn{c_i} : The sum of the \eqn{i} largest counts.
#' * \eqn{X_T} : Total of all counts.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' x <- sort(x[x > 0], decreasing = TRUE)
#' max(cumsum(x) / sum(x) - seq_along(x) / length(x))
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `strong()` is the same index.
#' 
#' @references
#' Strong, W. L. (2002). Assessing species abundance unevenness within and between plant communities. *Community Ecology*, 3, 237-246. \doi{10.1556/ComEc.3.2002.2.9}
#' 
#' @export
#' @examples
#'     strong(ex_counts)
strong <- function (counts, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'none'
  validate_args()
  
  .Call(C_alpha_div, ADIV_STRONG, counts, margin, norm, cpus, NULL)
}


#' Tsallis Entropy
#' 
#' The non-additive generalized entropy of order \eqn{q}.
#' 
#' @inherit documentation
#' @inherit adiv_percent_normalized
#' 
#' @family Diversity metrics
#' @seealso `alpha_div()`, `hill()`, `renyi()`, `vignette('adiv')`
#' 
#' @details
#' Tsallis's entropy of order \eqn{q} is defined as:
#' \deqn{{}^qS = \frac{1}{q - 1} \left(1 - \sum_{i = 1}^{n} P_i^q\right)}
#' 
#' with its limits \eqn{{}^1S = H} (Shannon's entropy) and
#' \eqn{{}^\infty S = 0}.
#' 
#' Where:
#' * \eqn{n} : The number of features.
#' * \eqn{P_i} : Proportional abundance of the \eqn{i}-th feature.
#' 
#' At \eqn{q = 0} it is the number of observed features minus 1, and at
#' \eqn{q = 2} the Gini-Simpson index, `simpson()`. A sample with one
#' feature has entropy 0. The default order, 2, is scikit-bio's.
#' 
#' **Base R Equivalent:**
#' ```r
#' x <- ex_counts[1,]
#' q <- 2
#' p <- x[x > 0] / sum(x)
#' (1 - sum(p ** q)) / (q - 1)
#' ```
#' 
#' @section Same names in other packages:
#' scikit-bio's `tsallis()` is the same (its `order` is `q`).
#' 
#' @references
#' Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics. *Journal of Statistical Physics*, 52, 479-487. \doi{10.1007/BF01016429}
#' 
#' @export
#' @examples
#'     tsallis(ex_counts)
tsallis <- function (counts, q = 2, margin = 1L, cpus = n_cpus()) {
  
  norm <- 'percent'
  validate_args()
  
  .Call(C_alpha_div, ADIV_TSALLIS, counts, margin, norm, cpus, as.double(q))
}
