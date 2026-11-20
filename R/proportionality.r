# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# vlr() follows vlr() and pairwise_vlr() from scikit-bio 0.7.4
# (skbio/stats/composition/_base.py), Copyright (c) 2013--, scikit-bio
# development team, under the BSD-3 license given in src/proportionality.c
# and inst/COPYRIGHTS. rho_p() and phi_s() are formulas (Erb & Notredame
# 2016; Quinn et al. 2017).



#' Variance of Log Ratios
#'
#' How far each pair of features is from being proportional: the variance,
#' across samples, of the log of their ratio. `0` means the two features
#' keep exactly the same ratio in every sample.
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of [replace_zeros()]
#'        are all accepted: each sample's scale cancels out of the log
#'        ratios. Pairs are taken between the **features** (the columns when
#'        `margin = 1`), and the variance is over the samples.
#'
#' @param pseudocount  Added to every value before taking logs. `NULL` (the
#'        default) adds nothing when `counts` has no zeros, and otherwise
#'        uses half the smallest non-zero value, with a warning. `0` adds
#'        nothing, and is an error when `counts` has zeros. See **Zeros**.
#'
#' @return A `dist` object **over the features** (the packed lower triangle,
#'   as [beta_div()] returns over samples), labelled with the feature names,
#'   of class `c("ecodive_dist", "dist")`, so that `plot()` draws it as a
#'   heatmap (see [plot.ecodive_dist()]). Besides `dist`'s `Size`, `Labels`,
#'   `Diag` and `Upper`, it has attributes `method` (`"vlr"`) and
#'   `display`, with its title (`"Variance of Log Ratios"`; see
#'   [ecodive_display]), as [beta_div()]'s results do, and with `clust`,
#'   `order`.
#'
#'   It is a dissimilarity, so it can be given as `dm` to [pcoa()],
#'   [dist_ord()] or [stats::hclust()] to ordinate or cluster the features.
#'   Use `as.matrix()` for the full \eqn{D \times D} matrix.
#'
#' @details
#' For features \eqn{i} and \eqn{j} measured in \eqn{n} samples,
#' \deqn{\mathrm{vlr}_{ij} = \mathrm{var}\left(\log\frac{x_i + c}{x_j +
#' c}\right) = \frac{1}{n - 1} \sum_{s=1}^{n} \left(\ell_{si} - \ell_{sj} -
#' \overline{\ell_i} + \overline{\ell_j}\right)^2}
#' with \eqn{\ell = \log(x + c)}, \eqn{c} the pseudocount, natural logs, and
#' the sample variance's divisor \eqn{n - 1} (scikit-bio's default `ddof =
#' 1`). It doesn't depend on the samples' totals, nor on any other feature,
#' and so equally the variance of the difference of the two features' clr
#' or alr coordinates. It is unbounded and grows with the features' own
#' variability; [rho_p()] and [phi_s()] rescale it.
#'
#' **Why not correlation?** Relative abundances are constrained to sum to a
#' constant, so their correlations are not those of the underlying
#' abundances: closure alone induces negative correlations, and the
#' correlation of two taxa changes when a third is dropped from the table
#' (Pearson 1897; Lovell et al. 2015). The ratio of two features is the same
#' whether their abundances are absolute or relative, so a constant ratio
#' (proportionality) is a property of the features themselves, and holds
#' after subsetting.
#'
#' Results match scikit-bio's `pairwise_vlr()` (with features as its
#' columns) and `vlr()`. scikit-bio computes the variance from a covariance
#' matrix (\eqn{c_{ii} + c_{jj} - 2c_{ij}}), which loses precision for
#' nearly proportional pairs; ecodive sums the differences directly. Its
#' `vlr(robust = TRUE)`, which leaves out the samples where either feature
#' is zero (so each pair uses different samples), is not offered.
#'
#' @section Zeros:
#'   Logs of zero are undefined, so zeros need replacing first. With
#'   `pseudocount = NULL` (the default), no pseudocount is used when
#'   `counts` has no zeros; otherwise half the smallest non-zero value is
#'   added to every value, with a warning. The recommended input is a table
#'   with its zeros replaced by [replace_zeros()], which has none and is
#'   used as it is.
#'
#'   propr (Quinn et al. 2017) replaces only the zeros, by 1, which is not
#'   `pseudocount = 1` (that adds 1 to every value). To reproduce it, write
#'   `x[x == 0] <- 1` and then call these functions on `x`.
#'
#'   Pseudocounts are not neutral: in sparse tables the many zeros, all
#'   given the same value, can dominate the variances. Filtering out rare
#'   features first helps.
#'
#' @section Size:
#'   The result has \eqn{D(D - 1)/2} values: 20,000 features give 2e8 pairs,
#'   1.6 GB. Filter out rare features first for large tables. The table
#'   itself stays sparse.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman & Hall.
#'
#' Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., & Bähler,
#' J. (2015). Proportionality: a valid alternative to correlation for
#' relative data. *PLoS Computational Biology*, 11, e1004075.
#'
#' Pearson, K. (1897). Mathematical contributions to the theory of
#' evolution. On a form of spurious correlation which may arise when indices
#' are used in the measurement of organs. *Proceedings of the Royal Society
#' of London*, 60, 489–498.
#'
#' @seealso [rho_p()] and [phi_s()], scaled versions; [log_ratio()] and
#'   [balance()] for one log ratio per sample; [replace_zeros()].
#'
#' @aliases proportionality
#' @export
#' @examples
#'     # Zeros replaced first
#'     x <- replace_zeros(ex_counts, 'mult')
#'     vlr(x)
#'
#'     # The same as the variance of each log ratio
#'     var(log(x[, 'Streptococcus'] / x[, 'Bacteroides']))
#'
#'     # Cluster the features
#'     plot(hclust(vlr(x)))
#'
vlr <- function (
    counts, margin = 1L, pseudocount = NULL, clust = NULL, cpus = n_cpus() ) {

  proportionality('vlr', 'the variance of log ratios', denominator = NULL)
}



#' Proportionality Correlation
#'
#' The proportionality coefficient \eqn{\rho_p} of Erb & Notredame (2016),
#' propr's default `"rho"`: a correlation-like measure of how close each
#' pair of features is to keeping a constant ratio. `1` means perfectly
#' proportional.
#'
#' @inheritParams vlr
#'
#' @param denominator  The reference that each feature's log abundance is
#'        taken relative to, for the variances \eqn{\mathrm{var}(z_i)}.
#'        `NULL` (the default) uses the centered log ratio (clr): each
#'        sample's geometric mean of all features, as propr's default `ivar
#'        = "clr"`. Otherwise features, as names, (1-based) indices, or a
#'        logical vector with one value per feature: their geometric mean.
#'        One feature gives the additive log ratio (alr, as
#'        `alr(denominator =)`), and that feature is then left out of the
#'        result (its alr is 0 in every sample); with several, every feature
#'        is kept.
#'
#' @return A `dist` object over the features, as [vlr()] returns, with
#'   `method` `"rho_p"`.
#'
#'   \eqn{\rho_p} is a **similarity** in a `dist` container: `as.matrix()`
#'   puts 0 on the diagonal, where \eqn{\rho_p} is 1, and it must not be
#'   given as `dm` to an ordination or test. Threshold it to find
#'   proportional pairs (`which(rho > 0.7)`, or `as.matrix(rho) > 0.7` for
#'   an adjacency matrix), or use `1 - rho_p(x)` (which is \eqn{\mathrm{vlr}_{ij}
#'   / (\mathrm{var}(z_i) + \mathrm{var}(z_j))}) for a dissimilarity.
#'
#' @details
#' With \eqn{z} the log-ratio-transformed abundances (clr or alr),
#' \deqn{\rho_{p,ij} = 1 - \frac{\mathrm{var}(z_i - z_j)}{\mathrm{var}(z_i)
#' + \mathrm{var}(z_j)} = \frac{2\,\mathrm{cov}(z_i, z_j)}{\mathrm{var}(z_i)
#' + \mathrm{var}(z_j)}}
#' where \eqn{\mathrm{var}(z_i - z_j)} is [vlr()]. It lies in \eqn{[-1,
#' 1]}: 1 for proportional features, 0 for none, -1 for reciprocal ones
#' (\eqn{z_i = -z_j}). It is Lin's concordance correlation coefficient of
#' \eqn{z_i} and \eqn{z_j}. Unlike [vlr()], it depends on the `denominator`.
#' A pair of features that are both constant relative to the reference has
#' no defined \eqn{\rho_p} (`NaN`).
#'
#' **Not a correlation of abundances.** \eqn{\rho_p} asks whether two
#' features keep a constant ratio: it is 1 only for a log ratio with slope 1,
#' and pairs related by a power law with another exponent, or with unequal
#' variances, score low. [sparcc()] starts from the same variance of log
#' ratios but divides by estimated variances of the log *absolute*
#' abundances instead, to estimate their Pearson correlation; that is
#' possible only under its assumption that most pairs are uncorrelated.
#'
#' This is the coefficient that propr computes with `propr(counts, metric =
#' "rho")` (see **Zeros** for its zero handling). propr also estimates false discovery
#' rates by permutation (`updateCutoffs()`), and offers `ivar = "iqlr"`, a
#' clr over the features of mid-range variance; neither is offered here.
#'
#' @inheritSection vlr Zeros
#' @inheritSection vlr Size
#' @inheritSection documentation Input Types
#'
#' @references
#' Erb, I., & Notredame, C. (2016). How should we measure proportionality on
#' relative gene expression data? *Theory in Biosciences*, 135, 21–36.
#'
#' Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., & Bähler,
#' J. (2015). Proportionality: a valid alternative to correlation for
#' relative data. *PLoS Computational Biology*, 11, e1004075.
#'
#' Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
#' propr: an R-package for identifying proportionally abundant features
#' using compositional data analysis. *Scientific Reports*, 7, 16252.
#'
#' @seealso [vlr()], [phi_s()], [sparcc()], [alr()], [normalize()]
#'   (`'clr'`), [replace_zeros()].
#'
#' @aliases propr
#' @export
#' @examples
#'     x   <- replace_zeros(ex_counts, 'mult')
#'     rho <- rho_p(x)
#'     rho
#'
#'     # Pairs above a cutoff
#'     m <- as.matrix(rho)
#'     which(m > 0.5 & upper.tri(m), arr.ind = TRUE)
#'
#'     # Relative to one feature (alr): that feature is left out
#'     rho_p(x, denominator = 'Streptococcus')
#'
rho_p <- function (
    counts, margin = 1L, denominator = NULL, pseudocount = NULL,
    clust = NULL, cpus = n_cpus() ) {

  proportionality('rho_p', 'proportionality', denominator)
}



#' Symmetric Proportionality
#'
#' The symmetric proportionality measure \eqn{\phi_s} of Quinn et al.
#' (2017), propr's `"phs"`: the variance of the log ratio of two features
#' relative to that of their log product. `0` means perfectly proportional.
#'
#' @inheritParams rho_p
#'
#' @return A `dist` object over the features, as [rho_p()] returns, with
#'   `method` `"phi_s"`. It is a dissimilarity, so it can be given as `dm`
#'   to [pcoa()], [dist_ord()] or [stats::hclust()].
#'
#' @details
#' With \eqn{z} the log-ratio-transformed abundances (clr or alr),
#' \deqn{\phi_{s,ij} = \frac{\mathrm{var}(z_i - z_j)}{\mathrm{var}(z_i +
#' z_j)} = \frac{1 - \rho_{p,ij}}{1 + \rho_{p,ij}}}
#' where \eqn{\mathrm{var}(z_i - z_j)} is [vlr()] and \eqn{\rho_p} is
#' [rho_p()]. It lies in \eqn{[0, \infty)}: 0 for proportional features, 1
#' for unrelated ones, and infinite for reciprocal ones (\eqn{z_i = -z_j}).
#' It is computed as \eqn{\mathrm{vlr}_{ij} / (2\mathrm{var}(z_i) +
#' 2\mathrm{var}(z_j) - \mathrm{vlr}_{ij})}, since \eqn{\mathrm{var}(a + b) =
#' 2\mathrm{var}(a) + 2\mathrm{var}(b) - \mathrm{var}(a - b)}. Unlike
#' [vlr()], it depends on the `denominator`. Near-reciprocal pairs (such as
#' the two features of a two-feature `denominator`, whose \eqn{z} are exact
#' opposites) have \eqn{\mathrm{var}(z_i + z_j)} near 0, where rounding
#' dominates: their value is `Inf` or merely very large.
#'
#' Lovell et al.'s (2015) original \eqn{\phi = \mathrm{vlr}_{ij} /
#' \mathrm{var}(z_i)} is not symmetric (\eqn{\phi_{ij} \ne \phi_{ji}}), so
#' it can't be a `dist`; Quinn et al. (2017) replaced it with \eqn{\phi_s}.
#' This is the measure that propr computes with `propr(counts, metric =
#' "phs")` (see **Zeros** for its zero handling).
#'
#' @inheritSection vlr Zeros
#' @inheritSection vlr Size
#' @inheritSection documentation Input Types
#'
#' @references
#' Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., & Bähler,
#' J. (2015). Proportionality: a valid alternative to correlation for
#' relative data. *PLoS Computational Biology*, 11, e1004075.
#'
#' Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
#' propr: an R-package for identifying proportionally abundant features
#' using compositional data analysis. *Scientific Reports*, 7, 16252.
#'
#' @seealso [vlr()], [rho_p()], [replace_zeros()].
#'
#' @export
#' @examples
#'     x <- replace_zeros(ex_counts, 'mult')
#'     phi_s(x)
#'
#'     # The same as (1 - rho_p) / (1 + rho_p)
#'     rho <- rho_p(x)
#'     all.equal(c(phi_s(x)), c((1 - rho) / (1 + rho)))
#'
phi_s <- function (
    counts, margin = 1L, denominator = NULL, pseudocount = NULL,
    clust = NULL, cpus = n_cpus() ) {

  proportionality('phi_s', 'proportionality', denominator)
}



# The shared body of vlr(), rho_p() and phi_s(), run in the caller's frame
# (for counts, margin, pseudocount and cpus).
proportionality <- function (method, what, denominator, env = parent.frame()) {

  validate_clust(env)
  info <- log_ratio_input(what, env = env, zero_ok = TRUE, dense_out = FALSE)
  counts      <- env$counts
  margin      <- env$margin
  pseudocount <- env$pseudocount
  cpus        <- env$cpus

  if (info$n < 2L)
    stop('`counts` must have at least 2 samples for ', what, '.', call. = FALSE)

  D    <- info$D
  keep <- seq_len(D)
  ref  <- NULL
  if (!is.null(denominator)) {
    ref <- feature_set(denominator, info$features, D, 'denominator')
    if (length(ref) == 1L) {
      keep <- keep[-ref]
      if (length(keep) < 2L)
        stop(
          '`counts` must have at least 2 features besides the `denominator` ',
          'for ', what, '.', call. = FALSE )
    }
  }

  labels <- info$features
  if (!is.null(labels))           labels <- labels[keep]
  else if (length(keep) < D)      labels <- as.character(keep)

  res <- .Call(
    C_proportionality, counts, margin,
    switch(method, vlr = 1L, rho_p = 2L, phi_s = 3L),
    if (!is.null(ref)) ref - 1L, keep - 1L, labels, pseudocount, cpus )

  attr(res, 'method') <- method
  res <- set_display(res, title = FEATURE_ASSOCS$name[match(method, FEATURE_ASSOCS$id)])
  dist_clust(res, env$clust, cpus)
}
