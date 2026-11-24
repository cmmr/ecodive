# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# batch_mmuphin(): an R and C port (with src/batch_mmuphin.c) of
# adjust_batch() and its helpers from MMUPHin 2.0.0 (R/adjust_batch.R,
# R/helpers_adjust_batch.R, R/helpers.R; identical in MMUPHin 1.16.0),
# restructured for ecodive, see plans/batch_correction.md. MMUPHin is
# Copyright (c) 2019 Harvard School of Public Health, under the MIT license
# given in src/batch_mmuphin.c and inst/COPYRIGHTS.



#' MMUPHin Batch Correction
#'
#' Removes batch effects (study, sequencing run, extraction kit) from a
#' count or relative abundance table, while keeping the effects of named
#' covariates: MMUPHin's zero-inflated empirical Bayes adjustment (Ma et al.
#' 2022), a port of `MMUPHin::adjust_batch()`. Zeros are left out of the
#' fits and stay zero.
#'
#' @inheritParams documentation
#'
#' @param counts  A table of non-negative abundances (samples \eqn{\times}
#'        features): counts, proportions or other non-negative values, or
#'        any of the objects listed under **Input Types**.
#'
#' @param batch  The batch of each sample: one value per sample (matched to
#'        the samples by name when both are named), or the name of a column
#'        of `data`. At least two batches, each of at least two samples.
#'
#' @param formula  A one-sided formula over the columns of `data` naming the
#'        covariates whose effects are **kept**, such as `~ disease + age`.
#'        Default `NULL`: none. Covariates left out have their differences
#'        between batches removed with the batch effects.
#'
#' @param data  A data.frame of sample metadata, with the columns that
#'        `batch` and `formula` name. Rows are matched to the samples by
#'        name, or else by position.
#'
#' @param zero_inflation  Leave zeros out of every fit (`TRUE`, the
#'        default), so that a feature is adjusted only in the batches where
#'        it is present, and its zeros stay zero. `FALSE` fits every value
#'        (MMUPHin's ComBat-like mode; see "Differences from ComBat").
#'
#' @param pseudocount  Added to every value before the log transform.
#'        Default `NULL`: half the smallest non-zero value (0.5 for counts),
#'        as MMUPHin. With `zero_inflation = TRUE`, `0` is allowed: zeros
#'        never enter the fits.
#'
#' @param tol  The convergence tolerance of the shrinkage iteration
#'        (MMUPHin's `conv`). Default `1e-4`.
#'
#' @param max_iters  The most iterations of the shrinkage (MMUPHin's
#'        `maxit`); more is an error. Default `1000`.
#'
#' @return The corrected table: a dense numeric `matrix` of class
#'   `ecodive_batch` in the orientation of `counts`, with its names. Each
#'   sample keeps its total, zeros stay zero, and a table of whole numbers
#'   (that is not all \eqn{\le 1}) gives whole numbers (rounded) back.
#'   Attributes:
#'   \describe{
#'     \item{`method`}{`"mmuphin"`.}
#'     \item{`batch_effects`}{A data frame with a row per feature and batch:
#'       `feature`, `batch`, `gamma_hat` and `delta_hat` (the feature's
#'       location and scale effects in the batch, on the standardised log
#'       scale), `gamma_star` and `delta_star` (after shrinkage). `NA` where
#'       the feature was not adjusted in the batch.}
#'     \item{`adjusted`}{A logical per feature: adjusted in any batch.
#'       Features that aren't are only re-closed, so they move slightly.}
#'     \item{`pseudocount`}{The pseudocount chosen, when
#'       `pseudocount = NULL`.}
#'     \item{`n_iters`}{The shrinkage iterations run in each batch.}
#'     \item{`batch_means`, `batch_means_adjusted`}{Feature \eqn{\times}
#'       batch matrices of mean relative abundances, before and after.}
#'     \item{`overall_means`}{Each feature's mean relative abundance over all
#'       samples, before.}
#'   }
#'   The class keeps [print()][print_ecodive] and
#'   [plot()][plot.ecodive_batch] short; the result is otherwise an
#'   ordinary matrix, accepted by every function that takes `counts`, and
#'   `x[i, j]` or `unclass(x)` drop the class and attributes. A sparse
#'   input gives a dense result; `Matrix::Matrix(x, sparse = TRUE)` makes it
#'   sparse again.
#'
#' @details
#' **The model.** Each feature's values are log2 relative abundances, \eqn{y
#' = \log_2((x + c) / \sum (x + c))} per sample, with pseudocount \eqn{c}.
#' With `zero_inflation`, only the feature's non-zero samples enter, and it
#' is adjusted in the batches present among them (at least two, with a
#' full-rank design of those batches and the covariates, and more samples
#' than columns); a feature present in one batch only is left as it is.
#' The feature is fitted by least squares on the batch indicators and the
#' covariates, and standardised:
#' \eqn{s = (y - \bar\alpha - X\beta) / \hat\sigma}, where \eqn{\bar\alpha}
#' is the mean of the fitted batch part and \eqn{\hat\sigma^2} the
#' residual variance. Each batch's location and scale,
#' \eqn{\hat\gamma} and \eqn{\hat\delta} (the mean and standard deviation of
#' \eqn{s} in the batch), are then shrunk towards their average over
#' features by empirical Bayes (normal and inverse-gamma priors, estimated
#' by moments), iterating until the largest relative change is below `tol`.
#' The adjusted value is \eqn{(s - \gamma^*) / \sqrt{\delta^*} \cdot
#' \hat\sigma + \bar\alpha + X\beta}: the covariates' effects are kept.
#' Finally \eqn{2^y}, with the original zeros put back, is re-closed to each
#' sample's original total.
#'
#' **What to correct for, and what with.** Covariates that differ between
#' batches (a study of patients next to one of controls) must be named in
#' `formula`, or their effects are removed with the batches'. A covariate
#' constant within each batch can't be kept: it is an error. Corrected
#' tables are for exploration and pooled analyses (ordination, clustering,
#' prediction). For differential abundance across studies, MMUPHin's authors
#' and Wang & Lê Cao (2020) advise modelling batch as a covariate, or a
#' meta-analysis of per-study results, rather than testing corrected values.
#' The correction is of log relative abundances, re-closed: a compositional
#' adjustment, which can move features that had no batch effect.
#'
#' **Checking the correction.** Ordinate the table before and after,
#' coloured by batch, and compare PERMANOVA's \eqn{R^2} for batch (see the
#' examples); `plot(x)` shows how much the shrinkage moved each batch
#' effect.
#'
#' @section Differences from ComBat:
#'
#' MMUPHin's adjustment with `zero_inflation = FALSE` is ComBat with
#' parametric priors (Johnson et al. 2007; [batch_combat()]) on log2
#' relative abundances, except for three details of MMUPHin's code, which
#' `batch_mmuphin()` keeps so that it matches MMUPHin:
#'
#' * \eqn{\hat\delta} is each batch's **standard deviation** of the
#'   standardised values, where ComBat uses the variance, so the
#'   inverse-gamma prior is fitted to standard deviations.
#' * The pooled variance divides by \eqn{n - 1} (ComBat: \eqn{n}).
#' * The pseudocount is added before the relative abundances are taken, so
#'   non-zero values move slightly too.
#'
#' @section Differences from MMUPHin:
#'
#' `batch_mmuphin()` ports MMUPHin 2.0.0's `adjust_batch()` (MIT license;
#' Siyuan Ma, Harvard School of Public Health), and its corrected tables
#' agree with it to rounding, except:
#'
#' * `counts` is samples \eqn{\times} features by default (`margin`);
#'   MMUPHin takes features \eqn{\times} samples.
#' * `formula` names the covariates, instead of a vector of column names
#'   (each its own term); interactions and transformations are allowed.
#' * `batch` and `data` follow ecodive's per-sample rules: `batch` can be
#'   values instead of a column name, and `data`'s rows are matched to the
#'   samples by name. Batches are the levels of a factor, else the sorted
#'   values (MMUPHin converts with a warning).
#' * Every batch needs at least two samples.
#' * Tables of non-negative values that are neither proportions nor whole
#'   numbers are accepted (and not rounded); MMUPHin stops.
#' * `pseudocount = 0` is allowed with `zero_inflation = TRUE`.
#' * The least-squares fits use a QR decomposition where MMUPHin solves the
#'   normal equations, so very ill-conditioned covariates may agree less
#'   closely.
#' * The result is the table, with the estimates kept as attributes; the
#'   diagnostic plot is drawn by `plot()` in base graphics instead of being
#'   written to a PDF file, and nothing is printed while it runs.
#' * A sparse input gives a dense result.
#'
#' **Same names in other packages:** none. MMUPHin's own function is
#' `adjust_batch()`.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Ma, S., Shungin, D., Mallick, H., Schirmer, M., Nguyen, L. H., Kolde, R.,
#' Franzosa, E., Vlamakis, H., Xavier, R., & Huttenhower, C. (2022).
#' Population structure discovery in meta-analyzed microbial communities
#' and inflammatory bowel disease using MMUPHin. *Genome Biology*, 23, 208.
#' \doi{10.1186/s13059-022-02753-4}
#'
#' Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects
#' in microarray expression data using empirical Bayes methods.
#' *Biostatistics*, 8, 118–127. \doi{10.1093/biostatistics/kxj037}
#'
#' Wang, Y., & Lê Cao, K.-A. (2020). Managing batch effects in microbiome
#' data. *Briefings in Bioinformatics*, 21, 1954–1970.
#' \doi{10.1093/bib/bbz105}
#'
#' @seealso [batch_combat()]; [plot.ecodive_batch()]; [permanova()] and
#'   [pcoa()] to check a correction.
#'
#' @export
#' @examples
#'     # Two runs of ten samples; the second reads ten features 3x higher.
#'     set.seed(1)
#'     counts <- matrix(
#'       rpois(20 * 30, 20), 20, 30,
#'       dimnames = list(paste0('S', 1:20), paste0('F', 1:30)) )
#'     batch <- rep(c('run1', 'run2'), each = 10)
#'     counts[11:20, 1:10] <- counts[11:20, 1:10] * 3
#'
#'     x <- batch_mmuphin(counts, batch)
#'     x
#'
#'     # Batch's share of the variation, before and after.
#'     permanova(bray(counts), batch)$statistic
#'     permanova(bray(x), batch)$statistic
#'
#'     plot(x)
#'
#'     # Keep a covariate's effect: name it in `formula`.
#'     meta <- data.frame(run = batch, sex = rep(c('F', 'M'), 10))
#'     batch_mmuphin(counts, 'run', ~ sex, meta)
batch_mmuphin <- function (
    counts, batch, formula = NULL, data = NULL, margin = 1L,
    zero_inflation = TRUE, pseudocount = NULL, tol = 1e-4, max_iters = 1000L,
    cpus = n_cpus() ) {

  if (missing(batch)) batch <- NULL
  info <- batch_input()

  zero_inflation <- batch_flag(zero_inflation, 'zero_inflation')
  it   <- batch_iter_args(tol, max_iters)
  vals <- info$vals

  pos <- vals[vals > 0]
  if (length(pos) == 0)
    stop('`counts` contains only zeros.', call. = FALSE)

  # check_feature_abd(): proportions, counts, or (ecodive only) other
  # non-negative values.
  whole <- !all(vals <= 1) && all(vals == floor(vals))

  chosen <- is.null(pseudocount)
  if (chosen) {
    pseudocount <- min(pos) / 2
  }
  else {
    if (!is.numeric(pseudocount) || length(pseudocount) != 1 || is.na(pseudocount) ||
        !is.finite(pseudocount) || pseudocount < 0)
      stop('`pseudocount` must be a single non-negative number.', call. = FALSE)
    pseudocount <- as.double(pseudocount)
    has_zeros <- info$sparse || any(vals == 0)
    if (pseudocount == 0 && !zero_inflation && has_zeros)
      stop(
        '`pseudocount = 0` needs `zero_inflation = TRUE`: otherwise the zeros ',
        'are fitted, and their log is undefined.', call. = FALSE )
  }

  design <- batch_design(batch, formula, data, info, min_size = 2L)
  B      <- length(design$levels)

  res <- .Call(
    C_batch_mmuphin, info$counts, info$margin, design$codes, design$X, B,
    zero_inflation, pseudocount, whole, it$tol, it$max_iters, cpus )

  lev <- design$levels
  switch(
    as.character(res$status),
    '1' = stop(
      'No feature can be adjusted. With `zero_inflation = TRUE`, a feature ',
      'needs non-zero values in at least two batches, and more of those ',
      'samples than batches and covariate columns, with a full-rank design ',
      '(MMUPHin: "All features are single-batch-specific").', call. = FALSE ),
    '2' = stop(
      'These batches can\'t be adjusted: they have fewer than two features ',
      'that are present in them and adjustable (non-zero in at least one ',
      'other batch, with a full-rank design). Remove them before the ',
      'correction: ', name_list(lev[res$bad]), '.', call. = FALSE ),
    '3' = stop(
      'The shrinkage of batch ', name_list(lev[res$which[[1]]]), '\'s effects ',
      'did not converge in ', it$max_iters, ' iterations (`max_iters`).',
      call. = FALSE ),
    '4' = stop(
      'The shrinkage of batch ', name_list(lev[res$which[[1]]]), '\'s effects ',
      'failed for feature ', name_list(if (is.null(info$features)) res$which[[2]] else info$features[res$which[[2]]]),
      ': its standardised values in the batch all equal its shrunken location, ',
      'so its scale is undefined. MMUPHin stops here too.', call. = FALSE ) )

  adjusted <- rowSums(res$ind) > 0
  n_iters  <- res$n_iters
  names(n_iters) <- lev

  batch_result(
    res, info, design, 'mmuphin', adjusted,
    attrs = c(
      if (chosen) list(pseudocount = pseudocount),
      list(n_iters = n_iters) ))
}
