# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# batch_combat(): an R and C port (with src/batch_combat.c) of ComBat() and
# its helpers from sva 3.60.0 (R/ComBat.R, R/helper.R), restructured for
# ecodive, see plans/batch_correction.md.
#
# sva's notice:
#
#   sva: Copyright (c) Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker,
#   Elana J. Fertig, Andrew E. Jaffe, Yuqing Zhang, John D. Storey and
#   Leonardo Collado Torres (the authors named in sva's DESCRIPTION, which
#   has no other copyright line). Licensed under the Artistic License 2.0,
#   below.
#
#   This file is a modified version of that code (Artistic License 2.0,
#   section 4), distributed as part of the package ecodive, which is not
#   sva and doesn't prevent sva from being installed. How it differs: the
#   computation is done in C (src/batch_combat.c, which lists its own
#   differences); `counts` is samples x features by default (`margin`) and
#   may be any of ecodive's table types, optionally normalised first (`norm`
#   'clr' or 'log'); `batch` and the covariates (a `formula` over `data`,
#   instead of a model matrix `mod`) follow ecodive's per-sample rules;
#   batches with one sample switch on `mean_only` with a warning instead of
#   a message; features constant within a batch are left unchanged
#   silently (sva prints how many), and are marked in the result's
#   `adjusted` attribute; missing values are refused; `tol` and `max_iters`
#   are arguments; nothing is printed while it runs; the result is the
#   corrected table with the estimates as attributes, and the prior plots
#   are replaced by plot.ecodive_batch(). The user-visible differences are
#   listed in ?batch_combat ("Differences from sva").
#
#   Artistic License 2.0
#
#   Copyright (c) 2000-2006, The Perl Foundation.
#
#   Everyone is permitted to copy and distribute verbatim copies of this
#   license document, but changing it is not allowed.
#
#   Preamble
#   ********
#
#   This license establishes the terms under which a given free software
#   Package may be copied, modified, distributed, and/or redistributed.  The
#   intent is that the Copyright Holder maintains some artistic control over
#   the development of that Package while still keeping the Package
#   available as open source and free software.
#
#   You are always permitted to make arrangements wholly outside of this
#   license directly with the Copyright Holder of a given Package.  If the
#   terms of this license do not permit the full use that you propose to
#   make of the Package, you should contact the Copyright Holder and seek a
#   different licensing arrangement.
#
#   Definitions
#   ***********
#
#   "Copyright Holder" means the individual(s) or organization(s) named in
#   the copyright notice for the entire Package.
#
#   "Contributor" means any party that has contributed code or other
#   material to the Package, in accordance with the Copyright Holder's
#   procedures.
#
#   "You" and "your" means any person who would like to copy, distribute, or
#   modify the Package.
#
#   "Package" means the collection of files distributed by the Copyright
#   Holder, and derivatives of that collection and/or of those files.  A
#   given Package may consist of either the Standard Version, or a Modified
#   Version.
#
#   "Distribute" means providing a copy of the Package or making it
#   accessible to anyone else, or in the case of a company or organization,
#   to others outside of your company or organization.
#
#   "Distributor Fee" means any fee that you charge for Distributing this
#   Package or providing support for this Package to another party.  It does
#   not mean licensing fees.
#
#   "Standard Version" refers to the Package if it has not been modified, or
#   has been modified only in ways explicitly requested by the Copyright
#   Holder.
#
#   "Modified Version" means the Package, if it has been changed, and such
#   changes were not explicitly requested by the Copyright Holder.
#
#   "Original License" means this Artistic License as Distributed with the
#   Standard Version of the Package, in its current version or as it may be
#   modified by The Perl Foundation in the future.
#
#   "Source" form means the source code, documentation source, and
#   configuration files for the Package.
#
#   "Compiled" form means the compiled bytecode, object code, binary, or any
#   other form resulting from mechanical transformation or translation of
#   the Source form.
#
#   Permission for Use and Modification Without Distribution
#   ********************************************************
#
#   (1) You are permitted to use the Standard Version and create and use
#   Modified Versions for any purpose without restriction, provided that you
#   do not Distribute the Modified Version.
#
#   Permissions for Redistribution of the Standard Version
#   ******************************************************
#
#   (2) You may Distribute verbatim copies of the Source form of the
#   Standard Version of this Package in any medium without restriction,
#   either gratis or for a Distributor Fee, provided that you duplicate all
#   of the original copyright notices and associated disclaimers.  At your
#   discretion, such verbatim copies may or may not include a Compiled form
#   of the Package.
#
#   (3) You may apply any bug fixes, portability changes, and other
#   modifications made available from the Copyright Holder.  The resulting
#   Package will still be considered the Standard Version, and as such will
#   be subject to the Original License.
#
#   Distribution of Modified Versions of the Package as Source
#   **********************************************************
#
#   (4) You may Distribute your Modified Version as Source (either gratis or
#   for a Distributor Fee, and with or without a Compiled form of the
#   Modified Version) provided that you clearly document how it differs from
#   the Standard Version, including, but not limited to, documenting any
#   non-standard features, executables, or modules, and provided that you do
#   at least ONE of the following:
#
#   (a) make the Modified Version available to the Copyright Holder of the
#   Standard Version, under the Original License, so that the Copyright
#   Holder may include your modifications in the Standard Version.
#
#   (b) ensure that installation of your Modified Version does not prevent
#   the user installing or running the Standard Version.  In addition, the
#   Modified Version must bear a name that is different from the name of the
#   Standard Version.
#
#   (c) allow anyone who receives a copy of the Modified Version to make the
#   Source form of the Modified Version available to others under
#
#   (i) the Original License or
#
#   (ii) a license that permits the licensee to freely copy, modify and
#   redistribute the Modified Version using the same licensing terms that
#   apply to the copy that the licensee received, and requires that the
#   Source form of the Modified Version, and of any works derived from it,
#   be made freely available in that license fees are prohibited but
#   Distributor Fees are allowed.
#
#   Distribution of Compiled Forms of the Standard Version or Modified
#   ******************************************************************
#   Versions without the Source
#   ***************************
#
#   (5) You may Distribute Compiled forms of the Standard Version without
#   the Source, provided that you include complete instructions on how to
#   get the Source of the Standard Version.  Such instructions must be valid
#   at the time of your distribution.  If these instructions, at any time
#   while you are carrying out such distribution, become invalid, you must
#   provide new instructions on demand or cease further distribution.  If
#   you provide valid instructions or cease distribution within thirty days
#   after you become aware that the instructions are invalid, then you do
#   not forfeit any of your rights under this license.
#
#   (6) You may Distribute a Modified Version in Compiled form without the
#   Source, provided that you comply with Section 4 with respect to the
#   Source of the Modified Version.
#
#   Aggregating or Linking the Package
#   **********************************
#
#   (7) You may aggregate the Package (either the Standard Version or
#   Modified Version) with other packages and Distribute the resulting
#   aggregation provided that you do not charge a licensing fee for the
#   Package.  Distributor Fees are permitted, and licensing fees for other
#   components in the aggregation are permitted.  The terms of this license
#   apply to the use and Distribution of the Standard or Modified Versions
#   as included in the aggregation.
#
#   (8) You are permitted to link Modified and Standard Versions with other
#   works, to embed the Package in a larger work of your own, or to build
#   stand-alone binary or bytecode versions of applications that include the
#   Package, and Distribute the result without restriction, provided the
#   result does not expose a direct interface to the Package.
#
#   Items That are Not Considered Part of a Modified Version
#   ********************************************************
#
#   (9) Works (including, but not limited to, modules and scripts) that
#   merely extend or make use of the Package, do not, by themselves, cause
#   the Package to be a Modified Version.  In addition, such works are not
#   considered parts of the Package itself, and are not subject to the terms
#   of this license.
#
#   General Provisions
#   ******************
#
#   (10) Any use, modification, and distribution of the Standard or Modified
#   Versions is governed by this Artistic License.  By using, modifying or
#   distributing the Package, you accept this license.  Do not use, modify,
#   or distribute the Package, if you do not accept this license.
#
#   (11) If your Modified Version has been derived from a Modified Version
#   made by someone other than you, you are nevertheless required to ensure
#   that your Modified Version complies with the requirements of this
#   license.
#
#   (12) This license does not grant you the right to use any trademark,
#   service mark, tradename, or logo of the Copyright Holder.
#
#   (13) This license includes the non-exclusive, worldwide, free-of-charge
#   patent license to make, have made, use, offer to sell, sell, import and
#   otherwise transfer the Package with respect to any patent claims
#   licensable by the Copyright Holder that are necessarily infringed by the
#   Package.  If you institute patent litigation (including a cross-claim or
#   counterclaim) against any party alleging that the Package constitutes
#   direct or contributory patent infringement, then this Artistic License
#   to you shall terminate on the date that such litigation is filed.
#
#   (14) Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT
#   HOLDER AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED
#   WARRANTIES.  THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#   PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT
#   PERMITTED BY YOUR LOCAL LAW.  UNLESS REQUIRED BY LAW, NO COPYRIGHT
#   HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT,
#   INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE
#   OF THE PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' ComBat Batch Correction
#'
#' Removes batch effects from a table of continuous values (log-scale
#' metabolite intensities, clr values, log relative abundances), while
#' keeping the effects of named covariates: the empirical Bayes location
#' and scale adjustment of Johnson et al. (2007), a port of `sva::ComBat()`,
#' with parametric or non-parametric priors.
#'
#' @inheritParams batch_mmuphin
#'
#' @param counts  A numeric table (samples \eqn{\times} features), or any of
#'        the objects listed under **Input Types**. With `norm = 'none'`,
#'        the values to correct, used as given (any finite numbers); with
#'        `'clr'` or `'log'`, non-negative abundances to transform first.
#'
#' @param batch  The batch of each sample: one value per sample (matched to
#'        the samples by name when both are named), or the name of a column
#'        of `data`. At least two batches; a batch of one sample turns on
#'        `mean_only`.
#'
#' @param norm  `"none"` (the default): correct the values as given.
#'        `"clr"`: centred log ratios, or `"log"`: log relative abundances
#'        (natural logs), computed first, as [normalize()] does. ComBat
#'        assumes roughly normal values in each batch, so counts should be
#'        transformed.
#'
#' @param pseudocount  For `norm = 'clr'` or `'log'`: added to every value
#'        before the log. Default `NULL`: none if there are no zeros, else
#'        half the smallest non-zero value, with a warning (as
#'        [normalize()]). Ignored with `norm = 'none'`.
#'
#' @param parametric  Shrink the batch effects with parametric priors
#'        (normal locations, inverse-gamma scales; `TRUE`, the default) or
#'        non-parametric ones (each feature's effects are averaged over the
#'        other features' estimates, weighted by how well they fit the
#'        feature; `FALSE`). sva's `par.prior`.
#'
#' @param mean_only  Adjust only the batches' means, not their variances.
#'        Default `FALSE`. A batch with a single sample makes it `TRUE`, with
#'        a warning (as sva).
#'
#' @param ref_batch  `NULL` (the default), or one of the batches: its
#'        values are left unchanged, and the other batches are adjusted to
#'        it (sva's `ref.batch`).
#'
#' @param tol  The convergence tolerance of the parametric shrinkage
#'        iteration. Default `1e-4`, which sva fixes.
#'
#' @param max_iters  The most iterations of the parametric shrinkage; more
#'        is an error. Default `1000` (sva has no cap).
#'
#' @return The corrected table: a dense numeric `matrix` of class
#'   `ecodive_batch` in the orientation of `counts`, with its names, on the
#'   scale of the values corrected (after `norm`). Attributes:
#'   \describe{
#'     \item{`method`}{`"combat"`.}
#'     \item{`batch_effects`}{A data frame with a row per feature and batch:
#'       `feature`, `batch`, `gamma_hat` and `delta_hat` (the feature's
#'       location effect and variance in the batch, on the standardised
#'       scale), `gamma_star` and `delta_star` (after shrinkage; 0 and 1 for
#'       `ref_batch`). `NA` for features not adjusted.}
#'     \item{`adjusted`}{A logical per feature: `FALSE` for features
#'       constant within a batch (of two or more samples), which are
#'       returned unchanged.}
#'     \item{`mean_only`}{Whether only the means were adjusted (`TRUE` also
#'       when a batch has one sample).}
#'     \item{`pseudocount`}{The pseudocount chosen, when `norm` is `'clr'`
#'       or `'log'` with `pseudocount = NULL` and the table has zeros.}
#'     \item{`n_iters`}{With parametric priors and `mean_only = FALSE`, the
#'       shrinkage iterations run in each batch.}
#'     \item{`batch_means`, `batch_means_adjusted`}{Feature \eqn{\times}
#'       batch matrices of the mean values (after `norm`), before and
#'       after.}
#'     \item{`overall_means`}{Each feature's mean value over all samples,
#'       before.}
#'   }
#'   See [batch_mmuphin()] for how the class behaves.
#'
#' @details
#' **The model.** Each feature \eqn{y} is fitted by least squares on the
#' batch indicators and the covariates, and standardised:
#' \eqn{s = (y - \bar\alpha - X\beta) / \hat\sigma}, where \eqn{\bar\alpha}
#' is the batches' coefficients averaged with weights by batch size (or the
#' reference batch's) and \eqn{\hat\sigma^2} the residual variance (divisor
#' \eqn{n}; with `ref_batch`, over the reference batch's samples). In each
#' batch, the location \eqn{\hat\gamma} (mean of \eqn{s}) and scale
#' \eqn{\hat\delta} (variance) are shrunk over features: with
#' `parametric = TRUE`, to their posterior means under a normal prior on
#' \eqn{\gamma} and an inverse-gamma prior on \eqn{\delta} (hyperparameters
#' by the method of moments), iterating until the largest relative change
#' is below `tol`; with `parametric = FALSE`, to averages of the other
#' features' estimates, each weighted by the likelihood of this feature's
#' values under it. The adjusted value is \eqn{(s - \gamma^*) /
#' \sqrt{\delta^*} \cdot \hat\sigma + \bar\alpha + X\beta}. Features
#' constant within a batch are left unchanged.
#'
#' **Which method.** ComBat assumes continuous, roughly normal values: use
#' it on log-transformed metabolomics or proteomics, or on clr or log
#' relative abundances. For count tables with many zeros, [batch_mmuphin()]
#' (zero-inflated ComBat on log relative abundances, re-closed to counts)
#' is designed for microbiome data. As there, covariates that differ
#' between batches must be named in `formula`, and corrected values are for
#' exploration and pooled analyses rather than for testing.
#'
#' The non-parametric priors take \eqn{O(F^2 n)} time for \eqn{F} features
#' and \eqn{n} samples (threaded over features and batches).
#'
#' @section Differences from sva:
#'
#' `batch_combat()` ports sva 3.60.0's `ComBat()` (Artistic License 2.0;
#' Jeffrey T. Leek, W. Evan Johnson and the other sva authors). Its
#' corrected values agree with `ComBat()` to rounding, except:
#'
#' * `counts` is samples \eqn{\times} features by default (`margin`); sva
#'   takes features \eqn{\times} samples. `norm` can transform it first.
#' * `formula` and `data` give the covariates, instead of a model matrix
#'   `mod`; `batch` and `data` follow ecodive's per-sample rules. Batches
#'   are the levels of a factor, else the sorted values.
#' * The non-parametric priors weigh the other features' estimates by their
#'   likelihoods on the log scale. sva multiplies the densities, which
#'   underflow to zero for large batches and then give `NaN` for every
#'   feature; where sva's results are finite, the two agree to rounding.
#' * `tol` and `max_iters` are arguments: sva fixes the tolerance at
#'   `1e-4` and has no cap on the iterations. An iteration that meets a
#'   `NaN` stops with an error (sva: "missing value where TRUE/FALSE
#'   needed").
#' * Missing values (`NA`) are refused; sva fits around them.
#' * A batch of one sample turns on `mean_only` with a warning (sva: a
#'   message); features constant within a batch are left unchanged without
#'   a message, and marked in the `adjusted` attribute. Nothing else is
#'   printed.
#' * The least-squares fits use a QR decomposition where sva solves the
#'   normal equations, so very ill-conditioned covariates may agree less
#'   closely.
#' * A constant covariate is confounded with the batches (an error); sva
#'   drops a covariate column of all ones.
#' * `prior.plots` and `BPPARAM` are not arguments: `plot()` draws the
#'   diagnostics, and `cpus` sets the threads.
#'
#' **Same names in other packages:** none. sva's own function is
#' `ComBat()`.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects
#' in microarray expression data using empirical Bayes methods.
#' *Biostatistics*, 8, 118–127. \doi{10.1093/biostatistics/kxj037}
#'
#' Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E., & Storey, J. D.
#' (2012). The sva package for removing batch effects and other unwanted
#' variation in high-throughput experiments. *Bioinformatics*, 28,
#' 882–883. \doi{10.1093/bioinformatics/bts034}
#'
#' Wang, Y., & Lê Cao, K.-A. (2020). Managing batch effects in microbiome
#' data. *Briefings in Bioinformatics*, 21, 1954–1970.
#' \doi{10.1093/bib/bbz105}
#'
#' @seealso [batch_mmuphin()]; [plot.ecodive_batch()]; [normalize()].
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
#'     # ComBat on clr values.
#'     x <- batch_combat(counts, batch, norm = 'clr')
#'     x
#'
#'     # Batch's share of the variation, before and after.
#'     clr <- normalize(counts, 'clr')
#'     permanova(euclidean(clr), batch)$statistic
#'     permanova(euclidean(x), batch)$statistic
#'
#'     plot(x)
#'
#'     # Non-parametric priors, with run1 as the reference.
#'     batch_combat(clr, batch, parametric = FALSE, ref_batch = 'run1')
batch_combat <- function (
    counts, batch, formula = NULL, data = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, parametric = TRUE, mean_only = FALSE,
    ref_batch = NULL, tol = 1e-4, max_iters = 1000L, cpus = n_cpus() ) {

  if (missing(batch)) batch <- NULL

  if (!is.character(norm) || length(norm) != 1 || is.na(norm) ||
      !tolower(norm) %in% c('none', 'clr', 'log'))
    stop('`norm` must be "none", "clr" or "log".', call. = FALSE)
  norm <- tolower(norm)

  info <- batch_input(nonneg = norm != 'none')

  parametric <- batch_flag(parametric, 'parametric')
  mean_only  <- batch_flag(mean_only, 'mean_only')
  it         <- batch_iter_args(tol, max_iters)

  # The pseudocount, as normalize() chooses it.
  chosen <- FALSE
  pc     <- 0
  norm_code <- NORM_NONE
  if (norm != 'none') {
    norm_code <- if (norm == 'clr') NORM_CLR else NORM_LOG
    chosen <- is.null(pseudocount)
    pc <- log_pseudocount(
      info$counts, pseudocount,
      if (norm == 'clr') 'CLR normalization' else 'log normalization' )
    chosen <- chosen && pc > 0
  }

  design <- batch_design(batch, formula, data, info, min_size = 1L)
  lev    <- design$levels
  B      <- length(lev)

  ref <- -1L
  if (!is.null(ref_batch)) {
    if (!is.atomic(ref_batch) || length(ref_batch) != 1 || is.na(ref_batch) ||
        !as.character(ref_batch) %in% lev)
      stop('`ref_batch` must be one of the batches: ', name_list(lev), '.', call. = FALSE)
    ref <- match(as.character(ref_batch), lev) - 1L
  }

  if (!mean_only && any(design$sizes == 1L)) {
    warning(
      'A batch has a single sample, so only the batches\' means are ',
      'adjusted (`mean_only = TRUE`), as in sva: ',
      name_list(lev[design$sizes == 1L]), '.', call. = FALSE )
    mean_only <- TRUE
  }

  # The design: batch indicators (the reference batch's all ones), then
  # the covariates.
  X <- design$X
  if (ref >= 0L) X[, ref + 1L] <- 1

  res <- .Call(
    C_batch_combat, info$counts, info$margin, norm_code, pc, design$codes, X,
    B, ref, parametric, mean_only, it$tol, it$max_iters, cpus )

  feat <- function (f) name_list(if (is.null(info$features)) f else info$features[f])
  switch(
    as.character(res$status),
    '1' = stop(
      'Fewer than two features vary within every batch (of two or more ',
      'samples), so the priors can\'t be estimated: ComBat needs at least two.',
      call. = FALSE ),
    '2' = stop(
      'The shrinkage of batch ', name_list(lev[res$which[[1]]]), '\'s effects ',
      'did not converge in ', it$max_iters, ' iterations (`max_iters`).',
      call. = FALSE ),
    '3' = stop(
      'The shrinkage of batch ', name_list(lev[res$which[[1]]]), '\'s effects ',
      'gave NaN (sva stops here too, with "missing value where TRUE/FALSE ',
      'needed"). Try `parametric = FALSE` or `mean_only = TRUE`.', call. = FALSE ),
    '4' = stop(
      'The non-parametric estimate for feature ', feat(res$which[[2]]),
      ' in batch ', name_list(lev[res$which[[1]]]), ' has no finite weight.',
      call. = FALSE ),
    '5' = stop('The design is rank deficient.', call. = FALSE) )  # nocov

  n_iters <- res$n_iters
  names(n_iters) <- lev

  batch_result(
    res, info, design, 'combat', res$keep,
    attrs = c(
      list(mean_only = mean_only),
      if (chosen) list(pseudocount = pc),
      if (parametric && !mean_only) list(n_iters = n_iters) ))
}
