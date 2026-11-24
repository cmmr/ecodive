# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# decontam(), decontam_not() and decontam_frequency_plot(): an R and C port
# (with src/decontam.c) of isContaminant(), isNotContaminant(),
# isContaminantFrequency(), isContaminantPrevalence() and fish.combine()
# (R/decontam.R) and of plot_frequency()'s fitted lines (R/plotting.R) from
# decontam 1.32.0 (R/decontam.R is identical in 1.22.0), restructured for
# ecodive, see plans/decontam.md.
#
# decontam's notice:
#
#   decontam: Copyright (c) Benjamin Callahan and Nicole Marie Davis (the
#   authors named in decontam's DESCRIPTION, which has no other copyright
#   line). Licensed under the Artistic License 2.0, below.
#
#   This file is a modified version of that code (Artistic License 2.0,
#   section 4), distributed as part of the package ecodive, which is not
#   decontam and doesn't prevent decontam from being installed. How it
#   differs: the scores are computed from per-feature sums made in C
#   (src/decontam.c) instead of lm() fits and table()s; the prevalence
#   test's prop.test() and fisher.test() calls are replaced by their
#   formulas, vectorized over features (the two-proportion statistic with
#   its continuity correction, used when no expected count is below 5, and
#   the hypergeometric mid-p otherwise); per-sample values come from `data`
#   by ecodive's rules instead of phyloseq sample variables; `normalize` is
#   `norm`; `detailed` is dropped; the result keeps the control counts and
#   the per-batch scores; samples whose `neg` is NA are left out of both
#   tests; a batch's unused factor levels are dropped; and
#   plot_frequency()'s ggplot2 figure is drawn with base graphics. The
#   user-visible differences are listed in ?decontam ("Differences from
#   decontam").
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



#' Contaminant Identification (decontam)
#'
#' Find the features (ASVs, OTUs, taxa, ...) that come from contaminating
#' DNA, by decontam's two signatures of contamination (Davis et al. 2018):
#'
#' * **Frequency:** a contaminant adds a roughly constant amount of DNA to
#'   every sample, so its relative abundance falls as the sample's total
#'   DNA concentration rises: it goes as \eqn{1/C}. A true sequence's
#'   relative abundance doesn't depend on \eqn{C}. Needs `conc`.
#' * **Prevalence:** a contaminant is present in more of the negative
#'   controls than of the true samples. Needs `neg`.
#'
#' `decontam()` scores each feature by either or both (decontam's
#' `isContaminant()`), and `decontam_not()` reverses the prevalence test's
#' null hypothesis, to find the features that are **not** contaminants in
#' low-biomass samples, where most may be (decontam's
#' `isNotContaminant()`). Neither removes anything: keep the features that
#' aren't called, `counts[, !res$contaminant]` (or
#' `counts[!res$contaminant, ]` with `margin = 2`).
#'
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features), or any of
#'        the objects listed under **Input Types**. Each sample's values are
#'        divided by its total (`norm = 'percent'`). Samples with no counts
#'        are left out, with a warning.
#'
#' @param conc  Each sample's DNA concentration before sequencing (for
#'        instance a fluorometer reading of the amplified DNA, or a qPCR
#'        count), in any unit: a positive number, for every sample, controls
#'        included. Only its inverse relation with a contaminant's frequency
#'        matters. A numeric vector with one value per sample, matched to the
#'        samples by name when both are named, else in the samples' order;
#'        or, when `data` is given, the name of one of its columns. Needed
#'        for the frequency test. Default: `NULL`.
#'
#' @param neg  Which samples are negative controls: a logical vector,
#'        `TRUE` for a control and `FALSE` for a true sample, matched as
#'        `conc`, or the name of a column of `data`. Samples with `NA` are
#'        left out of both tests. Needed for the prevalence test; for the
#'        frequency test, the controls are left out. Extraction blanks give
#'        the best results. Default: `NULL`.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `conc`, `neg` or `batch` as columns. Rows are
#'        matched to the samples in `counts` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as the samples. Default: `NULL`.
#'
#' @param method  Which test calls the contaminants:
#'        * `'auto'` (the default): `'frequency'` when only `conc` is given,
#'          `'prevalence'` when only `neg` is, and `'combined'` when both
#'          are.
#'        * `'frequency'`, `'prevalence'`: one test.
#'        * `'combined'`: the two scores combined by Fisher's method,
#'          \eqn{P(\chi^2_4 > -2 \log(s_f s_p))}.
#'        * `'minimum'`: the smaller of the two scores.
#'        * `'either'`, `'both'`: a contaminant by either test, or by both
#'          (`p` is then `NA`).
#'
#' @param batch  Optional batches (sequencing runs, extraction plates, ...):
#'        one value per sample, matched as `conc`, or the name of a column
#'        of `data`. Each batch is scored on its own, and the scores are
#'        combined by `batch_combine`. Each batch needs at least two
#'        samples, controls included (a warning at four or fewer). Default:
#'        `NULL`, one batch.
#'
#' @param batch_combine  How a feature's scores in the batches are combined
#'        into one: `'minimum'` (the default; the lowest score of the
#'        batches that have one), `'product'` (their product; 1 when no
#'        batch has a score, as decontam), or `'fisher'` (Fisher's method,
#'        with a missing score counted as 0.5). Applied to each test's
#'        scores, even with one batch.
#'
#' @param threshold  A feature is called a contaminant when its score is
#'        **strictly less** than `threshold`, a number from 0 to 1. Default:
#'        `0.1`, as decontam; `0.5` calls every feature whose prevalence
#'        score says it is more prevalent in the controls than in the
#'        samples. With `method = 'either'` or `'both'`, two values, the
#'        frequency and the prevalence thresholds (one value is used for
#'        both). For `decontam_not()`, a feature is called a
#'        **non**-contaminant when its score is less than `threshold`
#'        (default `0.5`).
#'
#' @param norm  `'percent'` (the default): each sample's values are divided
#'        by its total, as decontam's `normalize = TRUE`. `'none'`: the
#'        values are used as given (relative abundances already, or counts
#'        from samples of equal depth), as `normalize = FALSE`. Both tests
#'        need non-negative relative abundances, so no other norm is
#'        allowed.
#'
#' @return A data frame with one row per feature, of class
#'   `c("ecodive_decontam", "data.frame")` (`decontam()`) or
#'   `c("ecodive_decontam_not", "data.frame")` (`decontam_not()`), with
#'   columns:
#'   \describe{
#'     \item{`feature`}{The feature's name (its index if it has none).}
#'     \item{`freq`}{Its mean relative abundance over all the samples,
#'       controls included (after `norm`).}
#'     \item{`prev`}{The number of samples it is present in.}
#'     \item{`prev_neg`, `prev_pos`}{The number of negative controls and of
#'       true samples it is present in (`NA` without `neg`).}
#'     \item{`p_freq`}{The frequency score (`decontam()` only; `NA` when the
#'       frequency test wasn't run, or the feature is in fewer than two
#'       true samples of every batch).}
#'     \item{`p_prev`}{The prevalence score (`NA` when the prevalence test
#'       wasn't run, or the feature is in fewer than two samples of every
#'       batch). For `decontam_not()`, one minus the prevalence score.}
#'     \item{`p`}{The score compared with `threshold`, by `method` (for
#'       `decontam_not()`, `p_prev`).}
#'     \item{`contaminant`, or `not_contaminant`}{The call: `p <
#'       threshold`, with `NA` scores never called.}
#'   }
#'   Attributes: `method` (`"decontam"` or `"decontam_not"`), `test` (the
#'   `method` used, after `'auto'`), `threshold`, `batch_combine`
#'   (`decontam()` only), and `n_neg` and `n_pos`, the number of negative
#'   controls and of true samples tested (without `neg`, 0 and every
#'   sample). With `batch`, `p_freq_batch` and `p_prev_batch` are batch
#'   \eqn{\times} feature matrices of the scores in each batch (for the tests
#'   run).
#'
#' @details
#' **The frequency score.** In each batch, over the true samples (not
#' controls) where the feature is present, \eqn{k} of them, with frequency
#' \eqn{f} and concentration \eqn{c}: \eqn{SS_0} is the residual sum of
#' squares of \eqn{\log f} about its mean (the non-contaminant model,
#' frequency independent of \eqn{c}), and \eqn{SS_1} that of \eqn{\log f
#' + \log c} about its mean (the contaminant model, \eqn{\log f} falling
#' with slope \eqn{-1} in \eqn{\log c}). The score is
#' \eqn{F_{k-1,k-1}(SS_1 / SS_0)}, the lower tail of an F distribution:
#' small when the contaminant model fits much better. It is a heuristic, a
#' ratio of two sums of squares with equal degrees of freedom, **not a
#' p-value**, and Davis et al. call it a score. Features in fewer than two
#' true samples get `NA`.
#'
#' **The prevalence score.** In each batch, the feature's presence in the
#' controls and in the true samples is a 2 \eqn{\times} 2 table. The score is
#' the one-sided p-value that the feature is more prevalent in the
#' controls: a chi-squared test of two proportions with a continuity
#' correction (as `prop.test(alternative = "greater")`), or, when an
#' expected count is below 5, Fisher's exact test's mid-p
#' (\eqn{P(X > x) + P(X = x) / 2} under the hypergeometric distribution of
#' the table's margins). A feature present in every sample of the batch
#' scores 0.5; one in fewer than two samples, or a batch without both
#' controls and true samples, gives `NA`.
#'
#' **Choosing `threshold`.** Plot the scores (`plot(res)`): in Davis et
#' al.'s data they are bimodal, low for contaminants and high for the rest,
#' and the threshold goes in the trough. decontam's default is 0.1; 0.5
#' (with the prevalence test) removes every feature more prevalent in the
#' controls than in the samples, an aggressive choice suited to
#' low-biomass samples, where `decontam_not()` is the alternative.
#'
#' **Batches.** Contaminants differ between reagent lots and runs, so
#' samples processed together should be scored together, as `batch`.
#'
#' @section Controls and reporting:
#'
#' The guidelines of Fierer et al. (2025) for low-biomass studies ask for
#' negative controls at every stage (sampling, extraction, library
#' preparation), carried through sequencing like the samples, and for the
#' contaminant-removal workflow to be reported. decontam's own advice is
#' that extraction blanks give the best results. Statistical removal finds
#' contaminants that come from reagents and the laboratory; it doesn't undo
#' cross-contamination between samples ("well-to-well" leakage), which
#' looks like a true sequence. Report the method, the test, the threshold
#' and the number of features and reads removed.
#'
#' @section Differences from decontam:
#'
#' `decontam()` and `decontam_not()` port decontam 1.32.0's
#' `isContaminant()` and `isNotContaminant()` (Artistic License 2.0;
#' Benjamin Callahan and Nicole Marie Davis), and their scores and calls
#' agree with it to rounding, except:
#'
#' * **`NA` in `neg`** leaves those samples out of both tests. decontam
#'   leaves them out of the prevalence table, but still counts them towards
#'   a feature's two samples, and its frequency test subsets with the `NA`s
#'   and returns `NA` for every feature.
#' * `norm = 'percent'` or `'none'` replaces `normalize = TRUE` or `FALSE`.
#' * `conc`, `neg` and `batch` can name columns of `data`, matched to the
#'   samples by name, instead of phyloseq sample variables; a phyloseq
#'   object is read as any other count table.
#' * The result is always the table (decontam's `detailed = TRUE`), in
#'   snake case, with `prev_neg` and `prev_pos` added and the per-batch
#'   scores kept as attributes; `decontam_not()` drops the empty `p.freq`
#'   column.
#' * A batch needs at least two samples; levels of a factor `batch` with no
#'   samples are dropped (decontam counts them and stops), and an `NA`
#'   batch is an error. A batch with fewer than two true samples gives
#'   `NA` frequency scores, where decontam stops.
#' * One `threshold` for `'either'` or `'both'` is used for both tests
#'   without a message.
#' * decontam's `plot_frequency()` is [decontam_frequency_plot()], drawn
#'   with base graphics, and the score histogram and the prevalence plot
#'   are `plot()` methods.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Davis, N. M., Proctor, D. M., Holmes, S. P., Relman, D. A., & Callahan,
#' B. J. (2018). Simple statistical identification and removal of
#' contaminant sequences in marker-gene and metagenomics data.
#' *Microbiome*, 6, 226. \doi{10.1186/s40168-018-0605-2}
#'
#' Fierer, N., et al. (2025). Guidelines for preventing and reporting
#' contamination in low-biomass microbiome studies. *Nature Microbiology*,
#' 10, 1570–1580. \doi{10.1038/s41564-025-02035-2}
#'
#' Salter, S. J., et al. (2014). Reagent and laboratory contamination can
#' critically impact sequence-based microbiome analyses. *BMC Biology*, 12,
#' 87. \doi{10.1186/s12915-014-0087-z}
#'
#' Eisenhofer, R., et al. (2019). Contamination in low microbial biomass
#' microbiome studies: issues and recommendations. *Trends in
#' Microbiology*, 27, 105–117. \doi{10.1016/j.tim.2018.11.003}
#'
#' @seealso [plot.ecodive_decontam()] for the score histogram and the
#'   prevalence plot, and [decontam_frequency_plot()] for a feature's
#'   frequency against `conc`.
#'
#' @export
#' @examples
#'     # Ten samples and two extraction blanks, with DNA concentrations.
#'     # F1 and F2 are contaminants: a fixed amount of DNA, so their share
#'     # falls as the concentration rises, and they fill the blanks.
#'     set.seed(1)
#'     conc   <- c(exp(seq(1, 4, length.out = 10)), 0.5, 0.6)
#'     neg    <- rep(c(FALSE, TRUE), c(10, 2))
#'     counts <- matrix(rpois(12 * 8, 200), 12, 8, dimnames = list(
#'       c(paste0('S', 1:10), 'Blank1', 'Blank2'), paste0('F', 1:8) ))
#'     counts[, 1:2] <- round(counts[, 1:2] * 5 / conc)
#'     counts[11:12, 3:8] <- 0
#'
#'     # The frequency test
#'     res <- decontam(counts, conc)
#'     res
#'
#'     # Frequency and prevalence, combined
#'     res <- decontam(counts, conc, neg)
#'     res$feature[res$contaminant]
#'
#'     # Remove them
#'     clean <- counts[!neg, !res$contaminant]
#'
#'     # The values as columns of sample metadata
#'     meta <- data.frame(dna = conc, blank = neg, row.names = rownames(counts))
#'     res  <- decontam(counts, 'dna', 'blank', meta, method = 'prevalence')
#'
#'     # Low-biomass samples: the features that are not contaminants
#'     decontam_not(counts, neg)
#'
#'     # Plots
#'     plot(res)
#'     plot(res, 'prevalence')
#'     decontam_frequency_plot(counts, c('F1', 'F5'), conc, neg)
#'
decontam <- function (
    counts, conc = NULL, neg = NULL, data = NULL, margin = 1L,
    method = c('auto', 'frequency', 'prevalence', 'combined', 'minimum', 'either', 'both'),
    batch = NULL, batch_combine = c('minimum', 'product', 'fisher'),
    threshold = 0.1, norm = 'percent', cpus = n_cpus() ) {

  method <- decontam_choice(
    method, c('auto', 'frequency', 'prevalence', 'combined', 'minimum', 'either', 'both'),
    'method' )
  batch_combine <- decontam_choice(batch_combine, c('minimum', 'product', 'fisher'), 'batch_combine')

  info <- decontam_input()

  if (method == 'auto')
    method <- if (!is.null(conc) && is.null(neg)) 'frequency'
              else if (is.null(conc) && !is.null(neg)) 'prevalence'
              else 'combined'

  do_freq <- method %in% c('frequency', 'combined', 'minimum', 'either', 'both')
  do_prev <- method %in% c('prevalence', 'combined', 'minimum', 'either', 'both')

  if (do_prev && is.null(neg))
    stop('`neg` is required for the prevalence test (method "', method, '").', call. = FALSE)
  if (do_freq && is.null(conc))
    stop('`conc` is required for the frequency test (method "', method, '").', call. = FALSE)

  two <- method %in% c('either', 'both')
  threshold <- decontam_threshold(threshold, two)

  core <- decontam_core(info, conc, neg, batch, do_freq, do_prev, batch_combine, cpus)
  p_freq <- core$p_freq
  p_prev <- core$p_prev

  p <- switch(
    method,
    frequency  = p_freq,
    prevalence = p_prev,
    minimum    = pmin(p_freq, p_prev),
    combined   = stats::pchisq(-2 * log(p_freq * p_prev), df = 4, lower.tail = FALSE),
    rep(NA_real_, info$D) )

  call <- switch(
    method,
    either = (p_freq < threshold[[1]]) | (p_prev < threshold[[2]]),
    both   = (p_freq < threshold[[1]]) & (p_prev < threshold[[2]]),
    p < threshold )
  call[is.na(call)] <- FALSE

  res <- data.frame(
    core$table,
    p_freq      = p_freq,
    p_prev      = p_prev,
    p           = p,
    contaminant = call,
    stringsAsFactors = FALSE )

  attr(res, 'method')        <- 'decontam'
  attr(res, 'test')          <- method
  attr(res, 'threshold')     <- threshold
  attr(res, 'batch_combine') <- batch_combine
  attr(res, 'n_neg')         <- core$n_neg
  attr(res, 'n_pos')         <- core$n_pos
  if (core$batched) {
    if (do_freq) attr(res, 'p_freq_batch') <- core$p_freq_batch
    if (do_prev) attr(res, 'p_prev_batch') <- core$p_prev_batch
  }
  class(res) <- c('ecodive_decontam', 'data.frame')
  res
}


#' @rdname decontam
#' @export
decontam_not <- function (
    counts, neg, data = NULL, margin = 1L, threshold = 0.5, norm = 'percent',
    cpus = n_cpus() ) {

  if (missing(neg)) neg <- NULL
  info <- decontam_input()
  if (is.null(neg))
    stop('`neg` is required: which samples are negative controls.', call. = FALSE)
  threshold <- decontam_threshold(threshold, FALSE)

  core   <- decontam_core(info, NULL, neg, NULL, FALSE, TRUE, 'minimum', cpus)
  p_prev <- 1 - core$p_prev
  call   <- p_prev < threshold
  call[is.na(call)] <- FALSE

  res <- data.frame(
    core$table,
    p_prev          = p_prev,
    p               = p_prev,
    not_contaminant = call,
    stringsAsFactors = FALSE )

  attr(res, 'method')    <- 'decontam_not'
  attr(res, 'test')      <- 'prevalence'
  attr(res, 'threshold') <- threshold
  attr(res, 'n_neg')     <- core$n_neg
  attr(res, 'n_pos')     <- core$n_pos
  class(res) <- c('ecodive_decontam_not', 'data.frame')
  res
}



#------------------------------------------------------------------------------
# Checks decontam()'s and decontam_not()'s `counts`, `margin`, `cpus`, `norm`
# and `data` (in the caller's environment, as validate_*() do). Returns
# list(n, D, labels, features, idx (data's row for each sample), keep
# (samples with counts), norm (0 or 1)).
#------------------------------------------------------------------------------
decontam_input <- function (env = parent.frame()) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- env$counts
  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  remove('vals')

  norm <- env$norm
  if (!is.character(norm) || length(norm) != 1 || is.na(norm) ||
      !tolower(norm) %in% c('percent', 'none'))
    stop(
      '`norm` must be "percent" or "none": both tests need non-negative ',
      'relative abundances.', call. = FALSE )
  norm <- if (tolower(norm) == 'percent') NORM_PERCENT else NORM_NONE

  margin <- env$margin
  smp    <- count_samples(counts, margin)
  dn     <- switch(
    get_matrix_package(counts),
    base = dimnames(counts), slam = counts$dimnames, Matrix = counts@Dimnames )
  dims   <- switch(
    get_matrix_package(counts),
    base = dim(counts), slam = c(counts$nrow, counts$ncol), Matrix = counts@Dim )
  D <- as.integer(dims[[3L - margin]])
  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])

  data <- env$data
  idx  <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(smp$labels, smp$n, rn, nrow(data), ref = '`counts`')
  }

  keep <- sample_totals(counts, margin) > 0

  list(
    n = smp$n, D = D, labels = smp$labels, features = features,
    data = data, idx = idx, keep = keep, norm = norm,
    counts = counts, margin = margin )
}


# One of `choices`, by match.arg()'s rules (a unique prefix), for `arg`.
decontam_choice <- function (v, choices, arg) {
  if (identical(v, choices)) return (choices[[1]])
  i <- if (is.character(v) && length(v) == 1 && !is.na(v)) pmatch(tolower(v), choices) else NA
  if (is.na(i))
    stop('`', arg, '` must be one of ', paste0('"', choices, '"', collapse = ', '), '.', call. = FALSE)
  choices[[i]]
}


# `threshold`: numbers from 0 to 1, one (or, for 'either' and 'both', two:
# one is used for both).
decontam_threshold <- function (threshold, two) {
  if (!is.numeric(threshold) || length(threshold) == 0 || anyNA(threshold) ||
      any(threshold < 0 | threshold > 1))
    stop('`threshold` must be a number from 0 to 1.', call. = FALSE)
  if (two) {
    if (length(threshold) > 2)
      stop('`threshold` must be one or two numbers (frequency, prevalence).', call. = FALSE)
    return (rep_len(as.double(threshold), 2L))
  }
  if (length(threshold) != 1)
    stop(
      '`threshold` must be a single number; two (frequency, prevalence) ',
      'only with method "either" or "both".', call. = FALSE )
  as.double(threshold)
}


# A per-sample argument by per_sample()'s rule, and what to call it in
# messages: '`conc`', or '`data` column "x"'.
decontam_values <- function (v, arg, info) {
  src <- if (!is.null(info$data) && is.character(v) && length(v) == 1 &&
             !is.na(v) && v %in% names(info$data)) v
  val <- per_sample(v, arg, info$data, info$idx, info$labels, info$n, ref = '`counts`')
  if (is.factor(val) && arg != 'batch') val <- as.character(val)
  list(
    val  = val,
    what = if (is.null(src)) paste0('`', arg, '`') else paste0('`data` column ', name_list(src)) )
}


# Samples named in a message: their names, or else their positions.
decontam_samples <- function (bad, labels) {
  sprintf(
    '%d sample%s: %s', length(bad), if (length(bad) == 1) '' else 's',
    name_list(if (is.null(labels)) bad else labels[bad]) )
}


#------------------------------------------------------------------------------
# The two tests, per batch, and the batches combined. Returns list(table
# (feature, freq, prev, prev_neg, prev_pos), p_freq, p_prev, p_freq_batch,
# p_prev_batch (B x D), batched, n_neg, n_pos).
#------------------------------------------------------------------------------
decontam_core <- function (info, conc, neg, batch, do_freq, do_prev, batch_combine, cpus) {

  n    <- info$n
  D    <- info$D
  keep <- info$keep

  # Samples with no counts are left out, as decontam does.
  if (!any(keep))
    stop('Every sample in `counts` is empty (no counts).', call. = FALSE)
  if (!all(keep))
    warning(
      'Removed ', sum(!keep), ' sample', if (sum(!keep) == 1) '' else 's',
      ' with zero total counts (or frequency).', call. = FALSE )

  # Negative controls: 1, true samples: 0, unknown (NA): -1.
  neg_code <- rep(0L, n)
  if (!is.null(neg)) {
    v <- decontam_values(neg, 'neg', info)
    if (!is.logical(v$val))
      stop(v$what, ' must be logical: TRUE for a negative control, FALSE for a true sample.', call. = FALSE)
    neg_code <- ifelse(is.na(v$val), -1L, as.integer(v$val))
  }

  # log(conc), every sample > 0 (controls too, as decontam checks).
  logc <- NULL
  if (do_freq) {
    v <- decontam_values(conc, 'conc', info)
    if (!is.numeric(v$val))
      stop(v$what, ' must be numeric: each sample\'s DNA concentration.', call. = FALSE)
    bad <- which(keep & (is.na(v$val) | !is.finite(v$val) | v$val <= 0))
    if (length(bad))
      stop(
        v$what, ' must be finite and greater than 0, but is not for ',
        decontam_samples(bad, info$labels), call. = FALSE )
    logc <- ifelse(keep, log(as.double(v$val)), 0)
  }

  # Batches, of the samples kept: each with at least 2 samples.
  batched <- !is.null(batch)
  if (batched) {
    v <- decontam_values(batch, 'batch', info)
    if (anyNA(v$val[keep]))
      stop(v$what, ' must not contain NA: ', decontam_samples(which(keep & is.na(v$val)), info$labels),
           call. = FALSE)
    f    <- factor(v$val[keep])
    tab  <- table(f)
    if (any(tab < 2))
      stop(
        'Every batch needs at least 2 samples, but ', name_list(names(tab)[tab < 2]),
        if (sum(tab < 2) == 1) ' has' else ' have', ' fewer.', call. = FALSE )
    if (any(tab <= 4))
      warning(
        'Some batches have very few (<= 4) samples: ', name_list(names(tab)[tab <= 4]), '.',
        call. = FALSE )
    lev  <- levels(f)
    code <- rep(-1L, n)
    code[keep] <- as.integer(f) - 1L
  }
  else {
    lev  <- '1'
    code <- ifelse(keep, 0L, -1L)
  }
  B <- length(lev)

  res <- .Call(
    C_decontam, info$counts, info$margin, info$norm, code, B, neg_code, logc, cpus )

  features <- if (is.null(info$features)) as.character(seq_len(D)) else info$features
  used     <- code >= 0L
  n_neg_b  <- tabulate(code[used & neg_code == 1L] + 1L, B)
  n_pos_b  <- tabulate(code[used & neg_code == 0L] + 1L, B)

  # The frequency score per batch: F(SS1 / SS0; k - 1, k - 1), for k > 1.
  p_freq_batch <- matrix(NA_real_, B, D, dimnames = list(lev, features))
  if (do_freq) {
    k  <- res$k
    ok <- k > 1L
    p_freq_batch[ok] <- stats::pf(res$ss1[ok] / res$ss0[ok], k[ok] - 1, k[ok] - 1)
    p_freq_batch[is.nan(p_freq_batch)] <- NA_real_
  }

  # The prevalence score per batch, for features in 2 or more samples of a
  # batch with controls and true samples.
  p_prev_batch <- matrix(NA_real_, B, D, dimnames = list(lev, features))
  if (do_prev) {
    a  <- res$prev_neg
    c  <- res$prev_pos
    n1 <- matrix(n_neg_b, B, D)
    n2 <- matrix(n_pos_b, B, D)
    ok <- (a + c) > 1L & n1 > 0L & n2 > 0L
    p_prev_batch[ok] <- decontam_prev_score(a[ok], n1[ok], c[ok], n2[ok])
  }

  table <- data.frame(
    feature  = features,
    freq     = res$freq / sum(used),
    prev     = res$prev,
    prev_neg = if (is.null(neg)) NA_integer_ else as.integer(colSums(res$prev_neg)),
    prev_pos = if (is.null(neg)) NA_integer_ else as.integer(colSums(res$prev_pos)),
    stringsAsFactors = FALSE )

  list(
    table        = table,
    p_freq       = if (do_freq) decontam_combine(p_freq_batch, batch_combine) else rep(NA_real_, D),
    p_prev       = if (do_prev) decontam_combine(p_prev_batch, batch_combine) else rep(NA_real_, D),
    p_freq_batch = p_freq_batch,
    p_prev_batch = p_prev_batch,
    batched      = batched,
    n_neg        = sum(n_neg_b),
    n_pos        = sum(n_pos_b) )
}


#------------------------------------------------------------------------------
# isContaminantPrevalence()'s score from a 2 x 2 table of presence, given as
# vectors: `a` of the `n1` negative controls and `c` of the `n2` true
# samples have the feature.
#
#   * Present in every sample: 0.5.
#   * Else prop.test(table, alternative = "greater")$p.value: the one-sided
#     two-proportion test, with the continuity correction capped at the
#     difference between the proportions; when an expected count is below
#     5, prop.test() warns, and decontam takes instead
#   * fisher.test()'s one-sided p-value less half the probability of the
#     table (its mid-p), as decontam's fisher.pval() computes it.
#------------------------------------------------------------------------------
decontam_prev_score <- function (a, n1, c, n2) {

  a  <- as.double(a); c <- as.double(c); n1 <- as.double(n1); n2 <- as.double(n2)
  b  <- n1 - a
  d  <- n2 - c
  p  <- rep(0.5, length(a))
  go <- (b + d) > 0
  if (!any(go)) return (p)
  a <- a[go]; b <- b[go]; c <- c[go]; d <- d[go]; n1 <- n1[go]; n2 <- n2[go]

  # Expected counts under one pooled proportion.
  pp <- (a + c) / (n1 + n2)
  e1 <- n1 * pp; e2 <- n2 * pp; e3 <- n1 * (1 - pp); e4 <- n2 * (1 - pp)
  small <- e1 < 5 | e2 < 5 | e3 < 5 | e4 < 5

  # The two-proportion statistic with the continuity correction.
  delta <- a / n1 - c / n2
  yates <- pmin(0.5, abs(delta) / (1 / n1 + 1 / n2))
  stat  <- (abs(a - e1) - yates)^2 / e1 + (abs(c - e2) - yates)^2 / e2 +
           (abs(b - e3) - yates)^2 / e3 + (abs(d - e4) - yates)^2 / e4
  chi   <- stats::pnorm(sign(delta) * sqrt(stat), lower.tail = FALSE)

  # Fisher's mid-p: P(X >= a) - (P(X >= a) + P(X <= a) - 1) / 2, X the
  # controls with the feature, hypergeometric with the table's margins.
  m  <- a + c
  nn <- b + d
  pg <- stats::phyper(a - 1, m, nn, n1, lower.tail = FALSE)
  pl <- stats::phyper(a, m, nn, n1)
  fis <- pg - (pg + pl - 1) / 2

  p[go] <- ifelse(small, fis, chi)
  p
}


# A feature's scores in the batches (a B x D matrix, NA for none) combined
# into one, as isContaminant()'s batch.combine.
decontam_combine <- function (m, how) {
  B <- nrow(m)
  if (how == 'minimum') {
    r <- m[1L, ]
    for (b in seq_len(B)[-1L]) r <- pmin(r, m[b, ], na.rm = TRUE)
    return (unname(r))
  }
  if (how == 'product') {
    r <- rep(1, ncol(m))
    for (b in seq_len(B)) r <- r * ifelse(is.na(m[b, ]), 1, m[b, ])
    return (unname(r))
  }
  # fish.combine(na.replace = 0.5).
  r <- rep(1, ncol(m))
  for (b in seq_len(B)) r <- r * ifelse(is.na(m[b, ]), 0.5, m[b, ])
  unname(stats::pchisq(-2 * log(r), df = 2 * B, lower.tail = FALSE))
}



#' Plot decontam() results
#'
#' * `type = 'scores'` (the default): a histogram of the scores, with the
#'   `threshold` as a dashed line. Davis et al. (2018) choose the threshold
#'   in the trough between the low scores of the contaminants and the high
#'   scores of the rest. With `method = 'either'` or `'both'`, the
#'   frequency and prevalence scores are drawn over each other, each with
#'   its threshold.
#' * `type = 'prevalence'`: each feature's prevalence in the negative
#'   controls against its prevalence in the true samples (the share of each
#'   it is present in), with the features called highlighted. Contaminants
#'   sit above the diagonal. Needs a result made with `neg`.
#'
#' @param x  A [decontam()] or [decontam_not()] result.
#'
#' @param type  `'scores'` or `'prevalence'`.
#'
#' @param legend  Where the legend goes, as a [graphics::legend()] keyword
#'        (`"topleft"`, ...), or `NULL` for none.
#'
#' @param ...  Passed on to [graphics::hist()] (`'scores'`; such as
#'        `breaks`, `main` or `col`) or [graphics::plot()] (`'prevalence'`;
#'        such as `main`, `pch` or `cex`).
#'
#' @return `x`, invisibly.
#'
#' @seealso [decontam_frequency_plot()]
#'
#' @export
#' @examples
#'     set.seed(1)
#'     conc   <- c(exp(seq(1, 4, length.out = 10)), 0.5, 0.6)
#'     neg    <- rep(c(FALSE, TRUE), c(10, 2))
#'     counts <- matrix(rpois(12 * 8, 200), 12, 8)
#'     counts[, 1:2] <- round(counts[, 1:2] * 5 / conc)
#'     counts[11:12, 3:8] <- 0
#'
#'     res <- decontam(counts, conc, neg)
#'     plot(res)
#'     plot(res, 'prevalence')
#'     plot(decontam(counts, conc, neg, method = 'either'))
plot.ecodive_decontam <- function (x, type = c('scores', 'prevalence'), legend = 'topleft', ...) {
  type <- test_plot_type(type, c('scores', 'prevalence'))
  if (type == 'scores') decontam_plot_scores(x, legend, ...)
  else                  decontam_plot_prevalence(x, legend, ...)
  invisible(x)
}


#' @rdname plot.ecodive_decontam
#' @export
plot.ecodive_decontam_not <- plot.ecodive_decontam


# Whether x is a decontam_not() result, the name of its call column, and
# what its title calls it.
decontam_plot_info <- function (x) {
  not  <- inherits(x, 'ecodive_decontam_not')
  test <- attr(x, 'test', exact = TRUE)
  if (is.null(test)) test <- 'prevalence'
  list(
    not   = not,
    test  = test,
    call  = if (not) x$not_contaminant else x$contaminant,
    what  = if (not) 'non-contaminants' else 'contaminants',
    label = if (not) 'Not a contaminant' else 'Contaminant',
    name  = sprintf('%s, %s', if (not) 'decontam_not' else 'decontam', test) )
}


decontam_plot_scores <- function (x, legend, ...) {

  info <- decontam_plot_info(x)
  thr  <- attr(x, 'threshold', exact = TRUE)
  args <- list(...)
  pos  <- plot_legend_pos(legend)

  two <- info$test %in% c('either', 'both')
  vals <- if (two) list(x$p_freq, x$p_prev) else list(x$p)
  vals <- lapply(vals, function (v) v[is.finite(v)])
  if (all(lengths(vals) == 0))
    stop('There are no scores to plot: every score is NA.', call. = FALSE)

  xlab <- if (info$not) 'Non-contaminant score (1 - prevalence score)'
          else switch(
            info$test,
            frequency  = 'Frequency score',
            prevalence = 'Prevalence score',
            combined   = "Combined score (Fisher's method)",
            minimum    = 'Lower of the frequency and prevalence scores',
            'Score' )
  sub <- sprintf(
    '%d of %d features called %s (threshold %s)', sum(info$call), nrow(x), info$what,
    paste(format(unique(thr)), collapse = ' and ') )
  main <- sprintf(
    '%s scores (%s)', if (info$not) 'Non-contaminant' else 'Contaminant', info$name )

  brk <- plot_arg(args, 'breaks', seq(0, 1, length.out = 51))
  args <- plot_drop(args, 'breaks')

  if (!two) {
    plot_do(
      graphics::hist, c(list(vals[[1]], breaks = brk), args),
      list(main = main, sub = sub, xlab = xlab, col = 'grey85', border = 'grey45', las = 1) )
    graphics::abline(v = thr, col = PLOT_HIGHLIGHT, lwd = 2, lty = 2)
    return (invisible(NULL))
  }

  cols <- plot_palette(2L)
  h    <- lapply(vals, function (v) graphics::hist(v, breaks = brk, plot = FALSE))
  ymax <- max(1, unlist(lapply(h, `[[`, 'counts')))
  plot_do(
    graphics::hist, c(list(vals[[1]], breaks = brk), args),
    list(
      main = main, sub = sub, xlab = xlab, ylim = c(0, ymax), las = 1,
      col = grDevices::adjustcolor(cols[[1]], 0.5), border = cols[[1]] ))
  graphics::hist(
    vals[[2]], breaks = brk, add = TRUE,
    col = grDevices::adjustcolor(cols[[2]], 0.5), border = cols[[2]] )
  graphics::abline(v = thr, col = cols, lwd = 2, lty = 2)
  if (!is.null(pos))
    plot_legend(pos, c('Frequency', 'Prevalence'), cols, 15)
  invisible(NULL)
}


decontam_plot_prevalence <- function (x, legend, ...) {

  info  <- decontam_plot_info(x)
  n_neg <- attr(x, 'n_neg', exact = TRUE)
  n_pos <- attr(x, 'n_pos', exact = TRUE)
  if (all(is.na(x$prev_neg)) || !isTRUE(n_neg > 0) || !isTRUE(n_pos > 0))
    stop(
      'The prevalence plot needs negative controls and true samples: run ',
      'decontam() with `neg`.', call. = FALSE )

  args <- list(...)
  pos  <- plot_legend_pos(legend)
  col  <- ifelse(info$call, PLOT_HIGHLIGHT, PLOT_MUTED)

  plot_do(
    graphics::plot, c(list(x$prev_pos / n_pos, x$prev_neg / n_neg), args),
    list(
      main = sprintf('Prevalence in controls and samples (%s)', info$name),
      xlab = sprintf('Prevalence in true samples (of %d)', n_pos),
      ylab = sprintf('Prevalence in negative controls (of %d)', n_neg),
      xlim = c(0, 1), ylim = c(0, 1), col = col, pch = 16, las = 1 ))
  graphics::abline(0, 1, lty = 2, col = 'grey60')
  if (!is.null(pos))
    plot_legend(
      pos, c(info$label, if (info$not) 'Contaminant' else 'Not a contaminant'),
      c(PLOT_HIGHLIGHT, PLOT_MUTED), 16 )
  invisible(NULL)
}



#' Plot Feature Frequency against DNA Concentration (decontam)
#'
#' decontam's frequency test in pictures (its `plot_frequency()`): each
#' sample's relative abundance of a feature against its DNA concentration,
#' on log scales, with the two models that [decontam()]'s frequency score
#' compares. A contaminant falls along the line of slope \eqn{-1}
#' (frequency \eqn{\propto 1/C}); a true sequence along the horizontal line
#' (frequency independent of \eqn{C}).
#'
#' @inheritParams decontam
#'
#' @param features  The features to plot: their names, their 1-based
#'        indices, or a logical vector with one value per feature. Each gets
#'        its own panel.
#'
#' @param conc  Each sample's DNA concentration, as for [decontam()]:
#'        required.
#'
#' @param neg  Optional: which samples are negative controls, as for
#'        [decontam()]. Controls are drawn in their own colour and left out
#'        of the fitted lines, as they are out of the frequency test.
#'        Default: `NULL`.
#'
#' @param models  Whether to draw the two fitted lines. Default: `TRUE`.
#'
#' @param log  Whether both axes are on a log scale. Default: `TRUE`.
#'
#' @param legend  Where the legend goes, as a [graphics::legend()] keyword,
#'        or `NULL` for none. Default: `"topright"`.
#'
#' @param ...  Passed on to [graphics::plot()], such as `main`, `pch` or
#'        `cex`.
#'
#' @details
#' Over the true samples (not controls, and not `NA` in `neg`) where the
#' feature is present, the non-contaminant model is the horizontal line at
#' the geometric mean frequency, \eqn{\exp(\overline{\log f})}, and the
#' contaminant model is \eqn{f = \exp(\overline{\log f + \log c}) / c}.
#' Samples where the feature is absent can't be drawn on a log scale, and
#' are left out. Samples with no counts are left out, with a warning, as by
#' [decontam()]. With batches, plot each batch's samples separately.
#'
#' Differences from decontam's `plot_frequency()`: base graphics instead
#' of ggplot2, one panel per feature (`par(mfrow)`, restored on exit); the
#' lines leave out samples whose `neg` is `NA`, as the test does; and the
#' contaminant line is red and solid and the non-contaminant line black and
#' dashed, as decontam draws them.
#'
#' @return A data frame of what was drawn, invisibly: `feature`, `sample`
#'   (the sample's name, or index), `conc`, `neg` and `abundance` (after
#'   `norm`), one row per feature and sample.
#'
#' @seealso [decontam()], [plot.ecodive_decontam()]
#'
#' @export
#' @examples
#'     set.seed(1)
#'     conc   <- c(exp(seq(1, 4, length.out = 10)), 0.5, 0.6)
#'     neg    <- rep(c(FALSE, TRUE), c(10, 2))
#'     counts <- matrix(rpois(12 * 8, 200), 12, 8,
#'       dimnames = list(NULL, paste0('F', 1:8)) )
#'     counts[, 1:2] <- round(counts[, 1:2] * 5 / conc)
#'
#'     # F1 is a contaminant, F5 is not.
#'     decontam_frequency_plot(counts, c('F1', 'F5'), conc, neg)
#'
decontam_frequency_plot <- function (
    counts, features, conc, neg = NULL, data = NULL, margin = 1L,
    norm = 'percent', models = TRUE, log = TRUE, legend = 'topright', ... ) {

  cpus <- 1L
  if (missing(conc)) conc <- NULL
  info <- decontam_input()
  pos  <- plot_legend_pos(legend)
  if (!isTRUE(models) && !isFALSE(models))
    stop('`models` must be TRUE or FALSE.', call. = FALSE)
  if (!isTRUE(log) && !isFALSE(log))
    stop('`log` must be TRUE or FALSE.', call. = FALSE)
  if (is.null(conc))
    stop('`conc` is required: each sample\'s DNA concentration.', call. = FALSE)

  idx  <- feature_set(features, info$features, info$D, 'features')
  keep <- info$keep
  if (!any(keep))
    stop('Every sample in `counts` is empty (no counts).', call. = FALSE)
  if (!all(keep))
    warning(
      'Removed ', sum(!keep), ' sample', if (sum(!keep) == 1) '' else 's',
      ' with zero total counts (or frequency).', call. = FALSE )

  v <- decontam_values(conc, 'conc', info)
  if (!is.numeric(v$val))
    stop(v$what, ' must be numeric: each sample\'s DNA concentration.', call. = FALSE)
  bad <- which(keep & (is.na(v$val) | !is.finite(v$val) | v$val <= 0))
  if (length(bad))
    stop(v$what, ' must be finite and greater than 0, but is not for ',
         decontam_samples(bad, info$labels), call. = FALSE)
  cval <- as.double(v$val)

  negv <- rep(FALSE, info$n)
  if (!is.null(neg)) {
    v <- decontam_values(neg, 'neg', info)
    if (!is.logical(v$val))
      stop(v$what, ' must be logical: TRUE for a negative control, FALSE for a true sample.', call. = FALSE)
    negv <- v$val
  }

  # The selected features, dense, samples in rows, after `norm`.
  counts <- info$counts
  sub <- if (info$margin == 1L) counts[, idx, drop = FALSE] else counts[idx, , drop = FALSE]
  sub <- as.matrix(sub)
  if (info$margin == 2L) sub <- t(sub)
  storage.mode(sub) <- 'double'
  if (info$norm == NORM_PERCENT)
    sub <- sub / sample_totals(counts, info$margin)
  sub <- sub[keep, , drop = FALSE]
  cval <- cval[keep]
  negv <- negv[keep]

  fnames  <- if (is.null(info$features)) as.character(idx) else info$features[idx]
  samples <- if (is.null(info$labels)) which(keep) else info$labels[keep]

  args <- list(...)
  k    <- length(idx)
  if (k > 1L) {
    old <- graphics::par(mfrow = rev(grDevices::n2mfrow(k)))
    on.exit(graphics::par(old), add = TRUE)
  }

  ctrl_col <- plot_palette(2L)[[2]]
  ylab     <- if (info$norm == NORM_PERCENT) 'Frequency' else 'Abundance'
  any_neg  <- any(negv %in% TRUE)

  for (j in seq_len(k)) {

    f    <- sub[, j]
    show <- if (log) f > 0 else rep(TRUE, length(f))
    col  <- ifelse(negv %in% TRUE, ctrl_col, 'black')

    fit_rows <- negv %in% FALSE & f > 0
    lines_ok <- models && sum(fit_rows) > 1L
    if (lines_ok) {
      lf  <- log(f[fit_rows])
      lc  <- log(cval[fit_rows])
      a0  <- mean(lf)
      a1  <- mean(lf + lc)
      cc  <- exp(seq(log(min(cval)), log(max(cval)), length.out = 200))
      yc  <- exp(a1) / cc
    }

    yr <- range(c(f[show], if (lines_ok) c(exp(a0), yc)), finite = TRUE)
    if (!all(is.finite(yr))) yr <- c(0.1, 1)

    plot_do(
      graphics::plot, c(list(cval[show], f[show]), args),
      list(
        main = fnames[[j]], xlab = 'DNA concentration', ylab = ylab,
        log = if (log) 'xy' else '', ylim = yr, col = col[show], pch = 16, las = 1 ))

    if (lines_ok) {
      graphics::lines(cc, yc, col = PLOT_HIGHLIGHT, lwd = 2)
      graphics::abline(h = exp(a0), col = 'black', lwd = 2, lty = 2)
    }

    if (!is.null(pos) && j == 1L && (lines_ok || any_neg)) {
      lab <- c(if (any_neg) c('True sample', 'Negative control'),
               if (lines_ok) c('Contaminant model', 'Non-contaminant model'))
      cl  <- c(if (any_neg) c('black', ctrl_col), if (lines_ok) c(PLOT_HIGHLIGHT, 'black'))
      graphics::legend(
        pos, legend = lab, col = cl,
        pch = c(if (any_neg) c(16, 16), if (lines_ok) c(NA, NA)),
        lty = c(if (any_neg) c(NA, NA), if (lines_ok) c(1, 2)),
        lwd = c(if (any_neg) c(NA, NA), if (lines_ok) c(2, 2)),
        bg = 'white', cex = 0.8, inset = 0.01 )
    }
  }

  out <- data.frame(
    feature   = rep(fnames, each = length(samples)),
    sample    = rep(samples, k),
    conc      = rep(cval, k),
    neg       = rep(negv, k),
    abundance = as.vector(sub),
    stringsAsFactors = FALSE )
  invisible(out)
}
