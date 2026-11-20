# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# qmp_rarefy() ports the steps of rarefy_even_sampling_depth() in QMP.R
# (https://github.com/raeslab/QMP), which is:
#
#   BSD 3-Clause License
#
#   Copyright (c) 2017, Raes Lab
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   * Neither the name of the copyright holder nor the names of its
#     contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Modified port: the script's steps (sampling depths, the minimum, each
# sample's target, rarefying, scaling by the rarefied total and the load) are
# kept; ceiling() is not applied (counts must be integers), targets are
# floored with a tolerance for rounding, ecodive's rarefy() C code replaces
# phyloseq::rarefy_even_depth(), and any table class is accepted. See
# plans/qmp.md.



#' Quantitative Microbiome Profiling: Absolute Abundances from a Measured Load
#'
#' Turn relative abundances into estimates of absolute abundance, using each
#' sample's measured total microbial **load** (flow cytometry cell counts,
#' total 16S qPCR, a spike-in, DNA yield, ...).
#'
#' * `qmp_scale()`: "relative abundance \eqn{\times} load" (Jian et al.
#'   2020; most qPCR and spike-in studies). Each sample's counts are divided
#'   by its total and multiplied by its load:
#'   \deqn{q_{ij} = x_{ij} / n_i \times L_i}
#' * `qmp_rarefy()`: quantitative microbiome profiling as published by
#'   Vandeputte et al. (2017). Each sample is first rarefied to the **same
#'   sampling depth relative to its load** (reads per cell), so that every
#'   sample has the same detection limit in absolute units, then scaled to
#'   its load. The steps follow the authors' script, `QMP.R`:
#'   1. Each sample's sampling depth is \eqn{s_i = n_i / L_i}, its total
#'      number of reads over its load.
#'   2. \eqn{s_{min} = \min_i s_i}.
#'   3. Sample \eqn{i} is rarefied, without replacement, to
#'      \eqn{t_i = \lfloor s_{min} L_i \rfloor} reads (so \eqn{t_i \le n_i};
#'      the sample with the smallest \eqn{s_i} keeps all its reads).
#'   4. The rarefied counts \eqn{r_{ij}} are scaled to the load:
#'      \eqn{q_{ij} = r_{ij} / t_i \times L_i}.
#'
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features); see
#'        **Input Types**. `qmp_rarefy()` needs whole numbers (integer
#'        counts).
#'
#' @param load  The total microbial load of each sample, in any absolute
#'        unit (cells per gram, 16S copies per µl, ...): finite and greater
#'        than 0. A numeric vector with one value per sample, matched to the
#'        samples by name when both are named, else in the samples' order.
#'        Or, when `data` is given, the name of one of its columns (or a
#'        vector with one value per row of `data`).
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `load` as a column. Rows are matched to the samples
#'        in `counts` by row name, and rows for other samples are ignored;
#'        without row names, they must be in the same order as the samples.
#'        Default: `NULL`.
#'
#' @param seed  An integer seed for the random rarefaction. The same seed
#'        gives the same result, for any `cpus`. Default: `0`.
#'
#' @return A table of absolute abundances (doubles), in the form of `counts`,
#'   as [normalize()] returns it: the same orientation and dimnames, a base
#'   `matrix` for a base `matrix`, the same \pkg{Matrix} class (`dgCMatrix`,
#'   `dgTMatrix`, `dgeMatrix`) or a \pkg{slam} `simple_triplet_matrix` for
#'   those, with zeros not stored; and a `phyloseq`, `rbiom`, or
#'   `SummarizedExperiment` object gives its count matrix. Each sample's
#'   values sum to its `load`.
#'
#'   `qmp_rarefy()`'s result also has attributes `sampling_depth` (each
#'   sample's \eqn{s_i}) and `target` (its \eqn{t_i}, the number of reads
#'   it was rarefied to), named by sample. They are attached to every class,
#'   sparse ones too, but subsetting or converting a \pkg{Matrix} or
#'   \pkg{slam} object drops them, so read them first.
#'
#' @section Before and after:
#'
#' * **Neither method makes the data less noisy.** The error of the load
#'   measurement carries straight into every taxon (Vandeputte et al. 2017;
#'   Props et al. 2017): measure load in replicate.
#' * **16S copy-number correction** (dividing each taxon's counts by its 16S
#'   gene copy number) is a separate, earlier step in Vandeputte's pipeline,
#'   and isn't part of either function. Do it with `sweep()` (see the
#'   examples). The authors' script then rounds the corrected table **up**
#'   with `ceiling()`; `qmp_rarefy()` requires whole numbers instead, so call
#'   `qmp_rarefy(ceiling(x), load)` to do as the script does.
#' * **Spike-ins:** the load is usually
#'   `spike_added * total_reads / spike_reads`. Drop the spike-in feature
#'   from `counts` before calling `qmp_*()`.
#' * **The results are absolute abundances, not compositions.** Don't follow
#'   them with a centred log ratio (`norm = 'clr'`), which removes the load
#'   again. Bray-Curtis on QMP data is the usual next step: [bray()] or
#'   [beta_div()] with `norm = 'none'`.
#'
#' @section Predicted loads:
#'
#' Without a measured load, one can be predicted from a fecal profile with
#' the MLP tool of Nishijima et al. (2025;
#' \url{https://microbiome-tools.embl.de/mlp/}; an R package not on CRAN):
#' `pred <- MLP(profile, profiler, training_data, "load")` gives cells per
#' gram, to pass as `load`: `qmp_scale(counts, pred$load)`, with `counts`'
#' samples in the same order as `profile`'s. MLP's own `"qmp"` output is
#' the same calculation as `qmp_scale()`. But:
#'
#' * **A prediction is not a measurement.** It is computed from the same
#'   relative profile, so it adds no independent information, its errors
#'   are tied to composition, and it can't confirm a compositional result
#'   (ecodive's reasoning). Predicted and measured loads correlated at
#'   \eqn{r \approx} 0.67–0.68 within, and 0.56 across, the training
#'   cohorts.
#' * **Human stool only**, profiled with the supported tools (mOTUs,
#'   MetaPhlAn, DADA2 + RDP), and trained on European cohorts.
#' * That load is "the major determinant of gut microbiome variation" is a
#'   reason to measure or adjust for load, not a validation of any one
#'   prediction.
#'
#' @section Differences from QMP.R:
#'
#' `qmp_rarefy()` follows `rarefy_even_sampling_depth()` in the authors'
#' `QMP.R` (BSD-3, © 2017 Raes Lab), except:
#'
#' * It doesn't apply `ceiling()` to `counts`; non-integer counts are an
#'   error.
#' * Targets \eqn{s_{min} L_i} are floored (as the script's call to
#'   `phyloseq::rarefy_even_depth()` does), but a value within
#'   \eqn{10^{-9}} (relative) of a whole number counts as that number, so
#'   that floating-point rounding never costs a read; the sample with the
#'   smallest sampling depth keeps all its reads.
#' * A sample whose target is 0 reads (a load so small that
#'   \eqn{s_{min} L_i < 1}) is an error, naming it, where the script returns
#'   `NaN`s. So is an empty sample.
#' * Rarefaction uses ecodive's own [rarefy()] code and `seed`, so the
#'   random draws differ from the script's (`rngseed = 711`), while the
#'   sampling depths, targets, and totals are the same.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Vandeputte, D., Kathagen, G., D'hoe, K., et al. (2017). Quantitative
#' microbiome profiling links gut community variation to microbial load.
#' *Nature*, 551, 507–511. \doi{10.1038/nature24460}
#'
#' Jian, C., Luukkonen, P., Yki-Järvinen, H., Salonen, A., & Korpela, K.
#' (2020). Quantitative PCR provides a simple and accessible method for
#' quantitative microbiota profiling. *PLoS ONE*, 15, e0227285.
#' \doi{10.1371/journal.pone.0227285}
#'
#' Nishijima, S., Stankevic, E., Aasmets, O., et al. (2025). Fecal microbial
#' load is a major determinant of gut microbiome variation and a confounder
#' for disease associations. *Cell*, 188, 222–236.
#' \doi{10.1016/j.cell.2024.10.022}
#'
#' Props, R., Kerckhof, F.-M., Rubbens, P., et al. (2017). Absolute
#' quantification of microbial taxon abundances. *The ISME Journal*, 11,
#' 584–587. \doi{10.1038/ismej.2016.117}
#'
#' @seealso [normalize()], [rarefy()], and the compositional methods, which
#'   recover only ratios when no load was measured: [alr()], [log_ratio()],
#'   [rho_p()], [replace_zeros()].
#'
#' @name qmp
#' @export
#' @examples
#'     # Cell counts per gram for each sample of ex_counts.
#'     load <- c(Saliva = 2.0e10, Gums = 4.1e10, Nose = 1.0e10, Stool = 3.2e10)
#'
#'     # Relative abundance x load
#'     q <- qmp_scale(ex_counts, load)
#'     q[1:3, 1:3]
#'     all.equal(rowSums(q), load)
#'
#'     # Vandeputte et al.'s QMP: rarefy to an even sampling depth first
#'     q <- qmp_rarefy(ex_counts, load)
#'     q[1:3, 1:3]
#'     attr(q, 'sampling_depth')[1:3]  # reads per cell
#'     attr(q, 'target')[1:3]          # reads kept
#'
#'     # The load as a column of sample metadata
#'     meta <- data.frame(cells = load, row.names = names(load))
#'     q <- qmp_scale(ex_counts, 'cells', data = meta)
#'
#'     # 16S copy-number correction comes first, as in the authors' script
#'     copies <- seq(1, 4, length.out = ncol(ex_counts))  # per taxon
#'     x <- sweep(ex_counts, 2, copies, '/')
#'     q <- qmp_rarefy(ceiling(x), load)
#'
#'     # Usually followed by Bray-Curtis distances
#'     dm <- bray(qmp_rarefy(ex_counts, load))
#'
qmp_scale <- function (
    counts, load, data = NULL, margin = 1L, cpus = n_cpus() ) {

  if (missing(load)) load <- NULL
  info <- qmp_input()
  res  <- .Call(C_qmp_scale, counts, margin, info$load, info$dense, cpus)
  qmp_empty(attr(res, 'totals', exact = TRUE), info$labels)

  table_output(res, counts, info$dims, info$dense)
}


#' @rdname qmp
#' @export
qmp_rarefy <- function (
    counts, load, data = NULL, margin = 1L, seed = 0L, cpus = n_cpus() ) {

  if (missing(load)) load <- NULL
  info <- qmp_input()
  validate_seed()
  assert_integer_counts()

  if (is.matrix(counts) && !is.double(counts))
    storage.mode(counts) <- 'double'
  else if (inherits(counts, 'simple_triplet_matrix') && !is.double(counts$v))
    counts$v <- as.double(counts$v)

  # Steps 1-3: sampling depths and the number of reads to keep.
  totals <- sample_totals(counts, margin)
  qmp_empty(totals, info$labels)

  load   <- info$load
  depths <- totals / load
  target <- min(depths) * load

  # Floored, but a hair's breadth below a whole number is that number: the
  # sample with the smallest depth gives (n / L) * L, which can be n - 1e-12.
  near   <- round(target)
  whole  <- abs(target - near) <= 1e-9 * pmax(1, near)
  target <- ifelse(whole, near, floor(target))
  target <- pmin(target, totals)

  if (any(target < 1)) {
    bad <- which(target < 1)
    stop(
      sprintf(
        paste0(
          '%d sample%s would be rarefied to 0 reads (%s): %s load is too ',
          'small for the smallest sampling depth (%g reads per unit of load). ',
          'Check %s load, or remove %s.' ),
        length(bad), if (length(bad) == 1) '' else 's',
        name_list(if (is.null(info$labels)) bad else info$labels[bad]),
        if (length(bad) == 1) 'its' else 'their', min(depths),
        if (length(bad) == 1) 'its' else 'their',
        if (length(bad) == 1) 'it' else 'them' ),
      call. = FALSE )
  }

  # Step 3: rarefy each sample to its own target (as rarefy(depth =)).
  target <- as.integer(target)
  rare   <- .Call(C_rarefy, counts, target, seed, margin, cpus)

  # Step 4: scale by the rarefied total (= target) and the load.
  res <- .Call(C_qmp_scale, rare, margin, load, info$dense, cpus)
  res <- table_output(res, rare, info$dims, info$dense)

  attr(res, 'sampling_depth') <- stats::setNames(depths, info$labels)
  attr(res, 'target')         <- stats::setNames(target, info$labels)
  res
}



#------------------------------------------------------------------------------
# Checks qmp_scale()'s and qmp_rarefy()'s `counts`, `margin`, `cpus`, `load`
# and `data` (in the caller's environment, as validate_*() do). Returns
# list(load (doubles in the samples' order), labels, n, dims, dense).
#------------------------------------------------------------------------------
qmp_input <- function (env = parent.frame()) {

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

  smp  <- count_samples(counts, env$margin)
  dims <- switch(
    get_matrix_package(counts),
    base = dim(counts), slam = c(counts$nrow, counts$ncol), Matrix = counts@Dim )

  list(
    load   = qmp_load(env$load, env$data, smp$labels, smp$n),
    labels = smp$labels,
    n      = smp$n,
    dims   = dims,
    dense  = is.matrix(counts) || inherits(counts, 'dgeMatrix') )
}


# `load` (and `data`) by per_sample()'s rule, as `groups` in the tests: one
# finite value > 0 per sample, as doubles in the samples' order.
qmp_load <- function (load, data, labels, n) {

  if (is.null(load))
    stop('`load` is required: one value per sample, or a column name of `data`.',
         call. = FALSE)

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, dn, nrow(data), ref = '`counts`')
  }

  src <- if (!is.null(data) && is.character(load) && length(load) == 1 &&
             load %in% names(data)) load
  val <- per_sample(load, 'load', data, idx, labels, n, ref = '`counts`')
  what <- if (is.null(src)) '`load`' else paste0('`data` column ', name_list(src))

  if (!is.numeric(val) || is.factor(val))
    stop(what, ' must be numeric: each sample\'s total microbial load.', call. = FALSE)

  bad <- which(is.na(val) | !is.finite(val) | val <= 0)
  if (length(bad))
    stop(
      sprintf(
        '%s must be finite and greater than 0, but is not for %d sample%s: %s',
        what, length(bad), if (length(bad) == 1) '' else 's',
        name_list(if (is.null(labels)) bad else labels[bad]) ),
      call. = FALSE )

  stats::setNames(as.double(val), labels)
}


# Error for samples with no counts, which have no composition to scale.
qmp_empty <- function (totals, labels) {
  bad <- which(totals <= 0)
  if (length(bad))
    stop(
      sprintf(
        '%d sample%s in `counts` %s empty (no counts): %s',
        length(bad), if (length(bad) == 1) '' else 's',
        if (length(bad) == 1) 'is' else 'are',
        name_list(if (is.null(labels)) bad else labels[bad]) ),
      call. = FALSE )
}


# Each sample's total, for a validated count table.
sample_totals <- function (counts, margin) {
  as.vector(switch(
    get_matrix_package(counts),
    base   = if (margin == 1L) rowSums(counts) else colSums(counts),
    slam   = if (margin == 1L) slam::row_sums(counts) else slam::col_sums(counts),
    Matrix = if (margin == 1L) Matrix::rowSums(counts) else Matrix::colSums(counts) ))
}
