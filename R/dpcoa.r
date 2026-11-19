# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Double Principal Coordinate Analysis
#'
#' Ordinates samples and taxa in one space, using the phylogeny (or any
#' Euclidean dissimilarity) among the taxa: DPCoA (Pavoine, Dufour & Chessel
#' 2004). The taxa are placed by a principal coordinate analysis of their
#' dissimilarities, each sample sits at the abundance-weighted centre of its
#' taxa, and the axes are the samples' principal axes. Distances between
#' samples are those of Rao's quadratic entropy, a phylogenetic diversity.
#' This is ade4's `dpcoa()` and phyloseq's `ordinate(method = "DPCoA")`.
#'
#' @inheritParams documentation
#'
#' @param counts   A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} taxa, see `margin`), or any of the other objects
#'        listed in "Input Types". Values need not be integers. Every sample
#'        must have a positive total. The taxa must have names, to match
#'        them to the tips of `tree` or the labels of `taxa_dm`.
#'
#' @param tree   A `phylo` object with branch lengths, whose tips include
#'        every taxon of `counts` (see [read_tree()]). Tips that aren't taxa
#'        of `counts` are dropped. Can be omitted if a tree is embedded with
#'        the `counts` object or as `attr(counts, 'tree')`. The taxa
#'        dissimilarities are the square roots of the distances along the
#'        tree (patristic distances), which are always Euclidean.
#'
#' @param taxa_dm   Instead of `tree`: the dissimilarities among the taxa, as
#'        a `dist` object or symmetric matrix, used as given. They must be
#'        Euclidean (for instance, the square roots of patristic distances,
#'        or Euclidean distances between trait values). Labels are matched to
#'        the taxa of `counts`, and taxa not in `counts` are dropped; without
#'        labels, rows are taken in the order of `counts`' taxa. An embedded
#'        tree is ignored when `taxa_dm` is given.
#'
#' @param n_axes   The number of axes to return. Default: `3`
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   named `"DPCoA1"`, `"DPCoA2"`, etc., of class `"ecodive_ord"` (see
#'   "Ordination result" below), with attributes:
#'
#'   * `method`: `"dpcoa"`.
#'   * `eig`: The eigenvalues: each axis's weighted variance of the samples,
#'     as ade4's `eig`.
#'   * `prop`: Each axis's share of the total inertia, `trace`.
#'   * `taxa`: The taxa coordinates on the same axes (taxa \eqn{\times} axes,
#'     ade4's `dls`), named by taxon. A taxon absent from every sample
#'     doesn't shape the axes (its weight is 0), but is placed on them too.
#'   * `rao`: Each sample's Rao quadratic entropy (ade4's `RaoDiv`), named
#'     by sample.
#'   * `trace`: The total inertia: the weighted sum of the samples' squared
#'     distances from their centre. It is the among-sample part of Rao's
#'     apportionment of diversity: the Rao entropy of all samples pooled,
#'     less the weighted mean of `rao`.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{n_i} the total of sample \eqn{i}, the samples have weights
#' \eqn{w_i = n_i / \sum n} and relative abundances \eqn{p_{ik}}, and the
#' taxa have weights \eqn{\bar{p}_k = \sum_i w_i p_{ik}}, their shares of
#' the grand total. With \eqn{\delta_{kl}} the dissimilarity between taxa
#' \eqn{k} and \eqn{l}:
#'
#' 1. The taxa are placed at points \eqn{x_k} with
#'    \eqn{|x_k - x_l| = \delta_{kl}}, centred on their \eqn{\bar{p}}-weighted
#'    mean: a PCoA of \eqn{\delta} weighted by \eqn{\bar{p}}, with every
#'    axis kept.
#' 2. Sample \eqn{i} is placed at \eqn{y_i = \sum_k p_{ik} x_k}.
#' 3. The axes are the principal axes of the \eqn{y_i}, weighted by
#'    \eqn{w}: the eigenvectors of \eqn{\sum_i w_i y_i y_i'}, whose
#'    eigenvalues are `eig`. The samples' coordinates are their
#'    \eqn{y_i}, and the taxa's their \eqn{x_k}, projected on the axes.
#'
#' Rao's quadratic entropy of sample \eqn{i} is
#' \eqn{\frac{1}{2} \sum_{k,l} p_{ik} p_{il} \delta_{kl}^2}, and the
#' squared distance between samples \eqn{i} and \eqn{j} (over all axes) is
#' \eqn{\sum_{k,l} p_{ik} p_{jl} \delta_{kl}^2 - \mathrm{rao}_i -
#' \mathrm{rao}_j}.
#'
#' With a tree, \eqn{\delta_{kl}^2} is the length of the path between tips
#' \eqn{k} and \eqn{l}. So Rao's entropy is half the expected path length
#' between two reads drawn at random from the sample, and is computed over
#' the branches: for a branch of length \eqn{L} with a share \eqn{a} of
#' the sample below it, \eqn{L a (1 - a)}, summed over the branches.
#' Polytomies, the root's position and nodes with one child don't change
#' path lengths, so the tree may be unrooted or not binary.
#'
#' `dpcoa()` can also be run as `count_ord(counts, 'dpcoa', tree = tree)`.
#'
#'
#' **Speed and memory**
#'
#' No taxa \eqn{\times} taxa matrix is formed, and with a tree, nor are the
#' taxa's points: block Lanczos finds the top `n_axes` eigenpairs of the
#' samples' (weighted) inner products from products computed from the
#' sparse table. With a tree, each product takes time proportional to the
#' table's non-zero values plus the tree's size, from sums up and down the
#' tree. With `taxa_dm`, each product passes over the packed distances, and
#' they are first checked to be Euclidean in the same way (their smallest
#' PCoA eigenvalue may not fall below \eqn{-\sqrt{\epsilon}} times their
#' trace). Results are identical from run to run and for any number of
#' `cpus`.
#'
#'
#' **Differences from ade4**
#'
#' The sample and taxa coordinates, eigenvalues and Rao entropies agree with
#' ade4 1.7-24's `dpcoa(df, dis)` (with `dis` the square roots of the
#' patristic distances, for a tree) to about \eqn{10^{-12}}, up to each
#' axis's sign. Differences:
#'
#' * ade4 stops with an error when a taxon is absent from every sample;
#'   `dpcoa()` keeps it, with weight 0 (the results for the other taxa are
#'   ade4's without it), and gives its coordinates too.
#' * Axes with eigenvalues at rounding-error level, below
#'   `sqrt(.Machine$double.eps)` times the first, are dropped with a warning.
#'   `n_axes` is reduced to at most one less than the number of samples, or
#'   of taxa.
#' * ade4's taxa scores on its taxa PCoA axes (`c1`), its dissimilarities
#'   between samples (`RaoDis`) and its diversity table (`RaoDecodiv`) are
#'   not returned; `trace` is the table's between-samples diversity, and
#'   the within-samples diversity is `weighted.mean(rao, w)`.
#' * phyloseq's DPCoA applies a Cailliez correction to the distances when
#'   `ade4::is.euclid()` finds them non-Euclidean. The square roots of
#'   patristic distances are always Euclidean, so that only happens through
#'   rounding error; `dpcoa()` never corrects them.
#'
#'
#' **Same names in other packages**
#'
#' `ade4::dpcoa()` shares its name with `dpcoa()`; whichever package was
#' attached last masks the other's. Write `ecodive::dpcoa()`, or use
#' [count_ord()], to be sure.
#'
#' @references
#' Pavoine, S., Dufour, A.-B., & Chessel, D. (2004). From dissimilarities
#' among species to dissimilarities among communities: a double principal
#' coordinate analysis. *Journal of Theoretical Biology*, 228(4), 523-537.
#' \doi{10.1016/j.jtbi.2004.02.014}
#'
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
#' approach. *Theoretical Population Biology*, 21(1), 24-43.
#' \doi{10.1016/0040-5809(82)90004-1}
#'
#' Purdom, E. (2011). Analysis of a data matrix and a graph: metagenomic
#' data and the phylogenetic tree. *The Annals of Applied Statistics*,
#' 5(4), 2326-2358. \doi{10.1214/10-AOAS402}
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [count_ord()]; [pcoa()] of a UniFrac distance
#'   ([weighted_unifrac()]) for another phylogeny-aware ordination;
#'   [rpca_phylo()]; [read_tree()].
#'
#' @export
#' @examples
#'     ord <- dpcoa(ex_counts, ex_tree, n_axes = 2)
#'     plot(ord)
#'     str(ord)
#'
#'     # Taxa on the same axes, and each sample's Rao entropy.
#'     head(attr(ord, 'taxa'))
#'     attr(ord, 'rao')
#'
#'     # A Euclidean dissimilarity among the taxa instead of a tree:
#'     # distances between (made-up) trait values.
#'     traits <- data.frame(
#'       row.names = colnames(ex_counts),
#'       gram_pos  = c(1, 0, 1, 0, 1, 1),
#'       aerobic   = c(0, 0, 1, 1, 0, 1),
#'       size      = c(0.8, 1.5, 0.7, 1.2, 0.9, 1.0) )
#'     ord <- dpcoa(ex_counts, taxa_dm = dist(traits), n_axes = 2)
#'
dpcoa <- function (
    counts,
    tree    = NULL,
    taxa_dm = NULL,
    margin  = 1L,
    n_axes  = 3L,
    cpus    = n_cpus() ) {

  if (!is.null(tree) && !is.null(taxa_dm))
    stop('Give either `tree` or `taxa_dm`, not both.', call. = FALSE)
  use_dm <- !is.null(taxa_dm)

  validate_counts()
  validate_margin()
  validate_cpus()
  k <- ord_k(n_axes)

  if (use_dm) tree <- NULL  # An embedded tree isn't used.


  # Dimensions and names, with samples first.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  remove('vals')

  n        <- as.integer(dims[[1]])
  D        <- as.integer(dims[[2]])
  samples  <- if (length(dn[[1]])) as.character(dn[[1]])
  features <- if (length(dn[[2]])) as.character(dn[[2]])

  if (D < 2L)
    stop('`counts` must have at least 2 taxa for DPCoA.', call. = FALSE)

  k <- min(k, n - 1L, D - 1L)


  # The tree's edges, or taxa_dm's packed distances, in the taxa's order.
  epar <- echd <- elen <- dm <- NULL

  if (use_dm) {
    dm <- dpcoa_taxa_dm(taxa_dm, features, D)
  }
  else {
    if (is.null(tree))
      stop(
        '`tree` or `taxa_dm` is required: a phylo object whose tips are the ',
        'taxa of `counts` (or a phyloseq, rbiom or TreeSummarizedExperiment ',
        'object that holds one), or the dissimilarities among the taxa.',
        call. = FALSE )
    tr <- philr_tree(tree, features)
    if (is.null(tree$edge.length) || anyNA(tr$edge_len))
      stop('`tree` must have branch lengths, with none missing.', call. = FALSE)
    epar <- tr$edge_par
    echd <- tr$edge_child
    elen <- tr$edge_len
  }


  res <- .Call(C_dpcoa, counts, margin, epar, echd, elen, dm, k, cpus)

  if (!is.null(res$empty))
    stop(cca_empty_msg(res$empty, samples, 'DPCoA'), call. = FALSE)

  if (!is.null(res$neg))
    stop(
      '`taxa_dm` is not Euclidean: its PCoA has a negative eigenvalue (',
      signif(res$neg, 3), ' times its trace), and DPCoA needs Euclidean ',
      'dissimilarities. Square roots of distances are often Euclidean; or ',
      'apply a Lingoes or Cailliez correction first.', call. = FALSE )

  if (!res$converged) # nocov start
    warning('DPCoA did not converge; results may be inexact.', call. = FALSE)
  # nocov end

  keep  <- if (k > 0) ord_keep(res$eig, k) else logical(0)
  eig   <- res$eig[keep]
  trace <- res$trace

  axes <- if (sum(keep) > 0) paste0('DPCoA', seq_len(sum(keep)))
  taxa <- res$taxa[, keep, drop = FALSE]
  dimnames(taxa) <- list(features, axes)

  ord_result(
    res$points[, keep, drop = FALSE], 'dpcoa', samples, prefix = 'DPCoA',
    eig   = eig,
    prop  = if (trace > 0) eig / trace else eig,
    taxa  = taxa,
    rao   = stats::setNames(res$rao, samples),
    trace = trace )
}



# taxa_dm as a packed lower triangle (double) in the order of `features`
# (the count table's taxa, D of them).
dpcoa_taxa_dm <- function (taxa_dm, features, D) {

  v <- tryCatch(
    dist_vector(taxa_dm),
    error = function (e)
      stop(gsub('`dm`', '`taxa_dm`', conditionMessage(e), fixed = TRUE), call. = FALSE) )

  labels <- if (!is.null(v$labels)) as.character(v$labels)

  if (is.null(labels) || is.null(features)) {
    if (v$n != D)
      stop(
        '`taxa_dm` must have one row per taxon of `counts` (', D, '), or ',
        'labels to match them by.', call. = FALSE )
    return (v$dm)
  }

  if (anyDuplicated(labels))
    stop('`taxa_dm` has duplicated labels: ',
      dQuote(labels[anyDuplicated(labels)], FALSE), '.', call. = FALSE)
  if (anyDuplicated(features))
    stop('`counts` has duplicated taxa names: ',
      dQuote(features[anyDuplicated(features)], FALSE), '.', call. = FALSE)

  missing <- setdiff(features, labels)
  if (length(missing))
    stop(
      length(missing), ' tax', if (length(missing) > 1) 'a' else 'on', ' of `counts` ',
      if (length(missing) > 1) 'are' else 'is', ' not in `taxa_dm`: ',
      paste(dQuote(utils::head(missing, 5), FALSE), collapse = ', '),
      if (length(missing) > 5) ', ...', '.', call. = FALSE )

  i <- match(features, labels)
  if (identical(i, seq_len(v$n))) return (v$dm)

  # Reorder (and subset): element (r, c), r > c, of the packed triangle.
  n   <- v$n
  pos <- function (r, c) (c - 1) * n - (c - 1) * c / 2 + (r - c)
  rc  <- which(lower.tri(matrix(0, D, D)), arr.ind = TRUE)
  a   <- i[rc[, 1]]
  b   <- i[rc[, 2]]
  v$dm[pos(pmax(a, b), pmin(a, b))]
}
