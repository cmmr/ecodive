# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The phylogenetic isometric log-ratio transform (PhILR), written from
# Silverman et al. (2017, eLife, CC-BY) and checked against the outputs of
# the philr package (GPL-3), whose code was not read. See plans/philr.md.


# philr's part and ILR weightings.
PHILR_PART_WEIGHTS <- c(
  'uniform', 'gm.counts', 'anorm', 'anorm.x.gm.counts', 'enorm', 'enorm.x.gm.counts' )
PHILR_ILR_WEIGHTS <- c('uniform', 'blw', 'blw.sqrt', 'mean.descendants')



#' Phylogenetic Isometric Log-Ratio Transform (PhILR)
#'
#' Silverman et al.'s PhILR coordinates: an isometric log-ratio transform
#' whose balances follow a phylogeny. Each internal node of the tree gives
#' one coordinate, the log ratio of its two child clades' (weighted)
#' geometric means, so every coordinate has an evolutionary meaning: "this
#' clade against its sister clade".
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of counts or proportions (samples
#'        \eqn{\times} features, see `margin`), or an object of one of the
#'        types in **Input Types**. Its features must have names, matched to
#'        the tips of `tree`. Counts, proportions and the output of
#'        [replace_zeros()] are all accepted, with one caveat: the
#'        `'gm.counts'` part weights assume counts.
#'
#' @param tree  A `phylo` object whose tips include every feature of
#'        `counts`, or `NULL` to use the tree inside `counts` (a phyloseq,
#'        rbiom or TreeSummarizedExperiment object, or `attr(counts,
#'        'tree')`). See **The tree**.
#'
#' @param pseudocount  Added to every value before taking logs. `NULL` (the
#'        default) adds nothing when `counts` has no zeros, and otherwise
#'        uses half the smallest non-zero value, with a warning. See
#'        **Pseudocount**.
#'
#' @param part_weights  The weights \eqn{p_j} of the features (philr's
#'        `part.weights`), which set the reference measure of the simplex.
#'        One of
#'
#'   * `'uniform'` (the default): all 1, the usual ilr;
#'   * `'gm.counts'`: each feature's geometric mean over the samples (of
#'     `counts + pseudocount`);
#'   * `'enorm'`: the Euclidean norm of each feature's values over the
#'     samples, after closing them to sum to 1 (so that a feature found in
#'     only a few samples gets a weight near 1, and one spread evenly over
#'     \eqn{n} samples gets \eqn{1/\sqrt{n}});
#'   * `'anorm'`: the Aitchison norm of each feature's values over the
#'     samples, \eqn{\sqrt{\sum_k (\log x_{kj} - \overline{\log x_j})^2}};
#'   * `'enorm.x.gm.counts'` and `'anorm.x.gm.counts'`: products of the two;
#'   * a vector of positive numbers, one per feature, named by feature (or
#'     in the order of the features of `counts`).
#'
#'   Names may also be written with underscores (`'gm_counts'`). The weights
#'   computed from `counts` depend on every sample, so one sample's
#'   coordinates depend on the others.
#'
#' @param ilr_weights  The weights \eqn{f_i} of the coordinates (philr's
#'        `ilr.weights`), which scale each balance. One of
#'
#'   * `'uniform'` (the default): all 1;
#'   * `'blw'`: the sum of the node's two child branch lengths;
#'   * `'blw.sqrt'`: its square root;
#'   * `'mean.descendants'`: the sum of the two child branch lengths plus,
#'     for each child, its mean distance to its tips;
#'   * a vector of non-negative numbers, one per coordinate, named by
#'     coordinate (or in the order of the coordinates).
#'
#'   For the branch-length weights, zero-length branches to tips first get
#'   the smallest positive branch length of the tree, as in philr. They need
#'   `tree` to have branch lengths.
#'
#' @return A dense numeric `matrix` of the \eqn{D - 1} coordinates, one per
#'   internal node of the tree, in the orientation of `counts`: samples
#'   \eqn{\times} coordinates for `margin = 1`, coordinates \eqn{\times}
#'   samples for `margin = 2`. Sample names are kept. Coordinates are named
#'   by the tree's node labels when every node has a distinct one (nodes
#'   added to resolve a polytomy at node `"X"` get `"X.1"`, `"X.2"`, ...),
#'   and otherwise `n1`, `n2`, ..., in column order. Attributes:
#'
#'   * `method`: `'philr'`.
#'   * `part_weights`: the part weights used, named by feature, in the order
#'     of the features of `counts`.
#'   * `ilr_weights`: the ILR weights used, named by coordinate.
#'   * `tree`: the tree used, as a `phylo` object (see **The tree**). Its
#'     internal nodes are numbered and labelled in column order (node
#'     \eqn{D + k} is coordinate \eqn{k}), and for each node, its first
#'     child in `tree$edge` is the numerator ("+") clade.
#'   * `pseudocount`: the pseudocount used (`0` for none).
#'
#'   These are what philr's `return.all = TRUE` gives as `p`, `ilr.weights`
#'   and (through the tree) `sbp`. Its contrast matrix `V` is not kept: it
#'   has \eqn{D (D - 1)} entries.
#'
#' @details
#' For a sample \eqn{x} (with the pseudocount added) and part weights
#' \eqn{p}, write \eqn{y_j = x_j / p_j}. At an internal node with clades
#' \eqn{A} (its first child, the numerator) and \eqn{B} (its second), let
#' \eqn{n_A = \sum_{j \in A} p_j} and \eqn{\log g_A = \sum_{j \in A} p_j
#' \log y_j / n_A} (and likewise for \eqn{B}). The node's coordinate is
#' \deqn{f \sqrt{\frac{n_A n_B}{n_A + n_B}}\,\left(\log g_A - \log g_B\right),}
#' with \eqn{f} its ILR weight (Silverman et al. 2017, Methods). It is
#' computed in \eqn{O(D)} per sample, by summing over each clade from the
#' tips up, without forming the \eqn{D \times (D - 1)} basis.
#'
#' With uniform weights the coordinates are an ordinary ilr: they equal
#' `ilr(counts, basis = sbp)` for the tree's sequential binary partition, and
#' `euclidean(philr(x, tree))` is `aitchison(x)` (with the same
#' pseudocount), as for any orthonormal basis. It is the weights that make
#' distances phylogenetic: [philr_distance()] uses the paper's.
#'
#' The coordinates feed straight into [pca()], [rda()], [lda()], [envfit()]
#' and [euclidean()], which take negative values with `norm = 'none'`.
#'
#' `philr(x, tree, pseudocount = c, part_weights = pw, ilr_weights = iw)`
#' equals philr 1.28.0's `philr(x, tree, part.weights = pw, ilr.weights = iw,
#' pseudocount = c)` to rounding, for every weighting and for a binary tree
#' whose tips are the features of `x`.
#'
#' @section The tree:
#' PhILR needs a rooted binary tree whose tips are the features. `tree` is
#' made one in three steps (the result is the `tree` attribute):
#'
#' 1. **Tips that aren't features of `counts` are dropped**, and so are the
#'    internal nodes left without any. Every feature must be a tip.
#' 2. **Nodes left with one child are removed**, their branch length added
#'    to their child's (a root with one child is dropped), as
#'    `ape::drop.tip()` does.
#' 3. **Polytomies are resolved, with a warning:** a node with children
#'    \eqn{c_1, \dots, c_k}, \eqn{k > 2}, keeps \eqn{c_1} and gets a new
#'    child node, joined by a zero-length branch, that holds \eqn{c_2, \dots,
#'    c_k} (and so on), as `ape::multi2di(random = FALSE)` does. An unrooted
#'    tree (three children at the root) is resolved the same way. Results
#'    then depend on the order of the children: to choose the resolution,
#'    resolve the tree first.
#'
#' Coordinates are in the order of the internal nodes' numbers in `tree`
#' (the root first), or, when a step changed the tree, in preorder of the
#' new tree (which is also the order of `ape::multi2di()`'s and
#' `ape::drop.tip()`'s numbering). Branch lengths are needed only for the
#' branch-length ILR weights; they must be non-negative. A root edge is
#' ignored.
#'
#' **Same names in other packages**: `philr::philr()` computes the same
#' transform, with the same orientation and order. Whichever package was
#' attached last masks the other's `philr()`; write `ecodive::philr()` to be
#' sure.
#'
#' @inheritSection documentation Input Types
#' @inheritSection log_ratio_pseudocount_section Pseudocount
#'
#' @references
#' Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
#' (2017). A phylogenetic transform enhances analysis of compositional
#' microbiota data. *eLife*, 6, e21887. \doi{10.7554/eLife.21887}
#'
#' Egozcue, J. J., & Pawlowsky-Glahn, V. (2016). Changing the reference
#' measure in the simplex and its weighting effects. *Austrian Journal of
#' Statistics*, 45(4), 25–44. \doi{10.17713/ajs.v45i4.126}
#'
#' @seealso [philr_distance()], [ilr()], [aitchison()]
#'
#' @export
#' @examples
#'     x <- philr(ex_counts, tree = ex_tree, pseudocount = 1)
#'     dim(x)
#'     x[1:3, 1:4]
#'
#'     # With uniform weights, distances are Aitchison distances
#'     all.equal(
#'       as.vector(euclidean(x)),
#'       as.vector(aitchison(ex_counts, pseudocount = 1)) )
#'
#'     # The paper's weights
#'     x <- philr(
#'       ex_counts, tree = ex_tree, pseudocount = 1,
#'       part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt' )
#'     head(attr(x, 'ilr_weights'))
#'
philr <- function (
    counts, tree = NULL, margin = 1L, pseudocount = NULL,
    part_weights = 'uniform', ilr_weights = 'uniform', cpus = n_cpus() ) {

  philr_run(part_weights, ilr_weights)
}



#' PhILR Distance
#'
#' Euclidean distances between samples' PhILR coordinates, with the weights
#' that Silverman et al. (2017) recommend: a compositional, phylogenetic
#' beta diversity.
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of counts or proportions (samples
#'        \eqn{\times} features, see `margin`), or an object of one of the
#'        types in **Input Types**. Its features must have names, matched to
#'        the tips of `tree`. The part weights assume counts.
#'
#' @param tree  A `phylo` object with branch lengths, whose tips include
#'        every feature of `counts`, or `NULL` to use the tree inside
#'        `counts` (a phyloseq, rbiom or TreeSummarizedExperiment object, or
#'        `attr(counts, 'tree')`).
#'
#' @param pseudocount  Added to every value before taking logs. `NULL` (the
#'        default) adds nothing when `counts` has no zeros, and otherwise
#'        uses half the smallest non-zero value, with a warning. See
#'        **Pseudocount**.
#'
#' @return A distance matrix: a `dist` object of class `c("ecodive_dist",
#'   "dist")`, as every metric function returns (see [beta_div()]).
#'
#' @family Phylogenetic metrics
#' @seealso [philr()], `beta_div()`, [aitchison()], [weighted_unifrac()]
#'
#' @details
#' The distance between samples \eqn{x} and \eqn{y} is
#' \deqn{\sqrt{\sum_i \left(x^*_i - y^*_i\right)^2},}
#' where \eqn{x^*} are the [philr()] coordinates of \eqn{x} with part weights
#' `'enorm.x.gm.counts'` and ILR weights `'blw.sqrt'` (Silverman et al.
#' 2017, Methods), so that
#'
#' ```r
#' philr_distance(x, tree) ==
#'   euclidean(philr(x, tree, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt'))
#' ```
#'
#' For other weightings, use `euclidean(philr(x, tree, ...))`. With uniform
#' weights, PhILR distances are Aitchison distances ([aitchison()]): the
#' weights are what make them phylogenetic. The part weights are computed
#' from all the samples, so adding or removing a sample changes the
#' distances between the others. `pairs` limits which distances are
#' computed, not the samples used.
#'
#' The tree is handled as [philr()] handles it: tips that aren't features
#' are dropped, and polytomies are resolved, with a warning (see
#' **The tree** in [philr()]). It needs branch lengths.
#'
#' @inheritSection documentation Input Types
#' @inheritSection log_ratio_pseudocount_section Pseudocount
#'
#' @references
#' Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
#' (2017). A phylogenetic transform enhances analysis of compositional
#' microbiota data. *eLife*, 6, e21887. \doi{10.7554/eLife.21887}
#'
#' @export
#' @examples
#'     philr_distance(ex_counts, tree = ex_tree, pseudocount = 1)
#'
#'     # The same, from beta_div()
#'     d <- beta_div(ex_counts, 'philr', tree = ex_tree, pseudocount = 1)
philr_distance <- function (
    counts, tree = NULL, margin = 1L, pseudocount = NULL,
    pairs = NULL, clust = NULL, cpus = n_cpus() ) {

  x <- philr_run('enorm.x.gm.counts', 'blw.sqrt')

  # The coordinates replace the counts; the pseudocount used stays in this
  # frame for the `dist` record.
  pseudocount <- attr(x, 'pseudocount', exact = TRUE)
  counts      <- x
  validate_pairs()
  validate_clust()

  d <- .Call(C_beta_div, BDIV_EUCLIDEAN, counts, margin, NORM_NONE, pairs, cpus, 0, NULL)
  dist_result(d, 'philr_distance')
}



# philr() in the caller's frame, which holds counts, tree, margin,
# pseudocount and cpus. Validates them, then returns the coordinates with
# their attributes; the caller's pseudocount is set to the one used.
philr_run <- function (part_weights, ilr_weights, env = parent.frame()) {

  pw <- philr_weight_arg(part_weights, PHILR_PART_WEIGHTS, 'part_weights')
  iw <- philr_weight_arg(ilr_weights,  PHILR_ILR_WEIGHTS,  'ilr_weights')

  info   <- log_ratio_input(
    'the PhILR transform', env = env,
    check = function (features) philr_tree(env$tree, features) )
  counts <- env$counts
  margin <- env$margin
  pc     <- env$pseudocount
  tr     <- info$check

  if (tr$n_poly > 0L)
    warning(
      sprintf(
        '`tree` has %d %s with more than two children; resolved with zero-length branches, as ape::multi2di(random = FALSE) does. Resolve the tree first to choose how.',
        tr$n_poly, if (tr$n_poly == 1L) 'node' else 'nodes' ),
      call. = FALSE )

  p <- philr_part_weights(pw, counts, margin, pc, info)
  f <- philr_ilr_weights(iw, tr)  # Column order.

  res <- .Call(
    C_philr_coords, counts, margin, tr$c1, tr$c2, tr$col,
    unname(p), unname(f)[tr$col + 1L], pc, env$cpus )

  dimnames(res) <- if (margin == 1L) list(info$samples, tr$names)
                   else              list(tr$names, info$samples)

  attr(res, 'method')       <- 'philr'
  attr(res, 'part_weights') <- p
  attr(res, 'ilr_weights')  <- f
  attr(res, 'tree')         <- tr$phylo
  attr(res, 'pseudocount')  <- pc
  res
}



# A weighting: one of `choices` (with '_' or '.'), or a numeric vector,
# returned as given.
philr_weight_arg <- function (x, choices, arg) {

  if (is.numeric(x)) return (x)

  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    i <- match(gsub('_', '.', x, fixed = TRUE), choices)
    if (!is.na(i)) return (choices[[i]])
  }

  stop(
    '`', arg, '` must be one of ', paste0('"', choices, '"', collapse = ', '),
    ', or a numeric vector.', call. = FALSE )
}



# The tree as C_philr_tree() makes it (see src/philr.c), for the features
# `features` (in the count table's order), plus:
#   col:    each node's 0-based column;
#   names:  the coordinates' names, in column order;
#   phylo:  the tree used, tips in the features' order, internal nodes
#           numbered in column order.
philr_tree <- function (tree, features) {

  if (is.null(tree))
    stop(
      '`tree` is required: a phylo object whose tips are the features of ',
      '`counts` (or a phyloseq, rbiom or TreeSummarizedExperiment object ',
      'that holds one).', call. = FALSE )

  if (!inherits(tree, 'phylo') || !is.matrix(tree$edge) || ncol(tree$edge) != 2L ||
      is.null(tree$tip.label) || nrow(tree$edge) < 1L)
    stop('`tree` is not a valid phylo object.', call. = FALSE)

  if (is.null(features))
    stop(
      '`counts` must have feature names, to match them to the tips of `tree`.',
      call. = FALSE )

  tips <- as.character(tree$tip.label)
  if (anyDuplicated(features))
    stop('`counts` has duplicated feature names: ',
      dQuote(features[anyDuplicated(features)], FALSE), '.', call. = FALSE)
  if (anyDuplicated(tips))
    stop('`tree` has duplicated tip labels: ',
      dQuote(tips[anyDuplicated(tips)], FALSE), '.', call. = FALSE)

  missing <- setdiff(features, tips)
  if (length(missing))
    stop(
      length(missing), ' feature', if (length(missing) > 1) 's', ' of `counts` ',
      if (length(missing) > 1) 'are' else 'is', ' not in `tree`: ',
      paste(dQuote(utils::head(missing, 5), FALSE), collapse = ', '),
      if (length(missing) > 5) ', ...', '.', call. = FALSE )

  edge <- tree$edge
  if (!is.integer(edge)) {
    if (!is.numeric(edge) || any(edge %% 1 != 0, na.rm = TRUE))
      stop('`tree` is not a valid phylo object: malformed edge matrix.', call. = FALSE)
    edge <- matrix(as.integer(edge), ncol = 2L)
  }

  len <- tree$edge.length
  if (!is.null(len)) {
    if (!is.numeric(len))
      stop('`tree`\'s branch lengths must be numbers.', call. = FALSE)
    len <- as.double(len)
  }

  ntip <- length(tips)
  feat <- match(tips, features) - 1L
  feat[is.na(feat)] <- -1L

  tr <- .Call(C_philr_tree, edge, len, ntip, feat)

  D  <- length(features)
  ni <- D - 1L

  # Columns: the input's node numbers, unless the tree changed.
  tr$col <- if (tr$modified) seq_len(ni) - 1L else rank(tr$orig, ties.method = 'first') - 1L


  # Names: the node labels when every node has a distinct one.
  nm  <- NULL
  lab <- tree$node.label
  if (!is.null(lab) && length(lab) == as.integer(tree$Nnode)) {
    lab  <- as.character(lab)
    base <- lab[ifelse(tr$orig > 0L, tr$orig, tr$ladder) - ntip]
    nm   <- ifelse(tr$orig > 0L, base, paste0(base, '.', tr$rung))
    if (anyNA(base) || !all(nzchar(base)) || anyDuplicated(nm)) nm <- NULL
  }
  if (is.null(nm)) nm <- paste0('n', tr$col + 1L)

  tr$names <- character(ni)
  tr$names[tr$col + 1L] <- nm


  # The tree used: tip f + 1 is feature f; node D + k is column k.
  ape_id <- function (id) ifelse(id < D, id + 1L, D + 1L + tr$col[pmax(id - D, 0L) + 1L])
  phy <- list(
    edge      = cbind(ape_id(tr$edge_par), ape_id(tr$edge_child)),
    Nnode     = ni,
    tip.label = features )
  if (!is.null(len)) phy$edge.length <- tr$edge_len
  phy$node.label <- tr$names
  tr$phylo <- structure(phy, class = 'phylo', order = 'cladewise')

  tr
}



# The part weights, named by feature in the count table's order.
philr_part_weights <- function (pw, counts, margin, pc, info) {

  D <- info$D

  if (is.numeric(pw)) {
    p <- philr_named_weights(pw, info$features, D, 'part_weights', 'feature')
    if (!all(is.finite(p) & p > 0))
      stop('`part_weights` must be positive numbers.', call. = FALSE)
    return (stats::setNames(p, info$features))
  }

  if (pw == 'uniform')
    return (stats::setNames(rep(1, D), info$features))

  st <- .Call(C_philr_stats, counts, margin, pc)
  gm <- exp(st[, 3])
  en <- sqrt(st[, 2]) / st[, 1]
  an <- sqrt(st[, 4])

  p <- switch(
    EXPR                = pw,
    'gm.counts'         = gm,
    'enorm'             = en,
    'anorm'             = an,
    'enorm.x.gm.counts' = en * gm,
    'anorm.x.gm.counts' = an * gm )

  if (!all(is.finite(p) & p > 0)) {
    bad <- which(!(is.finite(p) & p > 0))
    stop(
      '`part_weights = "', pw, '"` gives ',
      if (startsWith(pw, 'anorm')) 'a weight of 0 to features that have the same value in every sample'
      else 'weights that are not positive numbers',
      ' (', if (is.null(info$features)) paste(utils::head(bad, 5), collapse = ', ')
            else paste(dQuote(utils::head(info$features[bad], 5), FALSE), collapse = ', '),
      if (length(bad) > 5) ', ...', '). Remove them, or choose other weights.', call. = FALSE )
  }

  stats::setNames(p, info$features)
}



# The ILR weights, named by coordinate in column order.
philr_ilr_weights <- function (iw, tr) {

  ni <- length(tr$c1)

  if (is.numeric(iw)) {
    f <- philr_named_weights(iw, tr$names, ni, 'ilr_weights', 'coordinate')
    if (!all(is.finite(f) & f >= 0))
      stop('`ilr_weights` must be non-negative numbers.', call. = FALSE)
    return (stats::setNames(f, tr$names))
  }

  if (iw == 'uniform') return (stats::setNames(rep(1, ni), tr$names))

  v <- if (iw == 'mean.descendants') tr$md else tr$blw
  if (anyNA(v))
    stop(
      '`ilr_weights = "', iw, '"` needs `tree` to have branch lengths, ',
      'with none missing and at least one positive.', call. = FALSE )
  if (iw == 'blw.sqrt') v <- sqrt(v)

  f <- numeric(ni)
  f[tr$col + 1L] <- v
  stats::setNames(f, tr$names)
}



# A numeric weight vector, matched to `nms` by name when it has names, else
# by position (length n).
philr_named_weights <- function (w, nms, n, arg, what) {

  w <- stats::setNames(as.double(w), names(w))

  if (!is.null(nms) && !is.null(names(w))) {
    miss <- setdiff(nms, names(w))
    if (length(miss))
      stop(
        '`', arg, '` has no value for ', length(miss), ' ', what,
        if (length(miss) > 1) 's', ': ',
        paste(dQuote(utils::head(miss, 5), FALSE), collapse = ', '),
        if (length(miss) > 5) ', ...', '.', call. = FALSE )
    if (anyDuplicated(names(w)))
      stop('`', arg, '` has duplicated names.', call. = FALSE)
    return (unname(w[nms]))
  }

  if (length(w) != n)
    stop(
      '`', arg, '` must have one value per ', what, ' (', n, '), or names.',
      call. = FALSE )
  unname(w)
}
