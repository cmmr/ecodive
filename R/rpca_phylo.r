# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Ported from gemelli 0.0.13 (gemelli/rpca.py, phylogenetic_rpca()), which
# carries this notice:
#
#   BSD 3-Clause License
#
#   Copyright (c) 2019, Cameron Martino
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
# The tree handling (src/rpca_phylo.c) also follows scikit-bio 0.7.4 and iow
# 1.0.8 (BSD-3); see that file and inst/COPYRIGHTS.



#' Phylogenetic Robust Aitchison PCA
#'
#' Robust Aitchison PCA ([rpca()]) of a count table expanded to every node of
#' a phylogenetic tree: each node's count is the total of the tips below it,
#' weighted by its branch length. Clades that change together then act as
#' one feature, and the result is a compositional counterpart of a UniFrac
#' PCoA, with loadings for every clade. This is gemelli's
#' `phylogenetic_rpca()` (Martino et al. 2022), QIIME 2's
#' `gemelli phylogenetic-rpca`.
#'
#' @inheritParams rpca
#' @inheritParams documentation
#'
#' @param counts   A numeric matrix of count data (samples \eqn{\times}
#'        features), or any of the other objects listed in "Input Types".
#'        Values must be finite and non-negative; they need not be integers.
#'        The features are the tips of `tree`.
#'
#' @param tree   A `phylo` object whose tips include every feature of
#'        `counts` (see [read_tree()]). Its branch lengths must be
#'        non-negative; missing ones (`NA`) count as 0. Its root length,
#'        `tree$root.edge`, only counts if the root is left with a single
#'        child (see "Method"). Can be omitted if a tree is embedded with
#'        the `counts` object or as `attr(counts, 'tree')`.
#'
#' @param n_axes   The rank of the low-rank model, which is also the number of
#'        axes returned: a whole number from 1 to 10. Default: `3`
#'
#' @param min_sample_count   Drop samples whose total count (after the
#'        feature filters) is not above this. Default: `0`
#'
#' @param min_feature_count   Drop tips whose total count is not above
#'        this. Default: `0`
#'
#' @param min_feature_frequency   Drop tips that are non-zero in no more
#'        than this percentage of samples (0 to 100). Default: `0`
#'
#' @param min_depth   Drop the internal nodes that have at most this many
#'        tips with counts below them (gemelli's `min_depth`). Tips are
#'        always kept. The default, `0`, keeps every node.
#'
#' @param node_counts   If `TRUE`, also return the count of every node in
#'        every sample (the total of its tips' counts, not weighted by
#'        branch length) in attribute `node_counts`. Default: `FALSE`
#'
#' @return A numeric matrix with one row per kept sample and one column per
#'   axis, named `"RPCA1"`, `"RPCA2"`, etc., with the same attributes as
#'   [rpca()]'s (see "Ordination result" below), except that:
#'
#'   * `method` is `"rpca_phylo"`.
#'   * `loadings` has one row per **node** of the tree that is used: the
#'     tips, named by their labels, and the internal nodes, named
#'     `"n<id>"` as gemelli names them (see "Method").
#'   * `tree`: the tree that was used, as a `phylo` object: `tree` cut down
#'     to the tips with counts, with nodes merged and removed as described in
#'     "Method", and with internal nodes labelled with the names used in
#'     `loadings`, so that loadings can be drawn on it. Its root length, if
#'     not 0, is in `root.edge`.
#'   * `node_counts`: Only with `node_counts = TRUE`. The node counts, with
#'     one row per kept sample and one column per node (in the order of
#'     `loadings`): a sparse `dgCMatrix` if the Matrix package is installed,
#'     else a base matrix.
#'
#' @details
#'
#' **Method**
#'
#' 1. **Filtering.** The tips (features) and samples are filtered exactly as
#'    by [rpca()], except that samples with a single non-zero tip are kept:
#'    they have counts on that tip's ancestors too.
#'
#' 2. **The tree.** The tree is cut down to its root, the tips with counts
#'    after filtering, and their ancestors. A node left with a single child
#'    is merged into it, adding its branch length to the child's. The root
#'    has no length (gemelli drops it), unless it is left with a single
#'    child: then the two become one node, the root, whose length is the
#'    child's (after the merges above) plus the root length,
#'    `tree$root.edge` (trees from [read_tree()] have none). Then, with
#'    `min_depth` above 0, each internal node with at most `min_depth` tips
#'    below it is removed, and its tips are attached, in order, after its
#'    parent's other children.
#'
#' 3. **Node table.** Each node's count in a sample is the total count of
#'    the tips below it (a tip's is its own count). Internal nodes are
#'    numbered as scikit-bio's `TreeNode.assign_ids()` numbers them: visiting
#'    the nodes in postorder, each node's children are numbered in turn, and
#'    the root last; they are named `n` followed by that number (from 0).
#'    The columns of the node table, and the rows of `loadings`, are in
#'    this order.
#'
#' 4. **Robust clr with branch lengths.** For each sample, the log of each
#'    node's count times its branch length, minus their mean over the
#'    sample's observed entries. An entry is observed when both the count
#'    and the branch length are positive. So a node with length 0 (such as
#'    the root, usually) is a feature with no
#'    observed values, which still counts in OptSpace's scaling and in the
#'    final centering, and has loadings. A sample whose observed values are
#'    all equal (for example one with a single tip, whose path to the root
#'    has branches of equal length) has a robust clr of 0 throughout, and
#'    has no observed values: as in gemelli, its position only reflects
#'    the other samples.
#'
#' 5. **Matrix completion and PCA**, exactly as [rpca()].
#'
#' `rpca_phylo()` can also be run as `count_ord(counts, 'rpca_phylo', tree =
#' tree)`.
#'
#'
#' **Speed and memory**
#'
#' The node table is never formed densely: each sample's node counts are
#' summed from its non-zero tips, in parallel, and the rest runs as for
#' [rpca()], whose "Speed and memory" notes apply with the nodes as
#' features. The tree is processed in time proportional to its size. Results
#' are identical from run to run and for any number of `cpus`.
#'
#'
#' **Differences from gemelli**
#'
#' `rpca_phylo()` is a port of `gemelli.rpca.phylogenetic_rpca()` (gemelli
#' 0.0.13), with the tree handling of iow 1.0.8 and scikit-bio 0.7.4 that it
#' calls. On the package's test tables (including gemelli's test tree),
#' the node counts are the same and the results agree with gemelli's to
#' within \eqn{2 \times 10^{-14}} (relative). The differences are those of [rpca()]
#' (its "Differences from gemelli" section: `max_iters` is used where gemelli
#' ignores `max_iterations` and always runs 5 iterations, `n_axes` may be 1,
#' axes of rounding-error variance are dropped, `eig` and the scores are
#' scaled as [pca()]'s), and these:
#'
#' * gemelli converts counts to whole numbers (64-bit integers) before
#'   summing them up the tree, so proportions become zeros. `rpca_phylo()`
#'   uses the counts as given.
#' * Samples with a single non-zero tip are kept, as gemelli keeps them
#'   (unlike [rpca()], which drops samples with a single non-zero feature).
#' * gemelli's OptSpace treats any robust clr value that is exactly 0 as
#'   missing, which depends on rounding. Here an entry is missing exactly
#'   when its count or its branch length is 0, or when all of its sample's
#'   values are equal.
#' * The root length is `tree$root.edge`. ape's `read.tree()` keeps a
#'   Newick root length there, but ecodive's [read_tree()] discards it. It
#'   only matters when the root is left with a single child: then, for a
#'   tree read by [read_tree()], the merged root's length lacks the root
#'   length that gemelli, reading the same Newick file, would add.
#' * Internal node labels in `tree` are replaced by the `n<id>` names, as
#'   gemelli does; gemelli's optional taxonomy output is not provided.
#'
#' @references
#' Martino, C., McDonald, D., Cantrell, K., Dilmore, A. H.,
#' Vázquez-Baeza, Y., Shenhav, L., Shaffer, J. P., Rahman, G., Armstrong,
#' G., Allaband, C., Song, S. J., & Knight, R. (2022). Compositionally aware
#' phylogenetic beta-diversity measures better resolve microbiomes associated
#' with phenotype. *mSystems*, 7(3), e00050-22.
#' \doi{10.1128/msystems.00050-22}
#'
#' Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
#' Knight, R., & Zengler, K. (2019). A novel sparse compositional technique
#' reveals microbial perturbations. *mSystems*, 4(1), e00016-19.
#' \doi{10.1128/mSystems.00016-19}
#'
#' Martino, C. (2019). *gemelli* (Python package, version 0.0.13). BSD
#' 3-Clause License. \url{https://pypi.org/project/gemelli/}
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [rpca()], the same without a tree; [count_ord()]; the UniFrac
#'   family ([unweighted_unifrac()], [weighted_unifrac()], ...), the
#'   phylogenetic distances whose compositional counterpart this is;
#'   [read_tree()].
#'
#' @export
#' @examples
#'     ord <- rpca_phylo(ex_counts, ex_tree, n_axes = 2)
#'     ord
#'     plot(ord)
#'     str(ord)
#'
#'     # Loadings for clades as well as tips.
#'     head(attr(ord, 'loadings'))
#'
rpca_phylo <- function (
    counts,
    tree                  = NULL,
    margin                = 1L,
    n_axes                = 3L,
    max_iters             = 5L,
    tol                   = 1e-5,
    min_sample_count      = 0,
    min_feature_count     = 0,
    min_feature_frequency = 0,
    min_depth             = 0L,
    distances             = FALSE,
    node_counts           = FALSE,
    cpus                  = n_cpus() ) {

  validate_counts()
  validate_margin()
  validate_tree()
  validate_cpus()

  k         <- ord_arg(n_axes, '`n_axes` must be a whole number from 1 to 10.', whole = TRUE, lo = 1, hi = 10)
  max_iters <- ord_arg(max_iters, '`max_iters` must be a whole number from 1 to 1000000.', whole = TRUE, lo = 1, hi = 1e6)
  tol       <- ord_arg(tol,       '`tol` must be a number of at least 0.')

  min_sample_count      <- ord_arg(min_sample_count,      '`min_sample_count` must be a number of at least 0.')
  min_feature_count     <- ord_arg(min_feature_count,     '`min_feature_count` must be a number of at least 0.')
  min_feature_frequency <- ord_arg(min_feature_frequency, '`min_feature_frequency` must be a number from 0 to 100.', hi = 100)
  min_depth             <- ord_arg(min_depth, '`min_depth` must be a whole number of at least 0.', whole = TRUE)

  if (!isTRUE(distances) && !isFALSE(distances))
    stop('`distances` must be TRUE or FALSE.', call. = FALSE)
  if (!isTRUE(node_counts) && !isFALSE(node_counts))
    stop('`node_counts` must be TRUE or FALSE.', call. = FALSE)


  # Sample names (validate_tree() has put the features in the order of
  # the tree's tips). Unnamed samples are numbered.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  sam_names <- if (is.null(dn[[1]])) as.character(seq_len(dims[[1]])) else as.character(dn[[1]])

  root_len <- tree[['root.edge']]
  root_len <- if (is.null(root_len)) 0 else as.double(root_len)[1]


  res <- .Call(
    C_rpca_phylo, counts, margin, tree$edge, tree$edge.length,
    as.integer(tree$Nnode), root_len, k, max_iters, tol,
    min_sample_count, min_feature_count, min_feature_frequency, min_depth,
    distances, node_counts, FALSE, 1e4, cpus )


  # Node names: tip labels, and n<id> (0-based ids) for internal nodes.
  K     <- length(res$node_tip)
  nodes <- ifelse(
    test = res$node_tip > 0L,
    yes  = as.character(tree$tip.label)[pmax(res$node_tip, 1L)],
    no   = paste0('n', seq_len(K) - 1L) )

  samples <- sam_names[res$samples]

  nc <- NULL
  if (node_counts) nc <- rpca_node_counts(res, samples, nodes)

  rpca_finish(
    res, k, 'rpca_phylo', samples, nodes, 'tree nodes', distances,
    tree        = rpca_tree(res, nodes),
    node_counts = nc )
}



# The tree that rpca_phylo() used, from C_rpca_phylo()'s list `res`, with
# nodes named `nodes` (in column order).
rpca_tree <- function (res, nodes) {

  K    <- length(nodes)
  col  <- res$phy_col          # Column of each ape node number.
  ntip <- sum(res$node_tip > 0L)

  phy <- list(
    edge        = res$phy_edge,
    edge.length = res$node_len[col[res$phy_edge[, 2]]],
    Nnode       = K - ntip,
    tip.label   = nodes[col[seq_len(ntip)]],
    node.label  = nodes[col[ntip + seq_len(K - ntip)]] )

  root_len <- res$node_len[[K]]   # The root is the last column.
  if (root_len > 0) phy$root.edge <- root_len

  structure(phy, class = 'phylo', order = 'cladewise')
}



# rpca_phylo()'s node counts, samples x nodes: a dgCMatrix if Matrix is
# installed, else a base matrix.
rpca_node_counts <- function (res, samples, nodes) {

  i  <- rep.int(seq_along(samples), diff(res$nc_pos))
  j  <- res$nc_idx + 1L
  dn <- list(samples, nodes)

  if (requireNamespace('Matrix', quietly = TRUE))
    return (Matrix::sparseMatrix(
      i = i, j = j, x = res$nc_val,
      dims = c(length(samples), length(nodes)), dimnames = dn ))

  m <- matrix(0, length(samples), length(nodes), dimnames = dn)
  m[cbind(i, j)] <- res$nc_val
  m
}
