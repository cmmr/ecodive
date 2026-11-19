#test_that("philr", {

  # philr() and philr_distance(): philr 1.28.0's outputs for every
  # weighting (data-raw/expected/philr.r), scikit-bio 0.7.4's tree_basis()
  # (data-raw/expected/philr.py), identities with ilr() and aitchison(), the
  # tree handling (pruning, polytomies), input types, threads, and argument
  # checks.

  max_cpus <- if (at_home()) 4L else 2L

  ref <- readRDS(file.path('expected', 'philr.rds'))

  # Relative difference, scaled by the largest value.
  rel   <- function (a, b) max(abs(a - b)) / max(1, max(abs(b)))
  plain <- function (m) { attributes(m) <- list(dim = dim(m)); m }
  quiet <- function (expr) suppressWarnings(expr)


  # -----------------------------------------------------------------------
  # philr 1.28.0, every weighting
  # -----------------------------------------------------------------------

  # philr's pseudocount 0 is ecodive's NULL (none, for a table without zeros).
  for (nm in names(ref))
    if (ref[[nm]]$pseudocount == 0) ref[[nm]]['pseudocount'] <- list(NULL)

  for (nm in setdiff(names(ref), c('skbio_small', 'skbio_random'))) {
    cs <- ref[[nm]]
    for (r in cs$results) {
      info <- paste(nm, r$part, r$ilr)
      run <- function () philr(
        cs$counts, cs$tree, pseudocount = cs$pseudocount,
        part_weights = r$part, ilr_weights = r$ilr, cpus = 1L )
      if (nm == 'polytomy') expect_warning(x <- run(), 'with more than two children')
      else                  x <- run()
      expect_true(rel(plain(x), plain(r$coords)) < 1e-10, info = info)
      expect_equal(rownames(x), rownames(cs$counts), info = info)
      expect_equal(attr(x, 'part_weights')[names(r$p)], r$p, tolerance = 1e-12, info = info)
      expect_equal(unname(attr(x, 'ilr_weights')), unname(r$ilr_weights), tolerance = 1e-12, info = info)
      if (!is.null(colnames(r$coords)))
        expect_identical(colnames(x), colnames(r$coords), info = info)
    }

    if (!is.null(cs$dist)) {
      d <- quiet(philr_distance(cs$counts, cs$tree, pseudocount = cs$pseudocount, cpus = 1L))
      expect_true(rel(as.vector(d), cs$dist) < 1e-10, info = nm)
      expect_identical(labels(d), rownames(cs$counts))
      expect_identical(quiet(beta_div(cs$counts, 'philr', tree = cs$tree, pseudocount = cs$pseudocount, cpus = 1L)), d)
    }
  }

  # Columns are named by the node labels, else n1, n2, ...
  cs <- ref$labelled
  x  <- philr(cs$counts, cs$tree, pseudocount = 0.5)
  expect_identical(colnames(x), cs$tree$node.label)
  expect_identical(colnames(philr(ref$unlabelled$counts, ref$unlabelled$tree)), paste0('n', 1:14))

  # No zeros, no pseudocount: none is used, without a warning.
  expect_silent(x <- philr(ref$unlabelled$counts, ref$unlabelled$tree))
  expect_identical(attr(x, 'pseudocount'), 0)


  # -----------------------------------------------------------------------
  # scikit-bio's tree_basis(): the same balances, with the opposite sign
  # (its numerator is the second child) and rows in level order. The trees
  # are read with read_tree(), whose node numbers must be ape's.
  # -----------------------------------------------------------------------

  for (nm in c('skbio_small', 'skbio_random')) {
    cs   <- ref[[nm]]
    tree <- read_tree(cs$newick)
    pc   <- cs$pseudocount

    for (r in cs$results) {
      x <- philr(cs$counts, tree, pseudocount = pc, part_weights = r$part, ilr_weights = r$ilr)
      expect_true(rel(plain(x), plain(r$coords)) < 1e-10, info = paste(nm, r$part))
      expect_identical(colnames(x), colnames(r$coords))
    }

    x   <- philr(cs$counts, tree, pseudocount = pc)
    clr <- normalize(cs$counts[, cs$tips] + pc, 'clr')
    sk  <- clr %*% t(cs$basis)
    colnames(sk) <- cs$nodes
    expect_true(rel(plain(x), -plain(sk[, colnames(x)])) < 1e-12, info = nm)

    d <- quiet(philr_distance(cs$counts, tree, pseudocount = pc))
    expect_true(rel(as.vector(d), cs$dist) < 1e-10, info = nm)
  }


  # -----------------------------------------------------------------------
  # Identities
  # -----------------------------------------------------------------------

  cs <- ref$labelled
  x  <- philr(cs$counts, cs$tree, pseudocount = 0.5)

  # Uniform weights: Aitchison distances, and ilr() with the tree's SBP.
  expect_equal(as.vector(euclidean(x)), as.vector(aitchison(cs$counts, pseudocount = 0.5)), tolerance = 1e-12)

  tr  <- attr(x, 'tree')
  D   <- length(tr$tip.label)
  kids <- split(tr$edge[, 2], tr$edge[, 1])
  tips_below <- function (v) if (v <= D) v else unlist(lapply(kids[[as.character(v)]], tips_below))
  sbp <- t(vapply(D + seq_len(tr$Nnode), function (v) {
    k <- kids[[as.character(v)]]
    s <- numeric(D)
    s[tips_below(k[1])] <-  1
    s[tips_below(k[2])] <- -1
    s }, numeric(D)))
  dimnames(sbp) <- list(tr$node.label, tr$tip.label)
  y <- ilr(cs$counts, basis = sbp, pseudocount = 0.5)
  expect_true(rel(plain(x), plain(y)) < 1e-12)
  expect_identical(colnames(x), colnames(y))
  expect_identical(tr$tip.label, colnames(cs$counts))

  # The tree attribute gives the same coordinates.
  for (w in list(c('uniform', 'uniform'), c('enorm.x.gm.counts', 'mean.descendants'))) {
    x1 <- philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = w[1], ilr_weights = w[2])
    x2 <- philr(cs$counts, attr(x1, 'tree'), pseudocount = 0.5, part_weights = w[1], ilr_weights = w[2])
    expect_identical(plain(x1), plain(x2))
    expect_identical(dimnames(x1), dimnames(x2))
  }

  # philr_distance() is euclidean() of the paper's coordinates.
  x <- philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt')
  d <- philr_distance(cs$counts, cs$tree, pseudocount = 0.5)
  expect_identical(as.vector(d), as.vector(euclidean(x)))
  expect_identical(attr(d, 'dist')$metric, 'philr_distance')
  expect_identical(attr(d, 'dist')$name, 'PhILR Distance')
  expect_identical(attr(d, 'dist')$pseudocount, 0.5)
  expect_true('philr_distance' %in% list_metrics('beta', 'id', phylo = TRUE))
  expect_identical(match_metric('philr')$id, 'philr_distance')
  expect_identical(match_metric('Phylogenetic ILR')$id, 'philr_distance')

  # Numeric weights: the computed ones, given back, by name or position.
  x  <- philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = 'anorm.x.gm.counts', ilr_weights = 'blw')
  pw <- attr(x, 'part_weights')
  iw <- attr(x, 'ilr_weights')
  expect_identical(plain(philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = rev(pw), ilr_weights = rev(iw))), plain(x))
  expect_identical(plain(philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = unname(pw), ilr_weights = unname(iw))), plain(x))
  expect_identical(plain(philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = 'anorm_x_gm_counts', ilr_weights = 'blw')), plain(x))
  expect_identical(plain(philr(cs$counts, cs$tree, pseudocount = 0.5, part_weights = c(pw, extra = 2), ilr_weights = iw)), plain(x))

  # Uniform weights are scale invariant, per sample.
  y <- ref$unlabelled$counts
  t15 <- ref$unlabelled$tree
  expect_true(rel(plain(philr(y * c(1, 10, 0.1, 3, 7, 2), t15)), plain(philr(y, t15))) < 1e-12)
  # enorm and anorm are invariant to scaling the whole table.
  for (pw in c('enorm', 'anorm'))
    expect_true(rel(plain(philr(y * 13, t15, part_weights = pw)), plain(philr(y, t15, part_weights = pw))) < 1e-12)

  # replace_zeros()'s output goes straight in.
  rz <- replace_zeros(ref$labelled$counts, 'mult')
  expect_silent(x <- philr(rz, ref$labelled$tree))
  expect_identical(attr(x, 'pseudocount'), 0)
  expect_true(all(is.finite(x)))


  # -----------------------------------------------------------------------
  # The tree: pruning, polytomies, numbering
  # -----------------------------------------------------------------------

  cs <- ref$polytomy

  # Only pruning (keep a binary part of the tree): no warning, the same
  # result as philr() on the pruned tree.
  rt <- cs$ref_tree
  sub <- rt$tip.label[1:12]
  expect_silent(x1 <- philr(cs$counts[, sub], rt, pseudocount = 1, ilr_weights = 'mean.descendants'))
  x2 <- philr(cs$counts[, sub], attr(x1, 'tree'), pseudocount = 1, ilr_weights = 'mean.descendants')
  expect_identical(plain(x1), plain(x2))
  expect_identical(ncol(x1), 11L)
  expect_identical(sort(attr(x1, 'tree')$tip.label), sort(sub))

  # Polytomies, with labels: the added nodes are named after theirs.
  tr <- list(
    edge = rbind(c(7L, 8L), c(8L, 1L), c(8L, 2L), c(8L, 3L), c(8L, 4L), c(7L, 9L), c(9L, 5L), c(9L, 6L)),
    edge.length = c(1, 1, 2, 3, 4, 2, 1, 1),
    tip.label = letters[1:6], Nnode = 3L, node.label = c('R', 'X', 'Y') )
  class(tr) <- 'phylo'
  m <- matrix(c(1:6, 6:1, 3, 1, 4, 1, 5, 9) + 0, 3, byrow = TRUE, dimnames = list(NULL, letters[1:6]))
  expect_warning(x <- philr(m, tr, ilr_weights = 'blw'), '1 node with more than two children')
  expect_identical(colnames(x), c('R', 'X', 'X.1', 'X.2', 'Y'))
  expect_equal(unname(attr(x, 'ilr_weights')), c(3, 1 + 0, 2 + 0, 3 + 4, 2))
  expect_identical(attr(x, 'tree')$edge, rbind(
    c(7L, 8L), c(8L, 1L), c(8L, 9L), c(9L, 2L), c(9L, 10L), c(10L, 3L), c(10L, 4L),
    c(7L, 11L), c(11L, 5L), c(11L, 6L) ))
  tr0 <- tr; tr0$node.label <- NULL
  expect_identical(colnames(quiet(philr(m, tr0))), paste0('n', 1:5))
  tr0$node.label <- c('R', '', 'Y')   # Incomplete labels: numbered.
  expect_identical(colnames(quiet(philr(m, tr0))), paste0('n', 1:5))

  # Single-child nodes collapse, adding their lengths; a single-child root is
  # dropped.
  tr <- list(
    edge = rbind(c(4L, 5L), c(5L, 6L), c(6L, 1L), c(6L, 7L), c(7L, 2L), c(7L, 3L)),
    edge.length = c(5, 1, 1, 2, 1, 1), tip.label = c('a', 'b', 'c'), Nnode = 4L )
  class(tr) <- 'phylo'
  expect_silent(x <- philr(m[, 1:3], tr, ilr_weights = 'blw'))
  expect_identical(attr(x, 'tree')$Nnode, 2L)
  expect_equal(unname(attr(x, 'ilr_weights')), c(1 + 2, 2))

  # Two features: one coordinate, log(a / b) / sqrt(2).
  tr2 <- structure(list(edge = rbind(c(3L, 1L), c(3L, 2L)), edge.length = c(1, 1), tip.label = c('a', 'b'), Nnode = 1L), class = 'phylo')
  x <- philr(m[, 1:2], tr2)
  expect_equal(as.vector(x), log(m[, 1] / m[, 2]) / sqrt(2))

  # Branch lengths are only needed for branch-length weights.
  tr3 <- tr2; tr3$edge.length <- NULL
  expect_identical(plain(philr(m[, 1:2], tr3)), plain(x))
  expect_error(philr(m[, 1:2], tr3, ilr_weights = 'blw'), 'needs `tree` to have branch lengths')
  tr3$edge.length <- c(0, 0)
  expect_error(philr(m[, 1:2], tr3, ilr_weights = 'blw.sqrt'), 'at least one positive')


  # -----------------------------------------------------------------------
  # Input types, orientation, pairs, threads
  # -----------------------------------------------------------------------

  cs <- ref$labelled
  pc <- 0.5
  x0 <- philr(cs$counts, cs$tree, pseudocount = pc, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt', cpus = 1L)
  d0 <- philr_distance(cs$counts, cs$tree, pseudocount = pc, cpus = 1L)

  xt <- philr(t(cs$counts), cs$tree, margin = 2L, pseudocount = pc, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt')
  expect_identical(t(plain(xt)), plain(x0))
  expect_identical(dimnames(xt), rev(dimnames(x0)))
  expect_identical(as.vector(philr_distance(t(cs$counts), cs$tree, margin = 2L, pseudocount = pc)), as.vector(d0))

  inputs <- list(slam = slam::as.simple_triplet_matrix(cs$counts))
  if (requireNamespace('Matrix', quietly = TRUE))
    inputs$dgC <- methods::as(methods::as(cs$counts, 'CsparseMatrix'), 'generalMatrix')
  for (nm in names(inputs)) {
    xs <- philr(inputs[[nm]], cs$tree, pseudocount = pc, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt')
    expect_true(rel(plain(xs), plain(x0)) < 1e-14, info = nm)
    expect_identical(dimnames(xs), dimnames(x0), info = nm)
    expect_true(rel(as.vector(philr_distance(inputs[[nm]], cs$tree, pseudocount = pc)), as.vector(d0)) < 1e-14, info = nm)
  }

  # A tree attached to the counts.
  y <- cs$counts
  attr(y, 'tree') <- cs$tree
  expect_identical(plain(philr(y, pseudocount = pc)), plain(philr(cs$counts, cs$tree, pseudocount = pc)))

  for (cpus in seq_len(max_cpus)[-1]) {
    expect_identical(philr(cs$counts, cs$tree, pseudocount = pc, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt', cpus = cpus), x0)
    expect_identical(philr_distance(cs$counts, cs$tree, pseudocount = pc, cpus = cpus), d0)
  }

  d <- philr_distance(cs$counts, cs$tree, pseudocount = pc, pairs = c(1, 5, 9))
  expect_identical(as.vector(d)[c(1, 5, 9)], as.vector(d0)[c(1, 5, 9)])
  expect_true(all(is.na(as.vector(d)[-c(1, 5, 9)])))
  expect_error(philr_distance(cs$counts, cs$tree, pseudocount = pc, pairs = 1:3, clust = 'average'), 'pairs')
  dc <- philr_distance(cs$counts, cs$tree, pseudocount = pc, clust = 'average')
  expect_identical(attr(dc, 'dist')$clust, 'average')

  # Zeros: the shared rule.
  expect_warning(x <- philr(cs$counts, cs$tree), 'Zeros detected')
  expect_identical(attr(x, 'pseudocount'), min(cs$counts[cs$counts > 0]) / 2)

  # Two samples.
  expect_identical(dim(philr(cs$counts[1:2, ], cs$tree, pseudocount = 1)), c(2L, 9L))
  expect_identical(attr(philr_distance(cs$counts[1:2, ], cs$tree, pseudocount = 1), 'Size'), 2L)


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  m <- cs$counts
  tr <- cs$tree
  expect_error(philr(m), '`tree` is required')
  expect_error(philr(m, tree = list()), 'not a valid phylo object')
  expect_error(philr(unname(m), tr), 'must have feature names')
  m2 <- m; colnames(m2)[1] <- 'nope'
  expect_error(philr(m2, tr), '1 feature of `counts` is not in `tree`: "nope"')
  m2 <- m; colnames(m2)[2] <- colnames(m2)[1]
  expect_error(philr(m2, tr), 'duplicated feature names')
  tr2 <- tr; tr2$tip.label[2] <- tr2$tip.label[1]
  expect_error(philr(m, tr2, pseudocount = 1), 'duplicated tip labels')
  tr2 <- tr; tr2$edge.length[1] <- -1
  expect_error(philr(m, tr2, pseudocount = 1), 'non-negative')
  tr2 <- tr; tr2$edge[1, 1] <- tr2$edge[1, 2]
  expect_error(philr(m, tr2, pseudocount = 1), 'malformed edge matrix')
  tr2 <- tr; tr2$edge[2, 2] <- tr2$edge[1, 2]
  expect_error(philr(m, tr2, pseudocount = 1), 'malformed edge matrix')

  expect_error(philr(m, tr, pseudocount = 1, part_weights = 'gm'), '`part_weights` must be one of')
  expect_error(philr(m, tr, pseudocount = 1, ilr_weights = 'sum.children'), '`ilr_weights` must be one of')
  expect_error(philr(m, tr, pseudocount = 1, part_weights = NA), '`part_weights` must be one of')
  expect_error(philr(m, tr, pseudocount = 1, part_weights = c(a = 1)), 'no value for 10 features')
  expect_error(philr(m, tr, pseudocount = 1, part_weights = 1:3), 'one value per feature \\(10\\)')
  expect_error(philr(m, tr, pseudocount = 1, part_weights = c(-1, rep(1, 9))), 'positive numbers')
  expect_error(philr(m, tr, pseudocount = 1, ilr_weights = c(node1 = 1)), 'no value for 8 coordinates')
  expect_error(philr(m, tr, pseudocount = 1, ilr_weights = c(NA, rep(1, 8))), 'non-negative numbers')
  expect_error(philr(m, tr, pseudocount = -1), '`pseudocount` must be')

  # anorm of a feature with the same value in every sample is 0.
  m2 <- m; m2[, 3] <- 7
  expect_error(philr(m2, tr, pseudocount = 1, part_weights = 'anorm'), 'a weight of 0 to features')
  expect_error(philr(m[1, , drop = FALSE], tr, pseudocount = 1, part_weights = 'anorm.x.gm.counts'), 'a weight of 0')

  expect_error(philr_distance(m[1, , drop = FALSE], tr, pseudocount = 1), 'at least two samples')

#})
