#test_that("rpca_phylo", {

  # Threads for cpus arguments: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  # Largest absolute difference relative to the largest value of `y`, after
  # matching each column's sign.
  rel <- function (x, y) {
    x <- unname(as.matrix(x)); y <- unname(as.matrix(y))
    attributes(x) <- list(dim = dim(x))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # rpca_phylo()'s C function, with its intermediate results.
  # The counts get a column for each of the tree's tips, in its order, as
  # validate_tree() gives them.
  run <- function (counts, tree, k = 3L, filters = c(0, 0, 0), min_depth = 0L,
                   nc = TRUE, cpus = 1L) {
    full <- matrix(0, nrow(counts), length(tree$tip.label),
                   dimnames = list(rownames(counts), tree$tip.label))
    full[, colnames(counts)] <- counts
    root <- tree$root.edge
    .Call(C_rpca_phylo, full, 1L, tree$edge, tree$edge.length,
          as.integer(tree$Nnode), if (is.null(root)) 0 else root,
          as.integer(k), 5L, 1e-5, filters[1], filters[2], filters[3],
          as.integer(min_depth), TRUE, nc, TRUE, 1e4, as.integer(cpus))
  }

  # A CSR list (pos, idx 0-based, val) as a dense matrix, `fill` elsewhere.
  dense <- function (pos, idx, val, n, K, fill = 0) {
    m <- matrix(fill, n, K)
    m[cbind(rep.int(seq_len(n), diff(pos)), idx + 1L)] <- val
    m
  }

  # The Newick root length, which read_tree() discards (ape keeps it).
  with_root <- function (newick) {
    tr <- read_tree(newick)
    rl <- regmatches(newick, regexpr(':[-0-9.eE+]+;$', newick))
    if (length(rl)) tr$root.edge <- as.numeric(sub(';', '', sub(':', '', rl)))
    tr
  }

  # Tables with planted groups, and a random tree.
  set.seed(1)
  group_counts <- function (sizes, p, depth = 2000, dropout = 0) {
    base <- rnorm(p, 0, 1.5)
    do.call(rbind, lapply(seq_along(sizes), function (g) {
      mu <- base + rnorm(p, 0, 1.5)
      t(replicate(sizes[g], {
        x <- rpois(p, depth * exp(rnorm(1, 0, 0.5)) * exp(mu + rnorm(p, 0, 0.5)) / sum(exp(mu)))
        x * (runif(p) >= dropout)
      }))
    }))
  }
  rand_newick <- function (tips) {
    build <- function (ts) {
      if (length(ts) == 1) return (paste0(ts, ':', round(runif(1, 0.05, 1), 3)))
      nc  <- if (length(ts) >= 3 && runif(1) < 0.3) 3L else 2L
      cut <- sort(sample(seq_len(length(ts) - 1), nc - 1))
      grp <- split(ts, findInterval(seq_along(ts), cut + 1))
      paste0('(', paste(vapply(grp, build, ''), collapse = ','), '):', round(runif(1, 0, 1), 3))
    }
    paste0(sub(':[0-9.]+$', '', build(tips)), ';')
  }

  m    <- group_counts(c(10, 10), 30, dropout = 0.4)
  dimnames(m) <- list(paste0('S', 1:20), paste0('t', 1:30))
  tree <- read_tree(rand_newick(paste0('t', sample(40))))



  # ====================================================================
  # gemelli 0.0.13 reference results (see data-raw/expected/rpca_phylo.r):
  # the node table (names, tree, lengths, counts, robust clr), then the
  # ordination, at 1e-10 (as rpca()'s; see test-rpca.r).
  # ====================================================================

  if (at_home()) {

    TOL      <- 1e-10
    expected <- readRDS('expected/rpca_phylo.rds')

    for (fx in expected) {

      name <- fx$name
      k    <- fx$k[1]
      tr   <- with_root(fx$newick)
      res  <- run(fx$counts, tr, k, fx$filters, fx$min_depth[1])
      n    <- length(res$samples)
      K    <- length(fx$fids)

      # Kept samples; the node table's columns, names, tree and lengths.
      expect_identical(rownames(fx$counts)[res$samples], fx$samples, info = name)
      expect_identical(length(res$node_tip), K, info = name)
      nodes <- ifelse(res$node_tip > 0, tr$tip.label[pmax(res$node_tip, 1L)],
                      paste0('n', seq_len(K) - 1L))
      expect_identical(nodes, fx$fids, info = name)
      expect_identical(res$node_par, as.integer(fx$parent), info = name)
      expect_equal(res$node_len, c(fx$lengths), tolerance = 1e-15, info = name)

      # Node counts (exact), and the robust clr values.
      nc <- dense(res$nc_pos, res$nc_idx, res$nc_val, n, K)
      expect_identical(nc, unname(fx$node_counts) + 0, info = name)
      rc  <- dense(res$rclr_pos, res$rclr_idx, res$rclr, n, K, fill = NA)
      ref <- unname(fx$rclr)
      ref[!is.finite(ref) | ref == 0] <- NA  # gemelli's missing values
      expect_identical(is.na(rc), is.na(ref), info = name)
      expect_true(max(abs(rc - ref), na.rm = TRUE) < TOL * max(abs(ref), na.rm = TRUE), info = name)
      expect_identical(which(rowSums(!is.na(rc)) == 0), which(fx$flat == 1), info = name)

      # The ordination. gemelli's sample loadings are u and its eigenvalues
      # s; the C code's scores are u * s and its eig s^2.
      expect_identical(res$k, as.integer(k), info = name)
      u <- sweep(res$scores, 2, sqrt(res$eig), '/')
      expect_true(rel(u, fx$u) < TOL, info = name)
      expect_true(rel(res$loadings, fx$v) < TOL, info = name)
      expect_equal(sqrt(res$eig), c(fx$eig), tolerance = TOL, info = name)
      expect_equal(res$eig / sum(res$eig), c(fx$prop), tolerance = TOL, info = name)
      expect_equal(res$dist, c(fx$dist), tolerance = TOL, info = name)

      # Through rpca_phylo() itself.
      ord <- rpca_phylo(
        fx$counts, tr, n_axes = k, min_sample_count = fx$filters[1],
        min_feature_count = fx$filters[2], min_feature_frequency = fx$filters[3],
        min_depth = fx$min_depth[1], distances = TRUE, node_counts = TRUE, cpus = 1 )
      expect_identical(rownames(ord), fx$samples, info = name)
      expect_identical(rownames(attr(ord, 'loadings')), fx$fids, info = name)
      expect_true(rel(ord, fx$u %*% diag(c(fx$eig), k)) < TOL, info = name)
      expect_true(rel(attr(ord, 'loadings'), fx$v) < TOL, info = name)
      expect_equal(sqrt(attr(ord, 'eig') * (n - 1)), c(fx$eig), tolerance = TOL, info = name)
      expect_equal(attr(ord, 'prop'), c(fx$prop), tolerance = TOL, info = name)
      expect_equal(c(attr(ord, 'dm')), c(fx$dist), tolerance = TOL, info = name)
      expect_identical(unname(as.matrix(attr(ord, 'node_counts'))), unname(fx$node_counts) + 0, info = name)
      expect_identical(dimnames(attr(ord, 'node_counts')), list(fx$samples, fx$fids), info = name)

      # The tree it used: the nodes named as the loadings, each with
      # gemelli's parent and length.
      phy <- attr(ord, 'tree')
      expect_true(inherits(phy, 'phylo'), info = name)
      lab <- c(phy$tip.label, phy$node.label)
      expect_identical(sort(lab), sort(fx$fids), info = name)
      col <- match(lab, fx$fids)
      par <- integer(length(lab))
      par[phy$edge[, 2]] <- col[phy$edge[, 1]]
      expect_identical(par[order(col)], as.integer(fx$parent), info = name)
      bl <- numeric(length(lab))
      bl[phy$edge[, 2]] <- phy$edge.length
      root <- length(phy$tip.label) + 1L
      bl[root] <- if (is.null(phy$root.edge)) 0 else phy$root.edge
      expect_equal(bl[order(col)], c(fx$lengths), tolerance = 1e-15, info = name)
      expect_identical(phy$node.label[1], fx$fids[K], info = name)

      # Through count_ord(), and with threads and sparse input.
      expect_identical(
        count_ord(fx$counts, 'rpca_phylo', tree = tr, n_axes = k,
                  min_sample_count = fx$filters[1], min_feature_count = fx$filters[2],
                  min_feature_frequency = fx$filters[3], min_depth = fx$min_depth[1],
                  distances = TRUE, node_counts = TRUE, cpus = 1),
        ord, info = name )
      for (cpus in 2:max_cpus)
        expect_identical(
          rpca_phylo(fx$counts, tr, n_axes = k, min_sample_count = fx$filters[1],
                     min_feature_count = fx$filters[2], min_feature_frequency = fx$filters[3],
                     min_depth = fx$min_depth[1], distances = TRUE, node_counts = TRUE,
                     cpus = cpus),
          ord, info = name )
      if (requireNamespace('Matrix', quietly = TRUE))
        expect_identical(
          rpca_phylo(Matrix::Matrix(fx$counts, sparse = TRUE), tr, n_axes = k,
                     min_sample_count = fx$filters[1], min_feature_count = fx$filters[2],
                     min_feature_frequency = fx$filters[3], min_depth = fx$min_depth[1],
                     distances = TRUE, node_counts = TRUE, cpus = max_cpus),
          ord, info = name )
    }


    # The cases cover what they are meant to.
    fx <- expected$single_root   # The root merged with a chain: 0.4 + 0.3 + 0.2 + 0.6.
    expect_equal(fx$lengths[length(fx$lengths)], 1.5)
    fx <- expected$star          # Tips, and a root with no observed values.
    expect_identical(fx$fids, c(paste0('t', 1:8), 'n8'))
    expect_true(all(is.nan(fx$rclr[, 9])))
    expect_identical(sum(expected$flat$flat), 2)
    expect_true(length(expected$random_depth$fids) < length(expected$random$fids))
    expect_true(any(is.na(read_tree(expected$random$newick)$edge.length)) ||
                  any(read_tree(expected$random$newick)$edge.length == 0))


    # Without a single-child root, the root length doesn't count.
    fx <- expected$random
    o1 <- rpca_phylo(fx$counts, with_root(fx$newick), cpus = 1)
    o2 <- rpca_phylo(fx$counts, read_tree(fx$newick), cpus = 1)
    expect_identical(unclass(o1)[], unclass(o2)[])
    fx <- expected$single_root
    o1 <- rpca_phylo(fx$counts, with_root(fx$newick), cpus = 1)
    o2 <- rpca_phylo(fx$counts, read_tree(fx$newick), cpus = 1)
    expect_true(max(abs(unclass(o1) - unclass(o2))) > 1e-3)
    expect_equal(attr(o1, 'tree')$root.edge, 1.5)
    expect_equal(attr(o2, 'tree')$root.edge, 0.9)
  }



  # ====================================================================
  # Result format.
  # ====================================================================

  ord <- rpca_phylo(m, tree, n_axes = 2, distances = TRUE, node_counts = TRUE, cpus = 1)

  expect_true(inherits(ord, 'ecodive_ord') && is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings',
      'cost', 'n_iters', 'dm', 'tree', 'node_counts', 'class') )
  expect_identical(dimnames(ord), list(rownames(m), c('RPCA1', 'RPCA2')))
  expect_identical(attr(ord, 'method'), 'rpca_phylo')
  expect_identical(attr(attr(ord, 'dm'), 'method'), 'rpca_phylo')
  expect_identical(attr(ord, 'n_iters'), 5L)
  expect_equal(attr(ord, 'prop'), attr(ord, 'eig') / attr(ord, 'trace'), tolerance = 1e-14)

  # Tips and internal nodes; the root last, with no length (and so no
  # observed values) here.
  load <- attr(ord, 'loadings')
  phy  <- attr(ord, 'tree')
  K    <- nrow(load)
  expect_identical(sort(rownames(load)[startsWith(rownames(load), 't')]), sort(colnames(m)))
  expect_identical(rownames(load)[K], paste0('n', K - 1L))
  expect_identical(sort(c(phy$tip.label, phy$node.label)), sort(rownames(load)))
  expect_identical(phy$Nnode, K - ncol(m))
  expect_null(phy$root.edge)
  expect_identical(dim(phy$edge), c(K - 1L, 2L))

  # Scores times loadings is the double-centered completed matrix.
  X <- unclass(ord) %*% t(load)
  expect_true(max(abs(rowSums(X)), abs(colSums(X))) < 1e-10 * max(abs(X)))
  expect_equal(unname(apply(ord, 2, var)), attr(ord, 'eig'), tolerance = 1e-12)

  # Node counts: each internal node's is the sum of its children's.
  nc <- as.matrix(attr(ord, 'node_counts'))
  expect_identical(dim(nc), c(nrow(m), K))
  expect_identical(nc[, colnames(m)], m + 0)
  ch <- tapply(phy$edge[, 2], phy$edge[, 1], c)
  lab <- c(phy$tip.label, phy$node.label)
  for (p in names(ch))
    expect_equal(nc[, lab[as.integer(p)]], rowSums(nc[, lab[ch[[p]]], drop = FALSE]))
  expect_true(all(nc[, K] == rowSums(m)))

  # Without node_counts and distances.
  o2 <- rpca_phylo(m, tree, n_axes = 2, cpus = 1)
  expect_identical(
    names(attributes(o2)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings',
      'cost', 'n_iters', 'tree', 'class') )
  expect_identical(c(o2), c(ord))
  expect_identical(attr(o2, 'loadings'), attr(ord, 'loadings'))

  # The tree from the counts' attribute.
  expect_identical(rpca_phylo(structure(m, tree = tree), n_axes = 2, cpus = 1), o2)



  # ====================================================================
  # Identical for any cpus, and for every input class.
  # ====================================================================

  big  <- group_counts(c(30, 30, 30), 120, depth = 5000, dropout = 0.5)
  dimnames(big) <- list(paste0('S', 1:90), paste0('t', 1:120))
  btr  <- read_tree(rand_newick(paste0('t', sample(150))))
  ref  <- rpca_phylo(big, btr, min_depth = 2, distances = TRUE, node_counts = TRUE, cpus = 1)
  for (cpus in 2:max_cpus)
    expect_identical(rpca_phylo(big, btr, min_depth = 2, distances = TRUE, node_counts = TRUE, cpus = cpus), ref)

  # validate_tree() pads tips missing from a table with samples as columns
  # (it gave the padding the wrong dimnames before rpca_phylo()).
  expect_identical(faith(t(m), tree, margin = 2, cpus = 1), faith(m, tree, cpus = 1))
  expect_identical(rpca_phylo(t(big), btr, margin = 2, min_depth = 2, distances = TRUE, node_counts = TRUE, cpus = 1), ref)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(
      rpca_phylo(Matrix::Matrix(big, sparse = TRUE), btr, min_depth = 2,
                 distances = TRUE, node_counts = TRUE, cpus = max_cpus),
      ref )
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(
      rpca_phylo(slam::as.simple_triplet_matrix(big), btr, min_depth = 2,
                 distances = TRUE, node_counts = TRUE, cpus = max_cpus),
      ref )

  # Compositional scale invariance.
  o2 <- rpca_phylo(big / rowSums(big), btr, min_depth = 2, cpus = 1)
  expect_true(rel(o2, ref) < 1e-9)
  expect_true(rel(attr(o2, 'loadings'), attr(ref, 'loadings')) < 1e-9)

  # Tips of the tree that aren't in the table, and the table's columns in
  # another order, change nothing.
  expect_identical(rpca_phylo(big[, 120:1], btr, min_depth = 2, distances = TRUE, node_counts = TRUE, cpus = 1), ref)



  # ====================================================================
  # A star tree: the tips, and the root (all missing) named last.
  # ====================================================================

  star <- read_tree(paste0('(', paste0('t', 1:30, ':1', collapse = ','), ');'))
  os   <- rpca_phylo(m, star, n_axes = 2, cpus = 1)
  expect_identical(rownames(attr(os, 'loadings')), c(colnames(m), 'n30'))
  expect_true(all(attr(os, 'loadings')[31, ] != 0))



  # ====================================================================
  # Samples whose node values are all equal: kept, as rows with no
  # observed values, as gemelli keeps them.
  # ====================================================================

  m2 <- m
  m2[1, ]  <- 0
  m2[1, 5] <- 7
  unit <- read_tree(paste0('((', paste0('t', 1:15, ':1', collapse = ','), '):1,(',
                           paste0('t', 16:30, ':1', collapse = ','), '):1);'))
  o <- rpca_phylo(m2, unit, cpus = 1)
  expect_identical(rownames(o), rownames(m2))
  r <- run(m2, unit)
  expect_identical(r$rclr_pos[2], 0L)



  # ====================================================================
  # Errors.
  # ====================================================================

  expect_error(rpca_phylo(m), 'phylo')
  expect_error(rpca_phylo(m, 'tree'), 'phylo')
  expect_error(rpca_phylo(`colnames<-`(m, NULL), tree), 'phylo')
  expect_error(rpca_phylo(`colnames<-`(m, paste0('x', 1:30)), tree), 'phylo')
  expect_error(rpca_phylo(m, tree, n_axes = 0),   '`n_axes`')
  expect_error(rpca_phylo(m, tree, n_axes = 11),  '`n_axes`')
  expect_error(rpca_phylo(m, tree, max_iters = 0), '`max_iters`')
  expect_error(rpca_phylo(m, tree, tol = -1),     '`tol`')
  expect_error(rpca_phylo(m, tree, min_depth = -1),  '`min_depth`')
  expect_error(rpca_phylo(m, tree, min_depth = 1.5), '`min_depth`')
  expect_error(rpca_phylo(m, tree, min_depth = 30),  '`min_depth` \\(30\\) must be less')
  expect_error(rpca_phylo(m, tree, min_feature_frequency = 101), '`min_feature_frequency`')
  expect_error(rpca_phylo(m, tree, distances = NA),   '`distances`')
  expect_error(rpca_phylo(m, tree, node_counts = 1),  '`node_counts`')
  expect_error(rpca_phylo(m, tree, margin = 3),  '`margin`')
  expect_error(rpca_phylo(m, tree, cpus = 0),    '`cpus`')
  expect_error(rpca_phylo(-m, tree),             'negative')
  expect_error(rpca_phylo(m * NA, tree),         'NA')
  expect_error(rpca_phylo(m[1, , drop = FALSE], tree), 'at least 2 samples')
  expect_error(rpca_phylo(m, tree, min_sample_count = 1e9), '0 samples')

  m1 <- m * 0
  m1[, 3] <- 1:20
  expect_error(rpca_phylo(m1, tree), 'only 1 tip')

  # All of each sample's values equal: nothing to ordinate.
  m3 <- m * 0
  m3[, 3] <- 1:20
  m3[, 4] <- 1:20
  t3 <- read_tree(paste0('((t3:1,t4:1):0,(', paste0('t', c(1:2, 5:30), ':1', collapse = ','), '):1);'))
  expect_error(rpca_phylo(m3, t3), "every node's count times branch length")
  # (Here the root, merged with the clade of t3 and t4, has length 0.) With
  # a length, the values differ, but are the same in every sample.
  t3$edge.length[t3$edge[, 2] == length(t3$tip.label) + 2L] <- 1
  expect_error(suppressWarnings(rpca_phylo(m3, t3)), 'all samples are the same')

  bad <- tree
  bad$edge.length[1] <- -1
  expect_error(rpca_phylo(m, bad), 'non-negative')
  bad <- tree
  bad$edge[1, 1] <- 1L
  expect_error(rpca_phylo(m, bad), 'phylo')
  bad <- tree
  bad$edge.length <- bad$edge.length[-1]
  expect_error(rpca_phylo(m, bad), 'edge lengths')

  # NA lengths count as 0.
  na <- tree
  na$edge.length[na$edge[, 2] == 1L] <- NA
  z  <- tree
  z$edge.length[z$edge[, 2] == 1L] <- 0
  expect_identical(rpca_phylo(m, na, cpus = 1), rpca_phylo(m, z, cpus = 1))

  # n_axes reduced with a warning.
  expect_warning(rpca_phylo(m[1:3, ], tree, n_axes = 5, cpus = 1), 'reduced to 2')

#})
