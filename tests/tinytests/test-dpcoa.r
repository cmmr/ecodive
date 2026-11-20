#test_that("dpcoa", {

  # Threads for cpus arguments: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  # Largest absolute difference relative to the largest value of `y`, after
  # matching each column's sign.
  rel <- function (x, y) {
    x <- unname(as.matrix(x)); y <- unname(as.matrix(y))
    attributes(x) <- list(dim = dim(x))
    attributes(y) <- list(dim = dim(y))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }
  relv <- function (x, y) max(abs(unname(x) - unname(y))) / max(abs(y))

  # Patristic distances between the tips of a phylo list, in plain R: each
  # tip's ancestors with their distances from it, then, for each pair, the
  # nearest common ancestor.
  patristic <- function (tr) {
    ntip <- length(tr$tip.label)
    par  <- integer(ntip + tr$Nnode); len <- numeric(ntip + tr$Nnode)
    par[tr$edge[, 2]] <- tr$edge[, 1]; len[tr$edge[, 2]] <- tr$edge.length
    up <- lapply(seq_len(ntip), function (t) {
      nodes <- t; dist <- 0
      while (par[nodes[length(nodes)]] > 0) {
        v <- nodes[length(nodes)]
        nodes <- c(nodes, par[v]); dist <- c(dist, dist[length(dist)] + len[v])
      }
      list(nodes = nodes, dist = dist)
    })
    d <- matrix(0, ntip, ntip, dimnames = list(tr$tip.label, tr$tip.label))
    for (i in seq_len(ntip)) for (j in seq_len(ntip)) {
      a <- match(up[[i]]$nodes, up[[j]]$nodes)
      f <- which(!is.na(a))[1]
      d[i, j] <- up[[i]]$dist[f] + up[[j]]$dist[a[f]]
    }
    d
  }

  # DPCoA in plain R, from the taxa dissimilarities `delta` (a matrix in
  # the order of counts' columns): eigen() of the samples' weighted Gram
  # matrix. Returns li, dls, eig, rao, trace.
  dpcoa_r <- function (counts, delta, k) {
    X <- as.matrix(counts)
    w <- rowSums(X) / sum(X); P <- X / rowSums(X); pbar <- colSums(X) / sum(X)
    A <- delta^2
    Q <- t(sweep(P, 2, pbar))                     # taxa x samples
    J <- diag(ncol(X)) - outer(rep(1, ncol(X)), pbar)
    B <- J %*% (-A / 2) %*% t(J)                  # taxa Gram, centred by pbar
    G <- sqrt(w) * t(sqrt(w) * (t(Q) %*% (-A / 2) %*% Q))
    e <- eigen(G, symmetric = TRUE)
    U <- e$vectors[, seq_len(k), drop = FALSE]
    s <- sqrt(e$values[seq_len(k)])
    list(
      li    = sweep(U, 2, s, '*') / sqrt(w),
      dls   = B %*% t(P) %*% (sqrt(w) * U) %*% diag(1 / s, k),
      eig   = e$values[seq_len(k)],
      rao   = rowSums((P %*% A) * P) / 2,
      trace = sum(e$values[e$values > 0]) )
  }


  # ade4 and phyloseq: the stored results (data-raw/expected/dpcoa.r). ====

  if (at_home()) {

    cases <- readRDS('expected/dpcoa.rds')

    for (nm in names(cases)) {

      cs <- cases[[nm]]
      o  <- if (is.null(cs$tree)) dpcoa(cs$counts, taxa_dm = cs$taxa_dm, n_axes = 4, cpus = 1)
            else                  dpcoa(cs$counts, cs$tree, n_axes = 4, cpus = 1)
      k  <- ncol(o)
      u  <- match(cs$used, colnames(cs$counts))

      expect_identical(k, min(4L, nrow(cs$counts) - 1L), info = nm)
      expect_identical(colnames(o), paste0('DPCoA', seq_len(k)), info = nm)
      expect_identical(rownames(o), rownames(cs$counts), info = nm)
      expect_true(rel(o, cs$li[, 1:k]) < 1e-12, info = nm)
      expect_true(rel(attr(o, 'taxa')[u, ], cs$dls[, 1:k]) < 1e-12, info = nm)
      expect_true(relv(attr(o, 'eig'), cs$eig[1:k]) < 1e-12, info = nm)
      expect_true(relv(attr(o, 'rao'), cs$rao) < 1e-12, info = nm)
      expect_equal(attr(o, 'trace'), cs$between, tolerance = 1e-12, info = nm)
      expect_equal(attr(o, 'prop'), cs$eig[1:k] / sum(cs$eig), tolerance = 1e-12, info = nm)
      expect_identical(rownames(attr(o, 'taxa')), colnames(cs$counts), info = nm)
      expect_identical(names(attr(o, 'rao')), rownames(cs$counts), info = nm)

      if (!is.null(cs$phyloseq)) {
        ps <- cs$phyloseq
        expect_true(rel(o, ps$li[, 1:k]) < 1e-12, info = nm)
        expect_true(rel(attr(o, 'taxa'), ps$dls[, 1:k]) < 1e-12, info = nm)
        expect_true(relv(attr(o, 'eig'), ps$eig[1:k]) < 1e-12, info = nm)
        expect_true(relv(attr(o, 'rao'), ps$rao) < 1e-12, info = nm)
      }

      # The tree cases again through taxa_dm, from plain-R patristic
      # distances; and against plain R.
      if (!is.null(cs$tree)) {
        pd    <- patristic(cs$tree)[colnames(cs$counts), colnames(cs$counts)]
        od    <- dpcoa(cs$counts, taxa_dm = sqrt(pd), n_axes = 4, cpus = 1)
        expect_true(rel(od, o) < 1e-12, info = nm)
        expect_true(rel(attr(od, 'taxa'), attr(o, 'taxa')) < 1e-12, info = nm)
        expect_true(relv(attr(od, 'rao'), attr(o, 'rao')) < 1e-12, info = nm)
        expect_equal(attr(od, 'trace'), attr(o, 'trace'), tolerance = 1e-12, info = nm)

        r <- dpcoa_r(cs$counts, sqrt(pd), k)
        expect_true(rel(o, r$li) < 1e-10, info = nm)
        expect_true(rel(attr(o, 'taxa'), r$dls) < 1e-10, info = nm)
        expect_true(relv(attr(o, 'eig'), r$eig) < 1e-10, info = nm)
        expect_true(relv(attr(o, 'rao'), r$rao) < 1e-12, info = nm)
        expect_equal(attr(o, 'trace'), r$trace, tolerance = 1e-10, info = nm)
      }
    }

    # Taxa absent from every sample: kept, with coordinates from the
    # formula (ade4 refuses them; plain R gives them too).
    cs <- cases$rtree
    expect_true(length(cs$used) < ncol(cs$counts))
    o  <- dpcoa(cs$counts, cs$tree, n_axes = 3, cpus = 1)
    pd <- patristic(cs$tree)[colnames(cs$counts), colnames(cs$counts)]
    r  <- dpcoa_r(cs$counts, sqrt(pd), 3)
    expect_true(rel(attr(o, 'taxa'), r$dls) < 1e-10)
    expect_true(all(is.finite(attr(o, 'taxa'))))

    # The tree with extra tips gives the same as the tree cut to the taxa.
    expect_true(length(cs$tree$tip.label) > ncol(cs$counts))
  }


  # Data for the rest. ====

  set.seed(3)
  newick <- function (tips) {
    build <- function (ts) {
      if (length(ts) == 1) return (paste0(ts, ':', round(runif(1, 0.05, 1), 3)))
      nc  <- if (length(ts) >= 3 && runif(1) < 0.25) 3L else 2L
      cut <- sort(sample(seq_len(length(ts) - 1), nc - 1))
      grp <- split(ts, findInterval(seq_along(ts), cut + 1))
      paste0('(', paste(vapply(grp, build, ''), collapse = ','), '):', round(runif(1, 0, 0.5), 3))
    }
    paste0(sub(':[0-9.]+$', '', build(tips)), ';')
  }
  D   <- 400
  tr  <- read_tree(newick(paste0('T', 1:D)))
  big <- matrix(rpois(300 * D, 2) * rbinom(300 * D, 1, 0.1), 300, D,
                dimnames = list(paste0('S', 1:300), sample(paste0('T', 1:D))))
  big[, 1] <- big[, 1] + 1L

  m  <- ecodive::ex_counts
  et <- ecodive::ex_tree


  # Threads, sparse input, margin, and count_ord(). ====

  o1 <- dpcoa(big, tr, n_axes = 5, cpus = 1)
  for (cpus in 2:max_cpus)
    expect_identical(dpcoa(big, tr, n_axes = 5, cpus = cpus), o1, info = cpus)

  pd_big <- patristic(tr)[colnames(big), colnames(big)]
  if (at_home()) {
    d1 <- dpcoa(big, taxa_dm = sqrt(pd_big), n_axes = 5, cpus = 1)
    for (cpus in 2:max_cpus)
      expect_identical(dpcoa(big, taxa_dm = sqrt(pd_big), n_axes = 5, cpus = cpus), d1, info = cpus)
    expect_true(rel(d1, o1) < 1e-10)
    expect_true(rel(attr(d1, 'taxa'), attr(o1, 'taxa')) < 1e-10)
    expect_true(relv(attr(d1, 'rao'), attr(o1, 'rao')) < 1e-12)
  }

  expect_identical(no_display(dpcoa(t(big), tr, margin = 2, n_axes = 5, cpus = 1)), no_display(o1))

  sp <- slam::as.simple_triplet_matrix(big)
  expect_true(rel(dpcoa(sp, tr, n_axes = 5, cpus = 1), o1) < 1e-12)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    dg <- methods::as(methods::as(big, 'CsparseMatrix'), 'generalMatrix')
    os <- dpcoa(dg, tr, n_axes = 5, cpus = 1)
    expect_true(rel(os, o1) < 1e-12)
    expect_true(rel(attr(os, 'taxa'), attr(o1, 'taxa')) < 1e-12)
    expect_true(relv(attr(os, 'rao'), attr(o1, 'rao')) < 1e-12)
  }

  expect_identical(count_ord(m, 'dpcoa', tree = et, cpus = 1), dpcoa(m, et, cpus = 1))
  expect_identical(count_ord(m, 'DPCoA', taxa_dm = sqrt(patristic(et)), cpus = 1),
                   dpcoa(m, taxa_dm = sqrt(patristic(et)), cpus = 1))


  # The result. ====

  o <- dpcoa(m, et, n_axes = 2, cpus = 1)
  expect_identical(names(attributes(o)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'taxa', 'rao', 'trace', 'class', 'display'))
  expect_identical(attr(o, 'method'), 'dpcoa')
  expect_identical(attr(o, 'display'), list(title = 'Double Principal Coordinate Analysis'))
  expect_identical(dim(attr(o, 'taxa')), c(ncol(m), 2L))
  expect_identical(colnames(attr(o, 'taxa')), c('DPCoA1', 'DPCoA2'))
  expect_true(all(diff(attr(o, 'eig')) <= 0))
  expect_true(sum(attr(o, 'prop')) <= 1 + 1e-12)
  expect_true(any(grepl('DPCoA1', capture.output(print(o)))))

  # Each sample's weighted variance on an axis is its eigenvalue, the
  # samples are centred, and the taxa project to the samples: a sample's
  # coordinates are its taxa's, weighted by its relative abundances.
  w <- rowSums(m) / sum(m); P <- m / rowSums(m)
  expect_equal(unname(colSums(w * unclass(o)[, 1:2]^2)), attr(o, 'eig'), tolerance = 1e-12)
  expect_equal(unname(colSums(w * unclass(o)[, 1:2])), c(0, 0), tolerance = 1e-12)
  expect_equal(unname(P %*% attr(o, 'taxa')), unname(unclass(o)[, 1:2]), tolerance = 1e-12)

  # Rao from the tree: L a (1 - a) summed over branches equals half the
  # expected patristic distance.
  pd <- patristic(et)[colnames(m), colnames(m)]
  expect_equal(unname(attr(o, 'rao')), unname(rowSums((P %*% pd) * P) / 2), tolerance = 1e-12)

  # The total inertia: Rao's diversity of the pooled samples, less the
  # weighted mean of the samples'.
  pbar <- colSums(m) / sum(m)
  expect_equal(attr(o, 'trace'), sum(pbar %*% pd %*% pbar) / 2 - sum(w * attr(o, 'rao')), tolerance = 1e-12)

  # An embedded tree; an embedded tree is ignored with taxa_dm.
  m2 <- m; attr(m2, 'tree') <- et
  expect_identical(dpcoa(m2, cpus = 1), dpcoa(m, et, cpus = 1))
  expect_identical(dpcoa(m2, taxa_dm = sqrt(pd), cpus = 1), dpcoa(m, taxa_dm = sqrt(pd), cpus = 1))


  # Trees: polytomies, unrooted, single-child nodes and extra tips change
  # nothing. ====

  t1 <- read_tree('((A:1,B:2):0.5,(C:1,(D:0.5,E:1.5):0.25):1);')
  t2 <- read_tree('(A:1,B:2,(C:1,(D:0.5,E:1.5):0.25):1.5);')          # unrooted
  t3 <- read_tree('((A:1,B:2):0.5,((C:1,((D:0.5):0,E:1.5):0.25):1,F:3):0);')  # single child, extra tip
  t4 <- read_tree('((A:1,B:2):0.5,(C:1.25,D:0.75,E:1.75):1);')         # polytomy
  x  <- matrix(c(5, 0, 3, 1, 2,  0, 4, 1, 0, 6,  2, 2, 2, 2, 2,  9, 1, 0, 0, 0), 4, 5,
               byrow = TRUE, dimnames = list(paste0('S', 1:4), c('E', 'D', 'C', 'B', 'A')))
  r1 <- dpcoa(x, t1, cpus = 1)
  expect_true(rel(dpcoa(x, t2, cpus = 1), r1) < 1e-12)
  expect_true(rel(dpcoa(x, t3, cpus = 1), r1) < 1e-12)
  expect_true(relv(attr(dpcoa(x, t3, cpus = 1), 'rao'), attr(r1, 'rao')) < 1e-12)
  expect_silent(r4 <- dpcoa(x, t4, cpus = 1))
  r <- dpcoa_r(x, sqrt(patristic(t4)[colnames(x), colnames(x)]), 3)
  expect_true(rel(r4, r$li) < 1e-10)
  expect_true(rel(attr(r4, 'taxa'), r$dls) < 1e-10)


  # Small and degenerate tables. ====

  # A sample with one taxon has a Rao entropy of 0.
  expect_equal(unname(attr(r1, 'rao')[4] > 0), TRUE)
  x1 <- x; x1[4, ] <- c(7, 0, 0, 0, 0)
  expect_equal(unname(attr(dpcoa(x1, t1, cpus = 1), 'rao')[4]), 0)

  # Two samples: one axis. One sample: none, but rao.
  o2 <- dpcoa(x[1:2, ], t1, cpus = 1)
  expect_identical(ncol(o2), 1L)
  expect_equal(unname(diff(unclass(o2)[, 1])^2), sum(patristic(t1)[colnames(x), colnames(x)] *
    outer(x[1, ] / sum(x[1, ]) - x[2, ] / sum(x[2, ]), x[1, ] / sum(x[1, ]) - x[2, ] / sum(x[2, ]))) / -2,
    tolerance = 1e-12)
  o1s <- dpcoa(x[1, , drop = FALSE], t1, cpus = 1)
  expect_identical(dim(o1s), c(1L, 0L))
  expect_identical(dim(attr(o1s, 'taxa')), c(5L, 0L))
  expect_equal(unname(attr(o1s, 'rao')), unname(attr(r1, 'rao')[1]), tolerance = 1e-12)

  # Two taxa: one axis at most.
  x2 <- x[, c('A', 'B')]; x2[3, ] <- c(1, 4); x2[4, ] <- c(3, 3)
  expect_identical(ncol(dpcoa(x2, t1, cpus = 1)), 1L)
  expect_identical(ncol(dpcoa(x2, taxa_dm = dist(c(A = 0, B = 1)), cpus = 1)), 1L)

  # Identical samples: no variation, no axes.
  xs <- x[c(1, 1, 1), ]; rownames(xs) <- paste0('S', 1:3)
  expect_warning(os <- dpcoa(xs, t1, cpus = 1), 'only 0 of the first 2 eigenvalues are > 0')
  expect_identical(ncol(os), 0L)
  expect_identical(attr(os, 'trace'), 0)
  expect_warning(os <- dpcoa(xs, taxa_dm = d_t1 <- sqrt(patristic(t1)), cpus = 1), 'only 0 of the first 2')
  expect_identical(attr(os, 'trace'), 0)
  xs[2, ] <- xs[2, ] * 7          # the same proportions
  expect_warning(os <- dpcoa(xs, t1, cpus = 1), 'only 0 of the first 2')
  expect_identical(ncol(os), 0L)

  # Proportions as well as counts; scaling a sample changes its weight only.
  expect_true(rel(dpcoa(x / 10, t1, cpus = 1), r1) < 1e-12)
  xw <- x; xw[1, ] <- xw[1, ] * 3
  expect_equal(attr(dpcoa(xw, t1, cpus = 1), 'rao'), attr(r1, 'rao'), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(attr(dpcoa(xw, t1, cpus = 1), 'eig'), attr(r1, 'eig'))))


  # taxa_dm: labels matched, extra taxa dropped, or by position. ====

  d1 <- sqrt(patristic(t1))
  rd <- dpcoa(x, taxa_dm = d1, cpus = 1)
  expect_true(rel(rd, r1) < 1e-12)
  expect_true(rel(dpcoa(x, taxa_dm = stats::as.dist(d1), cpus = 1), r1) < 1e-12)
  dx <- sqrt(patristic(t3))       # with F, in another order
  expect_true(rel(dpcoa(x, taxa_dm = dx[6:1, 6:1], cpus = 1), r1) < 1e-12)
  expect_true(rel(dpcoa(x, taxa_dm = unname(d1[colnames(x), colnames(x)]), cpus = 1), r1) < 1e-12)
  expect_identical(rownames(attr(rd, 'taxa')), colnames(x))


  # Errors. ====

  expect_error(dpcoa(x, t1, taxa_dm = d1), 'either `tree` or `taxa_dm`, not both')
  expect_error(dpcoa(x), '`tree` or `taxa_dm` is required')
  expect_error(dpcoa(x, read_tree('((A:1,B:2):0.5,(C:1,D:0.5):1);')), 'feature of `counts` is not in `tree`: "E"')
  t0 <- t1; t0$edge.length <- NULL
  expect_error(dpcoa(x, t0), 'must have branch lengths')
  t0 <- t1; t0$edge.length[2] <- NA
  expect_error(dpcoa(x, t0), 'must have branch lengths')
  t0 <- t1; t0$edge.length[2] <- -1
  expect_error(dpcoa(x, t0), 'finite and non-negative')
  xe <- x; xe[2, ] <- 0
  expect_error(dpcoa(xe, t1), "1 sample with no observations \\(all zero\\), which DPCoA cannot place: S2")
  xe <- x; xe[2, 2] <- -1
  expect_error(dpcoa(xe, t1), 'must not contain negative values')
  xe <- x; xe[2, 2] <- NA
  expect_error(dpcoa(xe, t1), 'NA, NaN, or infinite')
  expect_error(dpcoa(x[, 1, drop = FALSE], t1), 'at least 2 taxa')
  expect_error(dpcoa(unname(x), t1), 'must have feature names')
  expect_error(dpcoa(x, t1, n_axes = 0), '`n_axes` must be a whole number')
  expect_error(dpcoa(x, taxa_dm = patristic(t1)[1:4, 1:4]), 'taxon of `counts` is not in `taxa_dm`: "E"')
  expect_error(dpcoa(x, taxa_dm = unname(d1[1:4, 1:4])), 'one row per taxon of `counts` \\(5\\)')
  expect_error(dpcoa(x, taxa_dm = matrix(1:4, 2)), '`taxa_dm` must be a dist object or a symmetric numeric matrix')
  expect_error(dpcoa(x, taxa_dm = -d1), '`taxa_dm` must not contain negative distances')

  # Not Euclidean: a triangle inequality broken, and patristic distances
  # without the square root.
  bad <- matrix(c(0, 1, 3, 3,  1, 0, 1, 1,  3, 1, 0, 1,  3, 1, 1, 0), 4, 4,
                dimnames = list(c('A', 'B', 'C', 'D'), c('A', 'B', 'C', 'D')))
  expect_error(dpcoa(x[, c('A', 'B', 'C', 'D')], taxa_dm = bad), '`taxa_dm` is not Euclidean')
  expect_error(dpcoa(big, taxa_dm = pd_big), '`taxa_dm` is not Euclidean')

#})
