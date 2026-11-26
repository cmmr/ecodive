#test_that("dbmem", {

  ex <- readRDS('expected/dbmem.rds')

  autocors <- c('positive', 'negative', 'non-null', 'all')

  # The largest difference between the columns of a and b (MEMs of length
  # n): column by column up to sign where an eigenvalue (`vals`) is alone,
  # and as the projection onto their span where eigenvalues are equal (any
  # orthogonal basis of their space is right).
  mem_diff <- function (a, b, vals) {
    a <- unclass(a); b <- unclass(b)
    attributes(a) <- list(dim = dim(a)); attributes(b) <- list(dim = dim(b))
    cl <- cumsum(c(TRUE, abs(diff(vals)) > 1e-9 * max(abs(vals))))
    max(0, vapply(unique(cl), function (g) {
      j <- which(cl == g)
      if (length(j) == 1) {
        s <- if (sum(a[, j] * b[, j]) < 0) -1 else 1
        max(abs(a[, j] - s * b[, j]))
      } else {
        max(abs(tcrossprod(a[, j]) - tcrossprod(b[, j]))) / nrow(a)
      }
    }, 0))
  }

  # The longest edge of the minimum spanning tree (Kruskal's algorithm).
  mst_longest <- function (d) {
    d  <- as.matrix(d)
    n  <- nrow(d)
    ij <- which(lower.tri(d), arr.ind = TRUE)
    ij <- ij[order(d[ij]), , drop = FALSE]
    comp <- seq_len(n); longest <- 0; joined <- 0
    for (e in seq_len(nrow(ij))) {
      a <- comp[ij[e, 1]]; b <- comp[ij[e, 2]]
      if (a == b) next
      comp[comp == b] <- a
      longest <- max(longest, d[ij[e, 1], ij[e, 2]])
      if ((joined <- joined + 1) == n - 1) break
    }
    longest
  }



  # Same MEMs, eigenvalues and Moran's I as adespatial 0.3.30's dbmem()
  # (Moran's I from spdep 1.3.1's moran() with adespatial's weights), and the
  # same threshold as vegan's spantree(): stored results, so at home
  # only. ====

  if (at_home()) {
    for (cs in names(ex)) {

      k <- ex[[cs]]
      n <- if (inherits(k$space, 'dist')) attr(k$space, 'Size') else nrow(k$space)

      for (a in autocors) {

        info <- paste(cs, a)
        res  <- if (cs == 'mite_dist_t08') {
          expect_warning(
            dbmem(k$space, threshold = k$threshold, autocor = a, cpus = 1),
            'below the longest edge of the minimum spanning tree', info = info )
          suppressWarnings(dbmem(k$space, threshold = k$threshold, autocor = a, cpus = 1))
        } else {
          dbmem(k$space, threshold = k$threshold, autocor = a, cpus = 1)
        }

        v <- k$values[[a]]
        expect_identical(ncol(res), length(v), info = info)
        expect_equal(unname(attr(res, 'values')), v, tolerance = 1e-12, scale = 1, info = info)
        expect_equal(attr(res, 'threshold'), k$t, tolerance = 1e-15, info = info)
        expect_equal(attr(res, 'threshold'), if (is.null(k$threshold)) k$mst else k$threshold, info = info)

        # The stored MEMs are those of "all" (or of "positive", the first
        # columns of those of "all"); each autocor's MEMs are some of them,
        # in order.
        va  <- k$values$all
        tol <- sqrt(.Machine$double.eps) * max(abs(va))
        sel <- switch(a,
          positive   = which(va >  tol),
          negative   = which(va < -tol),
          `non-null` = which(abs(va) > tol),
          all        = seq_along(va) )
        cols <- which(sel <= ncol(k$vectors))
        j    <- sel[cols]
        if (length(j)) {
          expect_true(
            mem_diff(res[, cols, drop = FALSE], k$vectors[, j, drop = FALSE], va[j]) < 1e-10,
            info = info )
          expect_equal(
            unname(attr(res, 'moran_i'))[cols], k$moran[j],
            tolerance = 1e-12, scale = 1, info = info )
        }
      }

      # vegan's pcnm(): PCoA of the truncated distances, with eigenvalues
      # 8 t^2 (n * values + 1) (the positive ones are compared: with
      # duplicate positions, pcnm() leaves some others out), and (for mite)
      # its unit eigenvectors.
      if (!is.null(k$pcnm)) {
        all <- dbmem(k$space, threshold = k$threshold, autocor = 'all', cpus = 1)
        pv  <- 8 * k$t^2 * (n * unname(attr(all, 'values')) + 1)
        np  <- sum(k$pcnm$values > 1e-8 * max(k$pcnm$values))
        expect_true(np > 0 && all(pv[seq_len(np)] > 0) && !(pv[[np + 1]] > 1e-8 * max(pv)), info = cs)
        expect_equal(pv[seq_len(np)], k$pcnm$values[seq_len(np)], tolerance = 1e-12, scale = max(pv), info = cs)
        if (length(k$pcnm$values) == n - 1)
          expect_equal(pv, k$pcnm$values, tolerance = 1e-12, scale = max(pv), info = cs)
        expect_equal(k$pcnm$threshold, k$t, info = cs)
        if (!is.null(k$pcnm$vectors)) {
          expect_identical(ncol(k$pcnm$vectors), np, info = cs)
          expect_true(
            mem_diff(all[, seq_len(np)] / sqrt(n), k$pcnm$vectors, pv[seq_len(np)]) < 1e-10,
            info = cs )
        }
      }
    }
  }



  # Independent checks: the threshold is the MST's longest edge; the MEMs are
  # orthogonal, centred eigenvectors of the doubly centred weights; and their
  # Moran's I is that of moran_i() with the weights as given. ====

  set.seed(1)
  pts <- data.frame(east = runif(60, 0, 10), north = runif(60, 0, 5), row.names = paste0('P', 1:60))
  d   <- dist(pts)
  res <- dbmem(pts, autocor = 'all', cpus = 1)
  n   <- 60L
  t   <- mst_longest(d)

  expect_equal(attr(res, 'threshold'), t, tolerance = 1e-15)
  expect_identical(dim(res), c(60L, 59L))
  expect_identical(rownames(res), rownames(pts))
  expect_identical(colnames(res), paste0('MEM', 1:59))
  expect_true(max(abs(colMeans(unclass(res)))) < 1e-12)
  expect_equal(crossprod(unclass(res)), diag(n, 59), tolerance = 1e-12, scale = 1, check.attributes = FALSE)
  expect_true(all(diff(attr(res, 'values')) <= 0))

  D  <- as.matrix(d)
  W  <- ifelse(D <= t, 1 - (D / (4 * t))^2, 0); diag(W) <- 0
  C  <- diag(n) - 1 / n
  CW <- C %*% W %*% C
  M  <- unclass(res); attributes(M) <- list(dim = dim(M))
  lam <- n * unname(attr(res, 'values'))
  expect_true(max(abs(CW %*% M - sweep(M, 2, lam, '*'))) < 1e-12)
  expect_equal(unname(attr(res, 'moran_i')), n * lam / sum(W), tolerance = 1e-12, scale = 1)

  for (j in c(1, 10, 59)) {
    mi <- suppressWarnings(moran_i(M[, j], W, style = 'none', permutations = 0))
    expect_equal(mi$statistic[['I']], unname(attr(res, 'moran_i'))[[j]], tolerance = 1e-12, info = j)
  }

  # A larger set, against R's eigen() (BLAS-dependent, so at home only).
  if (at_home()) {
    set.seed(2)
    xy <- matrix(runif(600), ncol = 2)
    r  <- dbmem(xy, autocor = 'all', cpus = 2)
    Dx <- as.matrix(dist(xy)); tx <- attr(r, 'threshold')
    Wx <- ifelse(Dx <= tx, 1 - (Dx / (4 * tx))^2, 0); diag(Wx) <- 0
    Cx <- diag(300) - 1 / 300
    e  <- eigen(Cx %*% Wx %*% Cx, symmetric = TRUE)
    ev <- e$values[order(abs(e$values))][-1]       # without the ones vector's 0
    ev <- sort(ev, decreasing = TRUE)
    expect_equal(300 * unname(attr(r, 'values')), ev, tolerance = 1e-12, scale = 1)
    keep <- abs(e$values) > 1e-12
    ref  <- e$vectors[, keep] * sqrt(300)
    expect_true(mem_diff(r[, seq_len(sum(keep))], ref, e$values[keep]) < 1e-9)
  }



  # The forms of `space`, and `data`. ====

  r_d <- dbmem(d, cpus = 1)
  r_c <- dbmem(as.matrix(pts), cpus = 1)
  r_n <- dbmem(c('east', 'north'), pts, cpus = 1)
  expect_identical(unclass(r_d), unclass(r_c))
  expect_identical(r_c, r_n)
  expect_identical(ncol(r_d), sum(attr(res, 'values') > 0))
  expect_equal(unclass(r_d), unclass(res)[, seq_len(ncol(r_d))], check.attributes = FALSE)
  expect_identical(rownames(r_d), rownames(pts))
  expect_null(rownames(dbmem(unname(as.matrix(pts)), cpus = 1)))
  expect_identical(ncol(dbmem(d, autocor = 'pos', cpus = 1)), ncol(r_d))
  expect_identical(attr(r_d, 'method'), 'dbmem')
  expect_identical(class(r_d), c('ecodive_dbmem', 'matrix', 'array'))
  expect_identical(
    names(attributes(r_d)), c('dim', 'dimnames', 'values', 'moran_i', 'threshold', 'method', 'class') )
  expect_identical(names(attr(r_d, 'values')), colnames(r_d))

  # Parts taken with `[` are plain matrices; as.data.frame() works.
  expect_identical(class(r_d[, 1:2]), c('matrix', 'array'))
  expect_identical(names(as.data.frame(r_d)), colnames(r_d))



  # A custom threshold, and one below the MST's longest edge. ====

  r2 <- dbmem(d, threshold = 2 * t, cpus = 1)
  expect_identical(attr(r2, 'threshold'), 2 * t)
  expect_warning(
    dbmem(d, threshold = t * 0.6, cpus = 1),
    'below the longest edge of the minimum spanning tree' )

  # Two clusters, far apart: a sample with no neighbours.
  cl <- cbind(c(0, 1, 2, 50, 51, 60), 0)
  expect_warning(dbmem(cl, threshold = 1.5, cpus = 1), '1 sample has none')
  expect_error(dbmem(cl, threshold = 0.5, cpus = 1), 'No two samples are within `threshold` \\(0.5\\)')

  # Pairs at the threshold, to rounding, are neighbours.
  g  <- expand.grid(x = 0:4 * 0.1, y = 0:2 * 0.1)
  ra <- dbmem(g, threshold = 0.1, autocor = 'all', cpus = 1)
  rb <- dbmem(g, threshold = 0.1 * (1 - 1e-9), autocor = 'all', cpus = 1)
  expect_equal(unname(attr(ra, 'values')), unname(attr(rb, 'values')), tolerance = 1e-7)
  expect_identical(attr(dbmem(g, cpus = 1), 'threshold'), mst_longest(dist(g)))



  # Duplicates, collinear points, and small n. ====

  dup <- rbind(as.matrix(pts), Q1 = unlist(pts[3, ]), Q2 = unlist(pts[3, ]))
  rd  <- dbmem(dup, cpus = 1)
  expect_equal(unclass(rd)['P3', ], unclass(rd)['Q1', ])
  expect_equal(unclass(rd)['P3', ], unclass(rd)['Q2', ])

  tr <- dbmem(cbind(1:20, 0), autocor = 'all', cpus = 1)
  expect_identical(dim(tr), c(20L, 19L))
  expect_identical(attr(tr, 'threshold'), 1)
  # An evenly spaced line is symmetric, so every MEM is too, or antisymmetric.
  for (j in 1:19)
    expect_equal(abs(unclass(tr)[, j]), rev(abs(unclass(tr)[, j])), tolerance = 1e-10, info = j)

  r3 <- dbmem(data.frame(x = c(0, 1, 3)), autocor = 'all', cpus = 1)
  expect_identical(dim(r3), c(3L, 2L))
  expect_identical(attr(r3, 'threshold'), 2)
  expect_identical(
    ncol(dbmem(data.frame(x = c(0, 1, 3)), autocor = 'negative', cpus = 1)), 1L )

  # An equilateral triangle: every pair a neighbour, at the same distance, so
  # both eigenvalues are equal, and negative.
  tri <- data.frame(x = c(0, 1, 0.5), y = c(0, 0, sqrt(3) / 2))
  rt  <- dbmem(tri, autocor = 'all', cpus = 1)
  expect_identical(dim(rt), c(3L, 2L))
  expect_equal(attr(rt, 'values')[[1]], attr(rt, 'values')[[2]])
  expect_true(all(attr(rt, 'values') < 0))
  expect_warning(dbmem(tri, cpus = 1), 'No MEM has positive spatial autocorrelation')
  r0 <- suppressWarnings(dbmem(tri, cpus = 1))
  expect_identical(dim(r0), c(3L, 0L))



  # Results are identical for any `cpus`. ====

  set.seed(3)
  big <- matrix(runif(2 * 250), ncol = 2)
  cpus <- if (at_home()) c(1, 2, 3, 4, 8) else c(1, 2)
  for (a in c('positive', 'all')) {
    ref <- dbmem(big, autocor = a, cpus = 1)
    for (cp in cpus[-1])
      expect_identical(dbmem(big, autocor = a, cpus = cp), ref, info = paste(a, cp))
  }



  # Printing and plotting. ====

  out <- capture.output(print(r_d))
  expect_true(any(grepl("Moran's Eigenvector Maps", out, fixed = TRUE)))
  expect_true(any(grepl('60 samples', out, fixed = TRUE)))
  expect_true(any(grepl('Threshold:', out, fixed = TRUE)))
  expect_identical(capture.output(r_d), out)

  pdf(NULL)
  expect_identical(plot(r_d), r_d)
  expect_identical(plot(r_d, pts), r_d)
  expect_identical(plot(r_d, pts, mems = 'MEM2', main = 'MEM 2'), r_d)
  expect_identical(plot(r_d, pts, mems = 1:6), r_d)
  expect_error(plot(r_d, mems = 1), '`mems` needs `coords`')
  expect_error(plot(r_d, pts[1:5, ]), 'one row per sample')
  expect_error(plot(r_d, pts, mems = 99), '`mems` must be')
  expect_error(plot(r_d, pts, mems = 'MEM99'), '`mems` names MEMs')
  expect_error(plot(r0), 'no MEMs')
  x_lost <- r_d; attr(x_lost, 'moran_i') <- NULL
  expect_error(plot(x_lost), 'lost its Moran')
  dev.off()



  # Invalid arguments. ====

  expect_error(dbmem(), '`space` is required')
  expect_error(dbmem(d, autocor = 'up'), '`autocor` must be')
  expect_error(dbmem(d, autocor = c('positive', 'all')), '`autocor` must be')
  expect_error(dbmem(d, threshold = 0), '`threshold` must be a positive number')
  expect_error(dbmem(d, threshold = c(1, 2)), '`threshold` must be a positive number')
  expect_error(dbmem(d, threshold = NA), '`threshold` must be a positive number')
  expect_error(dbmem(d, threshold = Inf), '`threshold` must be a positive number')
  expect_error(dbmem(D), 'not a square matrix of spatial weights')
  expect_error(dbmem(dist(1:2)), 'at least 3 samples')
  expect_error(dbmem(data.frame(x = rep(1, 4))), 'at the same place')
  expect_error(dbmem(c('east', 'up'), pts), 'not in `data`')
  expect_error(dbmem(c('east', 'north')), '`data` is not given')
  expect_error(dbmem(data.frame(x = c(1, NA, 3))), 'must not contain NA')
  expect_error(dbmem(list(1, 2)), '`space` must be')
  expect_error(dbmem(d, cpus = 0))

#})
