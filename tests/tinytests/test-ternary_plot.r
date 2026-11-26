#test_that("ternary_plot", {

  # Draws on a null device; returns the result, with whether it was visible.
  tp <- function (...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    v <- withVisible(ternary_plot(...))
    attr(v$value, 'visible') <- v$visible
    v$value
  }

  # The title that a call draws, from the recorded display list.
  tp_title <- function (...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control('enable')
    ternary_plot(...)
    calls <- grDevices::recordPlot()[[1]]
    for (cl in calls) {
      fn <- cl[[2]][[1]]
      if (identical(fn$name, 'C_title')) return (cl[[2]][[2]])
    }
    NULL
  }

  V   <- rbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6))
  ilr <- function (P) log(as.matrix(P)) %*% t(V)
  cl  <- function (m) m / rowSums(m)

  set.seed(2)
  n <- 24; p <- 9
  grp <- rep(c('A', 'B', 'C'), each = n / 3)
  x <- matrix(stats::rpois(n * p, rep(exp(seq(4, 1, length.out = p)), each = n) *
                rep(c(1, 2, 0.5), each = n / 3)) + 1L, n, p,
              dimnames = list(paste0('S', 1:n), paste0('F', 1:p)))
  x3 <- x[, 1:3]
  P3 <- cl(x3)



  # Barycentric coordinates. ====

  bx <- matrix(c(5, 0, 0,  0, 5, 0,  0, 0, 5,  2, 2, 2), ncol = 3, byrow = TRUE,
               dimnames = list(c('a', 'b', 'c', 'd'), c('P1', 'P2', 'P3')))
  r <- tp(bx, parts = 1:3)
  expect_false(attr(r, 'visible'))
  expect_identical(names(r), c('P1', 'P2', 'P3', 'x', 'y'))
  expect_identical(rownames(r), c('a', 'b', 'c', 'd'))
  expect_equal(unname(as.matrix(r[1:3, c('x', 'y')])),
               rbind(c(0, 0), c(1, 0), c(1 / 2, sqrt(3) / 2)), tolerance = 1e-15)
  expect_equal(unname(unlist(r[4, c('x', 'y')])), c(1 / 2, sqrt(3) / 6), tolerance = 1e-15)
  expect_identical(attr(r, 'mode'), 'samples')
  expect_null(attr(r, 'center'))
  expect_identical(attr(r, 'pseudocount'), 0)
  expect_identical(attr(r, 'ellipses'), list())

  # Without logs, zeros are drawn on the edges, and no pseudocount is used.
  r <- tp(bx, parts = c('P3', 'P1', 'P2'), pseudocount = 1)
  expect_identical(names(r)[1:3], c('P3', 'P1', 'P2'))   # the order given
  expect_equal(unname(as.matrix(r[1:3, 1:3])), rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0)))
  expect_identical(attr(r, 'pseudocount'), 0)

  # Proportions, and x and y, of any table.
  r <- tp(x, parts = c('F1', 'F2', 'F3'))
  expect_equal(unname(as.matrix(r[, 1:3])), unname(P3), tolerance = 1e-15)
  expect_equal(r$x, unname(P3[, 2] + P3[, 3] / 2), tolerance = 1e-15)
  expect_equal(r$y, unname(sqrt(3) / 2 * P3[, 3]), tolerance = 1e-15)



  # Centring. ====

  r  <- tp(x, parts = 1:3, center = TRUE)
  Pc <- as.matrix(r[, 1:3])
  gm <- exp(colMeans(log(Pc)))
  expect_equal(unname(gm / sum(gm)), rep(1 / 3, 3), tolerance = 1e-14)
  expect_equal(unname(colMeans(ilr(Pc))), c(0, 0), tolerance = 1e-14)
  expect_equal(unname(rowSums(Pc)), rep(1, n), tolerance = 1e-15)

  g <- exp(colMeans(log(P3)))
  expect_equal(unname(attr(r, 'center')), unname(g / sum(g)), tolerance = 1e-14)
  expect_identical(names(attr(r, 'center')), c('F1', 'F2', 'F3'))

  # Differences between samples' ilr coordinates are unchanged.
  zu <- ilr(P3); zc <- ilr(Pc)
  expect_equal(unname(zc[-1, ] - zc[rep(1, n - 1), ]),
               unname(zu[-1, ] - zu[rep(1, n - 1), ]), tolerance = 1e-12)
  # ... and the centred points are the uncentred ones perturbed by 1/g.
  expect_equal(unname(Pc), unname(cl(sweep(P3, 2, g, '/'))), tolerance = 1e-14)



  # Ellipses. ====

  for (kind in c('confidence', 'prediction')) {
    for (lvl in c(0.95, 0.5)) {
      r  <- tp(x, parts = 1:3, groups = grp, ellipses = kind, conf_level = lvl)
      el <- attr(r, 'ellipses')
      expect_identical(names(el), c('A', 'B', 'C'))
      for (k in names(el)) {
        e  <- el[[k]]
        expect_identical(dim(e), c(200L, 5L))
        expect_identical(colnames(e), c('F1', 'F2', 'F3', 'x', 'y'))
        z  <- ilr(P3[grp == k, ])
        m  <- nrow(z)
        f  <- stats::qf(lvl, 2, m - 2)
        r2 <- if (kind == 'confidence') 2 * (m - 1) / (m * (m - 2)) * f
              else 2 * (m + 1) * (m - 1) / (m * (m - 2)) * f
        d  <- sweep(ilr(e[, 1:3]), 2, colMeans(z))
        q  <- rowSums((d %*% solve(stats::cov(z))) * d)
        expect_equal(q, rep(r2, 200), tolerance = 1e-10, info = paste(kind, lvl, k))
        expect_equal(unname(e[, 'x']), unname(e[, 2] + e[, 3] / 2), tolerance = 1e-15)
      }
    }
  }

  # Centred and uncentred curves differ by the same perturbation.
  ru <- tp(x, parts = 1:3, groups = grp, ellipses = 'prediction')
  rc <- tp(x, parts = 1:3, groups = grp, ellipses = 'prediction', center = TRUE)
  g  <- attr(rc, 'center')
  for (k in c('A', 'B', 'C'))
    expect_equal(attr(rc, 'ellipses')[[k]][, 1:3],
                 cl(sweep(attr(ru, 'ellipses')[[k]][, 1:3], 2, g, '/')), tolerance = 1e-13)

  # Without groups: one ellipse for all the samples.
  r <- tp(x, parts = 1:3, ellipses = 'confidence')
  expect_identical(names(attr(r, 'ellipses')), 'all')
  expect_false('group' %in% names(r))

  # Groups too small for an ellipse: one warning naming them.
  g2 <- c(rep('A', n - 3), 'B', 'B', 'C')
  expect_warning(
    r <- tp(x, parts = 1:3, groups = g2, ellipses = 'confidence'),
    "No ellipse for 2 groups with fewer than 3 samples or collinear log ratios: 'B', 'C'" )
  expect_identical(names(attr(r, 'ellipses')), 'A')
  expect_identical(levels(r$group), c('A', 'B', 'C'))

  # Collinear log ratios: a singular covariance.
  xs <- x3
  xs[, 3] <- xs[, 1]
  expect_warning(
    tp(xs, parts = 1:3, ellipses = 'confidence'),
    "No ellipse for 1 group with fewer than 3 samples or collinear log ratios: 'all'" )

  # Numeric groups: a gradient, no ellipses.
  num <- seq_len(n) / 3
  expect_warning(r <- tp(x, parts = 1:3, groups = num, ellipses = 'confidence'), 'gradient')
  expect_identical(r$group, num)
  expect_identical(attr(r, 'ellipses'), list())
  r <- tp(x, parts = 1:3, groups = num, grid = FALSE)
  expect_identical(r$group, num)

  # Groups by name, from `data`, matched by row name.
  meta <- data.frame(site = grp, row.names = rownames(x))[n:1, , drop = FALSE]
  r <- tp(x, parts = 1:3, groups = 'site', data = meta)
  expect_identical(as.character(r$group), grp)



  # Pseudocounts. ====

  xz <- x
  xz[1:3, 1] <- 0L
  xz[5, 6]   <- 0L   # outside the parts: still counts, by the whole-table rule
  expect_warning(
    r <- tp(xz, parts = list(1, 2, 3:4), center = TRUE),
    'Zeros detected in data. Using pseudocount = 0.5' )
  expect_identical(attr(r, 'pseudocount'), 0.5)
  y  <- cbind(xz[, 1], xz[, 2], xz[, 3] + xz[, 4]) + 0.5
  g  <- exp(colMeans(log(cl(y))))
  expect_equal(unname(as.matrix(r[, 1:3])), unname(cl(sweep(cl(y), 2, g, '/'))), tolerance = 1e-14)

  r <- tp(xz, parts = list(1, 2, 3:4), ellipses = 'prediction', pseudocount = 2)
  expect_identical(attr(r, 'pseudocount'), 2)
  expect_equal(unname(as.matrix(r[, 1:3])), unname(cl(y - 0.5 + 2)), tolerance = 1e-15)

  # Without logs, the pseudocount is not used.
  r <- tp(xz, parts = list(1, 2, 3:4), pseudocount = 2)
  expect_identical(attr(r, 'pseudocount'), 0)
  expect_equal(unname(r[1:3, 1]), c(0, 0, 0))

  expect_error(tp(x, parts = 1:3, center = TRUE, pseudocount = 0), 'single positive number')
  expect_error(tp(x, parts = 1:3, center = TRUE, pseudocount = -1), 'single positive number')

  # Samples with no counts in the three parts are left out.
  x0 <- x
  x0[c(2, 7), 1:3] <- 0L
  expect_warning(r <- tp(x0, parts = 1:3), "Left out 2 samples whose three parts are all zero: 'S2', 'S7'")
  expect_identical(rownames(r), rownames(x)[-c(2, 7)])
  expect_warning(r <- tp(x0, parts = 1:3, groups = grp), 'Left out 2 samples')
  expect_identical(as.character(r$group), grp[-c(2, 7)])



  # Amalgamation. ====

  r <- tp(x, parts = list(Low = c('F1', 'F2'), Mid = 4:6, High = NA))
  y <- cbind(x[, 1] + x[, 2], rowSums(x[, 4:6]), x[, 3] + rowSums(x[, 7:9]))
  expect_identical(names(r)[1:3], c('Low', 'Mid', 'High'))
  expect_equal(unname(as.matrix(r[, 1:3])), unname(cl(y)), tolerance = 1e-15)

  # Labels: list names, then a single feature's name, "Other" for NA,
  # "Part k" for an unnamed set.
  r <- tp(x, parts = list(c(TRUE, FALSE, TRUE, rep(FALSE, 6)), 'F2', NA))
  expect_identical(names(r)[1:3], c('Part 1', 'F2', 'Other'))
  expect_equal(unname(as.matrix(r[, 1:3])),
               unname(cl(cbind(x[, 1] + x[, 3], x[, 2], rowSums(x[, 4:9])))), tolerance = 1e-15)
  r <- tp(x, parts = c(a = 'F3', b = 'F1', 'F2'))
  expect_identical(names(r)[1:3], c('a', 'b', 'F2'))

  # Unnamed features.
  xu <- unname(x)
  r <- tp(xu, parts = list(1, 2, NA))
  expect_identical(names(r)[1:3], c('Feature 1', 'Feature 2', 'Other'))
  expect_identical(rownames(r), as.character(1:n))

  # Dense, sparse and margin = 2 agree.
  args <- list(parts = list(Low = c('F1', 'F2'), Mid = 4:6, High = NA),
               groups = grp, center = TRUE, ellipses = 'confidence')
  ref  <- do.call(tp, c(list(x), args))
  expect_identical(do.call(tp, c(list(t(x), margin = 2L), args)), ref)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(x, sparse = TRUE)
    expect_identical(do.call(tp, c(list(sp), args)), ref)
    expect_identical(do.call(tp, c(list(Matrix::t(sp), margin = 2L), args)), ref)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(do.call(tp, c(list(slam::as.simple_triplet_matrix(x)), args)), ref)



  # Features mode. ====

  r <- tp(x, groups = grp)
  expect_false(attr(r, 'visible'))
  expect_identical(attr(r, 'mode'), 'features')
  expect_identical(names(r), c('A', 'B', 'C', 'x', 'y', 'abundance'))
  expect_identical(rownames(r), colnames(x))
  rel <- cl(x)
  M <- rbind(colMeans(rel[grp == 'A', ]), colMeans(rel[grp == 'B', ]), colMeans(rel[grp == 'C', ]))
  M <- t(M)
  expect_equal(unname(as.matrix(r[, 1:3])), unname(M / rowSums(M)), tolerance = 1e-14)
  expect_equal(r$abundance, unname(colMeans(rel)), tolerance = 1e-14)
  expect_equal(r$x, r$B + r$C / 2, tolerance = 1e-15)
  expect_null(attr(r, 'center'))
  expect_identical(attr(r, 'ellipses'), list())

  # A factor's level order; NA in no group; features absent from all three.
  f  <- factor(grp, levels = c('C', 'A', 'B'))
  f[1] <- NA
  xf <- x
  xf[, 'F9'] <- 0L
  r <- tp(xf, groups = f)
  expect_identical(names(r)[1:3], c('C', 'A', 'B'))
  expect_identical(rownames(r), colnames(x)[1:8])
  rel <- cl(xf)
  expect_equal(r$abundance, unname(colMeans(rel[-1, 1:8])), tolerance = 1e-14)
  expect_equal(r$A, unname(colMeans(rel[2:8, 1:8]) / rowSums(cbind(
    colMeans(rel[grp == 'C', 1:8]), colMeans(rel[2:8, 1:8]), colMeans(rel[grp == 'B', 1:8])))),
    tolerance = 1e-14)

  # Samples with no counts are left out, with a warning.
  xf[2, ] <- 0L
  expect_warning(tp(xf, groups = grp), "Left out 1 sample whose counts are all zero: 'S2'")

  # Colours: one per feature, or by name; bad names are an error.
  r <- tp(x, groups = grp, col = c(F2 = 'red'), labels = Inf)
  expect_identical(nrow(r), as.integer(p))
  r <- tp(x, groups = grp, col = rep('red', p), pch = 1, labels = FALSE, legend = NULL)
  expect_error(tp(x, groups = grp, col = c(F99 = 'red')), "`col` names features not in `counts`: 'F99'")
  expect_error(tp(x, groups = grp, col = c('red', 'blue')), 'one per feature')

  # Numeric groups with three values.
  r <- tp(x, groups = rep(c(1, 2, 3), each = n / 3))
  expect_identical(names(r)[1:3], c('1', '2', '3'))



  # Titles and drawing. ====

  expect_identical(tp_title(x, parts = 1:3), 'Ternary plot of F1, F2, F3')
  expect_identical(tp_title(x, groups = grp), 'Mean relative abundance in A, B, C')
  expect_identical(tp_title(x, parts = 1:3, main = 'Mine'), 'Mine')
  r <- tp(x, parts = 1:3, grid = FALSE, cex = 2, sub = 'Mine', legend = NULL)
  expect_identical(nrow(r), as.integer(n))

  # R's random number generator is untouched.
  set.seed(5)
  seed <- .Random.seed
  tp(x, parts = 1:3, groups = grp, center = TRUE, ellipses = 'confidence')
  tp(x, groups = grp)
  expect_identical(.Random.seed, seed)



  # Errors. ====

  expect_error(tp(x, parts = 1:2), 'exactly 3 features, or be a list of 3 sets of them \\(it selects 2\\)')
  expect_error(tp(x, parts = 1:4), 'it selects 4')
  expect_error(tp(x, parts = list(1, 2)), 'A list `parts` must have 3 elements, one per part \\(it has 2\\)')
  expect_error(tp(x, parts = list(1, NA, NA)), 'Only one element of `parts` may be NA')
  expect_error(tp(x, parts = list(1:2, 2:3, 4)), "`parts` selects a feature more than once: 'F2'")
  expect_error(tp(x, parts = list(1:3, 4:6, 7:9, NA)), 'it has 4')
  expect_error(tp(x[, 1:3], parts = list(1, 2:3, NA)), 'there are none')
  expect_error(tp(x, parts = c('F1', 'F2', 'Nope')), "`parts` names a feature not in `counts`: 'Nope'")
  expect_error(tp(x, parts = c(1, 1, 2)), 'more than once')
  expect_error(tp(x, groups = rep(c('A', 'B'), each = n / 2)),
               'Without `parts`, `groups` must have exactly 3 levels \\(it has 2\\)')
  expect_error(tp(x, groups = rep(c('A', 'B', 'C', 'D'), each = n / 4)), 'it has 4')
  expect_error(tp(x), 'Give `parts`')
  expect_error(tp(x, groups = grp, center = TRUE), '`center` applies to samples mode only')
  expect_error(tp(x, groups = grp, ellipses = 'confidence'), '`ellipses` apply to samples mode only')
  expect_error(tp(x, parts = 1:3, conf_level = 0), 'strictly between 0 and 1')
  expect_error(tp(x, parts = 1:3, conf_level = 1), 'strictly between 0 and 1')
  expect_error(tp(x, parts = 1:3, ellipses = 'x'), '`ellipses` must be "none", "confidence" or "prediction"')
  expect_error(tp(x, parts = 1:3, center = NA), '`center` must be TRUE or FALSE')
  expect_error(tp(x, parts = 1:3, grid = 'yes'), '`grid` must be TRUE or FALSE')
  expect_error(tp(x, parts = 1:3, labels = -1), '`labels` must be a number')
  expect_error(tp(x, parts = 1:3, legend = 'nowhere'), '`legend` must be one of')
  expect_error(tp(x, parts = 1:3, groups = grp[-1]), 'has 23 values')
  expect_error(tp(-x, parts = 1:3), 'negative')
  x0 <- x
  x0[, 1:3] <- 0L
  expect_error(tp(x0, parts = 1:3), 'nothing to draw')
  ge <- factor(grp, levels = c('A', 'B', 'C'))
  expect_error(tp(x[grp != 'C', ], groups = ge[grp != 'C']), "Group has no samples with counts above 0: 'C'")

#})
