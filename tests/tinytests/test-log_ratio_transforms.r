#test_that("log_ratio_transforms", {

  # alr() and ilr(): scikit-bio 0.7.4's outputs (data-raw/expected/
  # log_ratio_transforms.py), compositions' as a cross-check, identities with
  # normalize(norm = 'clr') and aitchison(), and argument checks.

  max_cpus <- if (at_home()) 4L else 2L

  ref <- readRDS(file.path('expected', 'log_ratio_transforms.rds'))

  # Relative difference, scaled by the largest value.
  rel <- function (a, b) max(abs(a - b)) / max(1, max(abs(b)))

  # The values alone, as a matrix: no dimnames or other attributes.
  plain <- function (m) { attributes(m) <- list(dim = dim(m)); m }


  # -----------------------------------------------------------------------
  # scikit-bio's documented values
  # -----------------------------------------------------------------------

  x <- rbind(c(.1, .3, .4, .2))
  expect_equal(as.vector(alr(x, 1)), c(1.09861229, 1.38629436, 0.69314718), tolerance = 1e-8)
  expect_equal(as.vector(ilr(x)), c(-0.7768362, -0.68339802, 0.11704769), tolerance = 1e-7)
  e <- c(exp(1), 1, 1)
  expect_equal(as.vector(ilr(rbind(e / sum(e)))), c(0.70710678, 0.40824829), tolerance = 1e-8)

  # A vector is one sample, as a column: coordinates x 1.
  v <- ilr(c(.1, .3, .4, .2))
  expect_identical(dim(v), c(3L, 1L))
  expect_equal(as.vector(v), as.vector(ilr(x)))


  # -----------------------------------------------------------------------
  # Bases: scikit-bio's, to the last bit
  # -----------------------------------------------------------------------

  for (D in names(ref$helmert))
    expect_identical(helmert_basis(as.integer(D)), ref$helmert[[D]], info = D)

  for (nm in names(ref$sbp)) {
    expect_identical(unname(sbp_basis(ref$sbp[[nm]]$sbp)), ref$sbp[[nm]]$basis, info = nm)
    b <- ilr_basis(ref$sbp[[nm]]$sbp, ncol(ref$sbp[[nm]]$sbp))
    expect_identical(unname(b$basis), ref$sbp[[nm]]$basis, info = nm)
    expect_true(b$sbp)
  }

  # The default basis is the attribute, and the negated, normalized transpose
  # of contr.helmert().
  h <- stats::contr.helmert(5)
  expect_equal(plain(helmert_basis(5)), plain(-t(h) / sqrt(colSums(h^2))), tolerance = 1e-15)


  # -----------------------------------------------------------------------
  # scikit-bio's outputs
  # -----------------------------------------------------------------------

  for (i in seq_along(ref$skbio)) {
    cs  <- ref$skbio[[i]]
    tab <- ref$tables[[cs$table]]
    pc  <- if (cs$pc == 0) NULL else cs$pc
    info <- paste(i, cs$table, cs$fn, cs$ref_idx, cs$kind)

    res <- if (cs$fn == 'alr') {
      alr(tab, cs$ref_idx + 1L, pseudocount = pc, cpus = 1L)
    } else {
      ilr(tab, basis = cs$basis, pseudocount = pc, cpus = 1L)
    }

    if (cs$fn == 'alr') {
      # The same subtractions of the same logs: identical.
      expect_identical(plain(res), cs$result, info = info)
    } else {
      expect_true(rel(plain(res), cs$result) < 1e-12, info = info)
    }
    expect_identical(attr(res, 'pseudocount'), cs$pc, info = info)
  }


  # -----------------------------------------------------------------------
  # compositions (output only): its alr divides by the last part; its
  # default ilr basis is scikit-bio's with the opposite sign.
  # -----------------------------------------------------------------------

  if (!is.null(ref$compositions)) {
    for (nm in names(ref$tables)) {
      tab <- ref$tables[[nm]]
      pc  <- ref$skbio[[match(nm, vapply(ref$skbio, function (s) s$table, ''))]]$pc
      pc  <- if (pc == 0) NULL else pc
      a   <- alr(tab, ncol(tab), pseudocount = pc, cpus = 1L)
      b   <- ilr(tab, pseudocount = pc, cpus = 1L)
      expect_true(rel(as.vector(a),  as.vector(ref$compositions[[nm]]$alr)) < 1e-12, info = nm)
      expect_true(rel(as.vector(-b), as.vector(ref$compositions[[nm]]$ilr)) < 1e-12, info = nm)
    }
  }


  # -----------------------------------------------------------------------
  # Identities
  # -----------------------------------------------------------------------

  set.seed(1)
  n <- 130L; D <- 9L
  cm <- matrix(rpois(n * D, 3) * rep(c(1, 5, 20), length.out = D), n, D, byrow = TRUE,
    dimnames = list(paste0('S', 1:n), paste0('F', 1:D)) )
  sbp <- rbind(
    c( 1,  1,  1,  1, -1, -1, -1, -1, -1),
    c( 1,  1, -1, -1,  0,  0,  0,  0,  0),
    c( 1, -1,  0,  0,  0,  0,  0,  0,  0),
    c( 0,  0,  1, -1,  0,  0,  0,  0,  0),
    c( 0,  0,  0,  0,  1,  1, -1, -1, -1),
    c( 0,  0,  0,  0,  1, -1,  0,  0,  0),
    c( 0,  0,  0,  0,  0,  0,  1, -1, -1),
    c( 0,  0,  0,  0,  0,  0,  0,  1, -1) )
  colnames(sbp) <- colnames(cm)
  rownames(sbp) <- paste0('b', 1:8)
  qb <- qr.Q(qr(scale(matrix(rnorm(D * (D - 1)), D), scale = FALSE)))
  ob <- t(qb)

  clr <- normalize(cm, 'clr', pseudocount = 0.5)
  ait <- aitchison(cm, pseudocount = 0.5)

  for (basis in list(NULL, sbp, ob)) {
    r <- ilr(cm, basis = basis, pseudocount = 0.5)
    B <- attr(r, 'basis')
    expect_true(is.matrix(B) && identical(dim(B), c(D - 1L, D)))
    expect_true(max(abs(plain(r) %*% B - clr)) < 1e-12)
    expect_true(max(abs(euclidean(plain(r)) - ait)) < 1e-12)
    expect_identical(colnames(B), colnames(cm))

    # Same result for any cpus.
    for (cpus in 2:max_cpus)
      expect_identical(ilr(cm, basis = basis, pseudocount = 0.5, cpus = cpus), r)

    # margin = 2 gives the transpose.
    r2 <- ilr(t(cm), basis = basis, margin = 2L, pseudocount = 0.5)
    expect_identical(unname(t(plain(r2))), unname(plain(r)))
    expect_identical(dimnames(r2), rev(dimnames(r)))

    # Scale invariance: proportions (no zeros) give the same coordinates.
    p  <- cm + 0.5
    p2 <- p / rowSums(p) * c(1, 1000)
    expect_true(max(abs(ilr(p2, basis = basis) - plain(r))) < 1e-12)
  }

  # SBP: its sbp_basis() matrix, and row names as coordinate names.
  rs <- ilr(cm, basis = sbp, pseudocount = 0.5)
  expect_identical(colnames(rs), rownames(sbp))
  expect_equal(ilr(cm, basis = sbp_basis(sbp), pseudocount = 0.5), rs)
  expect_identical(attr(rs, 'basis'), sbp_basis(sbp))
  # The first balance, by hand: sqrt(rs / (r + s)) log(g(num) / g(den)).
  lg <- log(cm + 0.5)
  bal <- sqrt(4 * 5 / 9) * (rowMeans(lg[, 1:4]) - rowMeans(lg[, 5:9]))
  expect_true(max(abs(rs[, 'b1'] - bal)) < 1e-12)

  # Basis columns are matched to features by name.
  perm <- c(9, 1, 8, 2, 7, 3, 6, 4, 5)
  expect_equal(ilr(cm, basis = sbp[, perm], pseudocount = 0.5), rs)
  expect_equal(ilr(cm, basis = ob, pseudocount = 0.5)[, ],
               ilr(cm[, perm], basis = `colnames<-`(ob, colnames(cm))[, perm], pseudocount = 0.5)[, ])
  # ... or by position when either has no names.
  expect_equal(plain(ilr(unname(cm), basis = sbp, pseudocount = 0.5)),
               plain(rs))

  # alr: log ratios, names, attributes.
  for (k in c(1L, 5L, D)) {
    a <- alr(cm, k, pseudocount = 0.5)
    expect_equal(plain(a), plain(log((cm[, -k] + 0.5) / (cm[, k] + 0.5))), tolerance = 1e-14)
    expect_identical(colnames(a), colnames(cm)[-k])
    expect_identical(attr(a, 'denominator'), colnames(cm)[[k]])
    expect_identical(attr(a, 'method'), 'alr')
    expect_identical(alr(cm, colnames(cm)[[k]], pseudocount = 0.5), a)
    expect_identical(alr(cm, seq_len(D) == k, pseudocount = 0.5), a)
    for (cpus in 2:max_cpus)
      expect_identical(alr(cm, k, pseudocount = 0.5, cpus = cpus), a)
    a2 <- alr(t(cm), k, margin = 2L, pseudocount = 0.5)
    expect_identical(t(plain(a2)), plain(a))
  }
  au <- alr(unname(cm), 2, pseudocount = 0.5)
  expect_identical(colnames(au), paste0('alr', c(1, 3:D)))
  expect_identical(attr(au, 'denominator'), 2L)
  expect_null(rownames(au))
  iu <- ilr(unname(cm), pseudocount = 0.5)
  expect_identical(colnames(iu), paste0('ilr', 1:(D - 1)))
  expect_null(colnames(attr(iu, 'basis')))

  # Zeros: the shared pseudocount rule.
  expect_warning(z <- ilr(cm), 'pseudocount = 0.5')
  expect_identical(attr(z, 'pseudocount'), 0.5)
  expect_identical(z, ilr(cm, pseudocount = 0.5))
  expect_identical(attr(ilr(cm + 1), 'pseudocount'), 0)
  # A constant table is fine: every coordinate is 0.
  expect_true(all(alr(matrix(3, 2, 3), 1) == 0))

  # Sparse and dense agree.
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(cm, sparse = TRUE)
    expect_identical(ilr(sp, pseudocount = 0.5), ilr(cm, pseudocount = 0.5))
    expect_identical(ilr(sp, basis = sbp, pseudocount = 0.5), rs)
    expect_identical(alr(sp, 3, pseudocount = 0.5), alr(cm, 3, pseudocount = 0.5))
    expect_identical(ilr(Matrix::t(sp), margin = 2L, pseudocount = 0.5),
                     ilr(t(cm), margin = 2L, pseudocount = 0.5))
    expect_identical(ilr(methods::as(sp, 'TsparseMatrix'), pseudocount = 0.5),
                     ilr(cm, pseudocount = 0.5))
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(cm)
    expect_identical(alr(st, 3, pseudocount = 0.5), alr(cm, 3, pseudocount = 0.5))
    expect_identical(ilr(st, pseudocount = 0.5), ilr(cm, pseudocount = 0.5))
  }

  # Integer and logical tables.
  ci <- cm; storage.mode(ci) <- 'integer'
  expect_identical(ilr(ci, pseudocount = 0.5), ilr(cm, pseudocount = 0.5))

  # The default basis is not stored above 1000 features.
  big <- matrix(1:2002, 2, 1001)
  expect_identical(attr(ilr(big), 'basis'), 'helmert')
  expect_identical(dim(attr(ilr(big[, 1:1000]), 'basis')), c(999L, 1000L))
  # Many features: the recurrence agrees with the dense product.
  hb <- ilr(big[, 1:1000])
  expect_true(max(abs(plain(hb) %*% attr(hb, 'basis') -
                      normalize(big[, 1:1000], 'clr'))) < 1e-12)

  # Two features: one coordinate, log(x1 / x2) / sqrt(2).
  t2 <- matrix(c(1, 2, 3, 5, 7, 11), 3)
  expect_equal(as.vector(ilr(t2)), log(t2[, 1] / t2[, 2]) / sqrt(2), tolerance = 1e-15)
  expect_equal(as.vector(ilr(t2, basis = rbind(c(-1, 1)))), -log(t2[, 1] / t2[, 2]) / sqrt(2), tolerance = 1e-15)


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  expect_error(alr(cm), '`denominator` is required')
  expect_error(alr(cm, 'nope', pseudocount = 1), 'not a feature of `counts`')
  expect_error(alr(unname(cm), 'F1', pseudocount = 1), 'no feature names')
  expect_error(alr(cm, 0, pseudocount = 1), 'not a feature index')
  expect_error(alr(cm, D + 1, pseudocount = 1), 'not a feature index')
  expect_error(alr(cm, 1.5, pseudocount = 1), 'must be one feature')
  expect_error(alr(cm, c(1, 2), pseudocount = 1), 'must be one feature')
  expect_error(alr(cm, NA, pseudocount = 1), 'must be one feature')
  expect_error(alr(cm, c(TRUE, FALSE), pseudocount = 1), 'one value per feature')
  expect_error(alr(cm[, 1, drop = FALSE], 1), 'at least 2 features')
  expect_error(ilr(cm[, 1, drop = FALSE]), 'at least 2 features')

  neg <- cm; neg[1, 1] <- -1
  expect_error(alr(neg, 1, pseudocount = 1), 'negative values')
  expect_error(ilr(neg, pseudocount = 1), 'negative values')
  nas <- cm; nas[1, 1] <- NA
  expect_error(ilr(nas, pseudocount = 1), 'NA, NaN, or infinite')
  expect_error(ilr(cm * 0, pseudocount = 1), 'only zeros')
  expect_error(ilr(cm, pseudocount = 0), '`pseudocount` must be a single positive number')
  expect_error(ilr(cm, pseudocount = -1), '`pseudocount` must be a single positive number')
  expect_error(ilr(cm, margin = 3, pseudocount = 1), '`margin` must be 1 or 2')

  expect_error(ilr(cm, basis = sbp[-1, ], pseudocount = 1), 'must be a 8 x 9 matrix')
  expect_error(ilr(cm, basis = t(sbp), pseudocount = 1), 'Transposed')
  expect_error(ilr(cm, basis = 'helmert', pseudocount = 1), 'numeric matrix')
  bad <- sbp; bad[3, ] <- c(1, 1, 0, 0, 0, 0, 0, 0, 0)
  expect_error(ilr(cm, basis = bad, pseudocount = 1), 'row 3 has no -1')
  bad <- sbp; bad[3, ] <- c(1, 0, -1, 0, 0, 0, 0, 0, 0)
  expect_error(ilr(cm, basis = bad, pseudocount = 1), 'not a valid sequential binary partition')
  expect_error(ilr(cm, basis = ob * 1.1, pseudocount = 1), 'not orthonormal')
  # Orthonormal rows that don't sum to 0.
  qq <- t(qr.Q(qr(matrix(rnorm(D * D), D))))[1:(D - 1), ]
  expect_error(ilr(cm, basis = qq, pseudocount = 1), 'must sum to 0')
  wrong <- sbp; colnames(wrong)[1] <- 'X'
  expect_error(ilr(cm, basis = wrong, pseudocount = 1), "Not in `basis`: F1. Not in `counts`: X")
  dup <- sbp; colnames(dup)[2] <- 'F1'
  expect_error(ilr(cm, basis = dup, pseudocount = 1), 'duplicated column names')
  nab <- ob; nab[1, 1] <- NA
  expect_error(ilr(cm, basis = nab, pseudocount = 1), 'finite values')

#})
