#test_that("coda_diagnostics", {

  # variation_array() and zero_pattern(). References from
  # data-raw/expected/coda_diagnostics.r: compositions' variation() and
  # mvar(), and zCompositions' zPatterns() (outputs only).

  expected <- readRDS(file.path('expected', 'coda_diagnostics.rds'))

  # Draws on a null device; returns the result, with whether it was visible.
  drawn <- function (x, ...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    v <- withVisible(plot(x, ...))
    attr(v$value, 'visible') <- v$visible
    v$value
  }

  # The calls a plot records, as list(name, args).
  recorded <- function (x, ...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control('enable')
    plot(x, ...)
    lapply(grDevices::recordPlot()[[1]], function (cl) list(name = cl[[2]][[1]]$name, args = cl[[2]][-1]))
  }

  rel <- function (a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

  set.seed(4)
  n <- 30L; D <- 9L
  x <- matrix(
    stats::rpois(n * D, rep(exp(seq(4, 0.5, length.out = D)), each = n)) + 1L, n, D,
    dimnames = list(paste0('S', 1:n), paste0('F', 1:D)) )
  xz <- x - 1L                      # with zeros
  xz[1, 2] <- 0L
  up <- function (m) m[upper.tri(m)]
  lo <- function (m) m[lower.tri(m)]
  pm <- function (m) { a <- attributes(m); m <- unclass(m); attributes(m) <- a[c('dim', 'dimnames')]; m }



  # variation_array(): the identities. ====

  va  <- variation_array(x)
  tau <- as.matrix(vlr(x))
  L   <- log(x)
  mu  <- colMeans(L)
  clr <- L - rowMeans(L)

  expect_identical(class(va), c('ecodive_variation_array', 'matrix', 'array'))
  expect_identical(dim(va), c(D, D))
  expect_identical(dimnames(va), list(colnames(x), colnames(x)))
  expect_identical(up(pm(va)), up(tau))                               # vlr()'s worker
  expect_true(max(abs(lo(pm(va)) - lo(outer(mu, mu, '-')))) < 1e-14)  # row over column
  expect_true(all(is.na(diag(pm(va)))))
  expect_true(rel(attr(va, 'clr_var'), apply(clr, 2, stats::var)) < 1e-12)
  expect_identical(names(attr(va, 'clr_var')), colnames(x))
  expect_true(abs(sum(attr(va, 'clr_var')) - attr(va, 'total_var')) < 1e-12 * attr(va, 'total_var'))
  expect_true(abs(attr(va, 'total_var') - sum(tau) / (2 * D)) < 1e-12 * attr(va, 'total_var'))
  expect_identical(attr(va, 'n'), n)
  expect_identical(attr(va, 'pseudocount'), 0)
  expect_identical(attr(va, 'method'), 'variation_array')
  expect_identical(attr(va, 'display')$title, 'Variation Array')
  expect_null(attr(va, 'order'))
  expect_null(attr(va, 'display')$clust)

  # The lower triangle is log(row / column): positive for a more abundant row.
  expect_true(pm(va)['F9', 'F1'] < 0)
  expect_true(abs(pm(va)['F9', 'F1'] - mean(log(x[, 'F9'] / x[, 'F1']))) < 1e-13)

  # Against compositions' variation() and mvar(), at 1e-12.
  for (nm in names(expected$va)) {
    e <- expected$va[[nm]]
    v <- variation_array(e$x)
    expect_true(rel(up(pm(v)), up(e$variation)) < 1e-12, info = nm)
    expect_true(rel(attr(v, 'total_var'), e$mvar) < 1e-12, info = nm)
  }

  # Unnamed features are labelled with their positions.
  v <- variation_array(unname(x))
  expect_identical(dimnames(v), list(as.character(1:D), as.character(1:D)))
  expect_identical(up(pm(v)), up(pm(va)))

  # Two features, two samples.
  v <- variation_array(x[1:2, 1:2])
  expect_identical(dim(v), c(2L, 2L))
  expect_true(abs(attr(v, 'total_var') - pm(v)[1, 2] / 2) < 1e-15)



  # Zeros: vlr()'s pseudocount rule and warning. ====

  expect_warning(v <- variation_array(xz), 'pseudocount')
  pc <- min(xz[xz > 0]) / 2
  expect_identical(attr(v, 'pseudocount'), pc)
  expect_identical(up(pm(v)), up(as.matrix(suppressWarnings(vlr(xz)))))
  Lz <- log(xz + pc)
  expect_true(max(abs(lo(pm(v)) - lo(outer(colMeans(Lz), colMeans(Lz), '-')))) < 1e-14)
  expect_true(rel(attr(v, 'clr_var'), apply(Lz - rowMeans(Lz), 2, stats::var)) < 1e-12)
  v2 <- variation_array(xz, pseudocount = 1)
  expect_identical(attr(v2, 'pseudocount'), 1)
  expect_identical(up(pm(v2)), up(as.matrix(vlr(xz, pseudocount = 1))))
  expect_error(variation_array(xz, pseudocount = 0), 'pseudocount')



  # Dense, sparse, margin = 2 and cpus agree exactly. ====

  vz <- suppressWarnings(variation_array(xz))
  expect_identical(variation_array(t(x), margin = 2L), va)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(xz, sparse = TRUE)
    expect_identical(suppressWarnings(variation_array(sp)), vz)
    expect_identical(suppressWarnings(variation_array(Matrix::t(sp), margin = 2L)), vz)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(suppressWarnings(variation_array(slam::as.simple_triplet_matrix(xz))), vz)
  for (cpus in c(1L, 2L, 4L))
    expect_identical(suppressWarnings(variation_array(xz, cpus = cpus)), vz, info = cpus)



  # clust: reorders by sqrt(tau). ====

  for (cl in c('ward.D2', 'average', 'pcoa')) {
    v <- variation_array(x, clust = cl)
    o <- attr(v, 'order')
    ref <- if (cl == 'pcoa') order(pcoa(sqrt(vlr(x)), n_axes = 1L)[, 1])
           else stats::hclust(stats::as.dist(sqrt(tau)), method = cl)$order
    expect_identical(o, ref, info = cl)
    expect_identical(attr(v, 'display')$clust, cl, info = cl)
    expect_identical(rownames(v), colnames(x)[o], info = cl)
    expect_identical(up(pm(v)), up(tau[o, o]), info = cl)
    expect_true(max(abs(lo(pm(v)) - lo(outer(mu[o], mu[o], '-')))) < 1e-14, info = cl)
    expect_identical(attr(v, 'clr_var'), attr(va, 'clr_var')[o], info = cl)
  }
  for (cl in list(FALSE, 'none')) {
    v <- variation_array(x, clust = cl)
    expect_identical(attr(v, 'display')$clust, 'none')
    expect_null(attr(v, 'order'))
    expect_identical(pm(v), pm(va))
  }
  expect_identical(
    suppressWarnings(variation_array(xz, clust = 'ward.D2', cpus = 1L)),
    suppressWarnings(variation_array(xz, clust = 'ward.D2', cpus = 4L)) )



  # Errors. ====

  expect_error(variation_array(x, clust = 'x'), 'clust')
  expect_error(variation_array(x[, 1, drop = FALSE]), 'at least 2 features')
  expect_error(variation_array(x[1, , drop = FALSE]), 'at least 2 samples')
  expect_error(variation_array(-x), 'negative')



  # Arithmetic and subsetting give plain matrices; print(). ====

  expect_identical(class(va * 2), c('matrix', 'array'))
  expect_null(attr(va * 2, 'clr_var'))
  expect_identical(class(sqrt(abs(va))), c('matrix', 'array'))
  expect_identical(class(va[1:2, 1:2]), c('matrix', 'array'))

  out <- utils::capture.output(res <- withVisible(print(va)))
  expect_false(res$visible)
  expect_identical(res$value, va)
  expect_true(grepl('Variation array of 9 features \\(30 samples\\)', out[[1]]))
  expect_false(any(grepl('more features', out)))
  big <- matrix(stats::rpois(10 * 14, 20) + 1L, 10, 14)
  out <- utils::capture.output(print(variation_array(big)))
  expect_true(any(grepl('# ... and 4 more features', out, fixed = TRUE)))



  # plot(). ====

  r <- drawn(va)
  expect_false(attr(r, 'visible'))
  attr(r, 'visible') <- NULL
  expect_identical(r, va)

  # The title from `display`, and the cells' values written for D <= 15.
  calls <- recorded(va)
  nm    <- vapply(calls, function (cl) if (is.null(cl$name)) '' else cl$name, '')
  ttl   <- calls[nm == 'C_title']
  expect_identical(ttl[[1]]$args[[1]], 'Variation Array')
  txt   <- unlist(lapply(calls[nm == 'C_text'], function (cl) cl$args[[2]]))
  expect_true(formatC(round(tau['F1', 'F2'], 2), format = 'f', digits = 2) %in% txt)
  expect_true(formatC(round(mu[['F2']] - mu[['F1']], 2), format = 'f', digits = 2) %in% txt)
  v <- va
  attr(v, 'display')$title <- 'My array'
  calls <- recorded(v)
  nm    <- vapply(calls, function (cl) if (is.null(cl$name)) '' else cl$name, '')
  expect_identical(calls[nm == 'C_title'][[1]]$args[[1]], 'My array')

  # More than 15 features: no values written.
  wide <- matrix(stats::rpois(12 * 16, 20) + 1L, 12, 16)
  calls <- recorded(variation_array(wide), digits = 3)
  nm    <- vapply(calls, function (cl) if (is.null(cl$name)) '' else cl$name, '')
  txt   <- unlist(lapply(calls[nm == 'C_text'], function (cl) cl$args[[2]]))
  expect_false(any(grepl('^-?[0-9]+\\.[0-9]{3}$', txt)))

  expect_false(attr(drawn(variation_array(x, clust = 'ward.D2'), key = FALSE, digits = 0), 'visible'))
  expect_false(attr(drawn(variation_array(x[, 1:2])), 'visible'))
  expect_error(drawn(va, digits = -1), 'digits')
  expect_error(drawn(va, key = 'x'), 'key')

  # The plot changes neither the RNG nor par().
  seed <- .Random.seed
  grDevices::pdf(NULL)
  mai <- graphics::par('mai')
  plot(va)
  expect_identical(graphics::par('mai'), mai)
  grDevices::dev.off()
  expect_identical(.Random.seed, seed)



  # zero_pattern(): a hand-built table. ====

  hand <- rbind(
    c(1, 0, 2, 3), c(0, 0, 1, 1), c(1, 0, 2, 3),
    c(5, 1, 0, 1), c(0, 0, 1, 1), c(1, 0, 4, 9) )
  dimnames(hand) <- list(paste0('s', 1:6), paste0('f', 1:4))
  zp <- zero_pattern(hand)
  expect_identical(class(zp), 'ecodive_zero_pattern')
  expect_identical(as.vector(unclass(zp)), c(1L, 2L, 1L, 3L, 2L, 1L))
  expect_identical(names(zp), rownames(hand))
  expect_identical(attr(zp, 'n_samples'), c(3L, 2L, 1L))
  expect_identical(
    attr(zp, 'patterns'),
    matrix(
      c(FALSE, TRUE, FALSE, FALSE,
        TRUE,  TRUE, FALSE, FALSE,
        FALSE, FALSE, TRUE, FALSE), 3, 4, byrow = TRUE,
      dimnames = list(c('1', '2', '3'), colnames(hand)) ))
  expect_identical(attr(zp, 'zeros'), c(f1 = 2, f2 = 5, f3 = 1, f4 = 0) / 6)
  expect_identical(attr(zp, 'n'), 6L)
  expect_identical(attr(zp, 'D'), 4L)

  # Ties are numbered by their first sample.
  zt <- zero_pattern(hand[c(4, 2, 1, 6, 5, 3), ][c(1, 2, 3, 5), ])
  # Samples s4, s2, s1, s5: s2 and s5 share pattern 1; s4 and s1 have one
  # sample each, numbered in sample order.
  expect_identical(as.vector(unclass(zt)), c(2L, 1L, 3L, 1L))

  # Sparse = dense (with stored zeros too), margin = 2.
  expect_identical(zero_pattern(t(hand), margin = 2L), zp)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(hand, sparse = TRUE)
    expect_identical(zero_pattern(sp), zp)
    sp2 <- sp
    sp2@x[1] <- 0                       # an explicitly stored zero
    h2 <- hand
    h2[sp@i[1] + 1L, 1] <- 0
    expect_identical(zero_pattern(sp2), zero_pattern(h2))
    expect_identical(zero_pattern(Matrix::t(sp), margin = 2L), zp)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(zero_pattern(slam::as.simple_triplet_matrix(hand)), zp)

  # No zeros: one pattern.
  z1 <- zero_pattern(x)
  expect_identical(as.vector(unclass(z1)), rep(1L, n))
  expect_identical(attr(z1, 'n_samples'), n)
  expect_false(any(attr(z1, 'patterns')))

  # Unnamed input.
  zu <- zero_pattern(unname(hand))
  expect_null(names(zu))
  expect_null(colnames(attr(zu, 'patterns')))
  expect_identical(as.vector(unclass(zu)), as.vector(unclass(zp)))

  # Against zCompositions' zPatterns(): the same partition of the samples,
  # and the same pattern sizes.
  for (nm in names(expected$zp)) {
    e <- expected$zp[[nm]]
    z <- zero_pattern(e$x)
    expect_identical(length(attr(z, 'n_samples')), length(unique(e$id)), info = nm)
    expect_identical(attr(z, 'n_samples'), as.integer(sort(tabulate(e$id), decreasing = TRUE)), info = nm)
    expect_identical(sum(table(e$id, unclass(z)) > 0), length(unique(e$id)), info = nm)
    # Each pattern's zeros, from its samples.
    p <- attr(z, 'patterns')
    first <- match(seq_len(nrow(p)), unclass(z))
    expect_identical(unname(p), unname(e$x[first, , drop = FALSE] == 0), info = nm)
  }

  # print()
  out <- utils::capture.output(res <- withVisible(print(zp)))
  expect_false(res$visible)
  expect_identical(out[[1]], '# 3 zero patterns among 6 samples and 4 features')
  expect_identical(out[[2]], '# 1 feature never zero, 0 always zero')
  expect_identical(out[[3]], "# most frequent pattern (1): 3 samples, zeros in 'f2'")

  # Errors.
  expect_error(zero_pattern(-hand), 'negative')
  expect_error(zero_pattern(hand, margin = 3), 'margin')



  # plot() of zero patterns. ====

  r <- drawn(zp)
  expect_false(attr(r, 'visible'))
  attr(r, 'visible') <- NULL
  expect_identical(r, zp)
  expect_false(attr(drawn(zp, n_patterns = 1, n_features = 1), 'visible'))
  expect_false(attr(drawn(z1), 'visible'))           # no feature varies
  expect_false(attr(drawn(zp, n_patterns = Inf, n_features = TRUE, main = 'Z', sub = 'S'), 'visible'))
  expect_error(drawn(zp, n_patterns = 0), 'n_patterns')
  expect_error(drawn(zp, n_features = -1), 'n_features')

  calls <- recorded(zp)
  nm    <- vapply(calls, function (cl) if (is.null(cl$name)) '' else cl$name, '')
  mt    <- unlist(lapply(calls[nm == 'C_mtext'], function (cl) cl$args[[1]]))
  expect_true('Zero patterns' %in% mt)
  expect_true('3 of 3 patterns shown (100% of samples); 1 feature hidden (1 never zero)' %in% mt)

  seed <- .Random.seed
  grDevices::pdf(NULL)
  plot(zp)
  expect_identical(graphics::par('mfrow'), c(1L, 1L))
  grDevices::dev.off()
  expect_identical(.Random.seed, seed)

#})
