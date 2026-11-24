#test_that("batch correction", {

  # batch_mmuphin() and batch_combat(): against MMUPHin's adjust_batch()
  # (MIT) and sva's ComBat() (Artistic-2.0), stored by
  # data-raw/expected/batch_mmuphin.r and data-raw/expected/batch_combat.r;
  # invariants; threads and table classes; arguments; printing and
  # plotting.

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))
  exit_if_not(requireNamespace("slam",   quietly = TRUE))

  mm <- readRDS(file.path('expected', 'batch_mmuphin.rds'))
  cb <- readRDS(file.path('expected', 'batch_combat.rds'))

  # The values only, as a plain matrix.
  plain <- function (x) { m <- unclass(x); attributes(m) <- attributes(x)[c('dim', 'dimnames')]; m }

  # The largest relative difference, where both are numbers.
  rel_diff <- function (a, b) {
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) return (0)
    max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
  }
  effects <- function (x, col) matrix(attr(x, 'batch_effects')[[col]], ncol(plain(x)))


  #--------------------------------------------------
  # MMUPHin: the table (counts: identical after
  # rounding; proportions: 1e-10), and the estimates
  # (1e-10, with NA where MMUPHin has NA).
  #--------------------------------------------------
  for (nm in names(mm$cases)) {
    cs <- mm$cases[[nm]]
    r  <- mm$results[[nm]]
    x  <- mm$tables[[cs$table]]
    md <- mm$metas[[cs$table]]
    f  <- if (is.null(cs$cov)) NULL else stats::reformulate(cs$cov)
    res <- batch_mmuphin(x, mm$batches[[cs$table]], f, md, zero_inflation = cs$zi, pseudocount = cs$pc)

    if (cs$table == 'sim') expect_identical(plain(res), r$adjusted, info = nm)
    else                   expect_true(rel_diff(plain(res), r$adjusted) < 1e-10, info = nm)

    for (k in c('gamma_hat', 'delta_hat', 'gamma_star', 'delta_star')) {
      expect_identical(is.na(effects(res, k)), is.na(r[[k]]), info = paste(nm, k))
      expect_true(rel_diff(effects(res, k), r[[k]]) < 1e-10, info = paste(nm, k))
    }
    expect_identical(unname(attr(res, 'adjusted')), rowSums(r$ind_gamma) > 0, info = nm)
    if (is.null(cs$pc)) expect_equal(attr(res, 'pseudocount'), r$pseudocount, info = nm)
    else                expect_null(attr(res, 'pseudocount'))

    # Unrounded: the same counts scaled down (not whole numbers, so not
    # rounded), scaled back up.
    if (!is.null(r$unrounded)) {
      pc <- if (is.null(cs$pc)) r$pseudocount else cs$pc
      res2 <- batch_mmuphin(x / 1024, mm$batches[[cs$table]], f, md, zero_inflation = cs$zi, pseudocount = pc / 1024)
      expect_true(rel_diff(plain(res2) * 1024, r$unrounded) < 1e-10, info = nm)
    }
  }


  #--------------------------------------------------
  # ComBat: every case at 1e-10.
  #--------------------------------------------------
  for (nm in names(cb$cases)) {
    cs  <- cb$cases[[nm]]
    md  <- cb$metas[[cs$table]]
    f   <- if (is.null(cs$cov)) NULL else stats::reformulate(cs$cov)
    res <- if (nm == 'clr')
             batch_combat(cb$counts, 'batch', f, md, norm = 'clr', pseudocount = 1)
           else if (nm == 'single')
             expect_warning(
               batch_combat(cb$tables[[cs$table]], 'batch', f, md, parametric = cs$par),
               'single sample' )
           else
             batch_combat(
               cb$tables[[cs$table]], 'batch', f, md, parametric = cs$par,
               mean_only = isTRUE(cs$mo), ref_batch = cs$ref )
    if (nm == 'single') res <- suppressWarnings(batch_combat(cb$tables[[cs$table]], 'batch', f, md))
    expect_true(rel_diff(plain(res), cb$results[[nm]]) < 1e-10, info = nm)
    expect_identical(attr(res, 'mean_only'), isTRUE(cs$mo) || nm == 'single', info = nm)
  }

  # Features constant within a batch are unchanged, and marked.
  lx  <- cb$tables$logx
  res <- batch_combat(lx, cb$metas$logx$batch)
  expect_identical(unname(attr(res, 'adjusted')), c(rep(TRUE, 23), FALSE, FALSE))
  expect_identical(plain(res)[, 24:25], lx[, 24:25])
  expect_true(all(is.na(attr(res, 'batch_effects')$gamma_star[attr(res, 'batch_effects')$feature %in% c('F24', 'F25')])))

  # The reference batch is unchanged, with gamma_star 0 and delta_star 1.
  res <- batch_combat(lx, cb$metas$logx$batch, ref_batch = 'b2')
  b2  <- cb$metas$logx$batch == 'b2'
  expect_identical(plain(res)[b2, ], lx[b2, ])
  be  <- attr(res, 'batch_effects')
  expect_true(all(be$gamma_star[be$batch == 'b2' & !is.na(be$gamma_star)] == 0))
  expect_true(all(be$delta_star[be$batch == 'b2' & !is.na(be$delta_star)] == 1))
  expect_identical(names(attr(res, 'n_iters')), c('b1', 'b2', 'b3'))
  expect_null(attr(batch_combat(lx, cb$metas$logx$batch, parametric = FALSE), 'n_iters'))

  # The non-parametric priors with large batches, where sva's products of
  # densities underflow to 0 (and give NaN): finite.
  set.seed(1)
  big <- matrix(stats::rnorm(1600 * 12), 1600, 12)
  bb  <- rep(c('a', 'b'), each = 800)
  big[bb == 'b', ] <- big[bb == 'b', ] + 1
  res <- batch_combat(big, bb, parametric = FALSE, cpus = 2L)
  expect_true(all(is.finite(plain(res))))
  expect_true(abs(mean(plain(res)[bb == 'a', ]) - mean(plain(res)[bb == 'b', ])) < 0.1)


  #--------------------------------------------------
  # Invariants.
  #--------------------------------------------------
  x  <- mm$tables$sim
  md <- mm$metas$sim
  b  <- md$batch
  base <- batch_mmuphin(x, b)

  # Zeros stay zero; totals are kept (to rounding).
  expect_identical(plain(base) == 0, x == 0)
  expect_true(max(abs(rowSums(plain(base)) - rowSums(x))) <= ncol(x))

  # Permuting samples permutes the result; so does reordering features.
  p <- c(48:25, 1:24)
  expect_identical(plain(batch_mmuphin(x[p, ], b[p])), plain(base)[p, ])
  q <- c(30:1)
  expect_identical(plain(batch_mmuphin(x[, q], b)), plain(base)[, q])
  lxb <- cb$metas$logx$batch
  r0  <- batch_combat(lx, lxb)
  pp  <- c(45:1)
  expect_equal(plain(batch_combat(lx[pp, ], lxb[pp])), plain(r0)[pp, ], tolerance = 1e-12)

  # A feature present in one batch only is not adjusted (zero inflation):
  # F28 is only in b1.
  expect_false(attr(base, 'adjusted')[['F28']])
  be <- attr(base, 'batch_effects')
  expect_identical(is.na(be$gamma_hat[be$feature == 'F27']), c(FALSE, FALSE, TRUE))

  # A sample of all zeros stays all zeros.
  x0 <- x; x0[5, ] <- 0
  expect_true(all(plain(batch_mmuphin(x0, b))[5, ] == 0))

  # Batch's R^2 (of Bray-Curtis on relative abundances): with no batch
  # effect, low before and after; a planted effect is removed.
  r2 <- function (m, g) permanova(bray(normalize(m, 'percent')), g, permutations = 0)$r2
  set.seed(2)
  nb <- matrix(stats::rpois(40 * 30, 30), 40, 30)
  b2 <- rep(c('x', 'y'), 20)
  expect_true(r2(nb, b2) < 0.05)
  expect_true(r2(batch_mmuphin(nb, b2), b2) < 0.05)
  eb <- nb
  eb[b2 == 'y', 1:10] <- eb[b2 == 'y', 1:10] * 3
  expect_true(r2(eb, b2) > 0.5)
  expect_true(r2(batch_mmuphin(eb, b2), b2) < 0.05)
  expect_true(r2(base, b) < r2(x, b) / 2)
  r2e <- function (m, g) permanova(euclidean(m), g, permutations = 0)$r2
  expect_true(r2e(normalize(eb, 'clr', pseudocount = 1), b2) > 0.5)
  expect_true(r2e(batch_combat(eb, b2, norm = 'clr', pseudocount = 1), b2) < 0.05)

  # The result passes through ecodive's functions as a count table.
  expect_identical(bray(base), bray(plain(base)))
  expect_identical(
    unclass(pcoa(bray(base), n_axes = 2)),
    unclass(pcoa(bray(plain(base)), n_axes = 2)) )
  expect_identical(
    permanova(bray(base), b, seed = 1)$p_value,
    permanova(bray(plain(base)), b, seed = 1)$p_value )
  expect_identical(
    da_test(base, "da_wilcoxon", groups = b == "b1", pseudocount = 0.5)$p_value,
    da_test(plain(base), "da_wilcoxon", groups = b == "b1", pseudocount = 0.5)$p_value )
  expect_identical(normalize(base, 'percent'), normalize(plain(base), 'percent'))

  # `[` drops the class and attributes.
  expect_identical(class(base[1:3, ]), c('matrix', 'array'))
  expect_null(attr(base[1:3, ], 'batch_effects'))


  #--------------------------------------------------
  # Threads, table classes, orientation.
  #--------------------------------------------------
  if (packageVersion("Matrix") >= "1.5") {
    COMPRESSED <- "CsparseMatrix"; TRIPLET <- "TsparseMatrix"
  } else {
    COMPRESSED <- "dgCMatrix"; TRIPLET <- "dgTMatrix"
  }

  # Large enough for threads.
  set.seed(3)
  nbig <- 150; dbig <- 600
  bigb <- rep(c('p', 'q', 'r'), each = 50)
  mu   <- matrix(stats::rexp(dbig, 1 / 20), nbig, dbig, byrow = TRUE)
  mu[bigb == 'q', 1:200] <- mu[bigb == 'q', 1:200] * 2
  X    <- matrix(stats::rnbinom(nbig * dbig, mu = mu, size = 0.8), nbig, dbig)
  meta <- data.frame(g = rep(c('u', 'v'), 75), bat = bigb)

  for (zi in c(TRUE, FALSE)) {
    r1 <- batch_mmuphin(X, 'bat', ~ g, meta, zero_inflation = zi, cpus = 1L)
    for (cpus in c(2L, 4L))
      expect_identical(batch_mmuphin(X, 'bat', ~ g, meta, zero_inflation = zi, cpus = cpus), r1, info = cpus)
    expect_identical(batch_mmuphin(as(as(X, 'CsparseMatrix'), COMPRESSED), 'bat', ~ g, meta, zero_inflation = zi, cpus = 3L), r1)
    expect_identical(batch_mmuphin(as(as(X, 'CsparseMatrix'), TRIPLET), 'bat', ~ g, meta, zero_inflation = zi), r1)
    expect_identical(batch_mmuphin(slam::as.simple_triplet_matrix(X), 'bat', ~ g, meta, zero_inflation = zi, cpus = 2L), r1)
    rt <- batch_mmuphin(t(X), 'bat', ~ g, meta, margin = 2L, zero_inflation = zi, cpus = 2L)
    expect_identical(plain(rt), t(plain(r1)))
    expect_identical(attr(rt, 'batch_effects'), attr(r1, 'batch_effects'))
  }

  LX <- log(X + 1)
  for (par in c(TRUE, FALSE)) {
    r1 <- batch_combat(LX, bigb, parametric = par, cpus = 1L)
    for (cpus in c(2L, 4L))
      expect_identical(batch_combat(LX, bigb, parametric = par, cpus = cpus), r1, info = cpus)
    expect_identical(plain(batch_combat(t(LX), bigb, margin = 2L, parametric = par)), t(plain(r1)))
  }
  r1 <- batch_combat(X, bigb, norm = 'log', pseudocount = 1, cpus = 1L)
  expect_identical(batch_combat(as(as(X, 'CsparseMatrix'), COMPRESSED), bigb, norm = 'log', pseudocount = 1, cpus = 4L), r1)
  expect_identical(batch_combat(slam::as.simple_triplet_matrix(X), bigb, norm = 'log', pseudocount = 1, cpus = 2L), r1)
  expect_equal(
    plain(r1),
    plain(batch_combat(normalize(X, 'log', pseudocount = 1), bigb)),
    tolerance = 1e-14 )


  #--------------------------------------------------
  # Per-sample arguments, and their checks.
  #--------------------------------------------------
  # `batch` as values, a column name, or named values in another order.
  r_col <- batch_mmuphin(x, 'batch', data = md)
  expect_identical(r_col, base)
  nb_ <- stats::setNames(as.character(b), rownames(x))
  expect_identical(plain(batch_mmuphin(x, rev(nb_))), plain(base))
  expect_identical(plain(batch_mmuphin(x, 'batch', data = md[rev(rownames(md)), ])), plain(base))

  # A character formula; a two-sided one's response is ignored.
  r_f <- batch_mmuphin(x, 'batch', ~ cond + age, md)
  expect_identical(batch_mmuphin(x, 'batch', 'cond + age', md), r_f)
  expect_identical(batch_mmuphin(x, 'batch', y ~ cond + age, md), r_f)

  expect_error(batch_mmuphin(x), '`batch` is required')
  expect_error(batch_mmuphin(x, rep('a', 48)), 'at least two batches')
  expect_error(batch_mmuphin(x, c(rep('a', 47), 'b')), 'at least 2 samples')
  expect_error(batch_mmuphin(x, c(NA, as.character(b)[-1])), 'missing values')
  expect_error(batch_mmuphin(x, b[-1]), 'has 47 values')
  expect_error(batch_mmuphin(x, b, ~ cond), '`formula` needs `data`')
  expect_error(batch_mmuphin(x, b, ~ nope, md), "Not found in `data`: 'nope'")
  expect_error(batch_mmuphin(x, 'batch', ~ batch, md), 'Covariates and batch are confounded')
  md2 <- md; md2$cond2 <- md2$cond
  expect_error(batch_mmuphin(x, 'batch', ~ cond + cond2, md2), 'covariates in `formula` are confounded')
  md2$one <- 'k'
  expect_error(batch_mmuphin(x, 'batch', ~ one, md2), 'single value')
  md2$age[3] <- NA
  expect_error(batch_mmuphin(x, 'batch', ~ age, md2), 'Missing values')

  x_na <- x; x_na[2, 2] <- NA
  expect_error(batch_mmuphin(x_na, b), 'NA')
  expect_error(batch_combat(x_na, b), 'NA')
  x_neg <- x; x_neg[2, 2] <- -1
  expect_error(batch_mmuphin(x_neg, b), 'negative')
  expect_error(batch_combat(x_neg, b, norm = 'clr'), 'negative')
  expect_silent(batch_combat(log2(x + 1) - 3, b))
  expect_error(batch_mmuphin(x * 0, b), 'only zeros')

  expect_error(batch_mmuphin(x, b, zero_inflation = NA), '`zero_inflation`')
  expect_error(batch_mmuphin(x, b, tol = 0), '`tol`')
  expect_error(batch_mmuphin(x, b, max_iters = 1.5), '`max_iters`')
  expect_error(batch_mmuphin(x, b, pseudocount = -1), '`pseudocount`')
  expect_error(batch_mmuphin(x, b, zero_inflation = FALSE, pseudocount = 0), 'needs `zero_inflation = TRUE`')
  r_pc0 <- batch_mmuphin(x, b, pseudocount = 0)
  expect_identical(plain(r_pc0) == 0, x == 0)
  expect_error(batch_mmuphin(x, b, max_iters = 1), 'did not converge in 1 iterations')

  # A batch can't be adjusted: present in fewer than two adjustable
  # features.
  xb <- x; xb[b == 'b3', -1] <- 0
  expect_error(batch_mmuphin(xb, b), "can't be adjusted: .*'b3'")
  # No feature can be adjusted.
  xs <- matrix(0, 6, 3); xs[1:2, 1] <- 1; xs[3:4, 2] <- 1; xs[5:6, 3] <- 1
  expect_error(batch_mmuphin(xs, rep(1:3, each = 2)), 'No feature can be adjusted')

  expect_error(batch_combat(lx, lxb, norm = 'percent'), '`norm`')
  expect_error(batch_combat(lx, lxb, parametric = 'yes'), '`parametric`')
  expect_error(batch_combat(lx, lxb, mean_only = NA), '`mean_only`')
  expect_error(batch_combat(lx, lxb, ref_batch = 'zz'), '`ref_batch` must be one of')
  expect_error(batch_combat(lx, lxb, max_iters = 1), 'did not converge')
  expect_error(batch_combat(lx[, 24:25], lxb), 'Fewer than two features')
  expect_identical(
    batch_combat(lx, as.integer(lxb), ref_batch = 2),
    batch_combat(lx, as.integer(lxb), ref_batch = '2') )
  expect_warning(batch_combat(cb$counts, lxb, norm = 'clr'), 'pseudocount')
  r_clr <- suppressWarnings(batch_combat(cb$counts, lxb, norm = 'clr'))
  expect_equal(attr(r_clr, 'pseudocount'), 0.5)

  # Batch levels: a factor's order, else sorted.
  fb <- factor(as.character(b), levels = c('b3', 'b1', 'b2'))
  expect_identical(levels(attr(batch_mmuphin(x, fb), 'batch_effects')$batch), c('b3', 'b1', 'b2'))
  expect_identical(plain(batch_mmuphin(x, fb)), plain(base))


  #--------------------------------------------------
  # Attributes, printing, plotting.
  #--------------------------------------------------
  expect_identical(class(base), c('ecodive_batch', 'matrix', 'array'))
  expect_identical(attr(base, 'method'), 'mmuphin')
  expect_identical(
    names(attributes(base)),
    c('dim', 'dimnames', 'method', 'batch_effects', 'adjusted', 'pseudocount',
      'n_iters', 'batch_means', 'batch_means_adjusted', 'overall_means', 'class') )
  expect_identical(dim(attr(base, 'batch_means')), c(30L, 3L))
  expect_equal(
    unname(attr(base, 'batch_means')[, 'b2']),
    unname(colMeans((x / rowSums(x))[b == 'b2', ])), tolerance = 1e-14 )
  expect_equal(
    unname(attr(base, 'batch_means_adjusted')[, 'b1']),
    unname(colMeans((plain(base) / rowSums(plain(base)))[b == 'b1', ])), tolerance = 1e-14 )
  expect_equal(unname(attr(base, 'overall_means')), unname(colMeans(x / rowSums(x))), tolerance = 1e-14)
  expect_identical(attr(r0, 'method'), 'combat')
  expect_equal(
    unname(attr(r0, 'batch_means_adjusted')[, 'b3']),
    unname(colMeans(plain(r0)[lxb == 'b3', ])), tolerance = 1e-14 )

  out <- utils::capture.output(print(base))
  expect_true(any(grepl('48 samples', out)))
  expect_true(any(grepl('Batch-corrected (MMUPHin Batch Correction)', out, fixed = TRUE)))
  expect_true(any(grepl('features adjusted', out)))
  out <- utils::capture.output(print(rt))
  expect_true(any(grepl('600 features', out)))
  out <- utils::capture.output(print(r0))
  expect_true(any(grepl('(ComBat)', out, fixed = TRUE)))

  titles <- function (x, ...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control('enable')
    plot(x, ...)
    txt <- unlist(lapply(grDevices::recordPlot()[[1]], function (op) {
      a <- unlist(op[[2]], recursive = TRUE)
      a[vapply(a, is.character, logical(1))]
    }))
    gsub('\n', ' ', unique(unlist(txt)))
  }
  expect_true('MMUPHin Batch Correction' %in% titles(base))
  expect_true('MMUPHin Batch Correction' %in% titles(base, type = 'means'))
  expect_true('Batch mean relative abundance' %in% titles(base, type = 'means'))
  expect_true('ComBat' %in% titles(r0))
  expect_true('Batch mean value' %in% titles(r0, type = 'means'))
  expect_true('x' %in% titles(r0, main = 'x'))
  expect_true('Locations only (mean_only)' %in% titles(batch_combat(lx, lxb, mean_only = TRUE)))

  grDevices::pdf(NULL)
  expect_silent(plot(base, type = 'means', log = 'xy'))
  expect_silent(plot(base, legend = NULL))
  expect_error(plot(base, type = 'other'), '`type`')
  lost <- base
  attr(lost, 'batch_means') <- NULL
  expect_error(plot(lost), 'lost the attributes')
  grDevices::dev.off()

#})
