#test_that("ancom", {

  # ancom(): scikit-bio 0.7.4's ancom(), SciPy's f_oneway() and kruskal()
  # statistics and p-values, and SciPy's F and chi-squared tail
  # probabilities (data-raw/expected/ancom.py); base R's oneway.test(),
  # kruskal.test() and p.adjust(); threading, input types and argument
  # checks.

  ref <- readRDS(file.path('expected', 'ancom.rds'))

  details <- function (...) ancom_details(..., cpus = 1L)
  pc_of   <- function (t) if (t$pc > 0) t$pc
  upper   <- function (m) t(m)[lower.tri(m)]   # packed, row by row (i < j)
  maxrel  <- function (a, b) {
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), 1e-300)))
  }


  # -----------------------------------------------------------------------
  # scikit-bio: W and the decisions, exactly
  # -----------------------------------------------------------------------

  for (k in seq_len(nrow(ref$runs))) {
    r    <- ref$runs[k, ]
    t    <- ref$tables[[r$table]]
    info <- paste(r$table, r$sig_test, r$p_adjust, r$alpha, r$tau, r$theta)
    res  <- ancom(
      t$counts, t$groups, pseudocount = pc_of(t), alpha = r$alpha,
      tau = r$tau, theta = r$theta, p_adjust = r$p_adjust,
      sig_test = r$sig_test, cpus = 1L )
    expect_identical(res$W,      r$W[[1]],      info = info)
    expect_identical(res$reject, r$signif[[1]], info = info)
  }

  # scikit-bio's own test values.
  t <- ref$tables$skbio1
  expect_identical(ancom(t$counts, t$groups, p_adjust = 'none')$W, c(5L, 5L, 2L, 2L, 2L, 2L, 2L))
  expect_identical(ancom(t$counts, t$groups, p_adjust = 'none')$reject, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  t <- ref$tables$doc
  res <- ancom(t$counts, t$groups)
  expect_identical(res$W, c(0L, 4L, 0L, 1L, 1L, 0L, 1L))
  expect_identical(res$reject, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(attr(res, 'groups'), c('placebo', 'treatment'))


  # -----------------------------------------------------------------------
  # SciPy: each pair's statistic and p-value; scikit-bio's adjustment
  # -----------------------------------------------------------------------

  for (nm in names(ref$tables)) {
    t <- ref$tables[[nm]]
    for (st in c('anova', 'kruskal')) {
      info <- paste(nm, st)
      s    <- t[[st]]
      d    <- details(t$counts, t$groups, pseudocount = pc_of(t), sig_test = st, p_adjust = 'none')

      expect_identical(is.nan(d$statistic), is.nan(s$statistic), info = info)
      expect_identical(d$statistic[is.infinite(d$statistic)], s$statistic[is.infinite(s$statistic)], info = info)
      if (st == 'kruskal') {
        # The same operations as SciPy, so the same bits.
        expect_identical(d$statistic, s$statistic, info = info)
      } else {
        # SciPy's sums of squares cancel (to 3e-12); ecodive's two-pass ones
        # are within 2e-13 of the exact F of the same log ratios.
        expect_true(maxrel(d$statistic, s$statistic) < 1e-11, info = info)
      }

      p <- upper(d$p_value)
      expect_identical(is.nan(p), is.nan(s$p_value), info = info)
      expect_true(max(abs(p - s$p_value), na.rm = TRUE) < 1e-13, info = info)
      expect_equal(p, s$p_value, tolerance = 1e-10, info = info)

      for (a in names(s$adjusted)) {
        e <- unname(details(t$counts, t$groups, pseudocount = pc_of(t), sig_test = st, p_adjust = a)$p_adjusted)
        m <- s$adjusted[[a]]
        diag(e) <- NA
        expect_identical(is.na(e), is.na(m), info = paste(info, a))
        expect_true(max(abs(e - m), na.rm = TRUE) < 1e-13, info = paste(info, a))
      }
    }
  }

  # pf() and pchisq() against SciPy's fdtrc() and chdtrc().
  dd <- ref$dist
  mine <- numeric(nrow(dd))
  f    <- dd$dist == 'f'
  mine[f]  <- pf(dd$x[f], dd$df1[f], dd$df2[f], lower.tail = FALSE)
  mine[!f] <- pchisq(dd$x[!f], dd$df1[!f], lower.tail = FALSE)
  expect_true(max(abs(mine - dd$p)) < 1e-13)
  expect_equal(mine, dd$p, tolerance = 1e-12)


  # -----------------------------------------------------------------------
  # Base R: the tests and the adjustment
  # -----------------------------------------------------------------------

  set.seed(3)
  x   <- matrix(rlnorm(25 * 8), 25, 8)
  x[, 2] <- x[, 1] * 3                    # a constant log ratio: NaN
  grp <- rep(c('a', 'b', 'c'), c(8, 9, 8))
  x[grp == 'b', 3] <- x[grp == 'b', 3] * 5
  L   <- log(x)
  D   <- ncol(x)

  for (st in c('anova', 'kruskal')) {
    d <- details(x, grp, sig_test = st, p_adjust = 'none')
    for (i in 1:(D - 1)) for (j in (i + 1):D) {
      info <- paste(st, i, j)
      lr   <- L[, i] - L[, j]
      p    <- d$p_value[i, j]
      if (i == 1 && j == 2) next   # constant but for rounding in the logs
      bt <- if (st == 'anova') stats::oneway.test(lr ~ grp, var.equal = TRUE)
            else               stats::kruskal.test(lr, factor(grp))
      expect_equal(p, bt$p.value, tolerance = 1e-10, info = info)
      expect_identical(p, d$p_value[j, i], info = info)
    }
    expect_true(all(is.nan(diag(d$p_value))))

    # Each feature's family is c(0, its p-values), as p.adjust() adjusts it.
    for (a in c('holm', 'bonferroni', 'BH', 'BY', 'none')) {
      e <- details(x, grp, sig_test = st, p_adjust = a)
      for (i in seq_len(D)) {
        want <- stats::p.adjust(c(0, unname(d$p_value[i, -i])), a)[-1]
        expect_identical(unname(e$p_adjusted[i, -i]), want, info = paste(st, a, i))
      }
      expect_identical(unname(diag(e$p_adjusted)), rep(1, D))
      expect_identical(e$result$W, as.integer(rowSums(e$p_adjusted < 0.05, na.rm = TRUE) - (diag(e$p_adjusted) < 0.05)))
    }
  }

  # A log ratio the same in every sample has no p-value (x2 = 3 x1 isn't
  # exactly constant, by rounding in the logs; identical features are).
  y <- cbind(x, x[, 4])
  for (st in c('anova', 'kruskal')) {
    d <- details(y, grp, sig_test = st)
    expect_true(is.nan(d$p_value[4, 9]), info = st)
    expect_true(is.nan(d$statistic[3 * 9 - 6 + 5]), info = st)   # pair (4, 9)
    expect_identical(sum(is.nan(d$statistic)), 1L, info = st)
  }

  # Constant within every group, different between groups: F = Inf, p = 0.
  y <- cbind(f1 = rep(c(10, 20), c(4, 5)), f2 = 5, f3 = c(3, 9, 4, 8, 2, 7, 6, 5, 1))
  g <- rep(c('lo', 'hi'), c(4, 5))
  d <- details(y, g)
  expect_identical(d$statistic[1], Inf)
  expect_identical(d$p_value[1, 2], 0)
  expect_true(is.finite(d$statistic[2]) && is.finite(d$statistic[3]))


  # -----------------------------------------------------------------------
  # The cutoff
  # -----------------------------------------------------------------------

  # Nothing above theta: nothing rejected, cutoff NA.
  t   <- ref$tables$skbio3
  res <- ancom(t$counts, t$groups)
  expect_true(max(res$W) / nrow(res) < 0.1)
  expect_identical(res$reject, logical(nrow(res)))
  expect_identical(attr(res, 'cutoff'), NA_real_)

  # scikit-bio's rule on crafted W, including W on the boundary where
  # numpy.linspace(0.05, 0.25, 5)'s third value, one rounding step above
  # 0.15, decides.
  cf <- ref$cutoff
  expect_true(sum(cf$quirk) >= 1L)
  for (k in seq_len(nrow(cf))) {
    cut <- ancom_cutoff(cf$W[[k]], cf$D[k], cf$tau[k], cf$theta[k])
    expect_identical(cut$reject, cf$signif[[k]], info = k)
  }

  t   <- ref$tables$skbio1
  res <- ancom(t$counts, t$groups, p_adjust = 'none')
  expect_identical(attr(res, 'cutoff'), 5 / 7 - 0.1)
  expect_true(0x1.3333333333334p-3 != 0.15)
  expect_identical(attr(ancom(t$counts, t$groups, p_adjust = 'none', tau = 0.5), 'cutoff'), 5 / 7 - 0.1)


  # -----------------------------------------------------------------------
  # Invariance, threads, input types, orientation
  # -----------------------------------------------------------------------

  t <- ref$tables$three
  r0 <- details(t$counts, t$groups)
  y  <- t$counts * c(10, rep(1, 22), 1e-3)          # rescale two samples
  r1 <- details(y, t$groups)
  expect_identical(r1$result$W, r0$result$W)
  expect_identical(r1$result$reject, r0$result$reject)
  expect_equal(r1$p_value, r0$p_value, tolerance = 1e-10)

  set.seed(2)
  big <- matrix(rpois(60 * 150, 2), 60, 150,
    dimnames = list(paste0('s', 1:60), paste0('f', 1:150)) )
  bg  <- rep(c('x', 'y', 'z'), 20)
  big[bg == 'y', 1:10] <- big[bg == 'y', 1:10] * 4L
  sp  <- Matrix::Matrix(big, sparse = TRUE)
  sl  <- slam::as.simple_triplet_matrix(big)
  for (st in c('anova', 'kruskal')) {
    d1 <- ancom_details(big, bg, pseudocount = 0.5, sig_test = st, cpus = 1L)
    for (cpus in c(2L, 3L, 4L, 8L))
      expect_identical(ancom_details(big, bg, pseudocount = 0.5, sig_test = st, cpus = cpus), d1, info = paste(st, cpus))
    expect_identical(ancom_details(sp, bg, pseudocount = 0.5, sig_test = st, cpus = 4L), d1, info = st)
    expect_identical(ancom_details(sl, bg, pseudocount = 0.5, sig_test = st, cpus = 2L), d1, info = st)
    expect_identical(ancom_details(t(big), bg, margin = 2L, pseudocount = 0.5, sig_test = st, cpus = 2L), d1, info = st)
    expect_identical(ancom_details(Matrix::t(sp), bg, margin = 2L, pseudocount = 0.5, sig_test = st, cpus = 4L), d1, info = st)
    expect_true(any(d1$result$reject), info = st)
    expect_identical(
      ancom(big, bg, pseudocount = 0.5, sig_test = st, cpus = 2L),
      d1$result, info = st )
  }

  # Integer, double and data-frame tables.
  t <- ref$tables$doc
  expect_identical(ancom(t$counts, t$groups), ancom(t$counts * 1.0, t$groups))
  expect_identical(ancom(as.data.frame(t$counts), t$groups)$W, ancom(t$counts, t$groups)$W)


  # -----------------------------------------------------------------------
  # Zeros
  # -----------------------------------------------------------------------

  t <- ref$tables$zeros
  expect_warning(res <- ancom(t$counts, t$groups), 'Zeros detected')
  expect_identical(attr(res, 'pseudocount'), t$pc)
  expect_identical(res, ancom(t$counts, t$groups, pseudocount = t$pc))
  expect_error(ancom(t$counts, t$groups, pseudocount = 0), 'positive')
  expect_error(ancom(t$counts, t$groups, pseudocount = -1), 'positive')

  # replace_zeros() output goes straight in, with no pseudocount.
  z <- replace_zeros(t$counts, 'mult')
  expect_silent(res <- ancom(z, t$groups))
  expect_identical(attr(res, 'pseudocount'), 0)
  expect_identical(res$W, ancom(z, t$groups, pseudocount = 1e-300)$W)


  # -----------------------------------------------------------------------
  # Groups and the result's structure
  # -----------------------------------------------------------------------

  t   <- ref$tables$three
  x   <- t$counts
  dimnames(x) <- list(paste0('S', seq_len(nrow(x))), paste0('F', seq_len(ncol(x))))
  res <- ancom(x, t$groups, cpus = 1L)
  expect_inherits(res, 'ecodive_ancom')
  expect_inherits(res, 'data.frame')
  expect_identical(names(res), c('feature', 'W', 'reject'))
  expect_identical(res$feature, colnames(x))
  expect_true(is.integer(res$W) && is.logical(res$reject))
  expect_identical(attr(res, 'groups'), c('ctl', 'hi', 'lo'))   # sorted
  expect_identical(attr(res, 'sig_test'), 'anova')
  expect_identical(attr(res, 'p_adjust'), 'holm')
  expect_identical(attr(res, 'alpha'), 0.05)
  expect_identical(attr(res, 'tau'), 0.02)
  expect_identical(attr(res, 'theta'), 0.1)
  expect_identical(attr(res, 'pseudocount'), 0)
  expect_true(is.numeric(attr(res, 'cutoff')))

  # A factor keeps its levels' order; the results don't change.
  fg  <- factor(t$groups, levels = c('lo', 'hi', 'ctl', 'unused'))
  res2 <- ancom(x, fg, cpus = 1L)
  expect_identical(attr(res2, 'groups'), c('lo', 'hi', 'ctl'))
  expect_identical(res2$W, res$W)

  # Radix (C-locale) order, as numpy.unique.
  g2  <- ifelse(t$groups == 'ctl', 'B', ifelse(t$groups == 'hi', 'a', 'C'))
  expect_identical(attr(ancom(x, g2, cpus = 1L), 'groups'), c('B', 'C', 'a'))
  expect_identical(attr(ancom(x, match(t$groups, c('lo', 'hi', 'ctl')) * 10, cpus = 1L), 'groups'), c('10', '20', '30'))

  # `data`, and names matching samples.
  meta <- data.frame(row.names = rev(rownames(x)), grp = rev(t$groups))
  expect_identical(ancom(x, 'grp', meta, cpus = 1L)$W, res$W)
  expect_identical(ancom(x, setNames(rev(t$groups), rev(rownames(x))), cpus = 1L)$W, res$W)

  # Unnamed features.
  expect_identical(ancom(unname(x), t$groups, cpus = 1L)$feature, as.character(1:15))

  # Argument matching.
  expect_identical(attr(ancom(x, t$groups, sig_test = 'f_oneway'), 'sig_test'), 'anova')
  expect_identical(attr(ancom(x, t$groups, sig_test = 'kru'), 'sig_test'), 'kruskal')
  expect_identical(attr(ancom(x, t$groups, p_adjust = 'fdr'), 'p_adjust'), 'BH')
  expect_identical(attr(ancom(x, t$groups, p_adjust = 'bonf'), 'p_adjust'), 'bonferroni')
  expect_identical(ancom(x, t$groups, p_adjust = 'fdr'), ancom(x, t$groups, p_adjust = 'BH'))


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  expect_error(ancom(x), '`groups` is required')
  expect_error(ancom(x, rep('a', 24)), 'at least two different')
  expect_error(ancom(x, as.character(1:24)), 'group of its own')
  g <- t$groups; g[3] <- NA
  expect_error(ancom(x, g), 'NA')
  expect_error(ancom(x, t$groups[-1]), 'groups')
  for (arg in c('alpha', 'tau', 'theta'))
    for (v in list(0, 1, -0.1, 2, NA, c(0.1, 0.2), 'a'))
      expect_error(do.call(ancom, setNames(list(x, t$groups, v), c('counts', 'groups', arg))),
                   'strictly between 0 and 1', info = paste(arg, v))
  expect_error(ancom(x, t$groups, sig_test = 'ttest'), 'sig_test')
  expect_error(ancom(x, t$groups, sig_test = NA), 'sig_test')
  expect_error(ancom(x, t$groups, p_adjust = 'hochberg'), 'p_adjust')
  expect_error(ancom(x, t$groups, p_adjust = 'hommel'), 'p_adjust')
  expect_error(ancom(x, t$groups, p_adjust = 'nope'), 'p_adjust')
  expect_error(ancom(x, t$groups, p_adjust = NULL), 'p_adjust')
  expect_error(ancom(x[, 1, drop = FALSE], t$groups), 'at least 2 features')
  expect_error(ancom(x, t$groups, margin = 3), 'margin')
  expect_error(ancom(x, t$groups, cpus = 0), 'cpus')
  xn <- x; xn[1, 1] <- NA
  expect_error(ancom(xn, t$groups), 'NA')
  expect_error(ancom(-x, t$groups), 'negative')

#})
