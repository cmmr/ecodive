#test_that("mantel_regression", {

  ex <- readRDS('expected/mantel_regression.rds')
  x  <- ex$x
  n  <- attr(x, 'Size')

  geo   <- dist(ex$xy)
  dens  <- dist(ex$env$SubsDens)
  water <- dist(ex$env$WatrCont)
  topo  <- as.dist(outer(ex$env$Topo, ex$env$Topo, '!=') * 1)
  shrub <- as.dist(outer(ex$env$Shrub, ex$env$Shrub, '!=') * 1)
  terms <- list(geo = geo, dens = dens, water = water, topo = topo, shrub = shrub)



  # Same coefficients, R^2, F and p-values as ecodist::MRM(), given its
  # permutations: stored results, so at home only. ====

  if (at_home()) {
    for (cs in ex$cases) {

      info <- paste(cs$formula, cs$mrank)
      nms  <- all.vars(stats::as.formula(cs$formula))[-1]

      res <- do.call(mantel_regression, c(
        list(x), terms[nms],
        list(method = if (cs$mrank) 'spearman' else 'pearson',
             permutations = ex$perms, cpus = 1) ))
      cf <- res$coefficients

      expect_identical(cf$term, c('(Intercept)', nms), info = info)
      expect_equal(cf$estimate, cs$estimate, tolerance = 1e-12, info = info)
      # ecodist's R^2 has absolute rounding errors of about 1e-16.
      expect_true(abs(res$r_squared - cs$r_squared) < 1e-13, info = info)
      expect_equal(res$F, cs$F, tolerance = 1e-10, info = info)
      expect_identical(cf$p_value, cs$p_value, info = info)
      expect_identical(res$p_value, cs$F_p, info = info)
      expect_identical(res$p_value, cs$r2_p, info = info)
      expect_identical(res$n_perm, nrow(ex$perms), info = info)
      expect_identical(res$permutation, 'supplied matrix', info = info)
    }

    # The same terms from sample metadata: numbers' differences, a factor's
    # same (0) or different (1), coordinates' Euclidean distances.
    meta <- data.frame(ex$env, ex$xy)
    cs   <- ex$cases[[4]]  # x ~ water + topo + shrub
    res  <- mantel_regression(
      x ~ WatrCont + Topo + Shrub, data = meta, permutations = ex$perms, cpus = 1 )
    expect_equal(res$coefficients$estimate, cs$estimate, tolerance = 1e-12)
    expect_identical(res$coefficients$p_value, cs$p_value)
    expect_identical(res$coefficients$term, c('(Intercept)', 'WatrCont', 'Topo', 'Shrub'))

    cs  <- ex$cases[[5]]  # x ~ geo + water + topo, ranks
    res <- mantel_regression(
      x ~ cbind(x, y) + WatrCont + Topo, data = meta, method = 'spearman',
      permutations = ex$perms, cpus = 1 )
    expect_equal(res$coefficients$estimate, cs$estimate, tolerance = 1e-12)
    expect_identical(res$coefficients$p_value, cs$p_value)
    expect_identical(res$p_value, cs$F_p)
  }



  # Live references: lm() on the unfolded distances, and mantel(). ====

  set.seed(2)
  m    <- 25
  xy   <- matrix(runif(m * 2), m, 2, dimnames = list(paste0('S', 1:m), c('lat', 'lon')))
  meta <- data.frame(
    row.names = rownames(xy),
    ph   = round(runif(m, 5, 8), 1),
    temp = runif(m, 10, 20),
    soil = sample(c('clay', 'sand', 'loam'), m, replace = TRUE),
    wet  = sample(c(TRUE, FALSE), m, replace = TRUE),
    lat  = xy[, 1],
    lon  = xy[, 2],
    blk  = rep(1:5, 5),
    stringsAsFactors = FALSE )
  counts <- matrix(rpois(m * 15, 5), m, 15, dimnames = list(rownames(xy), NULL))
  counts[, 1:4] <- counts[, 1:4] + round(3 * meta$ph)
  dm <- bray(counts)

  vec  <- function (d) as.vector(d)
  dph  <- dist(stats::setNames(meta$ph, rownames(meta)))
  dtmp <- dist(stats::setNames(meta$temp, rownames(meta)))
  dgeo <- dist(xy)
  dsoil <- as.dist(outer(stats::setNames(meta$soil, rownames(meta)), meta$soil, '!=') * 1)

  for (meth in c('pearson', 'spearman')) {
    info <- meth
    tr <- if (meth == 'spearman') rank else identity
    fit <- summary(stats::lm(tr(vec(dm)) ~ tr(vec(dph)) + tr(vec(dtmp)) + tr(vec(dgeo)) + tr(vec(dsoil))))
    res <- mantel_regression(
      dm, ph = dph, temp = dtmp, geo = dgeo, soil = dsoil, method = meth,
      permutations = 99, cpus = 1 )
    cf <- res$coefficients

    expect_inherits(res, 'ecodive_mantel_regression', info = info)
    expect_equal(cf$estimate, unname(stats::coef(fit)[, 1]), tolerance = 1e-12, info = info)
    expect_equal(cf$t, unname(stats::coef(fit)[, 3]), tolerance = 1e-12, info = info)
    expect_equal(res$r_squared, fit$r.squared, tolerance = 1e-12, info = info)
    expect_equal(res$F, unname(fit$fstatistic[[1]]), tolerance = 1e-12, info = info)
    expect_identical(res$df, unname(c(fit$fstatistic[[2]], fit$fstatistic[[3]])), info = info)
    expect_identical(res$method, meth, info = info)
    expect_identical(res$n, as.integer(m), info = info)

    # Standardized coefficients: those of the standardized vectors.
    sc  <- function (v) (v - mean(v)) / stats::sd(v)
    fs  <- stats::lm(sc(tr(vec(dm))) ~ sc(tr(vec(dph))) + sc(tr(vec(dtmp))) + sc(tr(vec(dgeo))) + sc(tr(vec(dsoil))))
    expect_equal(cf$std_estimate[-1], unname(stats::coef(fs)[-1]), tolerance = 1e-12, info = info)
    expect_true(is.na(cf$std_estimate[[1]]), info = info)

    # Permuted statistics: the same lm() of the permuted response.
    perms <- t(replicate(5, sample(m)))
    res   <- mantel_regression(
      dm, ph = dph, temp = dtmp, geo = dgeo, soil = dsoil, method = meth,
      permutations = perms, cpus = 1 )
    dmm <- as.matrix(dm)
    for (k in 1:5) {
      yk  <- vec(stats::as.dist(dmm[perms[k, ], perms[k, ]]))
      fk  <- summary(stats::lm(tr(yk) ~ tr(vec(dph)) + tr(vec(dtmp)) + tr(vec(dgeo)) + tr(vec(dsoil))))
      expect_equal(unname(attr(res, 't_perm')[k, ]), unname(stats::coef(fk)[, 3]), tolerance = 1e-10, info = info)
      expect_equal(attr(res, 'f_perm')[[k]], unname(fk$fstatistic[[1]]), tolerance = 1e-10, info = info)
    }
    expect_identical(colnames(attr(res, 't_perm')), c('(Intercept)', 'ph', 'temp', 'geo', 'soil'), info = info)

    # One term: R^2 is mantel()'s r squared, and the test of F is its
    # two-sided test.
    r1 <- mantel_regression(dm, dph, method = meth, permutations = 199, cpus = 1)
    mt <- mantel(dm, dph, method = meth, alternative = 'two.sided', permutations = 199, cpus = 1)
    expect_equal(r1$r_squared, unname(mt$statistic)^2, tolerance = 1e-12, info = info)
    expect_identical(r1$p_value, mt$p_value, info = info)
    expect_equal(r1$coefficients$std_estimate[[2]], unname(mt$statistic), tolerance = 1e-12, info = info)
    expect_identical(r1$coefficients$term, c('(Intercept)', 'dph'), info = info)
  }



  # The ways of giving terms agree. ====

  ref <- mantel_regression(dm, ph = dph, temp = dtmp, soil = dsoil, geo = dgeo, permutations = 49)
  same <- function (res, info) {
    expect_identical(res$coefficients$term, ref$coefficients$term, info = info)
    expect_equal(res$coefficients[, -1], ref$coefficients[, -1], tolerance = 1e-13, info = info)
    expect_equal(res$F, ref$F, tolerance = 1e-13, info = info)
    expect_identical(res$p_value, ref$p_value, info = info)
  }

  geo_name <- function (res) { res$coefficients$term[[5]] <- 'geo'; res }
  same(geo_name(mantel_regression(dm ~ ph + temp + soil + cbind(lat, lon), data = meta, permutations = 49)), 'lhs formula')
  res <- mantel_regression(dm, data = meta, formula = ~ ph + temp + soil + cbind(lat, lon), permutations = 49)
  expect_identical(res$coefficients$term[[5]], 'cbind(lat, lon)')
  res$coefficients$term[[5]] <- 'geo'
  same(res, 'formula')
  res <- mantel_regression(dm, data = meta, formula = '~ ph + temp + soil', geo = xy, permutations = 49)
  same(res, 'character formula, coordinates in ...')
  same(mantel_regression(dm, ph = 'ph', temp = 'temp', soil = 'soil', geo = meta[, c('lat', 'lon')],
                         data = meta, permutations = 49), 'column names')
  same(mantel_regression(dm, ph = meta$ph, temp = meta$temp, soil = factor(meta$soil), geo = as.matrix(dgeo),
                         permutations = 49), 'vectors')
  same(mantel_regression(dm ~ ph + temp, soil = dsoil, geo = dgeo, data = meta, permutations = 49), 'formula and ...')

  # Named samples are matched: shuffled metadata, distances and vectors.
  sh <- sample(m)
  same(geo_name(mantel_regression(dm ~ ph + temp + soil + cbind(lat, lon), data = meta[sh, ], permutations = 49)), 'shuffled data')
  dsh <- stats::as.dist(as.matrix(dph)[sh, sh])
  same(mantel_regression(dm, ph = dsh, temp = stats::setNames(meta$temp, rownames(meta))[sh],
                         soil = dsoil, geo = xy[sh, ], permutations = 49), 'shuffled inputs')

  # A logical's distances are 0/1, as a factor's; an expression of a column.
  res <- mantel_regression(dm ~ wet + log(temp), data = meta, permutations = 0)
  dwet <- as.dist(outer(meta$wet, meta$wet, '!=') * 1)
  r2   <- mantel_regression(dm, wet = dwet, lt = dist(log(meta$temp)), permutations = 0)
  expect_equal(res$coefficients$estimate, r2$coefficients$estimate, tolerance = 1e-13)
  expect_identical(res$coefficients$term, c('(Intercept)', 'wet', 'log(temp)'))

  # Unnamed terms are named by what was written.
  res <- mantel_regression(dm, dph, dist(meta$temp), permutations = 0)
  expect_identical(res$coefficients$term, c('(Intercept)', 'dph', 'dist(meta$temp)'))

  # A matrix of distances.
  same(mantel_regression(as.matrix(dm), ph = dph, temp = dtmp, soil = dsoil, geo = dgeo, permutations = 49), 'matrix dm')



  # Permutations: seeds, strata, any number of cpus. ====

  a <- mantel_regression(dm, ph = dph, geo = dgeo, permutations = 99, seed = 1)
  b <- mantel_regression(dm, ph = dph, geo = dgeo, permutations = 99, seed = 1)
  d <- mantel_regression(dm, ph = dph, geo = dgeo, permutations = 99, seed = 2)
  expect_identical(a, b)
  expect_false(identical(attr(a, 'f_perm'), attr(d, 'f_perm')))
  expect_identical(a$permutation, 'free')
  expect_identical(a$n_perm, 99L)
  expect_identical(length(attr(a, 'f_perm')), 99L)
  expect_identical(dim(attr(a, 't_perm')), c(99L, 3L))
  expect_true(all(a$coefficients$p_value >= 1 / 100 & a$coefficients$p_value <= 1))
  expect_true(a$coefficients$p_value[[2]] <= 0.05)  # pH shifts the community

  # Within strata: the same as those permutations given as a matrix.
  s  <- mantel_regression(dm, ph = dph, geo = dgeo, data = meta, strata = 'blk', permutations = 99, seed = 3)
  ps <- t(.Call(ecodive:::C_perm_indices, m, 99L, as.integer(meta$blk) - 1L, 3L))
  expect_true(all(apply(ps, 1, function (p) all(meta$blk[p] == meta$blk))))
  sp <- mantel_regression(dm, ph = dph, geo = dgeo, permutations = ps)
  expect_identical(s$permutation, 'within strata')
  expect_identical(attr(s, 'f_perm'), attr(sp, 'f_perm'))
  expect_identical(attr(s, 't_perm'), attr(sp, 't_perm'))
  expect_identical(s$coefficients, sp$coefficients)
  s2 <- mantel_regression(dm, ph = dph, geo = dgeo, strata = stats::setNames(meta$blk, rownames(meta))[sh],
                          permutations = 99, seed = 3)
  expect_identical(s2$coefficients, s$coefficients)

  big <- mantel_regression(dm, ph = dph, temp = dtmp, geo = dgeo, method = 'spearman', permutations = 999, cpus = 1)
  for (cpus in if (at_home()) c(2, 3, 4, 8) else 2) {
    res <- mantel_regression(dm, ph = dph, temp = dtmp, geo = dgeo, method = 'spearman', permutations = 999, cpus = cpus)
    expect_identical(res, big, info = paste('cpus', cpus))
  }

  # No permutations: no p-values.
  res <- mantel_regression(dm, ph = dph, permutations = 0)
  expect_true(all(is.na(res$coefficients$p_value)) && is.na(res$p_value))
  expect_identical(res$n_perm, 0L)
  expect_identical(length(attr(res, 'f_perm')), 0L)
  expect_identical(dim(attr(res, 't_perm')), c(0L, 2L))



  # Small and tied inputs. ====

  # n = 3: 3 pairs, one term.
  d3  <- dist(c(a = 1, b = 2, c = 4))
  res <- mantel_regression(d3, z = c(a = 0, b = 1, c = 1), permutations = 5)
  expect_equal(res$r_squared, stats::cor(vec(d3), c(1, 1, 0))^2, tolerance = 1e-12)
  expect_identical(res$df, c(1, 1))
  expect_true(all(res$coefficients$p_value > 0 & res$coefficients$p_value <= 1))

  # Duplicate samples and heavily tied distances.
  dt  <- dist(rep(1:3, 4))
  res <- mantel_regression(dt, g = rep(c('a', 'b'), 6), h = rep(1:4, 3), method = 'spearman', permutations = 19)
  fit <- stats::lm(rank(vec(dt)) ~ rank(vec(dist(rep(c(1, 2), 6)))) + rank(vec(dist(rep(1:4, 3)))))
  expect_equal(res$coefficients$estimate, unname(stats::coef(fit)), tolerance = 1e-12)

  # A response explained exactly.
  res <- mantel_regression(dph, ph = meta$ph, permutations = 9)
  expect_equal(res$coefficients$estimate, c(0, 1), tolerance = 1e-12)
  expect_equal(res$r_squared, 1, tolerance = 1e-14)
  expect_identical(res$p_value, 0.1)



  # Results: title, print, as.data.frame, plot. ====

  res <- mantel_regression(dm, ph = dph, temp = dtmp, geo = dgeo, permutations = 99)
  expect_identical(attr(res, 'display')$title, 'Multiple Regression on Distance Matrices of Bray-Curtis Dissimilarity')
  expect_identical(
    attr(mantel_regression(dist(meta$temp), ph = dph, permutations = 9), 'display')$title,
    'Multiple Regression on Distance Matrices' )

  out <- capture.output(res)
  expect_identical(out[[1]], 'Multiple Regression on Distance Matrices of Bray-Curtis Dissimilarity')
  expect_true(any(grepl('^R2 = .*, F = .*, p = ', out)))
  expect_true(any(grepl('Permutations: 99 (free)', out, fixed = TRUE)))
  expect_true(any(grepl('Samples: 25 (300 pairs)', out, fixed = TRUE)))
  expect_true(any(grepl('(Intercept)', out, fixed = TRUE)))
  expect_true(any(grepl('geo', out, fixed = TRUE)))
  expect_true(any(grepl('ranks', capture.output(
    mantel_regression(dm, ph = dph, method = 'spearman', permutations = 0)))))

  expect_identical(as.data.frame(res), res$coefficients)
  expect_identical(names(res$coefficients), c('term', 'estimate', 'std_estimate', 't', 'p_value'))
  expect_identical(names(res), c('coefficients', 'r_squared', 'F', 'p_value', 'df', 'n_perm', 'method', 'n', 'permutation'))

  drawn <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  grDevices::pdf(NULL)
  expect_true(drawn(res))
  expect_true(drawn(res, alpha = 0.5, main = 'Mine', col = 'blue', sub = 'S'))
  expect_true(drawn(mantel_regression(dm, ph = dph, method = 'spearman', permutations = 0)))
  expect_error(plot(res, alpha = 2), '`alpha` must be a single number')
  expect_error(plot(structure(list(), class = 'ecodive_mantel_regression')), 'mantel_regression\\(\\) result')
  grDevices::dev.off()

  # plot() takes its title from the result, and R^2 for its subtitle.
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
  tt <- titles(res)
  expect_true('Multiple Regression on Distance Matrices of Bray-Curtis Dissimilarity' %in% tt)
  expect_true(any(grepl('^R2 = .*, p = ', tt)))
  expect_true(all(c('ph', 'temp', 'geo') %in% tt))
  expect_false('(Intercept)' %in% tt)
  attr(res, 'display')$title <- 'Custom'
  expect_true('Custom' %in% titles(res))
  expect_true('Mine' %in% titles(res, main = 'Mine'))



  # Errors. ====

  expect_error(mantel_regression(dm), 'at least one explanatory term')
  expect_error(mantel_regression(dm, ~ ph, data = meta), 'Name the formula argument')
  expect_error(mantel_regression(dm, data = ~ ph), 'Name the formula argument')
  expect_error(mantel_regression(dm, formula = ~ ph), '`data` is required')
  expect_error(mantel_regression(dm ~ ph, ph = dph, formula = ~ temp, data = meta), 'not both')
  expect_error(mantel_regression(dm ~ ph * temp, data = meta), 'interactions')
  expect_error(mantel_regression(dm ~ ph + Condition(temp), data = meta), 'Condition')
  expect_error(mantel_regression(dm ~ nope, data = meta), 'Not found in `data`')
  expect_error(mantel_regression(dm, formula = 'ph', data = meta), 'must be a formula')
  expect_error(mantel_regression(dm, ph = 'nope', data = meta), 'not a column of `data`')
  expect_error(mantel_regression(dm, ph = 'ph'), 'names a column, but `data` is not given')
  expect_error(mantel_regression(dm, ph = dph, ph = dtmp), 'Two terms are named')
  expect_error(mantel_regression(dm, a = dph, b = dph), 'Collinear terms')
  expect_error(mantel_regression(dm, a = dph, b = dist(rep(1, m))), "Collinear terms: 'b'")
  expect_error(mantel_regression(dm, a = dph, b = 2 * dph + dtmp, c = dtmp), 'Collinear terms')
  expect_error(mantel_regression(dist(rep(1, m)), a = dph), 'nothing to explain')
  expect_error(mantel_regression(d3, a = 1:3, b = c(1, 5, 2)), 'too few pairs')
  expect_error(mantel_regression(dist(1:2), a = 1:2), 'at least 3 samples')
  expect_error(mantel_regression(dm, a = 1:3), 'has 3 values')
  expect_error(mantel_regression(dm, a = list(1)), 'must be distances')
  expect_error(mantel_regression(dm, a = dist(1:4)), 'has 4 samples')
  expect_error(mantel_regression(dm, a = c(NA, 1:24)), 'NA')
  expect_error(mantel_regression(dm, a = c(NA, letters[1:24])), 'NAs')
  expect_error(mantel_regression(dm, a = matrix('a', m, 2)), 'numeric to be coordinates')
  expect_error(mantel_regression(dm, a = dph, method = 'kendall'), "'arg' should be one of")
  expect_error(mantel_regression(dm, a = dph, permutations = -1), '`permutations` must be')
  expect_error(mantel_regression(dm, a = dph, permutations = matrix(1L, 2, m)), 'not a permutation')
  expect_error(mantel_regression(dm, a = dph, data = meta, strata = 'blk', permutations = ps), 'not both')
  expect_error(mantel_regression(dm, a = dph, data = 1:3), 'must be a data.frame')
  expect_error(mantel_regression(dm, a = dph, seed = 'a'))
  expect_error(mantel_regression(dm, a = dph, cpus = 0))
  expect_error(mantel_regression(dm ~ log(nope2), data = meta), 'Not found in `data`')


  # Backward elimination. ====

  set.seed(5)
  noise <- lapply(1:3, function (i) dist(stats::setNames(rnorm(m), rownames(meta))))
  names(noise) <- paste0('noise', 1:3)
  all_terms <- c(list(ph = dph, temp = dtmp, soil = dsoil, geo = dgeo), noise)

  bw <- do.call(mantel_regression, c(list(dm), all_terms, list(permutations = 99, backward = TRUE, seed = 7)))
  el <- bw$elimination

  expect_identical(names(bw), c(
    'coefficients', 'r_squared', 'F', 'p_value', 'df', 'n_perm', 'method', 'n',
    'permutation', 'alpha', 'elimination' ))
  expect_identical(bw$alpha, 0.05)
  expect_identical(names(el), c('step', 'term', 't', 'p_value', 'r_squared_before', 'r_squared_after'))
  expect_true(nrow(el) >= 1)
  expect_identical(el$step, seq_len(nrow(el)))
  expect_true(all(bw$coefficients$p_value[-1] <= 0.05))
  expect_true('ph' %in% bw$coefficients$term)
  expect_identical(sort(c(bw$coefficients$term[-1], el$term)), sort(names(all_terms)))

  # Each step: refit the terms left, with the same permutations; the term
  # removed has the largest p-value (above alpha), and among ties the
  # smallest |t|.
  left <- names(all_terms)
  for (k in seq_len(nrow(el))) {
    fk <- do.call(mantel_regression, c(list(dm), all_terms[left], list(permutations = 99, seed = 7)))
    cf <- fk$coefficients[-1, ]
    info <- paste('step', k)
    expect_identical(el$term[[k]], cf$term[order(-cf$p_value, abs(cf$t))[[1]]], info = info)
    expect_true(el$p_value[[k]] > 0.05, info = info)
    expect_identical(el$p_value[[k]], cf$p_value[cf$term == el$term[[k]]], info = info)
    expect_identical(el$t[[k]], cf$t[cf$term == el$term[[k]]], info = info)
    expect_identical(el$r_squared_before[[k]], fk$r_squared, info = info)
    left <- setdiff(left, el$term[[k]])
  }
  expect_identical(el$r_squared_after[-nrow(el)], el$r_squared_before[-1])

  # The final model is the plain fit of the terms left.
  fin <- do.call(mantel_regression, c(list(dm), all_terms[left], list(permutations = 99, seed = 7)))
  expect_identical(bw$coefficients, fin$coefficients)
  expect_identical(bw$r_squared, fin$r_squared)
  expect_identical(utils::tail(el$r_squared_after, 1), fin$r_squared)
  expect_identical(attr(bw, 'f_perm'), attr(fin, 'f_perm'))
  expect_identical(attr(bw, 't_perm'), attr(fin, 't_perm'))
  expect_identical(bw[1:9], fin[1:9])

  # A supplied matrix and strata are used at every step; any cpus.
  pmat <- t(replicate(99, sample(m)))
  b1 <- do.call(mantel_regression, c(list(dm), all_terms, list(permutations = pmat, backward = TRUE, alpha = 0.2, cpus = 1)))
  for (cpus in if (at_home()) c(2, 3, 4, 8) else 2) {
    b2 <- do.call(mantel_regression, c(list(dm), all_terms, list(permutations = pmat, backward = TRUE, alpha = 0.2, cpus = cpus)))
    expect_identical(b2, b1, info = paste('cpus', cpus))
  }
  expect_identical(b1$alpha, 0.2)
  expect_true(all(b1$coefficients$p_value[-1] <= 0.2))
  bs <- mantel_regression(dm, ph = dph, temp = dtmp, soil = dsoil, data = meta, strata = 'blk',
                          permutations = 99, backward = TRUE, seed = 3)
  expect_identical(bs$permutation, 'within strata')

  # Nothing to remove: the full model, and an empty table.
  b0 <- mantel_regression(dm, ph = dph, permutations = 99, backward = TRUE)
  expect_identical(b0[1:9], mantel_regression(dm, ph = dph, permutations = 99)[1:9])
  expect_identical(nrow(b0$elimination), 0L)
  expect_identical(names(b0$elimination), names(el))
  expect_true(any(grepl('no term removed', capture.output(b0), fixed = TRUE)))

  # Everything removed: the intercept alone, the mean distance.
  ba <- mantel_regression(dm, n1 = noise$noise1, n2 = noise$noise2, permutations = 99,
                          backward = TRUE, alpha = 0)
  cf2 <- mantel_regression(dm, n1 = noise$noise1, n2 = noise$noise2, permutations = 99)$coefficients[-1, ]
  expect_identical(ba$elimination$term, cf2$term[order(-cf2$p_value, abs(cf2$t))])
  expect_identical(ba$coefficients$term, '(Intercept)')
  expect_equal(ba$coefficients$estimate, mean(vec(dm)), tolerance = 1e-14)
  fit0 <- summary(stats::lm(vec(dm) ~ 1))
  expect_equal(ba$coefficients$t, unname(stats::coef(fit0)[1, 3]), tolerance = 1e-12)
  expect_identical(ba$coefficients$p_value, 1)
  expect_identical(ba$r_squared, 0)
  expect_true(is.na(ba$F) && is.na(ba$p_value))
  expect_identical(ba$df, c(0, m * (m - 1) / 2 - 1))
  expect_identical(utils::tail(ba$elimination$r_squared_after, 1), 0)
  expect_identical(dim(attr(ba, 't_perm')), c(99L, 1L))
  expect_identical(length(attr(ba, 'f_perm')), 99L)
  expect_true(any(grepl('R2 = 0, F = NA, p = NA', capture.output(ba), fixed = TRUE)))
  grDevices::pdf(NULL)
  expect_error(plot(ba), 'no terms to plot')
  expect_true(drawn(bw))
  grDevices::dev.off()

  out <- capture.output(bw)
  expect_true(any(grepl(sprintf('Backward elimination (alpha = 0.05): %d term', nrow(el)), out, fixed = TRUE)))
  expect_true(any(grepl('r_squared_before', out, fixed = TRUE)))

  # Without backward elimination, no such fields.
  expect_false(any(c('alpha', 'elimination') %in% names(fin)))

  expect_error(mantel_regression(dm, ph = dph, backward = 'yes'), '`backward` must be TRUE or FALSE')
  expect_error(mantel_regression(dm, ph = dph, backward = NA), '`backward` must be TRUE or FALSE')
  expect_error(mantel_regression(dm, ph = dph, alpha = 2), '`alpha` must be a single number')
  expect_error(mantel_regression(dm, ph = dph, alpha = NA), '`alpha` must be a single number')
  expect_error(mantel_regression(dm, ph = dph, backward = TRUE, permutations = 0), 'needs permutations')

#})
