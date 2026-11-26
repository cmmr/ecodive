#test_that("principal_response_curves", {

  # Inputs and reference results, written by
  # data-raw/expected/principal_response_curves.r.
  ref <- readRDS(file.path('expected', 'principal_response_curves.rds'))
  dat <- ref$data

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  rel <- function (a, b) max(abs(a - b), na.rm = TRUE) / max(abs(b), na.rm = TRUE)

  fit_case <- function (cs, scaling = 'symmetric', permutations = 0L, cpus = 1L, ...) {
    rows <- if (is.null(cs$rows)) seq_len(nrow(dat$pyrifos)) else cs$rows
    y    <- dat$pyrifos[rows, ]
    args <- list(cpus = cpus, ...)
    if (cs$transform == 'hellinger') args$norm <- 'hellinger'
    if (cs$transform == 'counts') {
      y <- round((exp(y) - 1) / 10)
      args$norm <- 'clr'
      args$pseudocount <- 1
    }
    do.call(principal_response_curves, c(
      list(y, dat$dose[rows], dat$week[rows], control = cs$control, n_axes = 2L,
           scaling = scaling, permutations = permutations),
      args ))
  }

  # A treatments x times matrix of one axis's coefficients, without the
  # control, as vegan's summary(prc)$coefficients.
  coef_matrix <- function (x, axis) {
    cf <- x$coefficients
    cf <- cf[cf$treatment != x$control, , drop = FALSE]
    matrix(cf[[paste0("PRC", axis)]], ncol = nlevels(cf$time), byrow = TRUE)
  }



  # Same as vegan::prc(): stored results, so at home only. ====

  if (at_home()) {

    scalings <- c(sites = 1, species = 2, symmetric = 3)

    for (cs in ref$cases) {

      ex  <- cs$expected
      lbl <- cs$name

      for (sc in names(scalings)) {

        pr <- if (sc == 'symmetric') fit_case(cs, sc) else fit_case(cs, scalings[[sc]])
        expect_identical(pr$scaling, sc, info = lbl)

        for (a in 1:2) {
          info <- sprintf('%s, %s, axis %d', lbl, sc, a)
          cf <- coef_matrix(pr, a)
          vc <- ex$coef[[sc]][[a]]
          sg <- sign(sum(cf * vc, na.rm = TRUE))
          expect_identical(abs(sg), 1, info = info)
          expect_identical(is.na(cf), unname(is.na(vc)), info = info)
          expect_true(rel(sg * cf, vc) < 1e-10, info = info)
          expect_true(rel(sg * unname(pr$weights[, a]), unname(ex$sp[[sc]][[a]])) < 1e-10, info = info)
          expect_identical(names(ex$sp[[sc]][[a]]), rownames(pr$weights), info = info)
        }
        expect_equal(unname(pr$eig), ex$eig, tolerance = 1e-10, info = lbl)
      }

      # Signs: each axis's largest coefficient is positive.
      for (a in 1:2) {
        v <- pr$coefficients[[a + 2L]]
        expect_true(v[[which.max(abs(replace(v, is.na(v), 0)))]] > 0, info = lbl)
      }

      # The tests, with vegan's permutations.
      if (!is.null(cs$permutations)) {
        pr <- fit_case(cs, permutations = dat[[cs$permutations]])
        for (t in c('first', 'all')) {
          e   <- ex[[paste0('test_', t)]]
          row <- if (t == 'first') 'PRC1' else 'all'
          col <- if (t == 'first') 'PRC1' else 'all'
          expect_equal(pr$test[row, 'F'], e$F, tolerance = 1e-10, info = lbl)
          expect_identical(pr$test[row, 'p_value'], e$p, info = lbl)
          expect_identical(c(pr$test[row, 'df'], pr$df_residual), as.integer(e$df), info = lbl)
          expect_equal(pr$residual, e$residual, tolerance = 1e-10, info = lbl)
          expect_true(rel(pr$F_perm[, col], e$F_perm) < 1e-10, info = lbl)
        }
        expect_identical(pr$n_perm, 99L, info = lbl)
        expect_identical(pr$permutation, 'supplied matrix', info = lbl)
      }
    }


    # The overall test is PERMANOVA of the Euclidean distances, with the
    # indicators of the non-control cells as constraints (vegan's coding).
    # ====

    P  <- dat$perm_free
    pr <- fit_case(ref$cases[[1]], permutations = P)
    tc <- data.frame(time = dat$week, treatment = stats::relevel(dat$dose, '0'))
    d  <- data.frame(time = dat$week)
    d$X <- stats::model.matrix(~ time + time:treatment, tc)[, -(1:11)]
    pm <- permanova(euclidean(dat$pyrifos), data = d, formula = ~ X + Condition(time), permutations = P)
    expect_equal(pr$test['all', 'F'], unname(pm$statistic), tolerance = 1e-10)
    expect_identical(pr$test['all', 'p_value'], pm$p_value)
    expect_true(rel(pr$F_perm[, 'all'], as.vector(attr(pm, 'f_perm'))) < 1e-10)
  }



  # The underlying partial RDA. ====

  cs <- ref$cases[[1]]
  pr <- fit_case(cs)
  d  <- data.frame(treatment = stats::relevel(dat$dose, '0'), time = dat$week)
  o  <- suppressMessages(rda(dat$pyrifos, data = d, formula = ~ treatment:time + Condition(time), n_axes = 2))
  fl <- sign(colSums(unclass(pr$rda) * unclass(o)))
  expect_identical(abs(fl), c(RDA1 = 1, RDA2 = 1))
  expect_inherits(pr$rda, 'ecodive_ord')
  expect_equal(unclass(pr$rda) * rep(fl, each = nrow(o)), unclass(o), tolerance = 1e-12, check.attributes = FALSE)
  expect_equal(attr(pr$rda, 'loadings') * rep(fl, each = ncol(dat$pyrifos)), attr(o, 'loadings'), tolerance = 1e-12)
  expect_equal(attr(pr$rda, 'biplot') * rep(fl, each = nrow(attr(o, 'biplot'))), attr(o, 'biplot'), tolerance = 1e-12)
  expect_equal(attr(pr$rda, 'centroids') * rep(fl, each = nrow(attr(o, 'centroids'))), attr(o, 'centroids'), tolerance = 1e-12)
  for (nm in c('eig', 'prop', 'trace', 'inertia', 'constrained_eig', 'r_squared', 'adj_r_squared', 'method'))
    expect_identical(attr(pr$rda, nm), attr(o, nm), info = nm)
  expect_identical(attr(pr$rda, 'display')$title, attr(o, 'display')$title)
  expect_identical(unname(pr$eig), attr(o, 'eig')[1:2])
  expect_identical(unname(pr$prop), attr(o, 'prop')[1:2])

  # The coefficients are the cells' mean site scores, less the control's.
  w   <- unclass(pr$rda)[, 1] * sqrt(sqrt(pr$eig[[1]] / attr(pr$rda, 'trace')))
  mn  <- tapply(w, list(dat$dose, dat$week), mean)
  cf  <- sweep(mn, 2, mn['0', ])
  got <- pr$coefficients
  expect_equal(got$PRC1, cf[cbind(as.character(got$treatment), as.character(got$time))], tolerance = 1e-12)

  # The product of a coefficient and a weight is the same in every scaling.
  p1 <- fit_case(cs, 1)
  p2 <- fit_case(cs, 2)
  expect_equal(outer(p1$coefficients$PRC1, p1$weights[, 1]), outer(pr$coefficients$PRC1, pr$weights[, 1]), tolerance = 1e-12)
  expect_equal(outer(p2$coefficients$PRC2, p2$weights[, 2]), outer(pr$coefficients$PRC2, pr$weights[, 2]), tolerance = 1e-12)
  expect_identical(fit_case(cs, 'sites'), p1)
  expect_identical(fit_case(cs, 'form'), p1)
  expect_identical(fit_case(cs, 'covariance'), p2)
  expect_identical(fit_case(cs, 3), pr)



  # Layout of the result. ====

  expect_inherits(pr, 'ecodive_prc')
  expect_identical(
    names(pr),
    c('coefficients', 'weights', 'eig', 'prop', 'test', 'df_residual', 'residual',
      'F_perm', 'n_perm', 'permutation', 'scaling', 'control', 'rda') )
  cf <- pr$coefficients
  expect_identical(names(cf), c('treatment', 'time', 'PRC1', 'PRC2'))
  expect_identical(nrow(cf), 55L)
  expect_identical(levels(cf$treatment), levels(dat$dose))
  expect_identical(levels(cf$time), levels(dat$week))
  expect_identical(as.character(cf$treatment[1:11]), rep('0', 11))
  expect_identical(cf$time[1:11], factor(levels(dat$week), levels(dat$week)))
  expect_identical(as.character(unique(cf$treatment)), c('0', '0.1', '0.9', '6', '44'))
  expect_identical(cf$PRC1[1:11], rep(0, 11))
  expect_identical(dimnames(pr$weights), list(colnames(dat$pyrifos), c('PRC1', 'PRC2')))
  expect_identical(names(pr$eig), c('PRC1', 'PRC2'))
  expect_identical(pr$control, '0')
  expect_identical(pr$n_perm, 0L)
  expect_identical(pr$permutation, 'free')
  expect_identical(rownames(pr$test), c('PRC1', 'all'))
  expect_identical(names(pr$test), c('df', 'variance', 'F', 'z', 'p_value'))
  expect_identical(pr$test$df, c(1L, 44L))
  expect_identical(pr$df_residual, 77L)
  expect_true(all(is.na(pr$test$p_value)) && all(is.na(pr$test$z)))
  expect_true(all(is.finite(pr$test$F)))
  expect_identical(dim(pr$F_perm), c(0L, 2L))
  expect_equal(pr$test$variance, c(pr$eig[[1]], attr(pr$rda, 'inertia')[['constrained']]), tolerance = 1e-12, check.attributes = FALSE)
  expect_equal(pr$residual, attr(pr$rda, 'inertia')[['unconstrained']], tolerance = 1e-12)

  # One axis by default.
  p1 <- principal_response_curves(dat$pyrifos, dat$dose, dat$week, permutations = 0)
  expect_identical(names(p1$coefficients), c('treatment', 'time', 'PRC1'))
  expect_identical(p1$coefficients$PRC1, pr$coefficients$PRC1)
  expect_identical(p1$weights[, 1], pr$weights[, 1])
  expect_identical(ncol(p1$rda), 1L)



  # The same results for any cpus, and random permutations. ====

  sub <- seq_len(60)
  run <- function (cpus, permutations = 49L, ...)
    principal_response_curves(
      dat$pyrifos[sub, ], dat$dose[sub], dat$week[sub], n_axes = 2,
      permutations = permutations, cpus = cpus, ... )
  r1 <- run(1L)
  for (cp in c(2L, max_cpus))
    expect_identical(run(cp), r1, info = cp)
  expect_identical(r1$permutation, 'free')
  expect_identical(r1$n_perm, 49L)
  expect_identical(dim(r1$F_perm), c(49L, 2L))
  expect_true(all(r1$test$p_value > 0 & r1$test$p_value <= 1))
  expect_identical(run(1L, seed = 1L)$coefficients, r1$coefficients)
  expect_false(identical(run(1L, seed = 1L)$F_perm, r1$F_perm))
  expect_identical(run(1L, seed = 0L), r1)

  # Within strata: as a supplied matrix of the same permutations.
  rs <- run(1L, strata = dat$week[sub])
  expect_identical(rs$permutation, 'within strata')
  ps <- .Call(ecodive:::C_perm_indices, length(sub), 49L, as.integer(dat$week[sub]) - 1L, 0L)
  rm <- run(1L, permutations = t(ps))
  expect_identical(rm$F_perm, rs$F_perm)
  expect_identical(rm$test, rs$test)

  # Whole plots: each permutation moves whole series of samples.
  rp <- run(1L, plots = dat$ditch[sub])
  expect_identical(rp$permutation, 'whole plots')
  pp <- ecodive:::prc_plot_perms(dat$ditch[sub], droplevels(dat$week[sub]), 49L, 0L)
  for (k in 1:49) {
    p <- pp$perms[, k] + 1L
    expect_identical(sort(p), seq_along(sub))
    expect_identical(dat$week[sub][p], dat$week[sub])
    # Samples of one plot go to samples of one plot.
    expect_true(all(tapply(as.integer(dat$ditch[sub][p]), dat$ditch[sub], function (v) length(unique(v))) == 1))
  }
  expect_true(any(pp$perms + 1L != seq_along(sub)))
  rm <- run(1L, permutations = t(pp$perms + 1L))
  expect_identical(rm$F_perm, rp$F_perm)
  expect_identical(rm$test, rp$test)
  for (cp in c(2L, max_cpus))
    expect_identical(run(cp, plots = dat$ditch[sub]), rp, info = cp)



  # Inputs: columns of `data`, names, margin, sparse tables. ====

  meta <- data.frame(dose = dat$dose, week = dat$week, ditch = dat$ditch,
                     row.names = rownames(dat$pyrifos))
  shuf <- meta[rev(seq_len(nrow(meta))), ]
  a <- principal_response_curves(dat$pyrifos, 'dose', 'week', shuf, n_axes = 2, permutations = 0)
  expect_identical(a, pr)
  a <- principal_response_curves(dat$pyrifos, 'dose', 'week', shuf, n_axes = 2, permutations = 19, plots = 'ditch')
  b <- principal_response_curves(dat$pyrifos, dat$dose, dat$week, n_axes = 2, permutations = 19, plots = dat$ditch)
  expect_identical(a, b)
  a <- principal_response_curves(
    t(dat$pyrifos), stats::setNames(dat$dose, rownames(dat$pyrifos)), dat$week,
    margin = 2L, n_axes = 2, permutations = 0 )
  expect_identical(a, pr)
  a <- principal_response_curves(as.data.frame(dat$pyrifos), dat$dose, dat$week, n_axes = 2, permutations = 0)
  expect_identical(a, pr)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(dat$pyrifos, sparse = TRUE)
    sa <- principal_response_curves(sp, dat$dose, dat$week, n_axes = 2, permutations = 19)
    sb <- principal_response_curves(dat$pyrifos, dat$dose, dat$week, n_axes = 2, permutations = 19)
    expect_equal(sa$coefficients, sb$coefficients, tolerance = 1e-12)
    expect_equal(sa$weights, sb$weights, tolerance = 1e-12)
    expect_identical(sa$test$p_value, sb$test$p_value)
  }

  # A control that isn't the first level.
  a <- principal_response_curves(dat$pyrifos, dat$dose, dat$week, control = 44, permutations = 0)
  expect_identical(a$control, '44')
  expect_identical(as.character(a$coefficients$treatment[1]), '44')
  expect_identical(levels(a$coefficients$treatment), levels(dat$dose))

  # Time as numbers: levels in numeric order.
  tn <- as.numeric(as.character(dat$week))
  a  <- principal_response_curves(dat$pyrifos, dat$dose, tn, n_axes = 2, permutations = 0)
  expect_identical(levels(a$coefficients$time), as.character(sort(unique(tn))))
  expect_equal(a$coefficients[-2], pr$coefficients[-2])



  # Edge cases. ====

  # A time point without control samples: NA coefficients, with a warning.
  keep <- !(dat$dose == '0' & dat$week == '8')
  expect_warning(
    a <- principal_response_curves(dat$pyrifos[keep, ], dat$dose[keep], dat$week[keep], permutations = 0),
    "No control samples at time '8'" )
  expect_true(all(is.na(a$coefficients$PRC1[a$coefficients$time == '8'])))
  expect_true(all(is.finite(a$coefficients$PRC1[a$coefficients$time != '8'])))

  # One time point: no conditions, an RDA on the treatments.
  wk <- dat$week == '4'
  a  <- principal_response_curves(dat$pyrifos[wk, ], dat$dose[wk], dat$week[wk], permutations = 19)
  o  <- rda(dat$pyrifos[wk, ], dat$dose[wk], n_axes = 1)
  expect_equal(unname(a$eig), attr(o, 'eig'), tolerance = 1e-12)
  expect_identical(attr(a$rda, 'inertia')[['conditional']], 0)
  expect_identical(a$test$df, c(1L, 4L))
  pm <- permanova(euclidean(dat$pyrifos[wk, ]), dat$dose[wk], permutations = 19)
  expect_equal(a$test['all', 'F'], unname(pm$statistic), tolerance = 1e-10)
  expect_identical(a$test['all', 'p_value'], pm$p_value)

  # Saturated: one sample per cell leaves no residual degrees of freedom.
  one <- !duplicated(data.frame(dat$dose, dat$week))
  a   <- principal_response_curves(dat$pyrifos[one, ], dat$dose[one], dat$week[one], permutations = 9)
  expect_identical(a$df_residual, 0L)
  expect_true(all(is.na(a$test$F)) && all(is.na(a$test$p_value)))
  expect_true(all(is.finite(a$coefficients$PRC1)))

  # n_axes beyond the constrained axes.
  two <- dat$dose %in% c('0', '44') & dat$week %in% c('-4', '1', '4')
  expect_warning(
    a <- principal_response_curves(dat$pyrifos[two, ], droplevels(dat$dose[two]), dat$week[two], n_axes = 5, permutations = 0),
    'reduced to 3' )
  expect_identical(names(a$eig), c('PRC1', 'PRC2', 'PRC3'))
  expect_identical(ncol(a$rda), 3L)

  # Treatments that explain nothing beyond time.
  y <- matrix(rep(1:4, each = 3), 12, 3) + 0
  expect_error(
    principal_response_curves(y, rep(c('c', 't', 't'), 4), rep(1:4, each = 3), permutations = 0),
    'none of the variation' )



  # Printing and plotting. ====

  out <- capture.output(print(b))
  expect_true(any(grepl('Principal Response Curves: 132 samples, 4 treatments', out, fixed = TRUE)))
  expect_true(any(grepl('Permutation tests (19 permutations, whole plots)', out, fixed = TRUE)))
  expect_true(any(grepl('Coefficients on PRC2', out, fixed = TRUE)))
  expect_true(any(grepl('No permutation tests', capture.output(print(pr)), fixed = TRUE)))
  out <- capture.output(v <- withVisible(print(pr)))
  expect_identical(v$visible, FALSE)
  expect_identical(v$value, pr)

  grDevices::pdf(NULL)
  drawn <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  expect_true(drawn(pr))
  expect_true(drawn(b, axis = 2, species = FALSE, legend = NULL, main = 'PRC'))
  expect_true(drawn(pr, species = TRUE, col = 1:5, lwd = 2))
  expect_true(drawn(a))
  expect_true(drawn(principal_response_curves(
    dat$pyrifos, dat$dose, paste0('w', dat$week), permutations = 0 )))  # times not numbers
  expect_error(plot(pr, axis = 3), '`axis` must be an axis number from 1 to 2')
  expect_error(plot(pr, species = -1), '`species`')
  expect_error(plot(pr, legend = 'nowhere'), '`legend`')
  grDevices::dev.off()



  # Input checks. ====

  y <- dat$pyrifos
  expect_error(principal_response_curves(y, time = dat$week), '`treatment` is required')
  expect_error(principal_response_curves(y, dat$dose), '`time` is required')
  expect_error(principal_response_curves(y, NULL, dat$week), '`treatment` is NULL')
  expect_error(principal_response_curves(y, dat$dose[-1], dat$week), '`treatment`')
  expect_error(principal_response_curves(y, replace(dat$dose, 1, NA), dat$week), 'must not contain NAs')
  expect_error(principal_response_curves(y, rep('a', 132), dat$week), 'at least two values')
  expect_error(principal_response_curves(y, dat$dose, dat$week, control = 'x'), '`control` must be one of')
  expect_error(principal_response_curves(y, dat$dose, dat$week, control = c('0', '6')), '`control` must be one of')
  expect_error(principal_response_curves(y[1:2, ], dat$dose[1:2], dat$week[1:2]), 'at least 3 samples')
  expect_error(principal_response_curves(y, dat$dose, dat$week, scaling = 'contribution'), "scaling 1, 2 or 3")
  expect_error(principal_response_curves(y, dat$dose, dat$week, scaling = 4), '`scaling`')
  expect_error(principal_response_curves(y, dat$dose, dat$week, n_axes = 0), '`n_axes`')
  expect_error(principal_response_curves(y, dat$dose, dat$week, permutations = -1), '`permutations`')
  expect_error(principal_response_curves(y, dat$dose, dat$week, strata = dat$week, plots = dat$ditch), 'only one of')
  expect_error(principal_response_curves(y, dat$dose, dat$week, strata = dat$week, permutations = dat$perm_free), 'only one of')
  expect_error(principal_response_curves(y, dat$dose, dat$week, plots = rep(1:6, 22)), 'one sample at each time point')
  expect_error(principal_response_curves(y[1:24, ], dat$dose[1:24], dat$week[1:24], plots = rep(1, 24)), 'one sample at each time point')
  expect_error(principal_response_curves(y, 'dose', 'week', meta[, -1]), "`treatment` must be the name of a column of `data`")
  expect_error(principal_response_curves(y, 'dose', 'week', list(1)), '`data` must be a data.frame')
  expect_error(principal_response_curves(y, dat$dose, dat$week, seed = 'a'))
  expect_error(principal_response_curves(y, dat$dose, dat$week, cpus = 0), '`cpus`')
  bad <- y
  bad[1, 1] <- NA
  expect_error(principal_response_curves(bad, dat$dose, dat$week, permutations = 0), 'NA')

#})
