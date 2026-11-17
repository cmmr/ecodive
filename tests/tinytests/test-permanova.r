#test_that("permanova", {

  ex   <- readRDS('expected/permanova.rds')
  d    <- bray(ex$counts)
  meta <- ex$meta



  # Same tables as the reference, given the same permutations: stored
  # results, so at home only. Every permuted F, and so every p-value,
  # matches exactly. ====

  if (at_home()) {
    for (cs in ex$cases) {

      obj <- permanova(
        d, meta, stats::as.formula(cs$formula), by = cs$by,
        permutations = ex$perms, cpus = 1 )
      res <- obj$table

      info <- paste(cs$formula, cs$by)
      tbl  <- unname(as.matrix(res))
      rows <- rownames(res)

      if (cs$ref == 'dbrda') { # No Total row or R2 column.
        keep <- !rows %in% c('Conditional', 'Total')
        tbl  <- tbl[keep, -3, drop = FALSE]
        rows <- rows[keep]
      }

      expect_identical(rows, cs$rows, info = info)
      expect_equal(tbl[, -ncol(tbl)], cs$table[, -ncol(tbl)], tolerance = 1e-10, info = info)
      expect_identical(tbl[, ncol(tbl)], cs$table[, ncol(tbl)], info = info)

      if (!is.null(cs$f_perm))
        expect_equal(unname(attr(res, 'f_perm')), cs$f_perm, tolerance = 1e-10, info = info)

      # The model as a whole: the "Model" row of by = 'overall'.
      if (cs$by == 'overall') {
        expect_identical(obj$statistic[['F']], res[['F']][rownames(res) == 'Model'], info = info)
        expect_identical(obj$p_value, res[['Pr(>F)']][rownames(res) == 'Model'], info = info)
        expect_identical(attr(obj, 'f_perm'), unname(attr(res, 'f_perm')[, 'Model']), info = info)
      }
    }

    # Ordered factors, poly(), character and logical columns, unbalanced
    # and unused levels, permutations within blocks, and the one-factor
    # path with 500 samples.
    sets <- list(
      meta2  = list(d = d,                  data = ex$meta2,   perms = ex$perms),
      blocks = list(d = d,                  data = meta,       perms = ex$block_perms),
      n500   = list(d = bray(ex$counts500), data = ex$meta500, perms = ex$perms500) )

    for (cs in ex$more_cases) {

      s   <- sets[[cs$set]]
      res <- permanova(
        s$d, s$data, stats::as.formula(cs$formula), by = cs$by,
        permutations = s$perms, cpus = 2 )$table

      info <- paste(cs$set, cs$formula, cs$by)
      tbl  <- unname(as.matrix(res))

      expect_identical(rownames(res), cs$rows, info = info)
      expect_equal(tbl[, -ncol(tbl)], cs$table[, -ncol(tbl)], tolerance = 1e-10, info = info)
      expect_identical(tbl[, ncol(tbl)], cs$table[, ncol(tbl)], info = info)

      if (!is.null(cs$f_perm))
        expect_equal(unname(attr(res, 'f_perm')), cs$f_perm, tolerance = 1e-10, info = info)
    }
  }



  # Result format. ====

  obj <- permanova(d, meta, ~ trt + ph, by = 'terms', permutations = 99, cpus = 1)
  res <- obj$table

  expect_identical(class(obj), c('ecodive_permanova', 'ecodive_test'))
  expect_identical(
    names(obj),
    c('statistic', 'p_value', 'z', 'n_perm', 'r2', 'table', 'by', 'permutation') )
  expect_identical(names(obj$statistic), 'F')
  expect_identical(obj$n_perm, 99L)
  expect_identical(obj$by, 'terms')
  expect_identical(obj$permutation, 'free')
  expect_identical(length(attr(obj, 'f_perm')), 99L)
  expect_true(obj$z > 0)  # A real treatment effect.
  expect_equal(obj$z, (obj$statistic[['F']] - mean(attr(obj, 'f_perm'))) / sd(attr(obj, 'f_perm')))

  # The core is the model as a whole, whatever `by` is.
  ov <- permanova(d, meta, ~ trt + ph, permutations = 99, cpus = 1)
  for (by in c('terms', 'margin', 'onedf')) {
    o <- permanova(d, meta, ~ trt + ph, by = by, permutations = 99, cpus = 1)
    expect_identical(o[c('statistic', 'p_value', 'z', 'n_perm', 'r2')],
                     ov[c('statistic', 'p_value', 'z', 'n_perm', 'r2')], info = by)
    expect_identical(attr(o, 'f_perm'), attr(ov, 'f_perm'), info = by)
  }
  expect_identical(ov$statistic[['F']], ov$table$F[[1]])
  expect_identical(ov$p_value, ov$table[['Pr(>F)']][[1]])
  expect_identical(ov$r2, ov$table$R2[[1]])

  # One term: the core and the table's first row agree, for every `by`
  # with one row per term (the one-factor path and the general one).
  for (f in list(~ trt, ~ ph)) {
    for (by in c('overall', 'terms', 'margin')) {
      o <- permanova(d, meta, f, by = by, permutations = 99, cpus = 1)
      expect_identical(o$statistic[['F']], o$table$F[[1]], info = paste(deparse(f), by))
      expect_identical(o$p_value, o$table[['Pr(>F)']][[1]], info = paste(deparse(f), by))
      expect_identical(attr(o, 'f_perm'), unname(attr(o$table, 'f_perm')[, 1]))
    }
  }

  # One row of the shared core fields.
  df <- as.data.frame(obj)
  expect_identical(names(df), c('test', 'stat_name', 'statistic', 'p_value', 'z', 'n_perm'))
  expect_identical(df$test, 'permanova')
  expect_identical(df$stat_name, 'F')
  expect_identical(df$statistic, obj$statistic[['F']])

  out <- capture.output(print(obj))
  expect_identical(out[[1]], 'Permutational multivariate analysis of variance (PERMANOVA)')
  expect_true(any(grepl('^F = .*, p = .*, z = ', out)))
  expect_true(any(grepl('Terms added sequentially', out)))
  expect_true(any(grepl('^Residual', out)))

  expect_inherits(res, 'anova')
  expect_inherits(res, 'data.frame')
  expect_identical(names(res), c('Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'))
  expect_identical(rownames(res), c('trt', 'ph', 'Residual', 'Total'))
  expect_identical(dim(attr(res, 'f_perm')), c(99L, 2L))
  expect_identical(colnames(attr(res, 'f_perm')), c('trt', 'ph'))
  expect_equal(sum(res$SumOfSqs[1:3]), res$SumOfSqs[[4]], tolerance = 1e-12)
  expect_equal(res$R2[[4]], 1)
  expect_true(all(res[['Pr(>F)']][1:2] >= 1 / 100))
  expect_true(any(grepl('Number of permutations: 99', attr(res, 'heading'))))

  # Conditional row.
  obj <- permanova(d, meta, ~ trt + Condition(site), permutations = 9, cpus = 1)
  res <- obj$table
  expect_identical(obj$statistic[['F']], res$F[[2]])
  expect_identical(rownames(res), c('Conditional', 'Model', 'Residual', 'Total'))
  expect_equal(sum(res$SumOfSqs[1:3]), res$SumOfSqs[[4]], tolerance = 1e-12)
  expect_identical(res$Df, c(1, 2, 32, 35))

  # No permutations: no p-values.
  obj <- permanova(d, meta, ~ trt, permutations = 0, cpus = 1)
  res <- obj$table
  expect_true(is.na(res[['Pr(>F)']][[1]]))
  expect_identical(dim(attr(res, 'f_perm')), c(0L, 1L))
  expect_true(is.na(obj$p_value) && is.na(obj$z))
  expect_false(is.na(obj$statistic[['F']]))
  expect_identical(attr(obj, 'f_perm'), numeric(0))
  expect_false(any(grepl('p = ', capture.output(print(obj)), fixed = TRUE)))



  # Deterministic: the same for any number of cpus, and from run to run. ====

  for (f in list(~ trt, ~ trt * site + ph, ~ ph + Condition(site))) {
    r1 <- permanova(d, meta, f, by = 'terms', permutations = 199, seed = 7, cpus = 1)
    r3 <- permanova(d, meta, f, by = 'terms', permutations = 199, seed = 7, cpus = 3)
    expect_identical(r1, r3, info = deparse(f))
    expect_identical(r1, permanova(d, meta, f, by = 'terms', permutations = 199, seed = 7, cpus = 1))
  }

  # Large enough for run_parallel() to start threads.
  big  <- bray(big_mtx)
  bmeta <- data.frame(
    row.names = rownames(big_mtx),
    g = rep(c('a', 'b', 'c', 'd'), 26),
    x = sin(seq_len(104)) )
  r1 <- permanova(big, bmeta, ~ g + x, by = 'terms', permutations = 999, cpus = 1)
  r4 <- permanova(big, bmeta, ~ g + x, by = 'terms', permutations = 999, cpus = 4)
  expect_identical(r1, r4)
  r1 <- permanova(big, bmeta, ~ g, permutations = 999, cpus = 1)
  r4 <- permanova(big, bmeta, ~ g, permutations = 999, cpus = 4)
  expect_identical(r1, r4)

  # The seed changes the permutations, not the statistics.
  r2 <- permanova(d, meta, ~ trt, permutations = 99, seed = 2, cpus = 1)
  r3 <- permanova(d, meta, ~ trt, permutations = 99, seed = 3, cpus = 1)
  expect_identical(r2$table[, 1:4], r3$table[, 1:4])
  expect_false(identical(attr(r2, 'f_perm'), attr(r3, 'f_perm')))



  # The one-factor path gives the same sums of squares as the general one. ====

  dv  <- dist_vector(d)
  pm  <- perm_matrix(50L, dv$n)
  x   <- model_design(meta, ~ trt, dv$labels, dv$n)$x
  q   <- design_basis(x)$cols
  gen <- perm_ss(dv$dm, dv$n, cols = q, pm = pm, seed = 5L)
  fac <- perm_ss(dv$dm, dv$n, grp = as.integer(meta$trt) - 1L, pm = pm, seed = 5L)
  expect_equal(fac$ss[, 1], rowSums(gen$ss), tolerance = 1e-12)
  expect_identical(fac$total, gen$total)
  expect_equal(gen$total, sum(as.vector(d)^2) / attr(d, 'Size'), tolerance = 1e-14)



  # Strata. ====

  # A factor that is constant within strata can't change: every permuted F
  # equals the observed one.
  obj <- permanova(d, meta, ~ site, permutations = 99, strata = 'site', cpus = 1)
  res <- obj$table
  expect_identical(obj$permutation, 'within strata')
  expect_true(is.na(obj$z))
  expect_identical(res[['Pr(>F)']][[1]], 1)
  expect_equal(as.vector(attr(res, 'f_perm')), rep(res$F[[1]], 99), tolerance = 1e-10)

  # The same p-value as the reference, up to sampling error (stored, so at
  # home only).
  if (at_home()) {
    res <- permanova(d, meta, ~ ph, permutations = 9999, strata = meta$block, cpus = 1)$table
    expect_true(abs(res[['Pr(>F)']][[1]] - ex$strata_p) < 0.02)
  }

  # A strata column is not a term of the default formula.
  sb <- permanova(d, meta[c('ph', 'block')], strata = 'block', permutations = 19, cpus = 1)
  expect_identical(rownames(sb$table), c('Model', 'Residual', 'Total'))
  expect_identical(sb, permanova(d, meta, ~ ph, strata = 'block', permutations = 19, cpus = 1))

  expect_error(permanova(d, meta, ~ ph, permutations = ex$perms, strata = 'block'), 'not both')
  expect_error(permanova(d, meta, ~ ph, strata = 1:3), 'strata')
  expect_error(permanova(d, meta, ~ ph, strata = c(NA, meta$block[-1])), 'NA')



  # Metadata alignment. ====

  ref <- permanova(d, meta, ~ trt + ph, permutations = 19, cpus = 1)

  # Rows are matched by name.
  shuffled <- meta[rev(seq_len(nrow(meta))), ]
  expect_identical(permanova(d, shuffled, ~ trt + ph, permutations = 19, cpus = 1), ref)

  # ... or by position, without row names.
  unnamed <- meta
  rownames(unnamed) <- NULL
  expect_identical(permanova(d, unnamed, ~ trt + ph, permutations = 19, cpus = 1), ref)

  # A symmetric matrix works like a dist.
  expect_identical(permanova(as.matrix(d), meta, ~ trt + ph, permutations = 19, cpus = 1), ref)

  # A left-hand side is ignored; character columns are factors.
  expect_identical(permanova(d, meta, dm ~ trt + ph, permutations = 19, cpus = 1), ref)
  chr <- meta
  chr$trt <- as.character(chr$trt)
  expect_equal(
    permanova(d, chr, ~ trt + ph, permutations = 19, cpus = 1)$table$SumOfSqs,
    ref$table$SumOfSqs, tolerance = 1e-12 )

  # A single vector or factor, named or not.
  res <- permanova(d, meta$trt, permutations = 19, cpus = 1)
  expect_identical(rownames(res$table), c('Model', 'Residual', 'Total'))
  expect_equal(res$table$F[[1]], permanova(d, meta, ~ trt, permutations = 0)$table$F[[1]], tolerance = 1e-12)
  expect_identical(res, permanova(d, setNames(meta$trt, rownames(meta)), permutations = 19, cpus = 1))

  # All columns by default.
  expect_identical(
    permanova(d, meta[c('trt', 'ph')], by = 'terms', permutations = 19, cpus = 1),
    permanova(d, meta, ~ trt + ph, by = 'terms', permutations = 19, cpus = 1) )

  expect_error(permanova(d, meta[-1, ], ~ trt), "1 sample in `dm` has no row in `data`: 'S01'")
  expect_error(permanova(d, meta[-(1:2), ], ~ trt), '2 samples')
  extra <- rbind(meta, meta[1, ])
  rownames(extra)[nrow(extra)] <- 'X'
  # Rows for other samples are ignored, even with NAs.
  extra$ph[nrow(extra)] <- NA
  expect_identical(permanova(d, extra, ~ trt + ph, permutations = 19, cpus = 1), ref)

  # ... and factor levels found only in those rows are dropped.
  lv <- extra
  levels(lv$trt) <- c(levels(lv$trt), 'only_x')
  lv$trt[nrow(lv)] <- 'only_x'
  des <- model_design(lv, ~ trt, attr(d, 'Labels'), attr(d, 'Size'))
  expect_identical(levels(des$data$trt), levels(meta$trt))
  expect_identical(colnames(des$x), c('trtlow', 'trthigh'))
  expect_identical(
    permanova(d, lv, ~ trt, by = 'onedf', permutations = 19, cpus = 1),
    permanova(d, meta, ~ trt, by = 'onedf', permutations = 19, cpus = 1) )

  # ... including NAs in `strata`, as a column or as a vector.
  extra$block[nrow(extra)] <- NA
  sref <- permanova(d, meta, ~ ph, permutations = 19, strata = 'block', cpus = 1)
  expect_identical(permanova(d, extra, ~ ph, permutations = 19, strata = 'block', cpus = 1), sref)
  expect_identical(permanova(d, extra, ~ ph, permutations = 19, strata = extra$block, cpus = 1), sref)
  expect_error(permanova(d, extra[c(2:36, 1, 37), ], ~ ph, strata = c(NA, extra$block[-1])), 'NA')
  expect_error(permanova(d, unnamed[-1, ], ~ trt), 'has 35 rows')
  expect_error(permanova(d, meta, ~ nothere), "Not found in `data`: 'nothere'")
  na <- meta
  na$ph[3] <- NA
  expect_error(permanova(d, na, ~ ph), 'Missing values')
  expect_error(permanova(d), '`data` is required')
  expect_error(permanova(d, meta, ~ Condition(site)), 'no terms')
  expect_error(permanova(d, meta, 'trt'), 'formula')



  # Other arguments. ====

  expect_error(permanova(d, meta, ~ trt, permutations = -1), 'permutations')
  expect_error(permanova(d, meta, ~ trt, permutations = 1.5), 'permutations')
  expect_error(permanova(d, meta, ~ trt, permutations = ex$perms[, -1]), '36 columns')
  bad <- ex$perms
  bad[2, 1:2] <- 1L
  expect_error(permanova(d, meta, ~ trt, permutations = bad), 'Row 2')
  expect_error(permanova(d, meta, ~ trt, by = 'nope'))
  expect_error(permanova(d, meta, ~ trt, seed = 1.5), 'seed')
  expect_error(permanova(d, meta, ~ trt, cpus = 0), 'cpus')

  bad_d <- d
  bad_d[1] <- NA
  expect_error(permanova(bad_d, meta, ~ trt), 'NA')
  bad_d[1] <- -1
  expect_error(permanova(bad_d, meta, ~ trt), 'negative')

  # Saturated model: no residual degrees of freedom.
  sat <- data.frame(row.names = rownames(meta), id = rownames(meta))
  expect_warning(res <- permanova(d, sat, ~ id, permutations = 9), 'no residual')
  expect_true(is.na(res$table$F[[1]]))
  expect_true(is.na(res$statistic[['F']]) && is.na(res$p_value) && is.na(res$z))

  # Regression: values beyond the integer range became NA and slipped past
  # the per-row check into C as sample indices (a segfault).
  nd <- attr(d, 'Size')
  expect_error(permanova(d, meta, ~ trt, permutations = matrix(1e10, 1, nd)), 'one permutation of 1 to')
  expect_error(permanova(d, meta, ~ trt, permutations = matrix(0:(nd - 1), 1, nd)), 'one permutation of 1 to')

  # adonis2()'s style: the distances on the formula's left-hand side.
  expect_identical(permanova(d ~ trt + ph, meta, permutations = 19, cpus = 1), ref)
  expect_identical(permanova(d ~ trt + ph, data = meta, permutations = 19, cpus = 1), ref)
  expect_error(permanova(d ~ trt, meta, ~ ph), 'not both')

  # Degenerate tables.
  expect_warning(z <- permanova(d * 0, meta, ~ trt, permutations = 9, cpus = 1), 'All distances are 0')
  expect_true(is.na(z$statistic[['F']]) && is.na(z$r2) && all(is.na(attr(z, 'f_perm'))))
  z <- z$table
  expect_true(all(is.na(z$F)) && all(is.na(z$R2)) && all(is.na(attr(z, 'f_perm'))))
  const <- data.frame(row.names = rownames(meta), k = rep(1, nrow(meta)))
  expect_warning(z <- permanova(d, const, ~ k, permutations = 9, cpus = 1), 'No terms to test')
  expect_true(is.na(z$statistic[['F']]) && is.na(z$p_value))
  expect_identical(z$r2, 0)
  expect_identical(rownames(z$table), c('Residual', 'Total'))
  expect_true(any(grepl('Residual', capture.output(print(z)))))
  expect_warning(res <- permanova(d, sat, ~ id, permutations = 9), 'no residual')
  expect_true(all(is.na(attr(res, 'f_perm'))))
  al <- meta
  al$ph2 <- 2 * al$ph
  z <- suppressMessages(permanova(d, al, ~ trt + ph + ph2, by = 'margin', permutations = 9, cpus = 1))$table
  zero_df <- z$Df == 0 & !rownames(z) %in% c('Residual', 'Total')
  expect_true(all(is.na(z$F[zero_df])) && all(is.na(z[['Pr(>F)']][zero_df])))

  # Infinite values in a term.
  inf <- meta
  inf$ph[3] <- Inf
  expect_error(permanova(d, inf, ~ ph), 'Infinite values')

  # Duplicate sample names, with `data` matched by position.
  dd <- d
  attr(dd, 'Labels')[2] <- attr(dd, 'Labels')[1]
  expect_identical(
    permanova(dd, unnamed, ~ trt, permutations = 19, cpus = 1)$table[, 1:4],
    permanova(d,  unnamed, ~ trt, permutations = 19, cpus = 1)$table[, 1:4] )

  # A character formula is parsed, not evaluated (R 3.6 evaluated it).
  expect_error(permanova(d, meta, 'trt'),       '`formula` must be a formula')
  expect_error(permanova(d, meta, c('~ a', '~ b')), '`formula` must be a formula')
  expect_identical(
    permanova(d, meta, '~ trt + ph', permutations = 19, cpus = 1),
    permanova(d, meta,  ~ trt + ph,  permutations = 19, cpus = 1) )

#})
