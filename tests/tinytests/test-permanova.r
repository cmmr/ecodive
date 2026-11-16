#test_that("permanova", {

  ex   <- readRDS('expected/permanova.rds')
  d    <- bray(ex$counts)
  meta <- ex$meta



  # Same tables as the reference, given the same permutations. ====
  # Every permuted F, and so every p-value, matches exactly.

  for (cs in ex$cases) {

    res <- permanova(
      d, meta, stats::as.formula(cs$formula), by = cs$by,
      permutations = ex$perms, cpus = 1 )

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
  }



  # Result format. ====

  res <- permanova(d, meta, ~ trt + ph, by = 'terms', permutations = 99, cpus = 1)

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
  res <- permanova(d, meta, ~ trt + Condition(site), permutations = 9, cpus = 1)
  expect_identical(rownames(res), c('Conditional', 'Model', 'Residual', 'Total'))
  expect_equal(sum(res$SumOfSqs[1:3]), res$SumOfSqs[[4]], tolerance = 1e-12)
  expect_identical(res$Df, c(1, 2, 32, 35))

  # No permutations: no p-values.
  res <- permanova(d, meta, ~ trt, permutations = 0, cpus = 1)
  expect_true(is.na(res[['Pr(>F)']][[1]]))
  expect_identical(dim(attr(res, 'f_perm')), c(0L, 1L))



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
  expect_identical(r2[, 1:4], r3[, 1:4])
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
  res <- permanova(d, meta, ~ site, permutations = 99, strata = 'site', cpus = 1)
  expect_identical(res[['Pr(>F)']][[1]], 1)
  expect_equal(as.vector(attr(res, 'f_perm')), rep(res$F[[1]], 99), tolerance = 1e-10)

  # The same p-value as the reference, up to sampling error.
  res <- permanova(d, meta, ~ ph, permutations = 9999, strata = meta$block, cpus = 1)
  expect_true(abs(res[['Pr(>F)']][[1]] - ex$strata_p) < 0.02)

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
    permanova(d, chr, ~ trt + ph, permutations = 19, cpus = 1)$SumOfSqs,
    ref$SumOfSqs, tolerance = 1e-12 )

  # A single vector or factor, named or not.
  res <- permanova(d, meta$trt, permutations = 19, cpus = 1)
  expect_identical(rownames(res), c('Model', 'Residual', 'Total'))
  expect_equal(res$F[[1]], permanova(d, meta, ~ trt, permutations = 0)$F[[1]], tolerance = 1e-12)
  expect_identical(res, permanova(d, setNames(meta$trt, rownames(meta)), permutations = 19, cpus = 1))

  # All columns by default.
  expect_identical(
    permanova(d, meta[c('trt', 'ph')], by = 'terms', permutations = 19, cpus = 1),
    permanova(d, meta, ~ trt + ph, by = 'terms', permutations = 19, cpus = 1) )

  expect_error(permanova(d, meta[-1, ], ~ trt), "1 sample in `dm` has no row in `data`: 'S01'")
  expect_error(permanova(d, meta[-(1:2), ], ~ trt), '2 samples')
  extra <- rbind(meta, meta[1, ])
  rownames(extra)[nrow(extra)] <- 'X'
  expect_error(permanova(d, extra, ~ trt), "row of `data` is not in `dm`: 'X'")
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
  expect_true(is.na(res$F[[1]]))

#})
