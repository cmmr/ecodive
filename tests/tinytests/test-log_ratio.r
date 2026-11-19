#test_that("log_ratio", {

  # log_ratio() and balance(): scikit-bio 0.7.4's and Qurro 0.9.0's outputs
  # (data-raw/expected/log_ratio.py), compositions' balance() as a
  # cross-check, identities with ilr() and the base-R formulas, and argument
  # checks.

  ref <- readRDS(file.path('expected', 'log_ratio.rds'))

  # The values alone.
  plain <- function (v) { attributes(v) <- NULL; v }

  # A full SBP whose first row is num (+1) vs den (-1), the others chains.
  full_sbp <- function (num, den, D) {
    chain <- function (p) lapply(seq_len(length(p) - 1L), function (k) {
      r <- numeric(D); r[p[k]] <- 1; r[p[-seq_len(k)]] <- -1; r })
    first <- numeric(D); first[num] <- 1; first[den] <- -1
    rows  <- c(list(first), chain(num), chain(den))
    rest  <- setdiff(seq_len(D), c(num, den))
    if (length(rest)) {
      top <- numeric(D); top[c(num, den)] <- 1; top[rest] <- -1
      rows <- c(rows, list(top), chain(rest))
    }
    do.call(rbind, rows)
  }


  # -----------------------------------------------------------------------
  # References: scikit-bio, Qurro, compositions
  # -----------------------------------------------------------------------

  for (i in seq_along(ref$cases)) {
    cs <- ref$cases[[i]]
    x  <- ref$tables[[cs$table]]
    info <- paste(cs$table, i)

    if (cs$pc > 0) expect_warning(balance(x, cs$num, cs$den), 'Zeros detected')
    b <- balance(x, cs$num, cs$den, pseudocount = if (cs$pc > 0) cs$pc)
    expect_identical(attr(b, 'pseudocount'), cs$pc, info = info)
    expect_equal(plain(b), cs$skbio,  tolerance = 1e-12, info = info)
    expect_equal(plain(b), cs$closed, tolerance = 1e-12, info = info)
    if (!is.null(ref$compositions))
      expect_equal(plain(b), ref$compositions$values[[i]], tolerance = 1e-12, info = info)

    # The matching ilr() coordinate, to the last bit.
    sbp <- full_sbp(cs$num, cs$den, ncol(x))
    expect_identical(plain(b), unname(ilr(x, sbp, pseudocount = if (cs$pc > 0) cs$pc)[, 1]), info = info)

    # Qurro's log ratio of the raw counts: NA where it drops a sample.
    l <- suppressWarnings(log_ratio(x, cs$num, cs$den))
    if (cs$table != 'prop')
      expect_identical(plain(l), cs$qurro, info = info)
    base <- log(rowSums(x[, cs$num, drop = FALSE]) / rowSums(x[, cs$den, drop = FALSE]))
    base[!is.finite(base)] <- NA
    expect_equal(plain(l), base, tolerance = 1e-14, info = info)
    expect_equal(
      plain(log_ratio(x, cs$num, cs$den, pseudocount = 1)), cs$lr_pc1,
      tolerance = 1e-14, info = info )
  }


  # -----------------------------------------------------------------------
  # Selectors, names, attributes
  # -----------------------------------------------------------------------

  set.seed(1)
  x <- matrix(rpois(10 * 8, 4) + 1L, 10, 8,
    dimnames = list(paste0('s', 1:10), paste0('f', 1:8)) )
  lgl <- function (i) seq_len(8) %in% i

  for (fn in c(log_ratio, balance)) {
    a <- fn(x, c(2, 5), c(7, 1, 3))
    expect_identical(names(a), rownames(x))
    expect_identical(attr(a, 'numerator'),   c('f2', 'f5'))
    expect_identical(attr(a, 'denominator'), c('f1', 'f3', 'f7'))   # in feature order
    expect_identical(attr(a, 'pseudocount'), 0)
    expect_identical(fn(x, c('f5', 'f2'), c('f1', 'f3', 'f7')), a)
    expect_identical(fn(x, lgl(c(2, 5)), lgl(c(1, 3, 7))), a)
    expect_identical(fn(x, c('f2', 'f5'), c(1L, 3L, 7L)), a)          # mixed forms
    expect_identical(fn(x, factor(c('f2', 'f5')), lgl(c(1, 3, 7))), a)

    # margin = 2: the same named vector.
    expect_identical(fn(t(x), c(2, 5), c(7, 1, 3), margin = 2), a)

    # Other input types.
    expect_identical(fn(as.data.frame(x), c(2, 5), c(1, 3, 7)), a)
    expect_identical(fn(Matrix::Matrix(x, sparse = TRUE), c(2, 5), c(1, 3, 7)), a)
    expect_identical(fn(slam::as.simple_triplet_matrix(x), c(2, 5), c(1, 3, 7)), a)

    # Unnamed: indices as attributes, no names.
    u <- fn(unname(x), c(2, 5), c(7, 1, 3))
    expect_null(names(u))
    expect_identical(attr(u, 'numerator'),   c(2L, 5L))
    expect_identical(attr(u, 'denominator'), c(1L, 3L, 7L))
    expect_identical(plain(u), plain(a))

    # A vector is one sample.
    v <- fn(x[3, ], 1, 2)
    expect_identical(length(v), 1L)
    expect_identical(plain(v), plain(fn(x, 1, 2))[3])

    # Swapping the sides flips the sign; scale doesn't matter.
    expect_equal(plain(fn(x, c(1, 3, 7), c(2, 5))), -plain(a), tolerance = 1e-14)
    expect_equal(fn(x * 7.5, c(2, 5), c(1, 3, 7)), a, tolerance = 1e-14)
    expect_equal(fn(x / rowSums(x), c(2, 5), c(1, 3, 7)), a, tolerance = 1e-14)
  }

  # A balance's sign flip is exact (the same entries, negated).
  expect_identical(plain(balance(x, c(1, 3, 7), c(2, 5))), -plain(balance(x, c(2, 5), c(1, 3, 7))))

  # One feature on each side: the log ratio is the balance / sqrt(1/2).
  expect_equal(plain(log_ratio(x, 4, 6)), plain(balance(x, 4, 6)) / sqrt(1 / 2), tolerance = 1e-14)

  # Features in neither group don't change a balance (to rounding).
  expect_equal(plain(balance(x, 2:3, 5)), plain(balance(x[, c(2, 3, 5)], 1:2, 3)), tolerance = 1e-14)

  # The method attribute; two features is enough.
  expect_identical(attr(log_ratio(x, 1, 2), 'method'), 'log_ratio')
  expect_identical(attr(balance(x, 1, 2), 'method'), 'balance')
  expect_equal(plain(balance(x[, 1:2], 1, 2)), unname(sqrt(1 / 2) * log(x[, 1] / x[, 2])), tolerance = 1e-14)


  # -----------------------------------------------------------------------
  # Zeros
  # -----------------------------------------------------------------------

  z <- x
  z[1, 1:2] <- 0L
  z[4, 5]   <- 0L
  z[6, 3]   <- 0L

  # log_ratio(): zeros are fine while sums are positive; a zero sum is NA.
  expect_silent(l <- log_ratio(z, 1:3, 4:5))
  expect_true(all(is.finite(l)))
  expect_warning(
    l <- log_ratio(z, 1:2, c(3, 6)),
    "1 sample ('s1') has a numerator or denominator that sums to 0, so its log ratio is NA. A positive `pseudocount` gives it a value.",
    fixed = TRUE )
  expect_identical(is.na(l), c(s1 = TRUE, setNames(rep(FALSE, 9), paste0('s', 2:10))))
  expect_identical(plain(l)[-1], plain(log_ratio(z[-1, ], 1:2, c(3, 6))))
  z2 <- z; z2[2, 3] <- 0L; z2[2, 6] <- 0L
  expect_warning(log_ratio(z2, 1:2, c(3, 6)), "2 samples ('s1', 's2') have", fixed = TRUE)
  expect_warning(log_ratio(unname(z2), 1:2, c(3, 6)), "2 samples ('1', '2') have", fixed = TRUE)
  expect_silent(l <- log_ratio(z, 1:2, c(3, 6), pseudocount = 0.5))
  expect_equal(
    plain(l), unname(log(rowSums(z[, 1:2] + 0.5) / rowSums(z[, c(3, 6)] + 0.5))), tolerance = 1e-14 )
  expect_identical(attr(l, 'pseudocount'), 0.5)

  # log_ratio(pseudocount = NULL): the shared rule, on the whole table.
  expect_warning(l <- log_ratio(z, 1:2, 3:4, pseudocount = NULL), 'Zeros detected')
  expect_identical(attr(l, 'pseudocount'), 0.5)
  expect_identical(attr(log_ratio(x, 1:2, 3:4, pseudocount = NULL), 'pseudocount'), 0)

  # balance(): NULL is ilr()'s rule on the whole table, even when the zeros
  # are outside the selection.
  expect_warning(b <- balance(z, 6:7, 8), 'Zeros detected')
  expect_identical(attr(b, 'pseudocount'), 0.5)
  expect_identical(plain(b), unname(suppressWarnings(ilr(z, full_sbp(6:7, 8, 8)))[, 1]))
  expect_identical(attr(balance(z, 6:7, 8, pseudocount = 1), 'pseudocount'), 1)

  # balance(pseudocount = 0): NA where a selected feature is zero; zeros
  # elsewhere don't matter.
  expect_warning(
    b <- balance(z, 1:2, 5:6, pseudocount = 0),
    "2 samples ('s1', 's4') have a zero in a numerator or denominator feature, so their balances are NA. Use `pseudocount = NULL`, or replace_zeros() first, to give them a value.",
    fixed = TRUE )
  expect_identical(which(is.na(b)), c(s1 = 1L, s4 = 4L))
  expect_silent(b <- balance(z, 6:7, 8, pseudocount = 0))
  expect_equal(
    plain(b), unname(sqrt(2 / 3) * (rowMeans(log(z[, 6:7])) - log(z[, 8]))), tolerance = 1e-14 )
  expect_equal(plain(b), plain(balance(z[, 6:8], 1:2, 3)), tolerance = 1e-14)

  # The output of replace_zeros() needs no pseudocount.
  r <- replace_zeros(z, 'mult')
  expect_silent(b <- balance(r, 1:2, 5:6))
  expect_identical(attr(b, 'pseudocount'), 0)
  expect_true(all(is.finite(b)))
  expect_identical(plain(b), unname(ilr(r, full_sbp(1:2, 5:6, 8))[, 1]))
  expect_silent(log_ratio(r, 1:2, 5:6))

  # Sparse tables, with stored zeros.
  m <- Matrix::Matrix(z, sparse = TRUE)
  m@x[1] <- 0
  expect_identical(
    suppressWarnings(log_ratio(m, 1:2, c(3, 6))),
    suppressWarnings(log_ratio(as.matrix(m), 1:2, c(3, 6))) )
  expect_identical(
    suppressWarnings(balance(m, 1:3, 4:8)),
    suppressWarnings(balance(as.matrix(m), 1:3, 4:8)) )


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  for (fn in c(log_ratio, balance)) {
    expect_error(fn(x, denominator = 1), '`numerator` is required')
    expect_error(fn(x, 1), '`denominator` is required')
    expect_error(fn(x, integer(0), 1), '`numerator` selects no features.', fixed = TRUE)
    expect_error(fn(x, NULL, 1), '`numerator` must select features')
    expect_error(fn(x, list(1), 2), '`numerator` must select features')
    expect_error(fn(x, c(1, NA), 2), '`numerator` must not contain NA.', fixed = TRUE)
    expect_error(fn(x, rep(FALSE, 8), 2), 'all FALSE')
    expect_error(fn(x, c(TRUE, FALSE), 2),
      'A logical `numerator` needs one value per feature: `counts` has 8 features, not 2.', fixed = TRUE)
    expect_error(fn(x, 1, c('f2', 'nope', 'zip')),
      "`denominator` names 2 features not in `counts`: 'nope', 'zip'.", fixed = TRUE)
    expect_error(fn(x, 1, 'nope'), "`denominator` names a feature not in `counts`: 'nope'.", fixed = TRUE)
    expect_error(fn(x, c(0, 9, 2), 1),
      '`numerator` has indices outside 1 to 8 (`counts` has 8 features): 0, 9.', fixed = TRUE)
    expect_error(fn(x, 1.5, 2), '`numerator` must select features')
    expect_error(fn(x, -1, 2), 'an index outside 1 to 8')
    expect_error(fn(x, c(1, 2, 1), 3), '`numerator` selects a feature more than once: 1.', fixed = TRUE)
    expect_error(fn(x, c('f1', 'f1'), 3), "selects a feature more than once: 'f1'.", fixed = TRUE)
    expect_error(fn(unname(x), 'f1', 2), '`counts` has no feature names; give `numerator` as indices.', fixed = TRUE)
    expect_error(fn(x, 1:3, c(3, 2, 6)),
      "`numerator` and `denominator` must not share features, but they share 2 features: 'f2', 'f3'.", fixed = TRUE)
    expect_error(fn(unname(x), 1:3, 3:4), "they share 1 feature: '3'.", fixed = TRUE)
    expect_error(fn(x, 1, 2, pseudocount = -1), 'single non-negative number')
    expect_error(fn(x, 1, 2, pseudocount = c(1, 2)), 'single non-negative number')
    expect_error(fn(x, 1, 2, margin = 3), 'margin')
    expect_error(fn(x[, 1, drop = FALSE], 1, 1), 'at least 2 features')
    expect_error(fn(-x, 1, 2), 'negative')
    expect_error(fn(x * NA, 1, 2), 'NA')
    expect_error(fn(x * 0, 1, 2), 'only zeros')
  }

#})
