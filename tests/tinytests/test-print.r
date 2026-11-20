#test_that("print", {

  # ecodive's tables, ordinations and distance matrices print a preview, as
  # tibble does (?print_ecodive): a header, the first rows, the columns that
  # fit, and a footer.

  cpus <- 1L

  # The footer's ellipsis: one character in a UTF-8 locale, "..." otherwise
  # (print_sym()), as a regular expression.
  ell <- if (isTRUE(l10n_info()[['UTF-8']])) '\u2026' else '\\.\\.\\.'

  # How many rows, and how wide.
  expect_equal(print_n_rows(20, NULL), 20)
  expect_equal(print_n_rows(21, NULL), 10)
  expect_equal(print_n_rows(150, 3), 3)
  expect_equal(print_n_rows(150, Inf), 150)
  expect_equal(print_n_rows(0, NULL), 0)
  expect_error(print_n_rows(5, -1), '`n`')
  expect_error(print_n_rows(5, 'a'), '`n`')
  expect_error(print_width(5), '`width`')
  expect_identical(print_width(NULL), as.integer(getOption('width')))

  # Cells: three significant digits, lined up on the decimal point.
  expect_identical(print_cells(c(0.123456, 12.5, NA), 3)$s, c(' 0.123', '12.5  ', 'NA    '))
  expect_true(print_cells(1, 3)$right)
  expect_identical(print_cells(c(1e-10, 1), 3)$s, c('1e-10', '1'))
  expect_identical(print_cells(1:3, 3)$s, c('1', '2', '3'))
  expect_identical(print_cells(c(TRUE, NA), 3)$s, c('TRUE', 'NA'))
  expect_identical(print_cells(c('a', NA), 3), list(s = c('a', 'NA'), right = FALSE))
  expect_identical(print_cells(factor('b'), 3)$s, 'b')
  expect_identical(print_cells(list(1:3, 'a'), 3)$s, c('<integer [3]>', '<character [1]>'))
  expect_identical(nchar(print_cells(strrep('x', 50), 3)$s), 40L)
  expect_identical(print_cells(c(1, 2.5, 3), 3, blank = c(TRUE, FALSE, FALSE))$s, c('', '2.5', '3  '))
  expect_identical(
    vapply(list(1, 1L, 'a', factor('a'), ordered('a'), TRUE, list(), matrix(1)), print_type, ''),
    c('dbl', 'int', 'chr', 'fct', 'ord', 'lgl', 'list', 'mat') )


  # Data frames: a header, the column types, 10 of 30 rows, and a footer.
  x <- data.frame(a = 1:30, b = rep(letters, length.out = 30), c = seq(0.5, 15, 0.5))
  class(x) <- c('ecodive_simper', 'data.frame')
  attr(x, 'n_perm') <- 0L
  out <- capture.output(v <- withVisible(print(x)))
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_true(grepl('^# A data\\.frame: 30 . 3$', out[[1]]))
  expect_identical(out[[2]], '# SIMPER')
  expect_identical(out[[3]], '       a b         c')
  expect_identical(out[[4]], '   <int> <chr> <dbl>')
  expect_identical(out[[5]], ' 1     1 a       0.5')
  expect_identical(length(out), 2L + 2L + 10L + 3L)
  expect_true(grepl(paste0('^# ', ell, ' 20 more rows$'), out[[15]]))
  expect_identical(out[[16]], '# Use print(x, n = ...) to see more rows')
  expect_identical(out[[17]], '# Attributes: n_perm')
  expect_identical(length(capture.output(print(x, n = Inf))), 2L + 2L + 30L + 1L)
  expect_identical(length(capture.output(print(x, n = 0))), 2L + 2L + 3L)
  # The full data frame is still there.
  expect_identical(dim(as.data.frame(x)), c(30L, 3L))
  expect_identical(class(as.data.frame(x)), 'data.frame')

  # Row names are shown when the rows have them.
  rownames(x) <- paste0('row', 1:30)
  expect_true(startsWith(capture.output(print(x))[[5]], 'row1 '))

  # Columns that don't fit are named in the footer, which wraps.
  w <- as.data.frame(matrix(1:60, 2, 30))
  class(w) <- c('ecodive_simper', 'data.frame')
  out <- capture.output(print(w, width = 30))
  expect_true(all(nchar(out) <= 30))
  expect_true(any(grepl(paste0('^# ', ell, ' 26 more columns:'), out)))
  expect_true(any(grepl('^#   and 6 more$', out)))
  expect_false(any(grepl('V25', out, fixed = TRUE)))

  # No rows.
  out <- capture.output(print(x[0, ]))
  expect_identical(length(out), 2L + 2L + 1L)

  # Each data frame result has a header line naming it.
  expect_identical(print_df_label(structure(list(), class = c('ecodive_dirmult_ttest', 'data.frame'))), 'Dirichlet-multinomial t-test')
  expect_identical(print_df_label(structure(list(), groups = c('a', 'b'), class = c('ecodive_ancom', 'data.frame'))), 'ANCOM: a vs b')
  expect_identical(print_df_label(structure(list(), method = 'ancombc2', class = c('ecodive_ancombc_pairwise', 'data.frame'))), 'ANCOM-BC2 pairwise tests')
  expect_identical(print_df_label(structure(list(), method = 'ancombc', class = c('ecodive_ancombc', 'data.frame'))), 'ANCOM-BC')
  # Rows taken with `[` keep the class but not the attributes: no label line.
  sub_ <- structure(data.frame(a = 1:3), class = c('ecodive_ancombc', 'data.frame'))[2:3, , drop = FALSE]
  expect_null(print_df_label(sub_))
  out <- capture.output(print(sub_))
  expect_true(startsWith(out[[1]], '# A data.frame: 2 '))
  expect_false(any(out == '#'))
  expect_identical(print_df_label(structure(list(), partial = TRUE, class = c('ecodive_bioenv', 'data.frame'))), 'BIO-ENV (partial)')
  meta <- data.frame(Species = iris$Species)
  fit  <- envfit(pcoa(bray(iris[, 1:4], cpus = cpus), cpus = cpus), meta, permutations = 9, cpus = cpus)
  out  <- capture.output(print(fit))
  expect_identical(out[[2]], '# Fitted variables on Principal Coordinates Analysis of Bray-Curtis Dissimilarity')
  expect_true(any(out == '# Attributes: r2_perm, n_perm'))


  # Ordinations: size, title and % of variance, then 10 of 150 samples.
  o   <- pcoa(bray(iris[, 1:4], cpus = cpus), cpus = cpus)
  out <- capture.output(print(o))
  expect_true(grepl('^# A matrix: 150 samples . 3 axes$', out[[1]]))
  expect_identical(out[[2]], '# Principal Coordinates Analysis of Bray-Curtis Dissimilarity')
  expect_true(grepl('^# % of variance: [0-9.]+, [0-9.]+, [0-9.]+$', out[[3]]))
  expect_identical(strsplit(trimws(out[[4]]), ' +')[[1]], c('PCoA1', 'PCoA2', 'PCoA3'))
  expect_identical(length(out), 3L + 1L + 10L + 3L)
  expect_true(grepl(paste0('^# ', ell, ' 140 more samples$'), out[[15]]))
  expect_identical(out[[17]], '# Attributes: method, eig, prop, trace, correction, constant')
  expect_identical(length(capture.output(print(o, n = 2))), 3L + 1L + 2L + 3L)
  # No axes.
  z <- ord_result(matrix(numeric(0), 3, 0), 'pcoa', c('a', 'b', 'c'))
  expect_true(grepl('^# A matrix: 3 samples . 0 axes$', capture.output(print(z))[[1]]))


  # Distances: the top-left corner of the lower triangle.
  d   <- bray(ecodive::ex_counts, cpus = cpus)
  out <- capture.output(print(d))
  expect_identical(out[1:2], c('# A dist: 4 samples (lower triangle)', '# Bray-Curtis Dissimilarity'))
  expect_identical(strsplit(trimws(out[[3]]), ' +')[[1]], c('Saliva', 'Gums', 'Nose'))
  expect_true(grepl('^Gums +0\\.591$', out[[4]]))
  expect_identical(length(out), 2L + 1L + 3L + 1L)
  expect_identical(out[[7]], '# Attributes: method')
  # NA distances still show; the upper triangle is blank.
  dn <- d
  dn[2] <- NA
  expect_true(grepl('^Nose +NA ', capture.output(print(dn))[[5]]))
  # Clustered matrices say so; matrices over features count features.
  expect_true(any(capture.output(print(bray(ecodive::ex_counts, clust = 'average', cpus = cpus))) == '# Clustered (average linkage)'))
  expect_identical(capture.output(print(rho_p(ecodive::ex_counts + 1, cpus = cpus)))[[1]], '# A dist: 6 features (lower triangle)')
  # 10 rows of 29.
  d30 <- euclidean(matrix(1:60, 30), cpus = cpus)
  out <- capture.output(print(d30))
  expect_true(any(grepl(paste0('^# ', ell, ' 19 more samples$'), out)))
  expect_identical(length(capture.output(print(d30, n = Inf, width = 1000))), 2L + 1L + 29L + 1L)
  # Tiny matrices.
  one <- structure(numeric(0), Size = 1L, class = c('ecodive_dist', 'dist'))
  expect_identical(capture.output(print(one)), '# A dist: 1 sample (lower triangle)')
  expect_identical(capture.output(print(structure(numeric(0), Size = 0L, class = c('ecodive_dist', 'dist')))), '# A dist: 0 samples (lower triangle)')

#})
