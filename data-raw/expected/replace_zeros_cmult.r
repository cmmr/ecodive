# Writes tests/tinytests/expected/replace_zeros_cmult.rds: inputs and
# reference results for replace_zeros_czm(), replace_zeros_gbm(),
# replace_zeros_sq() and replace_zeros_bl(), so that the tests need no
# reference package.
#
# Run from the package root, with zCompositions 1.6.2 installed:
#
#     Rscript data-raw/expected/replace_zeros_cmult.r
#
# zCompositions (GPL; output only, never its source): cmultRepl() with
# output = "prop" and z.delete = FALSE. Its z.warning check (which stops
# with z.delete = FALSE when a row or column is at least 80% zeros) is
# turned off with z.warning = 1, as ecodive never drops a sample or
# feature. The formulas behind ecodive's implementation were established
# by data-raw/expected/zero_replacement_probe.r.

library(zCompositions)
stopifnot(as.character(packageVersion('zCompositions')) == '1.6.2')

cr <- function (x, method, frac, threshold, adjust) {
  r <- suppressMessages(cmultRepl(
    x, label = 0, method = method, output = 'prop', frac = frac,
    threshold = threshold, adjust = adjust, z.warning = 1, z.delete = FALSE,
    suppress.print = TRUE ))
  unname(as.matrix(r))
}


# The tables. ====

tables <- list(

  # One zero, three samples.
  probe_3x3 = rbind(c(5, 0, 3), c(2, 4, 1), c(10, 7, 9)),

  # One or two zeros per sample, unequal depths (20, 20, 40, 50).
  probe_4x5 = rbind(c(10, 0, 5, 3, 2), c(4, 6, 0, 0, 10), c(1, 20, 7, 8, 4), c(30, 2, 3, 0, 15)),

  # A shallow sample whose imputed values exceed the features' smallest
  # observed proportions (adjust).
  probe_adjust = rbind(c(1, 0, 0, 1), c(100, 1, 3, 50), c(3, 5, 2, 0), c(40, 2, 0, 9), c(7, 8, 9, 10)),

  # Feature 2 is non-zero in one sample only: GBM refuses it.
  probe_single = rbind(c(10, 0, 5, 3, 2), c(4, 0, 0, 1, 10), c(1, 20, 7, 8, 4), c(30, 0, 3, 0, 15)),

  # Two samples.
  probe_2x4 = rbind(c(3, 0, 5, 2), c(1, 4, 6, 0)),

  # Not whole numbers.
  probe_real = rbind(c(10, 0, 5, 3, 2), c(4, 6, 0, 0, 10), c(1, 20, 7, 8, 4), c(30, 2, 3, 0, 15)) / 7,

  ex_counts = unname(ecodive::ex_counts) )

# A sparse random table: 25 samples x 60 features, about 70% zeros, depths
# of a few hundred to a few thousand; features non-zero in at least 2
# samples (for GBM).
set.seed(55)
sp <- matrix(rnbinom(25 * 60, mu = 40, size = 0.4), 25, 60)
sp[runif(length(sp)) < 0.65] <- 0
sp <- sp[, colSums(sp > 0) >= 2]
stopifnot(all(rowSums(sp) > 0), mean(sp == 0) > 0.6)
tables$sparse <- sp


# The cases. ====

cases <- list()
add <- function (table, method, frac = 0.65, threshold = 0.5, adjust = TRUE) {
  x <- tables[[table]]
  id <- sprintf('%s-%s-%s-%s-%s', table, method, frac, threshold, adjust)
  cases[[id]] <<- list(
    id = id, table = table, method = tolower(method), frac = frac,
    threshold = threshold, adjust = adjust,
    result = cr(x, method, frac, threshold, adjust) )
}

for (tb in names(tables))
  for (m in c('CZM', 'GBM', 'SQ', 'BL'))
    for (adj in c(TRUE, FALSE)) {
      if (m == 'GBM' && tb %in% c('probe_single', 'probe_2x4')) next  # refused
      add(tb, m, adjust = adj)
    }

for (tb in c('probe_adjust', 'sparse')) {
  for (m in c('GBM', 'SQ', 'BL')) add(tb, m, frac = 0.3)
  add(tb, 'CZM', frac = 0.3, threshold = 0.8)
  add(tb, 'CZM', frac = 0.9, threshold = 0.2, adjust = FALSE)
}

# Each table's refusals, for the record.
gbm_single <- tryCatch(cr(tables$probe_single, 'GBM', 0.65, 0.5, TRUE), error = conditionMessage)
stopifnot(is.character(gbm_single))


saveRDS(
  list(tables = tables, cases = unname(cases)),
  file.path('tests', 'tinytests', 'expected', 'replace_zeros_cmult.rds'),
  version = 2 )
