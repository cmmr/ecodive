# Writes tests/tinytests/expected/rarefaction_curve.rds: reference values
# for test-rarefaction_curve.r, so that the tests need no reference package.
#
# Run from the package root, with ecodive installed:
#
#     Rscript data-raw/expected/rarefaction_curve.r    # writes Python's inputs
#     python data-raw/expected/rarefaction_curve.py    # NumPy's and pandas' outputs
#     Rscript data-raw/expected/rarefaction_curve.r    # folds them in
#
# vegan and iNEXT are GPL: only their outputs are used, never their source.
# Generated with R 4.3.3, vegan 2.7.8, iNEXT 3.0.2, NumPy 2.5.5, pandas
# 3.1.1 and mpmath 1.4.1 (versions stored in `versions`).
#
# Tables (`tables`, each samples x features):
#   ex_counts  ecodive's example data (4 samples, N = 345 to 1,011).
#   edge       small samples: one observation; all singletons; one feature;
#              no doubletons (f2 = 0); no singletons (f1 = 0); one
#              singleton and no doubletons; a mixed sample.
#   deep       two simulated samples with about 10^5 observations, to test
#              the running products' range.
#
# Fields:
#   vegan      For each table, a data frame (sample, depth, S, se):
#              vegan::rarefy(x[i, ], depth, se = TRUE), Hurlbert's expected
#              richness and Heck et al.'s standard error, at every depth
#              from 1 to the sample's total for ex_counts and edge, and at
#              `deep_depths` for deep.
#   inext      The same rows (sample, depth, qD, SC) from iNEXT::iNEXT(x[i,
#              ], q = 0, datatype = 'abundance', size = depths, nboot = 0):
#              qD is the expected richness and SC the expected sample
#              coverage (at the sample's total, Chao et al. 2014's Eq. 12
#              estimate). Samples that iNEXT() refuses are left out.
#   exact      The same rows (sample, depth, S, sd, coverage) computed by
#              data-raw/expected/rarefaction_curve.py with mpmath at 50
#              significant digits, then rounded to doubles: Hurlbert's
#              expected richness, Heck et al.'s standard deviation, and the
#              expected coverage (Chao & Jost 2012; at the sample's total,
#              Chao et al. 2014's Eq. 12, with Appendix G for f2 = 0), from
#              the published formulas, as a reference without rounding
#              error.
#   grid       NumPy's np.linspace(1, max_depth, num = steps, dtype = int)
#              for the cases in grid_cases (max_depth, steps): a list of
#              integer vectors, in that order.
#   q2_median  pandas' DataFrame.describe() 50% percentile (q2-diversity's
#              per-sample summary of its iterations) of `iters`, the values
#              of Shannon's index from 7 rarefactions of `q2_table` at
#              depths `q2_depths` (rows ordered by sample, then depth; one
#              column per rarefaction), each rarefaction t being
#              alpha_div(rarefy(q2_table, d, seed = t - 1), 'shannon').

library(ecodive)
library(vegan)
library(iNEXT)
stopifnot(packageVersion('vegan') >= '2.7.8')

dir <- file.path('data-raw', 'expected', 'rarefaction_curve')
dir.create(dir, showWarnings = FALSE)


# Tables -----------------------------------------------------------------

edge <- matrix(0L, 7, 6, dimnames = list(
  c('one_obs', 'singletons', 'one_feature', 'no_f2', 'no_f1', 'f1_no_f2', 'mixed'),
  sprintf('F%d', 1:6) ))
edge['one_obs',     1]   <- 1L
edge['singletons',  1:5] <- 1L
edge['one_feature', 3]   <- 50L
edge['no_f2',       1:4] <- c(1L, 1L, 5L, 9L)
edge['no_f1',       1:3] <- c(2L, 3L, 10L)
edge['f1_no_f2',    1:4] <- c(1L, 4L, 4L, 7L)
edge['mixed',       1:6] <- c(1L, 1L, 2L, 3L, 8L, 30L)

set.seed(38)
deep <- t(vapply(1:2, function (i) {
  lam <- exp(rnorm(400, 3, 2.2))
  x   <- rpois(400, lam / sum(lam) * 1e5)
  x[1:30] <- pmax(x[1:30], 1L)   # some rare features
  x[31:40] <- 1L                 # singletons
  x[41:45] <- 2L                 # doubletons
  as.integer(x)
}, integer(400)))
dimnames(deep) <- list(c('D1', 'D2'), sprintf('F%03d', 1:400))
deep_depths <- function (N) sort(unique(c(1, 2, 10, 100, 1000, 5000, 10000,
  round(N * c(0.25, 0.5, 0.75, 0.9, 0.99)), N - 100, N - 10, N - 2, N - 1, N )))

tables <- list(ex_counts = ex_counts, edge = edge, deep = deep)


# vegan and iNEXT --------------------------------------------------------

vegan_rows <- function (x, depths_of) {
  do.call(rbind, lapply(rownames(x), function (s) {
    v <- x[s, ]
    d <- depths_of(sum(v))
    r <- suppressWarnings(vegan::rarefy(v, sample = d, se = TRUE))
    data.frame(sample = s, depth = as.integer(d), S = r[1, ], se = r[2, ],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

inext_rows <- function (x, depths_of) {
  do.call(rbind, lapply(rownames(x), function (s) {
    v <- x[s, ]
    v <- v[v > 0]
    d <- depths_of(sum(v))
    o <- tryCatch(
      iNEXT::iNEXT(v, q = 0, datatype = 'abundance', size = d, nboot = 0),
      error = function (e) NULL )
    if (is.null(o)) return (NULL)
    sb <- o$iNextEst$size_based
    sb <- sb[sb$m %in% d, ]
    data.frame(sample = s, depth = as.integer(sb$m), qD = sb$qD, SC = sb$SC,
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

every <- function (N) seq_len(N)
depths_of <- list(ex_counts = every, edge = every, deep = deep_depths)

vegan_fx <- Map(vegan_rows, tables, depths_of)
inext_fx <- Map(inext_rows, tables, depths_of)


# Inputs for Python ------------------------------------------------------

grid_cases <- data.frame(
  max_depth = c(1011, 1011, 345, 100, 7, 2, 1, 100007, 12345, 1000, 50, 99999),
  steps     = c(  20,   10,  20,  10, 10, 20, 5,     20,    37, 1000,  3,    7) )
utils::write.csv(grid_cases, file.path(dir, 'grid_cases.csv'), row.names = FALSE)

q2_table  <- ex_counts
q2_depths <- c(1L, 50L, 200L, 345L, 600L)
iters <- do.call(cbind, lapply(0:6, function (seed) {
  unlist(lapply(rownames(q2_table), function (s) {
    vapply(q2_depths, function (d) {
      if (sum(q2_table[s, ]) < d) return (NA_real_)
      unname(ecodive::alpha_div(ecodive::rarefy(q2_table, d, seed = seed), "shannon")[s])
    }, 0)
  }))
}))
rownames(iters) <- paste(rep(rownames(q2_table), each = length(q2_depths)), q2_depths)
iters <- iters[stats::complete.cases(iters), , drop = FALSE]
utils::write.csv(
  data.frame(row = rownames(iters), apply(iters, 2, sprintf, fmt = '%a'),
             check.names = FALSE),
  file.path(dir, 'iters.csv'), row.names = FALSE )

# The tables (long: table, sample, count, non-zero counts only) and the
# rows (table, sample, depth) for the high-precision reference.
long <- do.call(rbind, lapply(names(tables), function (nm) {
  x <- tables[[nm]]
  do.call(rbind, lapply(rownames(x), function (s) {
    v <- x[s, ]
    data.frame(table = nm, sample = s, count = as.integer(v[v > 0]))
  }))
}))
utils::write.csv(long, file.path(dir, 'counts.csv'), row.names = FALSE)
rows <- do.call(rbind, Map(function (nm, v) data.frame(table = nm, v[, c('sample', 'depth')]),
                           names(vegan_fx), vegan_fx))
utils::write.csv(rows, file.path(dir, 'rows.csv'), row.names = FALSE)


# Python's outputs -------------------------------------------------------

py <- file.path(dir, 'numpy.csv')
if (!file.exists(py)) {
  message('Now run: python data-raw/expected/rarefaction_curve.py; then this script again.')
  quit(save = 'no')
}

np   <- utils::read.csv(py, stringsAsFactors = FALSE)
grid <- lapply(seq_len(nrow(grid_cases)), function (i) as.integer(np$depth[np$case == i]))

pd <- utils::read.csv(file.path(dir, 'pandas.csv'), colClasses = 'character')
stopifnot(identical(pd$row, rownames(iters)))
q2_median <- as.numeric(pd$median)  # as.numeric() reads hexadecimal floats exactly
names(q2_median) <- pd$row

mp <- utils::read.csv(file.path(dir, 'mpmath.csv'), colClasses = 'character')
stopifnot(identical(paste(mp$table, mp$sample, mp$depth), paste(rows$table, rows$sample, rows$depth)))
exact <- lapply(split(mp, factor(mp$table, names(tables))), function (d) data.frame(
  sample   = d$sample,
  depth    = as.integer(d$depth),
  S        = as.numeric(d$S),
  sd       = as.numeric(d$sd),
  coverage = as.numeric(d$coverage),
  stringsAsFactors = FALSE ))

versions <- c(
  R      = paste(R.version$major, R.version$minor, sep = '.'),
  vegan  = as.character(packageVersion('vegan')),
  iNEXT  = as.character(packageVersion('iNEXT')),
  numpy  = readLines(file.path(dir, 'versions.txt'))[[1]],
  pandas = readLines(file.path(dir, 'versions.txt'))[[2]],
  mpmath = readLines(file.path(dir, 'versions.txt'))[[3]] )

fx <- list(
  tables      = tables,
  deep_depths = lapply(rowSums(deep), deep_depths),
  vegan       = vegan_fx,
  inext       = inext_fx,
  exact       = exact,
  grid_cases  = grid_cases,
  grid        = grid,
  q2_table    = q2_table,
  q2_depths   = q2_depths,
  iters       = iters,
  q2_median   = q2_median,
  versions    = versions )

saveRDS(fx, file.path('tests', 'tinytests', 'expected', 'rarefaction_curve.rds'), version = 2)
