# Writes tests/tinytests/expected/alpha_metrics.rds: inputs and reference
# results for the alpha diversity metrics in test-alpha_div.r, so that the
# tests need no reference package.
#
# Run from the package root, with ecodive installed:
#
#     Rscript data-raw/expected/alpha_metrics.r inputs   # the input tables
#     python data-raw/expected/alpha_metrics.py          # scikit-bio's outputs
#     Rscript data-raw/expected/alpha_metrics.r          # folds them in
#
# References:
# * scikit-bio 0.7.4 (BSD-3): skbio.diversity.alpha on every sample, with
#   the orders q = 0, 0.5, 1, 2, 3, Inf for hill, renyi and tsallis.
#   Values are hexadecimal floats, read back exactly.
# * vegan 2.7.8 (GPL-2; output only): renyi(x, scales, hill = TRUE),
#   renyi(x, scales), diversity(x) / log(specnumber(x)) and
#   diversity(x, 'simpson'), on the same tables.

dir <- file.path('data-raw', 'expected', 'alpha_metrics')
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(TRUE)


# ---- Input tables ------------------------------------------------------

make_tables <- function () {

  # Edge cases: one feature; all singletons; an empty sample; ties (for
  # strong); singletons with doubletons; singletons without doubletons;
  # doubletons without singletons; all equal; one dominant feature.
  edge <- rbind(
    c(7, 0, 0, 0, 0, 0, 0, 0),
    c(1, 1, 1, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0),
    c(5, 5, 5, 2, 2, 0, 0, 0),
    c(1, 1, 2, 3, 8, 0, 4, 0),
    c(1, 1, 1, 3, 8, 0, 4, 9),
    c(2, 2, 3, 3, 0, 0, 0, 0),
    c(4, 4, 4, 4, 4, 4, 4, 4),
    c(1000, 1, 1, 0, 0, 0, 0, 0) )

  # The same samples, with 12 more (absent) features: the Gini index
  # changes, the other metrics don't.
  edge_wide <- cbind(edge, matrix(0, nrow(edge), 12))

  # Typical amplicon-like counts, many zeros.
  set.seed(1)
  rand <- matrix(rnbinom(30 * 40, mu = 20, size = 0.3), 30, 40)

  # Deep samples, large counts.
  set.seed(2)
  deep <- matrix(rnbinom(6 * 25, mu = 1e5, size = 0.5), 6, 25)

  tabs <- list(
    ex_counts = unname(ecodive::ex_counts),
    edge      = edge,
    edge_wide = edge_wide,
    rand      = rand,
    deep      = deep )

  lapply(tabs, function (m) { storage.mode(m) <- 'integer'; m })
}


if (identical(args, 'inputs')) {
  tabs <- make_tables()
  for (nm in names(tabs))
    write.table(
      tabs[[nm]], file.path(dir, sprintf('table-%s.csv', nm)),
      sep = ',', row.names = FALSE, col.names = FALSE )
  quit(save = 'no')
}


# ---- Fold in scikit-bio's results -----------------------------------------

read_table <- function (nm) {
  m <- as.matrix(read.csv(file.path(dir, sprintf('table-%s.csv', nm)), header = FALSE))
  dimnames(m) <- NULL
  storage.mode(m) <- 'integer'
  m
}

res <- read.csv(file.path(dir, 'results.csv'), colClasses = 'character')

tables <- lapply(setNames(nm = unique(res$table)), read_table)
stopifnot(identical(tables, make_tables()[names(tables)]))

skbio <- data.frame(
  table  = res$table,
  sample = as.integer(res$sample),
  metric = res$metric,
  q      = ifelse(res$q == '', NA_real_, ifelse(res$q == 'inf', Inf, as.numeric(res$q))),
  value  = as.numeric(res$value),
  stringsAsFactors = FALSE )
stopifnot(!anyNA(skbio$value[!is.nan(skbio$value)]))


# ---- vegan (output only) ------------------------------------------------

orders <- c(0, 0.5, 1, 2, 3, Inf)

vegan <- do.call(rbind, lapply(names(tables), function (nm) {
  x  <- tables[[nm]]
  ok <- rowSums(x) > 0
  x  <- x[ok, , drop = FALSE]
  smp <- which(ok)

  h <- vegan::renyi(x, scales = orders, hill = TRUE)
  r <- vegan::renyi(x, scales = orders)
  h <- matrix(as.matrix(h), nrow(x))
  r <- matrix(as.matrix(r), nrow(x))

  rbind(
    data.frame(table = nm, sample = rep(smp, length(orders)), metric = 'hill',
               q = rep(orders, each = nrow(x)), value = as.vector(h)),
    data.frame(table = nm, sample = rep(smp, length(orders)), metric = 'renyi',
               q = rep(orders, each = nrow(x)), value = as.vector(r)),
    data.frame(table = nm, sample = smp, metric = 'pielou', q = NA_real_,
               value = vegan::diversity(x) / log(vegan::specnumber(x))),
    data.frame(table = nm, sample = smp, metric = 'dominance', q = NA_real_,
               value = 1 - vegan::diversity(x, 'simpson')) )
}))
vegan$value <- unname(vegan$value)


saveRDS(
  file    = file.path('tests', 'tinytests', 'expected', 'alpha_metrics.rds'),
  version = 2,
  object  = list(
    tables   = tables,
    skbio    = skbio,
    vegan    = vegan,
    versions = c(skbio = '0.7.4', vegan = as.character(packageVersion('vegan'))) ))
