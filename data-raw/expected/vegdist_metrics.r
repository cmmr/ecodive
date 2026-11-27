# Writes tests/tinytests/expected/vegdist_metrics.rds: input tables and
# reference results for the vegdist() metrics in test-vegdist_metrics.r
# (kulczynski(), mountford(), raup_crick(), binomial_deviance(),
# chao_jaccard(), cao(), mahalanobis_distance(), alt_gower(),
# chisq_distance(), and soergel() as vegdist()'s "jaccard"), so that the
# tests need no reference package.
#
# Run from the package root, with vegan 2.7-8 and abdiv 0.2.0 installed,
# in two steps (the second step folds in SciPy's Mahalanobis distances):
#
#     Rscript data-raw/expected/vegdist_metrics.r inputs
#     python  data-raw/expected/vegdist_metrics.py
#     Rscript data-raw/expected/vegdist_metrics.r
#
# References (outputs only; vegan's source was not read):
# * vegan 2.7-8 (GPL-2): vegdist() with every method above, and with
#   binary = TRUE for kulczynski, binomial and altGower (the other
#   metrics take no `norm` in ecodive, or are presence/absence already);
#   on the protocol grids for Mountford (a, b, j) and Raup-Crick
#   (S, a, b, j).
# * abdiv 0.2.0 (MIT): weighted_kulczynski_second(), binomial_deviance()
#   and modified_mean_character_difference() on every pair of samples.
# * SciPy (BSD-3): scipy.spatial.distance.mahalanobis() with VI =
#   inv(cov) (data-raw/expected/vegdist_metrics.py).

stopifnot(
  packageVersion('vegan') == '2.7.8',
  packageVersion('abdiv') == '0.2.0' )

out_dir <- file.path('data-raw', 'expected', 'vegdist_metrics')
dir.create(out_dir, showWarnings = FALSE)


# ---- Input tables ------------------------------------------------------

make_tables <- function () {

  # Edge cases: an empty sample (2) and another (13); a sample nested in
  # another (4 in 3); identical samples (5, 6); a disjoint pair (7, 8);
  # singletons and doubletons for Chao's estimator (9, 10); one feature
  # (11, 12).
  edge <- rbind(
    c(3, 1, 0, 2, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(5, 4, 3, 2, 1, 0, 0, 0, 0, 0),
    c(2, 0, 3, 0, 1, 0, 0, 0, 0, 0),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
    c(4, 4, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 9, 1, 0, 0, 0),
    c(1, 2, 1, 2, 7, 0, 0, 0, 0, 0),
    c(2, 1, 1, 0, 1, 2, 0, 0, 0, 0),
    c(6, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 6, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )

  # Sparse 30 x 200 integer counts, with the same kinds of samples. At most
  # about 40 features per sample, where vegan's Mountford solution is
  # accurate.
  set.seed(1)
  sparse <- matrix(rnbinom(30 * 200, mu = 3, size = 0.3), 30, 200)
  sparse[matrix(runif(30 * 200), 30) < 0.88] <- 0
  sparse[2, ]  <- 0                                     # empty
  sparse[4, ]  <- sparse[3, ] * (runif(200) < 0.6)      # nested in 3
  sparse[6, ]  <- sparse[5, ]                           # identical
  sparse[7, 101:200] <- 0; sparse[8, 1:100] <- 0        # disjoint
  stopifnot(max(rowSums(sparse > 0)) <= 45)

  ints <- lapply(
    list(ex_counts = unname(ecodive::ex_counts), edge = edge, sparse = sparse),
    function (m) { storage.mode(m) <- 'double'; dimnames(m) <- NULL; m } )

  e <- new.env()
  utils::data('dune', 'varespec', package = 'vegan', envir = e)
  dune     <- unname(as.matrix(e$dune))
  varespec <- unname(as.matrix(e$varespec))
  storage.mode(dune) <- 'double'

  c(ints, list(
    dune     = dune,
    varespec = varespec,                                 # not integers
    percent  = ints$ex_counts / rowSums(ints$ex_counts) ))
}


# Tables for Mahalanobis: more samples than features, none collinear.
make_maha <- function () {
  e <- new.env()
  utils::data('varespec', 'varechem', package = 'vegan', envir = e)
  set.seed(2)
  gauss <- matrix(rnorm(30 * 6), 30, 6) %*% matrix(rnorm(36), 6, 6)
  list(
    varechem = unname(as.matrix(e$varechem)[, 1:8]),
    gauss    = gauss,
    counts   = matrix(as.double(rpois(25 * 5, 4)), 25, 5) )
}


tables <- make_tables()
maha   <- make_maha()


if (identical(commandArgs(TRUE), 'inputs')) {
  for (nm in names(maha))
    utils::write.table(
      matrix(sprintf('%.17g', maha[[nm]]), nrow(maha[[nm]])),
      file.path(out_dir, paste0('maha-', nm, '.csv')),
      sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE )
  quit(save = 'no')
}


# ---- vegan: every table and method -------------------------------------

v <- function (d) as.vector(d)

# vegdist() method of each ecodive metric, and whether vegan's binary =
# TRUE is compared with norm = 'binary'.
METHODS <- list(
  kulczynski        = list(m = 'kulczynski', binary = TRUE),
  mountford         = list(m = 'mountford',  binary = FALSE),
  raup_crick        = list(m = 'raup',       binary = FALSE),
  binomial_deviance = list(m = 'binomial',   binary = TRUE),
  chao_jaccard      = list(m = 'chao',       binary = FALSE),
  cao               = list(m = 'cao',        binary = FALSE),
  alt_gower         = list(m = 'altGower',   binary = TRUE),
  chisq_distance    = list(m = 'chisq',      binary = FALSE),
  soergel           = list(m = 'jaccard',    binary = FALSE) )

vegan_ref <- lapply(tables, function (x) {
  res <- list()
  for (f in names(METHODS)) {
    if (f == 'chao_jaccard' && any(x %% 1 != 0)) next
    res[[f]] <- suppressWarnings(v(vegan::vegdist(x, METHODS[[f]]$m)))
    if (METHODS[[f]]$binary)
      res[[paste0(f, '_binary')]] <- suppressWarnings(
        v(vegan::vegdist(x, METHODS[[f]]$m, binary = TRUE)) )
  }
  res
})

vegan_maha <- lapply(maha, function (x) v(vegan::vegdist(x, 'mahalanobis')))
vegan_maha$counts_binary <- v(vegan::vegdist(maha$counts > 2, 'mahalanobis', binary = TRUE))


# ---- Protocol grids ----------------------------------------------------

# A 2-row presence/absence table: richness a and b, j shared, in S
# features.
pa_pair <- function (a, b, j, S = a + b - j) {
  x <- matrix(0, 2, S)
  x[1, seq_len(a)] <- 1
  if (b > 0) x[2, (a - j) + seq_len(b)] <- 1
  x
}

# Mountford: every (a, b, j) with a, b in 1..12 and 0 <= j <= min(a, b),
# plus pairs with an empty sample.
mountford_grid <- expand.grid(a = 0:12, b = 0:12, j = 0:12)
mountford_grid <- subset(mountford_grid, j <= pmin(a, b))
mountford_grid$vegan <- mapply(
  function (a, b, j) suppressWarnings(v(vegan::vegdist(pa_pair(a, b, j, a + b - j + 1), 'mountford'))),
  mountford_grid$a, mountford_grid$b, mountford_grid$j )

# Raup-Crick: every (S, a, b, j) with S in {5, 10, 25, 60}, a, b <= 12
# (and <= S), and every possible j; plus far tails at S = 100.
raup_grid <- do.call(rbind, lapply(c(5, 10, 25, 60), function (S) {
  g <- expand.grid(S = S, a = 0:min(S, 12), b = 0:min(S, 12), j = 0:12)
  subset(g, j <= pmin(a, b) & a + b - j <= S)
}))
raup_grid <- rbind(raup_grid, data.frame(
  S = 100, a = c(40, 40, 40, 40, 70), b = c(40, 40, 40, 40, 60), j = c(0, 5, 30, 40, 42) ))
raup_grid$vegan <- mapply(
  function (S, a, b, j) suppressWarnings(v(vegan::vegdist(pa_pair(a, b, j, S), 'raup'))),
  raup_grid$S, raup_grid$a, raup_grid$b, raup_grid$j )
rownames(mountford_grid) <- rownames(raup_grid) <- NULL


# ---- abdiv: every pair -------------------------------------------------

abdiv_pairs <- function (x) {
  pr <- utils::combn(nrow(x), 2)
  f  <- function (fn) apply(pr, 2, function (p) fn(x[p[1], ], x[p[2], ]))
  list(
    kulczynski        = f(abdiv::weighted_kulczynski_second),
    binomial_deviance = f(abdiv::binomial_deviance),
    alt_gower         = f(abdiv::modified_mean_character_difference) )
}
abdiv <- lapply(tables[c('ex_counts', 'dune', 'varespec', 'sparse')], abdiv_pairs)


# ---- SciPy: Mahalanobis ------------------------------------------------

scipy_maha <- lapply(stats::setNames(nm = names(maha)), function (nm) {
  f <- file.path(out_dir, paste0('scipy-', nm, '.csv'))
  if (!file.exists(f)) stop('Run data-raw/expected/vegdist_metrics.py first.')
  as.numeric(readLines(f))
})


saveRDS(
  file    = file.path('tests', 'tinytests', 'expected', 'vegdist_metrics.rds'),
  version = 2,
  object  = list(
    tables         = tables,
    vegan          = vegan_ref,
    maha           = maha,
    vegan_maha     = vegan_maha,
    scipy_maha     = scipy_maha,
    mountford_grid = mountford_grid,
    raup_grid      = raup_grid,
    abdiv          = abdiv,
    versions       = c(
      vegan = as.character(packageVersion('vegan')),
      abdiv = as.character(packageVersion('abdiv')) ) ))
