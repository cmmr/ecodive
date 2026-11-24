# Writes tests/tinytests/expected/beta_partition.rds: input tables and
# reference results for the partition of beta diversity in
# test-beta_partition.r, so that the tests need no reference package.
#
# Run from the package root, with ecodive, betapart and abdiv installed:
#
#     Rscript data-raw/expected/beta_partition.r
#
# References (outputs only; betapart's source was not read):
# * betapart 1.6.1 (GPL-2 | GPL-3): beta.pair(), beta.pair.abund(),
#   beta.multi() and beta.multi.abund(), with every index family; and
#   beta.multi() / beta.multi.abund() of every subset of 4 of the first 9
#   samples of a table, which each of ecodive's random draws must match;
#   and beta.sample() / beta.sample.abund() means, which use R's RNG, for a
#   statistical comparison.
# * abdiv 0.2.0 (MIT): jaccard_turnover(), jaccard_nestedness(),
#   sorenson_turnover(), sorenson_nestedness(), bray_curtis_balanced(),
#   bray_curtis_gradient(), ruzicka_balanced() and ruzicka_gradient() on
#   every pair of samples.
# * betapart's own example data (ceram.s, ceram.n) are compared here but
#   not stored: the maximum differences are printed.

stopifnot(
  packageVersion('betapart') == '1.6.1',
  packageVersion('abdiv')    == '0.2.0' )


# ---- Input tables ------------------------------------------------------

make_tables <- function () {

  # Edge cases: an empty sample (2); a sample nested in another (4 in 3);
  # identical samples (5, 6); a disjoint pair (7, 8); equal richness with
  # turnover (9, 10); equal totals in different features (11, 12).
  edge <- rbind(
    c(3, 1, 0, 2, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(5, 4, 3, 2, 1, 0, 0, 0, 0, 0),
    c(2, 0, 3, 0, 1, 0, 0, 0, 0, 0),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    c(4, 4, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 9, 1, 0, 0, 0),
    c(0, 2, 2, 2, 0, 0, 0, 0, 0, 0),
    c(0, 0, 2, 2, 2, 0, 0, 0, 0, 0),
    c(6, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 6, 0, 0, 0, 0, 0, 0, 0, 0) )

  # Sparse 40 x 300 amplicon-like counts, with the same kinds of samples.
  set.seed(1)
  sparse <- matrix(rnbinom(40 * 300, mu = 3, size = 0.2), 40, 300)
  sparse[sparse > 0 & matrix(runif(40 * 300), 40) < 0.5] <- 0
  sparse[2, ]  <- 0                                     # empty
  sparse[4, ]  <- sparse[3, ] * (runif(300) < 0.6)      # nested in 3
  sparse[6, ]  <- sparse[5, ]                           # identical
  sparse[7, 151:300] <- 0; sparse[8, 1:150] <- 0        # disjoint
  sparse[10, ] <- sparse[9, ] * 3                       # abundance gradient

  lapply(
    list(ex_counts = unname(ecodive::ex_counts), edge = edge, sparse = sparse),
    function (m) { storage.mode(m) <- 'integer'; dimnames(m) <- NULL; m } )
}

tables <- make_tables()


# ---- betapart: pairwise ------------------------------------------------

pair_ref <- function (x) {
  pa <- (x > 0) * 1
  s  <- betapart::beta.pair(pa, 'sorensen')
  j  <- betapart::beta.pair(pa, 'jaccard')
  b  <- betapart::beta.pair.abund(x, 'bray')
  r  <- betapart::beta.pair.abund(x, 'ruzicka')
  v  <- function (d) as.vector(d)
  list(
    sorensen = cbind(v(s$beta.sim),     v(s$beta.sne),     v(s$beta.sor)),
    jaccard  = cbind(v(j$beta.jtu),     v(j$beta.jne),     v(j$beta.jac)),
    bray     = cbind(v(b$beta.bray.bal), v(b$beta.bray.gra), v(b$beta.bray)),
    ruzicka  = cbind(v(r$beta.ruz.bal),  v(r$beta.ruz.gra),  v(r$beta.ruz)) )
}

multi_ref <- function (x) {
  pa <- (x > 0) * 1
  c(
    sorensen = list(unname(unlist(betapart::beta.multi(pa, 'sorensen')))),
    jaccard  = list(unname(unlist(betapart::beta.multi(pa, 'jaccard')))),
    bray     = list(unname(unlist(betapart::beta.multi.abund(x, 'bray')))),
    ruzicka  = list(unname(unlist(betapart::beta.multi.abund(x, 'ruzicka')))) )
}

betapart_pair  <- lapply(tables, pair_ref)
betapart_multi <- lapply(tables, multi_ref)

# Every subset of 4 of the first 9 samples of `sparse`: the values each
# random draw of 4 of those samples must equal.
sub_x   <- tables$sparse[1:9, ]
subsets <- utils::combn(9, 4)
betapart_subsets <- lapply(
  stats::setNames(nm = c('sorensen', 'jaccard', 'bray', 'ruzicka')),
  function (idx) t(apply(subsets, 2L, function (s) multi_ref(sub_x[s, ])[[idx]])) )


# ---- betapart: resampled means (R's RNG) -------------------------------

set.seed(42)
pa <- (tables$sparse > 0) * 1
sample_means <- list(
  sorensen = unname(betapart::beta.sample(pa, 'sorensen', sites = 10, samples = 2000)$mean.values),
  bray     = unname(betapart::beta.sample.abund(tables$sparse, 'bray', sites = 10, samples = 2000)$mean.values) )


# ---- abdiv: pairwise ---------------------------------------------------

abdiv_pairs <- function (x) {
  m <- utils::combn(nrow(x), 2)
  f <- function (fn) apply(m, 2L, function (p) fn(x[p[1], ], x[p[2], ]))
  cbind(
    sorensen_turnover   = f(abdiv::sorenson_turnover),
    sorensen_nestedness = f(abdiv::sorenson_nestedness),
    jaccard_turnover    = f(abdiv::jaccard_turnover),
    jaccard_nestedness  = f(abdiv::jaccard_nestedness),
    bray_balanced       = f(abdiv::bray_curtis_balanced),
    bray_gradient       = f(abdiv::bray_curtis_gradient),
    soergel_balanced    = f(abdiv::ruzicka_balanced),
    soergel_gradient    = f(abdiv::ruzicka_gradient) )
}

# combn() order (1,2), (1,3), ... is the order of a dist's lower triangle.
abdiv <- lapply(tables, abdiv_pairs)


# ---- betapart's example data: compared, not stored ---------------------

ceram_dir <- system.file('data', package = 'betapart')
if (nzchar(ceram_dir) && requireNamespace('ecodive', quietly = TRUE)) {
  for (nm in c('ceram.s', 'ceram.n')) {
    e <- new.env()
    utils::data(list = nm, package = 'betapart', envir = e)
    x  <- as.matrix(e[[nm]])
    bp <- betapart::beta.pair(x, 'sorensen')
    ep <- ecodive::beta_partition(x, 'sorensen')
    bm <- unlist(betapart::beta.multi(x, 'sorensen'))
    em <- ecodive::beta_partition_multi(x, 'sorensen')$value
    cat(sprintf(
      '%s (%d x %d): pairwise max diff %.3g, multi max diff %.3g\n', nm, nrow(x), ncol(x),
      max(abs(c(as.vector(bp$beta.sim) - as.vector(ep$turnover),
                as.vector(bp$beta.sne) - as.vector(ep$nestedness),
                as.vector(bp$beta.sor) - as.vector(ep$total)))),
      max(abs(bm - em)) ))
  }
}


saveRDS(
  file    = file.path('tests', 'tinytests', 'expected', 'beta_partition.rds'),
  version = 2,
  object  = list(
    tables           = tables,
    betapart_pair    = betapart_pair,
    betapart_multi   = betapart_multi,
    subsets          = subsets,
    betapart_subsets = betapart_subsets,
    sample_means     = sample_means,
    abdiv            = abdiv,
    versions         = c(
      betapart = as.character(packageVersion('betapart')),
      abdiv    = as.character(packageVersion('abdiv')) ) ))
