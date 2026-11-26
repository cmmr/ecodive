# Expected results for tests/tinytests/test-mantel_correlogram.r, from
# vegan::mantel.correlog() on the oribatid mites (vegan's `mite` and
# `mite.xy`). Run from the package root:
#
#     Rscript data-raw/expected/mantel_correlogram.r
#
# Writes tests/tinytests/expected/mantel_correlogram.rds.
#
# vegan::mantel.correlog() has no argument for a permutation matrix: it runs
# vegan::mantel() once per class, and each run draws its own permutations
# through vegan's internal getPermuteMatrix(). To share one matrix with the
# tests, this script makes getPermuteMatrix() return that matrix (every
# class then uses it). mantel.correlog() passes its class matrix to
# mantel() first, so the permutations shuffle the classes, with D.eco held
# fixed, as mantel_correlogram() does with a supplied matrix. P-values (raw
# and adjusted) must then match exactly, not just statistically.

library(vegan)

data(mite)
data(mite.xy)

n     <- nrow(mite)
x     <- vegdist(decostand(mite, 'hellinger'), 'euclidean')  # community
xy    <- as.matrix(mite.xy)                                  # coordinates, m
geo   <- dist(xy)

set.seed(1)
perms <- shuffleSet(n, 199)
attributes(perms) <- list(dim = dim(perms))

assignInNamespace('getPermuteMatrix', function (perm, N, strata = NULL) perms, 'vegan')


# One case: vegan's classes and tests.
run <- function (...) {
  res <- mantel.correlog(x, geo, nperm = nrow(perms), ...)
  mr  <- res$mantel.res
  list(
    args       = list(...),
    breaks     = res$break.pts,
    class_mid  = unname(mr[, 'class.index']),
    n_dist     = unname(mr[, 'n.dist']),
    r          = unname(mr[, 'Mantel.cor']),
    p_value    = unname(mr[, 'Pr(Mantel)']),
    p_adjusted = if ('Pr(corrected)' %in% colnames(mr)) unname(mr[, 'Pr(corrected)']),
    n_tests    = res$n.tests )
}

cases <- list(
  run(),                                                       # Sturges, Pearson, progressive Holm
  run(r.type = 'spearman', break.pts = c(0.15, 1, 2, 3, 4, 6, 10)),
  run(r.type = 'kendall', n.class = 6, cutoff = FALSE, mult = 'BH', progressive = FALSE),
  run(n.class = 8, mult = 'sidak'),
  run(r.type = 'spearman', break.pts = seq(0, 10, by = 0.5), mult = 'bonferroni') )

saveRDS(
  object = list(
    x        = x,
    xy       = xy,
    perms    = perms,
    cases    = cases,
    versions = list(vegan = as.character(packageVersion('vegan'))) ),
  file     = 'tests/tinytests/expected/mantel_correlogram.rds',
  version  = 2 )
