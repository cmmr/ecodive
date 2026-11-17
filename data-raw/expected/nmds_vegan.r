# Reference solution quality for tests/tinytests/test-nmds.r, from
# vegan::metaMDS(), which also uses Kruskal's primary approach to ties. Run
# from the package root:
#
#     Rscript data-raw/expected/nmds_vegan.r
#
# Writes tests/tinytests/expected/nmds_vegan.rds: a named list of cases, each
# a list of
#   dm      the distances (no zeros, so vegan's stress and nmds()'s agree for
#           the same configuration);
#   stress  metaMDS()'s stress;
#   points  metaMDS()'s configuration (n x 2).
#
# Only vegan's outputs are stored; the test compares nmds() with them.

library(vegan)

set.seed(1)
n   <- 40
pts <- cbind(runif(n), runif(n) * 0.6, rnorm(n, sd = 0.02)) # near a plane
rownames(pts) <- paste0('P', seq_len(n))

data(dune, package = 'vegan', envir = environment())

dists <- list(
  plane = dist(pts),
  dune  = vegdist(dune, method = 'bray') )

cases <- lapply(dists, function (d) {
  stopifnot(all(d > 0))
  set.seed(2)
  vg <- metaMDS(d, k = 2, trace = 0)
  attr(d, 'call') <- NULL
  list(
    dm     = d,
    stress = vg$stress,
    points = unname(vg$points[, 1:2]) )
})

saveRDS(cases, file.path('tests', 'tinytests', 'expected', 'nmds_vegan.rds'), version = 2)
