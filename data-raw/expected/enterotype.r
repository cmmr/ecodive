# Expected results for tests/tinytests/test-enterotype.r, from
# cluster::pam(), cluster::silhouette() and clusterSim::index.G1(). Run from
# the package root:
#
#     Rscript data-raw/expected/enterotype.r
#
# Writes tests/tinytests/expected/enterotype.rds.
#
# cluster and clusterSim are GPL: they are only run here, to make fixtures
# (see plans/community_typing.md). Their source was not read.
#
# Cases:
#   dirichlet: Jensen-Shannon distances of 60 simulated genus profiles,
#              drawn as Dirichlet around three community types, 30 genera.
#   ex_counts: Jensen-Shannon distances of ecodive's ex_counts (4 samples).
#   varespec:  Jensen-Shannon distances of vegan's varespec (24 samples).
#   ties:      Manhattan distances of points on a small integer grid, with
#              duplicated points, so that BUILD and SWAP meet exact ties.
#   euclid:    Euclidean distances of iris, for the centroid form of the
#              CH index (index.G1() with centrotypes = "centroids" on the
#              data equals Calinski & Harabasz's index).
#
# For every case and number of types k: pam()'s id.med, clustering and
# objective; the silhouette widths and neighbours, in the samples' order;
# and index.G1() with centrotypes = "medoids", given the distances and, as
# data, the samples' proportions (as Arumugam et al.'s script does) or, for
# `ties`, the points. The inputs are stored (except ex_counts, which ships
# with ecodive, and iris); the tests recompute the distances.

library(ecodive)
library(cluster)
library(clusterSim)
data(varespec, package = 'vegan')

rdirichlet <- function (n, alpha) {
  g <- matrix(rgamma(n * length(alpha), alpha), n, byrow = TRUE)
  g / rowSums(g)
}

set.seed(45)
S     <- 30
base  <- rep(0.3, S)
types <- list(replace(base, 1, 12), replace(base, 2, 12), replace(base, 3:4, 6))
grp   <- rep(1:3, c(25, 20, 15))
props <- do.call(rbind, lapply(seq_along(grp), function (i) rdirichlet(1, types[[grp[i]]] * 3)))
dirichlet_counts <- round(props * 2000)
dirichlet_counts[dirichlet_counts == 0 & runif(length(dirichlet_counts)) < 0.5] <- 0
dimnames(dirichlet_counts) <- list(sprintf('S%02d', seq_along(grp)), sprintf('G%02d', 1:S))
storage.mode(dirichlet_counts) <- 'integer'
stopifnot(all(rowSums(dirichlet_counts) > 0))

set.seed(7)
ties_x <- matrix(sample(0:4, 40, replace = TRUE), 20, 2)
ties_x <- rbind(ties_x, ties_x[c(1, 1, 5, 9), ])
dimnames(ties_x) <- list(paste0('P', 1:nrow(ties_x)), c('x', 'y'))

varespec <- as.matrix(varespec)
props    <- function (m) m / rowSums(m)

# `x` is what index.G1() is given as data: the samples' proportions, as in
# Arumugam et al.'s script, or the points themselves.
cases <- list(
  dirichlet = list(d = jensen(dirichlet_counts), x = props(dirichlet_counts), k = 2:6),
  ex_counts = list(d = jensen(ecodive::ex_counts), x = props(ecodive::ex_counts), k = 2:3),
  varespec  = list(d = jensen(varespec), x = props(varespec), k = 2:6),
  ties      = list(d = dist(ties_x, 'manhattan'), x = ties_x, k = 2:6) )

ref <- function (d, x, k) {
  p <- pam(d, k, diss = TRUE)
  s <- silhouette(p$clustering, d)  # in the samples' order
  list(
    k          = k,
    id_med     = p$id.med,
    clustering = unname(p$clustering),
    objective  = unname(p$objective),
    silhouette = unname(s[, 'sil_width']),
    neighbor   = unname(s[, 'neighbor']),
    ch_medoids = index.G1(x, p$clustering, d = d, centrotypes = 'medoids') )
}

expected <- lapply(cases, function (cs) {
  list(fits = lapply(cs$k, function (k) ref(cs$d, cs$x, k)))
})

# The centroid form, on Euclidean data.
x  <- as.matrix(iris[, 1:4])
de <- dist(x)
expected$euclid <- list(
  fits = lapply(2:5, function (k) {
    p <- pam(de, k, diss = TRUE)
    list(
      k            = k,
      clustering   = unname(p$clustering),
      ch_centroids = index.G1(x, p$clustering) )
  }) )

expected$dirichlet_counts <- dirichlet_counts
expected$ties_x           <- ties_x
expected$varespec_counts  <- varespec

expected$versions <- c(
  cluster    = as.character(packageVersion('cluster')),
  clusterSim = as.character(packageVersion('clusterSim')),
  R          = paste(R.version$major, R.version$minor, sep = '.') )

saveRDS(expected, 'tests/tinytests/expected/enterotype.rds', version = 2)
