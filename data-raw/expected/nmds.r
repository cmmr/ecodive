# Expected results for tests/tinytests/test-nmds.r, from vegan::monoMDS(),
# vegan::metaMDS() and MASS::isoMDS(). Run from the package root:
#
#     Rscript data-raw/expected/nmds.r
#
# Writes tests/tinytests/expected/nmds.rds.

library(ecodive)

set.seed(1)
rand_counts <- function (n, p) {
  m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
  m[, 1] <- m[, 1] + 1
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
  storage.mode(m) <- 'integer'
  m
}

# Inputs. `plane` is points near a plane, with noise: a clear NMDS solution.
inputs <- list(
  bray  = rand_counts(60, 40),  # used as bray()
  jacc  = rand_counts(50, 8),   # used as jaccard(): many tied distances
  plane = cbind(runif(60), runif(60) * 0.6, rnorm(60, sd = 0.02)) )
rownames(inputs$plane) <- paste0('P', 1:60)

# A plain numeric matrix, without dimnames or other attributes.
plain <- function (m) matrix(as.vector(m), nrow(m))

dists <- list(
  bray  = bray(inputs$bray),
  jacc  = jaccard(inputs$jacc),
  plane = dist(inputs$plane) )


# vegan::monoMDS() from given starting configurations (uniform on [0, 1)),
# with up to 25 iterations, where results should agree closely on any
# build, and with the default 200.
monomds <- function (data, k, weakties = TRUE, stress = 1, starts = 2) {

  d    <- dists[[data]]
  n    <- attr(d, 'Size')
  init <- array(runif(n * k * starts), c(n, k, starts))

  run <- function (maxit) {
    res <- lapply(seq_len(starts), function (i) vegan::monoMDS(
      d, y = init[, , i], k = k, pc = FALSE, maxit = maxit,
      weakties = weakties, stress = stress ))
    list(
      points = array(sapply(res, function (r) r$points), c(n, k, starts)),
      stress = sapply(res, function (r) r$stress),
      iters  = sapply(res, function (r) r$iters),
      icause = sapply(res, function (r) r$icause) )
  }

  list(
    data = data, k = k, weakties = weakties, stress = stress, init = init,
    maxit25 = run(25), maxit200 = run(200)[c('stress', 'iters', 'icause')] )
}

monomds_cases <- list(
  monomds('bray',  k = 2),
  monomds('bray',  k = 3),
  monomds('bray',  k = 1),
  monomds('jacc',  k = 2),
  monomds('jacc',  k = 2, weakties = FALSE),
  monomds('bray',  k = 2, stress = 2),
  monomds('plane', k = 2) )

# Five iterations from one start.
init5 <- matrix(runif(60 * 2), 60, 2)
mono5 <- vegan::monoMDS(dists$bray, y = init5, k = 2, pc = FALSE, maxit = 5)


# metaMDS(), without and with half-change scaling (a `maxdist` attribute).
set.seed(2)
meta_plane <- vegan::metaMDS(dists$plane, trace = 0)

d_hc <- dists$plane / max(dists$plane)
attr(d_hc, 'maxdist') <- 1
set.seed(2)
meta_hc <- vegan::metaMDS(d_hc, trace = 0)

iso <- MASS::isoMDS(dists$plane, trace = FALSE)


expected <- list(

  versions = list(
    R     = R.version.string,
    vegan = as.character(packageVersion('vegan')),
    MASS  = as.character(packageVersion('MASS')) ),

  inputs = inputs,

  # Each case: data (name in `inputs`), monoMDS() arguments, starting
  # configurations (n x k x starts), and results.
  monomds = monomds_cases,

  monomds5 = list(
    init   = init5,
    points = plain(mono5$points),
    stress = mono5$stress,
    iters  = mono5$iters,
    icause = mono5$icause ),

  metamds_plane = list(points = plain(meta_plane$points), stress = meta_plane$stress),
  metamds_hc    = list(points = plain(meta_hc$points),    stress = meta_hc$stress),

  isomds_plane = list(points = plain(iso$points), stress = iso$stress) )

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/nmds.rds', version = 2)
