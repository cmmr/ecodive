# Expected results for tests/tinytests/test-tsne.r, from Rtsne::Rtsne() on a
# distance matrix. Run from the package root:
#
#     Rscript data-raw/expected/tsne.r
#
# Writes tests/tinytests/expected/tsne.rds.
#
# Rtsne(is_distance = TRUE) finds Barnes-Hut neighbors with a vantage-point
# tree, which is exact only for metric distances without ties, so the input
# is Euclidean distances between random points.

library(Rtsne)

set.seed(1)
n   <- 60
ctr <- matrix(rep(c(0, 3, 6), c(20, 20, 20)), n, 5)
pts <- ctr + matrix(rnorm(n * 5), n, 5)
d   <- dist(pts)

# Starting coordinates: small noise, as Rtsne draws them, and PCoA.
init2 <- matrix(rnorm(n * 2, sd = 1e-4), n, 2)
init3 <- matrix(rnorm(n * 3, sd = 1e-4), n, 3)
init1 <- matrix(rnorm(n * 1, sd = 1e-4), n, 1)
pco2  <- unname(cmdscale(d, k = 2))

# Arguments, named as in ecodive::tsne(); `init` names one of the above.
cases <- list(
  list(k = 2L, perplexity = 10, theta = 0,   maxit = 500L, init = 'init2', stop_lying_iter = 100L, mom_switch_iter = 100L),
  list(k = 2L, perplexity = 10, theta = 0.5, maxit = 500L, init = 'init2', stop_lying_iter = 100L, mom_switch_iter = 100L),
  list(k = 2L, perplexity =  5, theta = 0.5, maxit = 200L, init = 'pco2',  stop_lying_iter =   0L, mom_switch_iter =   0L),
  list(k = 3L, perplexity = 12, theta = 0.3, maxit = 300L, init = 'init3', stop_lying_iter =  50L, mom_switch_iter =  80L),
  list(k = 1L, perplexity =  8, theta = 0.5, maxit = 300L, init = 'init1', stop_lying_iter =  50L, mom_switch_iter =  50L)
)

inits <- list(init1 = init1, init2 = init2, init3 = init3, pco2 = pco2)

for (i in seq_along(cases)) {
  a <- cases[[i]]
  r <- Rtsne(
    d, dims = a$k, perplexity = a$perplexity, theta = a$theta,
    max_iter = a$maxit, Y_init = inits[[a$init]], pca = FALSE,
    stop_lying_iter = a$stop_lying_iter, mom_switch_iter = a$mom_switch_iter )
  cases[[i]]$Y    <- r$Y
  cases[[i]]$cost <- tail(r$itercosts, 1)

  # After 20 iterations only, before rounding differences between
  # platforms (e.g. fused multiply-adds) can grow.
  r <- Rtsne(
    d, dims = a$k, perplexity = a$perplexity, theta = a$theta,
    max_iter = 20, Y_init = inits[[a$init]], pca = FALSE,
    stop_lying_iter = a$stop_lying_iter, mom_switch_iter = a$mom_switch_iter )
  cases[[i]]$Y20    <- r$Y
  cases[[i]]$cost20 <- tail(r$itercosts, 1)
}

expected <- list(

  versions = list(
    R     = R.version.string,
    Rtsne = as.character(packageVersion('Rtsne')) ),

  # Input: Euclidean distances between the rows of `pts`.
  pts   = pts,
  inits = inits,

  # Each case's arguments, then Rtsne's coordinates and final KL divergence
  # (the last of its `itercosts`), after `maxit` iterations and after 20.
  cases = cases
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/tsne.rds', version = 2)
