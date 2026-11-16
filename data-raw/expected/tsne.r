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

# More cases, each with its own points and start: early exaggeration on 600
# points (a deep Barnes-Hut tree), 3 dimensions on 200 points (an octree),
# and the largest perplexity that 60 points allow, (n - 1) / 3.
set.seed(2)
blobs <- function (n, p, groups = 3) {
  ctr <- matrix(rep(seq(0, by = 3, length.out = groups), length.out = n), n, p)
  ctr + matrix(rnorm(n * p), n, p)
}

# `pts = NULL` means the points above.
more_cases <- list(
  list(pts = blobs(600, 2), k = 2L, perplexity = 30,     theta = 0.5, maxit = 150L, stop_lying_iter = 100L, mom_switch_iter = 100L),
  list(pts = blobs(200, 3), k = 3L, perplexity = 15,     theta = 0.5, maxit = 150L, stop_lying_iter =  50L, mom_switch_iter =  50L),
  list(pts = NULL,          k = 2L, perplexity = 59 / 3, theta = 0.5, maxit = 200L, stop_lying_iter =  50L, mom_switch_iter =  50L) )

for (i in seq_along(more_cases)) {
  a  <- more_cases[[i]]
  if (is.null(a$pts)) a$pts <- pts
  dd <- dist(a$pts)
  n  <- nrow(a$pts)
  more_cases[[i]]$init <- init <- matrix(rnorm(n * a$k, sd = 1e-4), n, a$k)
  for (it in c(20L, a$maxit)) {
    r <- Rtsne(
      dd, dims = a$k, perplexity = a$perplexity, theta = a$theta,
      max_iter = it, Y_init = init, pca = FALSE,
      stop_lying_iter = a$stop_lying_iter, mom_switch_iter = a$mom_switch_iter )
    if (it == 20L) {
      more_cases[[i]]$Y20    <- r$Y
      more_cases[[i]]$cost20 <- tail(r$itercosts, 1)
    } else {
      more_cases[[i]]$Y    <- r$Y
      more_cases[[i]]$cost <- tail(r$itercosts, 1)
    }
  }
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
  cases = cases,

  # The same, each case with its own start `init`, and points `pts` (NULL:
  # the points above).
  more_cases = more_cases
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/tsne.rds', version = 2)
