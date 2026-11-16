# Expected results for tests/tinytests/test-pcoa.r, from ape::pcoa() and
# vegan::wcmdscale(). Run from the package root:
#
#     Rscript data-raw/expected/pcoa.r
#
# Writes tests/tinytests/expected/pcoa.rds.

library(ecodive)

set.seed(1)
counts <- matrix(rpois(150 * 40, 3) * rbinom(150 * 40, 1, 0.4), 150, 40)
counts[, 1] <- counts[, 1] + 1
dimnames(counts) <- list(paste0('S', 1:150), paste0('OTU', 1:40))
storage.mode(counts) <- 'integer'

d  <- bray(counts)
ap <- ape::pcoa(d)
vg <- vegan::wcmdscale(d, k = 3, eig = TRUE)


# More cases: the input, either counts (for bray()) or points (for dist()),
# the axes requested, and ape's and vegan's results for the axes that
# pcoa() returns (those with positive eigenvalues, at most k).
set.seed(2)
rand_counts <- function (n, p = 20) {
  m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
  m[, 1] <- m[, 1] + 1L
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
  storage.mode(m) <- 'integer'
  m
}

dups <- rand_counts(30)
dups[c(6, 7), ] <- dups[c(5, 5), ]
dups[20, ]      <- dups[12, ]

pts <- matrix(rnorm(40 * 4), 40, 4, dimnames = list(paste0('P', 1:40), NULL))

case <- function (name, k, counts = NULL, pts = NULL) {

  d  <- if (is.null(pts)) bray(counts) else dist(pts)
  ap <- ape::pcoa(d)
  vg <- suppressWarnings(vegan::wcmdscale(d, k = k, eig = TRUE))
  j  <- seq_len(min(k, ncol(ap$vectors)))

  # pcoa() warns when it returns fewer than min(k, n - 1) axes.
  k1   <- min(k, attr(d, 'Size') - 1)
  warn <- if (length(j) < k1) sprintf('only %d of the first %d', length(j), k1)

  list(
    name = name, counts = counts, pts = pts, k = k, warn = warn,
    ape  = list(
      vectors      = matrix(as.vector(ap$vectors[, j]), nrow(ap$vectors)),
      eig          = ap$values$Eigenvalues[j],
      relative_eig = ap$values$Relative_eig[j],
      trace        = ap$trace ),
    vegan = list(
      points = matrix(as.vector(vg$points[, j]), nrow(vg$points)),
      eig    = vg$eig[j],
      prop   = unname(summary(vegan::eigenvals(vg))[2, j]) ))
}

cases <- list(
  # k = n - 1: more axes than positive eigenvalues.
  case('bray, k = n - 1',   k = 11, counts = rand_counts(12)),
  # Euclidean distances: no negative eigenvalues; exactly 4 positive.
  case('euclidean',         k = 4,  pts = pts),
  case('bray, n = 2',       k = 3,  counts = rand_counts(2)),
  case('bray, n = 3',       k = 3,  counts = rand_counts(3)),
  # Duplicate samples: distances of 0.
  case('bray, duplicates',  k = 3,  counts = dups),
  case('bray, n = 300',     k = 3,  counts = rand_counts(300)) )

stopifnot(!is.null(cases[[1]]$warn))

expected <- list(

  versions = list(
    R     = R.version.string,
    ape   = as.character(packageVersion('ape')),
    vegan = as.character(packageVersion('vegan')) ),

  # Input: Bray-Curtis distances between the rows of `counts`.
  counts = counts,

  # ape::pcoa(d): first three axes. Bray-Curtis has negative eigenvalues.
  ape = list(
    vectors      = matrix(as.vector(ap$vectors[, 1:3]), nrow(ap$vectors)),
    eig          = ap$values$Eigenvalues[1:3],
    relative_eig = ap$values$Relative_eig[1:3],
    min_eig      = min(ap$values$Eigenvalues),
    trace        = ap$trace ),

  # vegan::wcmdscale(d, k = 3, eig = TRUE), and its "Proportion Explained".
  vegan = list(
    points = matrix(as.vector(vg$points), nrow(vg$points)),
    eig    = vg$eig[1:3],
    prop   = unname(summary(vegan::eigenvals(vg))[2, 1:3]) ),

  # The cases above.
  cases = cases
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/pcoa.rds', version = 2)
