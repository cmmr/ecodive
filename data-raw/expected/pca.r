# Generates tests/tinytests/expected/pca.rds: stats::prcomp() results that
# pca() must match, so that the tests need no reference implementation.
#
# Run from the package root:  Rscript data-raw/expected/pca.r

set.seed(20261116)

rand_counts <- function (n, p) {
  m <- matrix(rpois(n * p, 4) * rbinom(n * p, 1, 0.35), n, p)
  m[, 1] <- m[, 1] + 1L
  storage.mode(m) <- 'integer'
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('F', seq_len(p)))
  m
}

# The normalizations that pca() offers, in base R.
transform <- function (m, norm, pseudocount) {
  m <- m + 0
  switch(
    EXPR    = norm,
    none    = m,
    percent = m / rowSums(m),
    chord   = m / sqrt(rowSums(m^2)),
    binary  = (m > 0) + 0,
    clr     = { l <- log(m + pseudocount); l - rowMeans(l) },
    rclr    = {
      for (i in seq_len(nrow(m))) {
        x <- m[i, ]
        m[i, x > 0] <- log(x[x > 0]) - mean(log(x[x > 0]))
      }
      m } )
}

wide <- rand_counts(25, 150) # n - 1 <= p: sample-side solver
tall <- rand_counts(120, 20) # p < n - 1:  feature-side solver
tall <- tall[, colSums(tall) > 0]

cases <- list(
  list(counts = wide, norm = 'none',    pseudocount = NULL, scale = FALSE, k = 3L),
  list(counts = wide, norm = 'clr',     pseudocount = 1,    scale = FALSE, k = 5L),
  list(counts = wide, norm = 'percent', pseudocount = NULL, scale = TRUE,  k = 3L),
  list(counts = wide, norm = 'rclr',    pseudocount = NULL, scale = FALSE, k = 3L),
  list(counts = tall, norm = 'none',    pseudocount = NULL, scale = FALSE, k = 3L),
  list(counts = tall, norm = 'clr',     pseudocount = 0.5,  scale = TRUE,  k = 4L),
  list(counts = tall, norm = 'chord',   pseudocount = NULL, scale = FALSE, k = 3L),
  list(counts = tall, norm = 'binary',  pseudocount = NULL, scale = FALSE, k = 2L) )


# More cases, each also with:
#   sparse  Given to pca() as a features x samples dgCMatrix, margin = 2.
#   k_ret   The axes that pca() returns: fewer than k when the centered
#           table's rank is less than k, with warning `warn`.

set.seed(20261117)

# Constant features: absent everywhere, and the same count everywhere.
with_const <- function (m) {
  m[, 2] <- 0L
  m[, 3] <- 5L
  m
}

# Integer counts whose centered table has rank 3: an offset plus 3 profiles.
rank3 <- function (n, p) {
  w <- cbind(1L, matrix(sample(0:3, n * 3, TRUE), n, 3))
  m <- w %*% matrix(sample(0:4, 4 * p, TRUE), 4, p)
  storage.mode(m) <- 'integer'
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('F', seq_len(p)))
  m
}

wide2 <- rand_counts(25, 150)
tall2 <- rand_counts(120, 20)

more_cases <- list(
  list(counts = wide2,             norm = 'clr',     pseudocount = 1,    scale = FALSE, k = 3L, sparse = TRUE),
  list(counts = tall2,             norm = 'percent', pseudocount = NULL, scale = FALSE, k = 3L, sparse = TRUE),
  list(counts = with_const(wide2), norm = 'none',    pseudocount = NULL, scale = TRUE,  k = 3L, sparse = FALSE),
  list(counts = with_const(tall2), norm = 'none',    pseudocount = NULL, scale = TRUE,  k = 4L, sparse = FALSE),
  list(counts = with_const(tall2), norm = 'chord',   pseudocount = NULL, scale = TRUE,  k = 3L, sparse = TRUE),
  list(counts = rank3(30, 40),     norm = 'none',    pseudocount = NULL, scale = FALSE, k = 5L, sparse = FALSE),
  list(counts = rank3(60, 10),     norm = 'none',    pseudocount = NULL, scale = FALSE, k = 5L, sparse = FALSE) )


expect_prcomp <- function (cs) {

  x <- transform(cs$counts, cs$norm, cs$pseudocount)

  # prcomp(scale. = TRUE) stops on constant features; pca() leaves them out.
  if (cs$scale) x <- x[, apply(x, 2L, function (v) any(v != v[1]))]

  pr <- stats::prcomp(x, scale. = cs$scale)

  # Axes with positive variances, as pca() keeps them.
  k_ret <- min(cs$k, sum(pr$sdev^2 > sqrt(.Machine$double.eps) * pr$sdev[1]^2))
  j     <- seq_len(k_ret)

  list(
    x        = pr$x[, j, drop = FALSE],
    rotation = pr$rotation[, j, drop = FALSE],
    sdev     = pr$sdev )
}

for (i in seq_along(cases))
  cases[[i]]$expected <- expect_prcomp(cases[[i]])

for (i in seq_along(more_cases)) {
  cs    <- more_cases[[i]]
  ex    <- expect_prcomp(cs)
  k_ret <- ncol(ex$x)
  more_cases[[i]]$k_ret    <- k_ret
  more_cases[[i]]$warn     <- if (k_ret < cs$k) sprintf('only %d of the first %d eigenvalues', k_ret, cs$k)
  more_cases[[i]]$expected <- ex
}

# Constant features are present, and two cases have rank 3.
stopifnot(
  sapply(more_cases[3:5], function (cs) nrow(cs$expected$rotation) < ncol(cs$counts)),
  sapply(more_cases[6:7], function (cs) cs$k_ret == 3L) )

saveRDS(
  file    = 'tests/tinytests/expected/pca.rds',
  version = 2,
  object  = list(
    cases      = cases,
    more_cases = more_cases,
    generated = list(
      reference = 'stats::prcomp',
      versions  = list(stats = as.character(packageVersion('stats'))),
      R         = R.version.string )))
