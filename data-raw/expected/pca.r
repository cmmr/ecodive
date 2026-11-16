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


for (i in seq_along(cases)) {

  cs <- cases[[i]]
  x  <- transform(cs$counts, cs$norm, cs$pseudocount)

  # prcomp(scale. = TRUE) stops on constant features; pca() leaves them out.
  if (cs$scale) x <- x[, apply(x, 2L, function (v) any(v != v[1]))]

  pr <- stats::prcomp(x, scale. = cs$scale)
  j  <- seq_len(cs$k)

  cases[[i]]$expected <- list(
    x        = pr$x[, j, drop = FALSE],
    rotation = pr$rotation[, j, drop = FALSE],
    sdev     = pr$sdev )
}

saveRDS(
  file    = 'tests/tinytests/expected/pca.rds',
  version = 2,
  object  = list(
    cases     = cases,
    generated = list(
      reference = 'stats::prcomp',
      versions  = list(stats = as.character(packageVersion('stats'))),
      R         = R.version.string )))
