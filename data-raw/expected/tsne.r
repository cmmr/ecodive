# Expected results for tests/tinytests/test-tsne.r. Run from the package
# root:
#
#     Rscript data-raw/expected/tsne.r
#
# Writes tests/tinytests/expected/tsne.rds, a list with one element per
# dataset ("blobs" and "bray"):
#
#   * scikit-learn 1.9.1's reference values, read from
#     data-raw/expected/tsne/<dataset>.txt. data-raw/expected/tsne.py wrote
#     those files (see its docstring for each value); there is no Python
#     here, so the committed text files are the source of truth.
#   * `rtsne`: the quality of Rtsne::Rtsne()'s layouts of the same
#     distances, for the same perplexity (`qperp`), with exact gradients
#     (theta = 0) and Rtsne's other defaults, from random starts drawn
#     with set.seed(0), ..., set.seed(9): the final KL divergence (`kl`),
#     trustworthiness with 5 neighbors (`trust`), and leave-one-out
#     5-nearest-neighbor label accuracy (`knn`). The functions for the last
#     two are the same as in test-tsne.r.
#
# Needs Rtsne and ecodive (for bray()).

library(ecodive)
stopifnot(requireNamespace('Rtsne'))


# One "name v1 v2 ..." line per value.
fx_read <- function (name) {
  x <- strsplit(readLines(file.path('data-raw', 'expected', 'tsne', paste0(name, '.txt'))), ' ', fixed = TRUE)
  structure(lapply(x, function (v) as.numeric(v[-1])), names = vapply(x, `[`, '', 1))
}

# Trustworthiness with m neighbors, as sklearn.manifold.trustworthiness().
trust <- function (d, y, m = 5) {
  D <- as.matrix(d); diag(D) <- Inf
  E <- as.matrix(dist(y)); diag(E) <- Inf
  n <- nrow(D)
  rk <- t(apply(D, 1, function (v) rank(v, ties.method = 'first')))
  t  <- 0
  for (i in seq_len(n)) {
    r <- rk[i, order(E[i, ])[seq_len(m)]] - m
    t <- t + sum(r[r > 0])
  }
  1 - t * 2 / (n * m * (2 * n - 3 * m - 1))
}

# Leave-one-out m-nearest-neighbor label accuracy on the plot.
knn_acc <- function (y, labels, m = 5) {
  E <- as.matrix(dist(y)); diag(E) <- Inf
  mean(vapply(seq_along(labels), function (i) {
    tb <- table(labels[order(E[i, ])[seq_len(m)]])
    as.numeric(names(tb)[which.max(tb)]) == labels[[i]]
  }, logical(1)))
}

dist_of <- function (v) {
  n <- (1 + sqrt(1 + 8 * length(v))) / 2
  structure(v, Size = as.integer(n), Diag = FALSE, Upper = FALSE, class = 'dist')
}


expected <- list(blobs = fx_read('blobs'), bray = fx_read('bray'))

for (name in names(expected)) {

  fx <- expected[[name]]

  # The distances that test-tsne.r gives tsne() for its quality checks.
  d <- if (name == 'bray')
    bray(matrix(fx$counts, nrow = fx$counts_dim[[1]], byrow = TRUE))
  else
    dist_of(fx$dist)

  runs <- lapply(0:9, function (seed) {
    set.seed(seed)
    r <- Rtsne::Rtsne(
      as.matrix(d), is_distance = TRUE, pca = FALSE, dims = 2,
      perplexity = fx$qperp, theta = 0, max_iter = 1000 )
    c(kl    = r$itercosts[[length(r$itercosts)]],
      trust = trust(d, r$Y),
      knn   = knn_acc(r$Y, fx$labels))
  })
  runs <- do.call(rbind, runs)

  expected[[name]]$rtsne <- list(
    kl    = unname(runs[, 'kl']),
    trust = unname(runs[, 'trust']),
    knn   = unname(runs[, 'knn']) )
}

str(lapply(expected, `[[`, 'rtsne'))

saveRDS(expected, 'tests/tinytests/expected/tsne.rds', version = 2)
