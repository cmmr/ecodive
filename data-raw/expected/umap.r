# Generates tests/tinytests/expected/umap.rds: uwot::umap()'s results for
# test-umap.r, so that the tests don't need uwot.
#
# Run from the package root:  Rscript data-raw/expected/umap.r
#
# Needs uwot, Matrix, RSpectra and vegan. Every stage that uwot computes
# deterministically is recorded; so are full layouts by uwot's batch mode
# (batch = TRUE, which ecodive's umap() always uses) from given starting
# coordinates and seeds, and quality measures of a default-settings layout.
# In batch mode, uwot draws two 32-bit seeds per vertex per epoch from R's
# RNG, each as static_cast<uint32_t>(runif(1) * UINT32_MAX); test-umap.r
# regenerates them from the seed recorded here.

library(uwot)

stopifnot(requireNamespace('Matrix'), requireNamespace('RSpectra'))


# A plain matrix, without dimnames or scale()'s attributes.
strip <- function (x) matrix(as.double(x), nrow(x))

# uwot's fuzzy graph, as 0-based column pointers, row indices and values.
fgraph <- function (u) {
  V <- methods::as(u$fgraph, 'CsparseMatrix')
  list(p = as.integer(V@p), i = as.integer(V@i), x = V@x)
}

# Unit eigenvectors of the normalized Laplacian for its 2nd to (k+1)th
# smallest eigenvalues, as uwot's spectral init before scaling, but to a
# tight tolerance. Each vector's largest absolute coordinate is positive.
spectral <- function (u, k) {
  L   <- uwot:::form_normalized_laplacian(u$fgraph)
  ev  <- RSpectra::eigs_sym(L, k + 1, which = 'SM', opts = list(tol = 1e-14))
  vec <- uwot:::sort_eigenvectors(ev, k)
  for (j in seq_len(k)) {
    v <- vec[, j]
    if (v[[which.max(abs(v))]] < 0) vec[, j] <- -v
  }
  list(vectors = vec, values = sort(ev$values)[seq_len(k) + 1])
}


# Layout quality (the same functions are in test-umap.r).

# Fraction of each sample's m nearest neighbours in `d` that are among its
# m nearest in the embedding `x`.
nn_preservation <- function (d, x, m = 10) {
  dm <- as.matrix(d)
  de <- as.matrix(dist(x))
  mean(sapply(seq_len(nrow(dm)), function (i)
    length(intersect(order(dm[, i])[2:(m + 1)], order(de[, i])[2:(m + 1)])) / m ))
}

# Trustworthiness (Venna & Kaski 2001) of the embedding's m nearest
# neighbours.
trustworthiness <- function (d, x, m = 10) {
  dm <- as.matrix(d)
  de <- as.matrix(dist(x))
  n  <- nrow(dm)
  s  <- sum(sapply(seq_len(n), function (i) {
    rk <- integer(n)
    rk[order(dm[, i])] <- 0:(n - 1)          # 0 = self
    nn <- order(de[, i])[2:(m + 1)]
    sum(pmax(rk[nn] - m, 0))
  }))
  1 - 2 * s / (n * m * (2 * n - 3 * m - 1))
}

# UMAP's fuzzy cross-entropy between the graph (0-based column pointers p,
# rows i, values x) and the embedding's similarities 1 / (1 + a d^(2b)),
# over all pairs.
cross_entropy <- function (graph, x, ab) {
  n <- nrow(x)
  v <- matrix(0, n, n)
  v[cbind(graph$i + 1, rep(seq_len(n), diff(graph$p)))] <- graph$x
  q <- 1 / (1 + ab[[1]] * as.matrix(dist(x))^(2 * ab[[2]]))
  q <- pmin(pmax(q, 1e-12), 1 - 1e-12)
  ce <- -(v * log(q) + (1 - v) * log(1 - q))
  sum(ce[lower.tri(ce)])
}


#--------------------------------------------------
# Case 1: three overlapping Gaussian clusters, Euclidean,
# with a connected graph (for the spectral start).
#--------------------------------------------------
set.seed(1)
n  <- 100
X  <- rbind(
  matrix(rnorm(40 * 5),           ncol = 5),
  matrix(rnorm(30 * 5, mean = 1.5), ncol = 5),
  matrix(rnorm(30 * 5, mean = c(0, 3)), ncol = 5, byrow = TRUE) )
d1 <- dist(X)
attr(d1, 'call') <- NULL

u1 <- umap(
  d1, n_neighbors = 15, n_epochs = 0, init = 'random',
  ret_extra = c('sigma', 'nn', 'fgraph'), verbose = FALSE )
ab1 <- uwot:::find_ab_params(spread = 1, min_dist = 0.01)

# Batch-mode layouts from given starting coordinates and seed.
set.seed(2)
init1 <- matrix(runif(n * 2, -10, 10), n)
sgd1  <- list(
  init = init1, seed = 3, n_epochs = 100,
  embedding = strip(umap(
    d1, n_neighbors = 15, init = init1, seed = 3, n_epochs = 100,
    batch = TRUE, n_sgd_threads = 0, verbose = FALSE )))

# The same, but only 5 epochs: too few for rounding differences between
# platforms to grow, so it can be compared tightly everywhere.
sgd1$short <- strip(umap(
  d1, n_neighbors = 15, init = init1, seed = 3, n_epochs = 5,
  batch = TRUE, n_sgd_threads = 0, verbose = FALSE ))

# Non-default settings, 3 dimensions.
set.seed(4)
init2 <- matrix(runif(n * 3, -10, 10), n)
args2 <- list(
  n_neighbors = 8, min_dist = 0.3, spread = 1.5, learning_rate = 0.5,
  negative_sample_rate = 3, repulsion_strength = 2, n_epochs = 60 )
sgd2 <- list(
  init = init2, seed = 5, args = args2,
  embedding = strip(do.call(umap, c(list(
    X = d1, n_components = 3, init = init2, seed = 5, batch = TRUE,
    n_sgd_threads = 0, verbose = FALSE), args2))) )
sgd2$ab <- uwot:::find_ab_params(spread = 1.5, min_dist = 0.3)

# Quality of default-settings layouts from uwot's spectral start: batch mode
# (what ecodive does) and uwot's default in-place SGD, for comparison.
set.seed(6)
init_q <- strip(uwot:::spectral_init(u1$fgraph, ndim = 2))
quality <- function (x) c(
  preservation    = nn_preservation(d1, x),
  trustworthiness = trustworthiness(d1, x),
  cross_entropy   = cross_entropy(fgraph(u1), x, ab1) )
u_batch  <- umap(d1, init = init_q, seed = 7, batch = TRUE, verbose = FALSE)
u_serial <- umap(d1, init = init_q, seed = 7, verbose = FALSE)

case1 <- list(
  d           = d1,
  n_neighbors = 15,
  nn_idx      = u1$nn[[1]]$idx,
  nn_dist     = u1$nn[[1]]$dist,
  sigma       = u1$sigma,
  rho         = u1$rho,
  graph       = fgraph(u1),
  ab          = ab1,
  spectral    = spectral(u1, 2),
  sgd1        = sgd1,
  sgd2        = sgd2,
  init_q      = init_q,
  quality     = rbind(batch = quality(u_batch), serial = quality(u_serial)) )

print(case1$quality)


#--------------------------------------------------
# Case 2: Bray-Curtis of small integer counts, with
# duplicate samples (zero distances) and many ties,
# fuzzy union mixed with intersection, and
# local_connectivity between neighbours.
#--------------------------------------------------
set.seed(7)
counts <- matrix(rpois(40 * 6, 3), 40)
counts[c(5, 17, 30), ] <- counts[rep(3, 3), ]
counts[22, ] <- counts[9, ]
d2 <- vegan::vegdist(counts)
attr(d2, 'call')    <- NULL
attr(d2, 'maxdist') <- NULL

u2 <- umap(
  d2, n_neighbors = 10, n_epochs = 0, init = 'random',
  set_op_mix_ratio = 0.5, local_connectivity = 1.5,
  ret_extra = c('sigma', 'nn', 'fgraph'), verbose = FALSE )

case2 <- list(
  d                  = d2,
  n_neighbors        = 10,
  set_op_mix_ratio   = 0.5,
  local_connectivity = 1.5,
  nn_idx             = u2$nn[[1]]$idx,
  nn_dist            = u2$nn[[1]]$dist,
  sigma              = u2$sigma,
  rho                = u2$rho,
  graph              = fgraph(u2) )


expected <- list(
  case1    = case1,
  case2    = case2,
  versions = list(
    uwot     = as.character(packageVersion('uwot')),
    RSpectra = as.character(packageVersion('RSpectra')),
    Matrix   = as.character(packageVersion('Matrix')),
    vegan    = as.character(packageVersion('vegan')),
    R        = R.version.string ))

saveRDS(expected, 'tests/tinytests/expected/umap.rds', version = 2)

cat('Wrote tests/tinytests/expected/umap.rds:',
    file.size('tests/tinytests/expected/umap.rds'), 'bytes\n')
