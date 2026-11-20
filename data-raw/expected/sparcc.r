# Writes tests/tinytests/expected/sparcc.rds: inputs and reference results
# for test-sparcc.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/sparcc.r
#
# SparCC's random part is the Dirichlet draws. Given the same log
# fractions (one n x D matrix per draw), the rest is deterministic, and is
# computed here with a plain base-R implementation of Friedman & Alm's
# (2012) steps, written from the paper and plans/sparcc.md: the variation
# matrix from var() of each log ratio, the basis variances from a dense
# solve() of the system (not ecodive's closed form and Woodbury updates),
# the exclusion loop, and the medians over draws from median(). The tests
# match ecodive's sparcc_details() to it at 1e-12.
#
# FastSpar and SpiecEasi (GPL-3; output only) could not be run when these
# fixtures were made (see plans/sparcc.md), so no outputs of theirs are
# stored.


# One draw: L is n x D log abundances. Returns the packed rho (NA for pairs
# of a dropped component, NaN for a non-positive variance), the variances
# (NA when dropped), and the excluded pairs' packed indices in order.
sparcc_draw_ref <- function (L, max_exclusions, threshold) {

  D <- ncol(L)
  T <- matrix(0, D, D)
  for (i in seq_len(D)) for (j in seq_len(D)) if (i != j) T[i, j] <- var(L[, i] - L[, j])

  # Pairs (i < j) in the order of a dist's packed lower triangle.
  pi <- unlist(lapply(seq_len(D - 1L), function (i) rep(i, D - i)))
  pj <- unlist(lapply(seq_len(D - 1L), function (i) (i + 1L):D))
  P  <- length(pi)

  active   <- rep(TRUE, D)
  excluded <- rep(FALSE, P)
  order_ex <- integer(0)

  variances <- function () {
    keep <- matrix(TRUE, D, D)
    keep[cbind(pi[excluded], pj[excluded])] <- FALSE
    keep[cbind(pj[excluded], pi[excluded])] <- FALSE
    diag(keep) <- FALSE
    idx <- which(active)
    K   <- keep[idx, idx, drop = FALSE]
    # t_i = sum over kept j of t_ij = deg_i w_i + sum over kept j of w_j
    M   <- diag(rowSums(K), length(idx)) + K * 1
    ti  <- rowSums(T[idx, idx, drop = FALSE] * K)
    if (rcond(M) < 1e-12) return (NULL)
    w <- rep(NA_real_, D)
    w[idx] <- solve(M, ti)
    w
  }

  rho_of <- function (w) {
    r <- rep(NA_real_, P)
    for (k in seq_len(P)) {
      i <- pi[k]; j <- pj[k]
      if (!active[i] || !active[j]) next
      r[k] <- if (w[i] > 0 && w[j] > 0)
        (w[i] + w[j] - T[i, j]) / (2 * sqrt(w[i]) * sqrt(w[j])) else NaN
    }
    r
  }

  w    <- variances()
  rho  <- rho_of(w)
  stop <- 'max_exclusions'

  for (it in seq_len(max_exclusions)) {

    a <- abs(rho)
    a[excluded | is.na(a)] <- -Inf
    if (!any(a > threshold)) { stop <- 'threshold'; break }
    k <- which.max(a)                     # the first of any ties

    excluded[k] <- TRUE
    order_ex    <- c(order_ex, k)
    for (f in c(pi[k], pj[k])) {
      partners <- setdiff(which(active), f)
      done <- vapply(partners, function (g) {
        excluded[which((pi == min(f, g)) & (pj == max(f, g)))]
      }, logical(1))
      if (all(done)) active[f] <- FALSE
    }
    if (sum(active) < 4) { stop <- 'three'; break }

    w_new <- variances()
    if (is.null(w_new)) {                 # singular: undo, and stop
      excluded[k] <- FALSE
      order_ex    <- order_ex[-length(order_ex)]
      active[c(pi[k], pj[k])] <- TRUE
      stop <- 'singular'
      break
    }
    w   <- w_new
    rho <- rho_of(w)
  }

  w[!active] <- NA_real_
  rho[!(active[pi] & active[pj])] <- NA_real_
  list(rho = rho, variances = w, exclusions = order_ex, stop = stop)
}


sparcc_ref <- function (draws, max_exclusions, threshold) {
  res <- lapply(draws, sparcc_draw_ref, max_exclusions, threshold)
  rho <- do.call(rbind, lapply(res, `[[`, 'rho'))
  w   <- do.call(rbind, lapply(res, `[[`, 'variances'))
  list(
    draw_rho       = rho,
    draw_variances = w,
    exclusions     = lapply(res, `[[`, 'exclusions'),
    stop           = vapply(res, `[[`, '', 'stop'),
    rho            = apply(rho, 2, stats::median, na.rm = TRUE),
    variances      = apply(w,   2, stats::median, na.rm = TRUE) )
}


# Log fractions of n_draws Dirichlet(counts + 1) draws, with R's rgamma().
dirichlet_logs <- function (counts, n_draws) {
  lapply(seq_len(n_draws), function (d)
    matrix(log(stats::rgamma(length(counts), counts + 1)), nrow(counts)) )
}


set.seed(20261120)
cases <- list()

# Correlated features among uncorrelated ones; the paper's settings.
n <- 25; D <- 12
base <- matrix(rnorm(n * D), n, D)
base[, 2] <- base[, 1] + rnorm(n, sd = 0.4)
base[, 5] <- -base[, 4] + rnorm(n, sd = 0.4)
counts <- matrix(rpois(n * D, exp(base + 4)), n, D)
cases$correlated <- list(
  counts = counts, draws = dirichlet_logs(counts, 5),
  max_exclusions = 20L, threshold = 0.1 )

# No exclusions: the closed-form variances.
cases$basic <- list(
  counts = counts, draws = cases$correlated$draws[1:3],
  max_exclusions = 0L, threshold = 0.1 )

# Pairs excluded down to an exact fit among 5 features: components are
# dropped, and the loop stops at 3 components or a singular system. The
# threshold is just above rounding error: once the pairs still included
# determine the variances exactly, their correlations are 0 up to rounding,
# which mustn't decide what is excluded next.
n <- 15; D <- 5
counts <- matrix(rpois(n * D, exp(matrix(rnorm(n * D, sd = 1.5), n, D) + 3)), n, D)
cases$drop <- list(
  counts = counts, draws = dirichlet_logs(counts, 6),
  max_exclusions = 20L, threshold = 1e-6 )

# Six features, likewise: exercises the 4-feature case.
n <- 12; D <- 6
counts <- matrix(rpois(n * D, exp(matrix(rnorm(n * D), n, D) + 3)), n, D)
cases$small <- list(
  counts = counts, draws = dirichlet_logs(counts, 6),
  max_exclusions = 20L, threshold = 1e-6 )

# 20 exclusions among 30 features: Woodbury against dense solves.
n <- 20; D <- 30
counts <- matrix(rpois(n * D, exp(matrix(rnorm(n * D), n, D) + 3)), n, D)
cases$many <- list(
  counts = counts, draws = dirichlet_logs(counts, 3),
  max_exclusions = 20L, threshold = 1e-6 )

# Ties: three identical features, whose three pairs share the largest
# |rho| exactly (the same arithmetic gives the same value), and then two
# of them; the first in packed order wins. (Ties between two separate
# groups of identical features would depend on rounding.)
n <- 10; D <- 8
L <- matrix(rnorm(n * D), n, D)
L[, 2] <- L[, 1]; L[, 3] <- L[, 1]
cases$ties <- list(
  counts = matrix(1L, n, D), draws = list(L, L + 0.5),
  max_exclusions = 5L, threshold = 0.1 )

# The reference results (6 draws in some cases: medians of the middle two).
for (nm in names(cases)) {
  cc <- cases[[nm]]
  cases[[nm]]$expected <- sparcc_ref(cc$draws, cc$max_exclusions, cc$threshold)
}



# The reference implementations' outputs (output only): FastSpar 1.0.0 and
# SpiecEasi 1.1.1, from data-raw/expected/sparcc/, written by
# sparcc_refs.r and sparcc_fastspar.sh. Correlations only (FastSpar prints
# 4 decimals), with FastSpar's covariance diagonals (the variances), and
# FastSpar's p-values for one run.
dir <- file.path('data-raw', 'expected', 'sparcc')
read_tsv <- function (file) as.matrix(utils::read.table(
  file.path(dir, file), sep = '\t', header = TRUE, row.names = 1,
  comment.char = '', check.names = FALSE ))
refs <- list(tables = lapply(c(tight = 'tight', moderate = 'moderate', d8 = 'd8'),
  function (nm) t(read_tsv(paste0(nm, '.tsv'))) ))   # samples x features
refs$fastspar <- list()
for (f in list.files(dir, '^fastspar_.*_cor[.]tsv$')) {
  key <- sub('^fastspar_(.*)_cor[.]tsv$', '\\1', f)
  refs$fastspar[[key]] <- list(
    cor       = read_tsv(f),
    variances = diag(read_tsv(sub('_cor', '_cov', f))) )
  pv <- sub('_cor', '_pvalues', f)
  if (file.exists(file.path(dir, pv))) refs$fastspar[[key]]$p_value <- read_tsv(pv)
}
refs$spieceasi <- list()
for (f in list.files(dir, '^spieceasi_.*[.]tsv$'))
  refs$spieceasi[[sub('^spieceasi_(.*)[.]tsv$', '\\1', f)]] <- read_tsv(f)

cases$refs <- refs
saveRDS(cases, file.path('tests', 'tinytests', 'expected', 'sparcc.rds'), version = 2)
