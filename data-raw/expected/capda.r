# Writes tests/tinytests/expected/capda.rds: inputs and reference results for
# test-capda.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/capda.r
#
# References (GPL; outputs only, never their source):
#
#   BiodiversityR 2.18-1's CAPdiscrim() (with vegan 2.7-8 and MASS
#   7.3-60.0.1): for m = 0 (chosen), mmax = 10, the chosen m, `percent`,
#   `percent.level`, the leave-one-out classes `CV`, `F` (MASS's svd^2, from
#   which delta^2 = F (g - 1) / (F (g - 1) + n - g)), Pillai's trace from
#   `manova`, and the scores `x`; and `percent` for every fixed m = 1, 2, ...
#   (the success for each m, and the choice of m beyond mmax = 10).
#
#   The leave-one-out classes for every m by refitting without each sample
#   (exact_loo() below), in the log domain: the success for every m, where
#   MASS's posteriors underflow.
#
#   MASS::lda(CV = TRUE) on cmdscale()'s axes, for a grouping with a single
#   sample in one group, which CAPdiscrim() rejects: the leave-one-out
#   classes for each m, NA for the singleton.
#
#   The permutation test of the trace statistic (sum of delta^2, Pillai's
#   trace for a MANOVA of the first m PCoA axes): with a shared matrix of
#   permutations, Pillai's trace from summary(manova()) of cmdscale()'s
#   first m axes for each permutation of the groups, and vegan::adonis2() of
#   the Euclidean distances of those axes scaled to unit length (whose F is
#   an increasing function of the trace, so its p-value is the same, given
#   the inverse permutations). The groups of permutation p are groups[p], as
#   in vegan's mrpp() and ecodive's capda(). And,
#   for a statistical comparison, adonis2()'s p-value from 9999 random
#   permutations within strata.

suppressMessages({
  library(BiodiversityR)
  library(vegan)
  library(MASS)
  library(permute)
})

data(dune, dune.env)
data(mite, mite.env)

bc <- function (x) vegdist(x, 'bray')

cap <- function (dm, groups, m, mmax = 10) {
  # CAPdiscrim() looks up the formula's left-hand side in the global
  # environment.
  env <- data.frame(G = factor(groups, ordered = FALSE))
  assign('cap_dm', dm, envir = globalenv())
  suppressWarnings(CAPdiscrim(cap_dm ~ G, data = env, axes = 2, m = m, mmax = mmax))
}

# Leave-one-out classes by refitting without each sample in turn, for the
# first m columns of `pts`, for every m: the group means and pooled
# covariance of the other samples (divisor n - 1 - g), with the priors of
# all n, and the class with the largest log score. CAPdiscrim() (that is,
# MASS::lda(CV = TRUE)) agrees, but where a sample is far from every group:
# for mite by Shrub, sample 67's scores are near -800 for m >= 20, its
# posterior underflows to 0 / 0 = NaN in MASS, and its class there is
# arbitrary (16 of 37 values of m differ by that one sample).
exact_loo <- function (pts, g) {
  n   <- nrow(pts)
  lev <- levels(g)
  pr  <- as.vector(table(g)) / n
  sapply(seq_len(ncol(pts)), function (m) vapply(seq_len(n), function (i) {
    x  <- pts[-i, seq_len(m), drop = FALSE]
    gi <- g[-i]
    nk <- as.vector(table(gi))
    if (any(nk == 0)) return (NA_character_)
    mu <- rowsum(x, gi) / nk
    S  <- crossprod(x - mu[as.integer(gi), , drop = FALSE]) / (n - 1 - length(lev))
    sc <- vapply(seq_along(lev), function (k) {
      d <- pts[i, seq_len(m)] - mu[k, ]
      log(pr[k]) - sum(d * solve(S, d)) / 2
    }, numeric(1))
    lev[[which.max(sc)]]
  }, character(1)))
}

# CAPdiscrim()'s results for the chosen m (mmax = 10), and its success for
# every fixed m up to `max_m`.
cap_case <- function (name, dm, groups, max_m) {

  groups <- factor(groups, ordered = FALSE)
  g      <- nlevels(groups)
  n      <- length(groups)
  r      <- cap(dm, groups, 0)
  k      <- min(g - 1, r$m)
  F      <- r$F[seq_len(k)]

  by_m <- sapply(seq_len(max_m), function (m) cap(dm, groups, m)$percent)

  # CAPdiscrim()'s MANOVA of one axis has no Pillai's trace: it is then the
  # R^2 of the axis on the groups.
  pillai <- if (r$m > 1) r$manova$stats[1, 2] else
    summary(lm(cmdscale(dm, k = 1) ~ groups))$r.squared

  pts   <- suppressWarnings(cmdscale(dm, k = max_m))
  exact <- exact_loo(pts, groups)

  list(
    name   = name,
    dm     = dm,
    groups = groups,
    expected = list(
      m             = r$m,
      percent       = r$percent,
      percent_level = r$percent.level,
      cv            = as.character(r$CV),
      F             = F,
      delta2        = F * (g - 1) / (F * (g - 1) + n - g),
      pillai        = pillai,
      x             = unname(r$x[, seq_len(k), drop = FALSE]),
      by_m          = by_m,
      by_m_exact    = colMeans(!is.na(exact) & exact == as.character(groups)) * 100 ) )
}

dune_dm <- bc(dune)
mite_dm <- bc(mite)

# Positive PCoA axes: 14 for dune, 37 for mite.
n_pos <- function (dm) sum(cmdscale(dm, k = attr(dm, 'Size') - 1, eig = TRUE)$eig > 1e-8)

cases <- list(
  cap_case('dune_management', dune_dm, dune.env$Management, 14),
  cap_case('dune_use',        dune_dm, dune.env$Use,        14),
  cap_case('dune_moisture',   dune_dm, dune.env$Moisture,   14),
  cap_case('mite_shrub',      mite_dm, mite.env$Shrub,      min(50, n_pos(mite_dm))) )


# A singleton group: MASS::lda(CV = TRUE) on cmdscale()'s first m axes.
single_groups <- factor(c('a', rep(c('b', 'c'), length.out = 19)))
pco <- suppressWarnings(cmdscale(dune_dm, k = 19, eig = TRUE))$points
single_cv <- sapply(1:14, function (m) {
  cv <- lda(pco[, 1:m, drop = FALSE], single_groups, CV = TRUE)$class
  as.character(cv)
})


# The permutation test of the trace, dune by Use at m = 2 (the m chosen),
# where the p-value is moderate.
set.seed(20261119)
perms <- shuffleSet(20, 199)
perm_m <- 2L
pts    <- cmdscale(dune_dm, k = perm_m)
grp    <- factor(dune.env$Use, ordered = FALSE)

pillai <- function (grp) summary(manova(pts ~ grp))$stats[1, 2]

trace_obs  <- pillai(grp)
trace_perm <- apply(perms, 1, function (p) pillai(grp[p]))
trace_p    <- (sum(trace_perm >= trace_obs - sqrt(.Machine$double.eps)) + 1) / (nrow(perms) + 1)

q   <- sweep(pts, 2, sqrt(colSums(pts^2)), '/')
# adonis2() permutes the rows of the response, the same as giving the
# groups the inverse permutation.
# (Kept as a shuffleSet() object: adonis2() of vegan 2.7-8 with permute
# 0.10-1 rejects a plain matrix.)
inv    <- perms
inv[]  <- t(apply(perms, 1, order))
ado <- adonis2(dist(q) ~ grp, permutations = inv)

# Within strata (Management), 9999 random permutations.
set.seed(2)
strata_p <- adonis2(
  dist(q) ~ grp,
  permutations = how(blocks = dune.env$Management, nperm = 9999) )[['Pr(>F)']][[1]]


expected <- list(

  versions = list(
    R             = R.version.string,
    BiodiversityR = as.character(packageVersion('BiodiversityR')),
    vegan         = as.character(packageVersion('vegan')),
    MASS          = as.character(packageVersion('MASS')),
    permute       = as.character(packageVersion('permute')) ),

  cases = cases,

  single = list(
    dm     = dune_dm,
    groups = single_groups,
    cv     = single_cv ),

  perm = list(
    dm         = dune_dm,
    groups     = grp,
    strata     = dune.env$Management,
    m          = perm_m,
    perms      = matrix(as.integer(perms), nrow(perms)),
    trace_obs  = trace_obs,
    trace_perm = trace_perm,
    p_value    = trace_p,
    adonis_p   = ado[['Pr(>F)']][[1]],
    strata_p   = strata_p ) )

saveRDS(expected, file.path('tests', 'tinytests', 'expected', 'capda.rds'), version = 2)
