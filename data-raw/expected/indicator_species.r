# Expected results for tests/tinytests/test-indicator_species.r, from
# indicspecies 1.8.0 (multipatt(), strassoc()) and labdsv 2.3.1 (indval()).
# Both are GPL: only their outputs are used (clean room; see
# plans/indicator_species.md). Run from the package root, with ecodive
# installed (for its permutations):
#
#     Rscript data-raw/expected/indicator_species.r
#
# Writes tests/tinytests/expected/indicator_species.rds.
#
# The permutations are ecodive's own (permutation k of seed s, as
# indicator_species(permutations = 199, seed = s) draws them), given to
# multipatt() as a matrix (one permutation per row), so the tests can check
# that ecodive's p-values are identical with random and with supplied
# permutations. multipatt() breaks ties between combinations (within a
# relative 1e-5) at random, so `near_tie` flags the features whose two best
# combinations are that close: the tests compare their values, not their
# combinations or p-values.

library(indicspecies)
library(labdsv)

stopifnot(packageVersion('indicspecies') == '1.8.0')

set.seed(50)


#--------------------------------------------------------
# Tables: 40 samples x 30 features in three unequal groups
# (sorted names, so that the groups' order is the same
# everywhere), with features that live in one group, in
# two, everywhere, nowhere, or are the same everywhere.
#--------------------------------------------------------
n   <- 40L
p   <- 30L
grp <- sample(rep(c('forest', 'marsh', 'meadow'), c(8, 13, 19)))
lam <- rep(runif(p, 0.3, 4), each = n) * matrix(runif(n * p, 0, 2), n, p)
x   <- matrix(rpois(n * p, lam) * rbinom(n * p, 1, 0.6), n, p)
dimnames(x) <- list(sprintf('S%02d', seq_len(n)), sprintf('F%02d', seq_len(p)))
x[grp == 'forest', 1:3]  <- x[grp == 'forest', 1:3] + rpois(8 * 3, 4)
x[grp != 'forest', 1]    <- 0
x[grp == 'marsh',  4:5]  <- x[grp == 'marsh', 4:5] + rpois(13 * 2, 3)
x[grp != 'meadow', 6:8]  <- x[grp != 'meadow', 6:8] + rpois(21 * 3, 3)
x[, 9]  <- 0                       # absent everywhere
x[, 10] <- 3                       # the same in every sample
x[, 11] <- x[, 11] + 1             # present everywhere
storage.mode(x) <- 'double'
xpa  <- (x > 0) * 1
xrel <- x / rowSums(x)

# Four groups, for combination ids beyond three.
grp4 <- sample(rep(c('a', 'b', 'c', 'd'), c(7, 9, 11, 13)))
x4   <- x
x4[grp4 %in% c('a', 'c'), 12:13] <- x4[grp4 %in% c('a', 'c'), 12:13] + 4


#--------------------------------------------------------
# ecodive's permutations: free, and within strata.
#--------------------------------------------------------
n_perm <- 199L
seed   <- 7L
perms  <- t(.Call(ecodive:::C_perm_indices, n, n_perm, NULL, seed))
strata <- rep(1:4, each = 10)
sperms <- t(.Call(ecodive:::C_perm_indices, n, n_perm, strata - 1L, seed))
stopifnot(all(apply(sperms, 1, function (r) all(strata[r] == strata))))


#--------------------------------------------------------
# One multipatt() run, and what the tests need from it.
#--------------------------------------------------------
run <- function (data, cluster, func, ..., permutations = perms) {

  mp <- multipatt(data, cluster, func = func, control = NULL, permutations = permutations, ...)

  str  <- mp$str
  sc   <- str
  if (isTRUE(list(...)$allow.negative)) {
    # A combination and its complement tie in absolute value, and
    # multipatt() then picks the one of fewer groups, as documented: with
    # three groups, the single groups. Near ties are among those.
    sc <- abs(str)[, !grepl('+', colnames(str), fixed = TRUE), drop = FALSE]
  }
  near <- apply(sc, 1, function (s) {
    s <- sort(s[is.finite(s)], decreasing = TRUE)
    length(s) > 1 && s[[1]] - s[[2]] <= 1e-5 * abs(s[[1]])
  })

  idx <- mp$sign$index
  ab  <- function (m) if (is.null(m)) NULL else m[cbind(seq_len(nrow(m)), idx)]
  list(
    str      = unclass(str),
    comb     = colnames(mp$comb),
    best     = colnames(mp$comb)[idx],
    stat     = mp$sign$stat,
    p_value  = mp$sign$p.value,
    a        = ab(mp$A),
    b        = ab(mp$B),
    near_tie = unname(near) )
}

fixtures <- list(
  x = x, xpa = xpa, xrel = xrel, grp = grp, x4 = x4, grp4 = grp4,
  n_perm = n_perm, seed = seed, strata = strata )

cases <- list()
for (func in c('IndVal', 'IndVal.g', 'r', 'r.g')) {
  cases[[paste0(func, '/abund')]] <- run(x,    grp, func)
  cases[[paste0(func, '/pa')]]    <- run(xpa,  grp, func)
  cases[[paste0(func, '/rel')]]   <- run(xrel, grp, func)
  cases[[paste0(func, '/duleg')]] <- run(x,    grp, func, duleg = TRUE)
  cases[[paste0(func, '/strata')]] <- run(x,   grp, func, permutations = sperms)
}
# multipatt(allow.negative = TRUE) stops on a feature whose correlation is
# undefined (the same value in every sample: F09, F10, and in xpa those
# present everywhere), so they are left out; ecodive gives them NA.
vary <- function (m) m[, apply(m, 2, function (v) any(v != v[[1]]))]
fixtures$x_vary   <- vary(x)
fixtures$xpa_vary <- vary(xpa)
for (func in c('r', 'r.g')) {
  cases[[paste0(func, '/negative')]]    <- run(vary(x),   grp, func, allow.negative = TRUE)
  cases[[paste0(func, '/negative_pa')]] <- run(vary(xpa), grp, func, allow.negative = TRUE)
}
cases[['IndVal.g/max2']]     <- run(x, grp, 'IndVal.g', max.order = 2)
cases[['IndVal.g/restcomb']] <- run(x, grp, 'IndVal.g', restcomb = c(1, 3, 5, 6))
cases[['IndVal.g/restcomb_all']]  <- run(x, grp, 'IndVal.g', restcomb = c(1, 2, 7))
cases[['IndVal.g/restcomb_max2']] <- run(x, grp, 'IndVal.g', restcomb = c(1, 2, 4, 7), max.order = 2)
cases[['IndVal.g/four']]     <- run(x4, grp4, 'IndVal.g')
cases[['IndVal.g/four_min2']] <- run(x4, grp4, 'IndVal.g', min.order = 2, max.order = 3)
cases[['r.g/four_restcomb']] <- run(x4, grp4, 'r.g', restcomb = c(2, 5, 7, 11, 14))
fixtures$cases <- cases


#--------------------------------------------------------
# strassoc(): each index, and A and B, for single groups.
#--------------------------------------------------------
fixtures$strassoc <- lapply(
  c(IndVal = 'IndVal', IndVal.g = 'IndVal.g', r = 'r', r.g = 'r.g', A = 'A', A.g = 'A.g', B = 'B'),
  function (f) unclass(suppressWarnings(strassoc(x, grp, func = f))) )


#--------------------------------------------------------
# labdsv::indval(): the IndVal of Dufrene & Legendre (A B,
# group-equalized), single groups. Its permutations use R's
# random numbers, so only the values are compared.
#--------------------------------------------------------
lv <- labdsv::indval(x[, -c(9, 10)], match(grp, sort(unique(grp))), numitr = 1)
fixtures$labdsv <- list(
  features = colnames(x)[-c(9, 10)],
  indval   = unname(as.matrix(lv$indval)),
  maxcls   = unname(lv$maxcls),
  indcls   = unname(lv$indcls) )


#--------------------------------------------------------
# Probes that settled the clean-room questions (see the
# plan's "Decisions made while implementing").
#--------------------------------------------------------

# Ties: a permuted value below the observed one by a few ulps is not a tie
# for multipatt() (p = 0.5 with this one permutation, not 1).
pg <- rep(1:2, c(4, 4))
px <- sapply(10^-(15:2), function (d) c(5, 4, 3, 2, 2 - d, 0, 0, 0))
pp <- matrix(c(1, 2, 3, 5, 4, 6, 7, 8), 1)
fixtures$probe_ties <- sapply(c('IndVal', 'IndVal.g', 'r', 'r.g'), function (f)
  multipatt(px, pg, func = f, duleg = TRUE, control = NULL, permutations = pp)$sign$p.value)

dir.create('tests/tinytests/expected', showWarnings = FALSE, recursive = TRUE)
saveRDS(fixtures, 'tests/tinytests/expected/indicator_species.rds', version = 2)
