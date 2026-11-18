# Writes tests/tinytests/expected/plsda.rds: inputs and reference results
# for test-plsda.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/plsda.r      # writes the inputs
#     python data-raw/expected/plsda.py      # scikit-learn's results
#     Rscript data-raw/expected/plsda.r      # folds them in
#
# References:
#
#   sklearn   scikit-learn 1.9.1's PLSRegression (BSD-3), which plsda() ports,
#             fitted to the normalized table and the 0/1 dummy matrix of the
#             groups, from the CSV files that data-raw/expected/plsda.py
#             writes to data-raw/expected/plsda/. Run with tol = 1e-28 (so
#             that NIPALS converges fully) and also, for the two-group case,
#             with its default tol, where NIPALS converges at once.
#   mixomics  mixOmics::plsda() (GPL; outputs only, never its source),
#             fitted to the normalized table with tol = 1e-15: variates$X,
#             loadings$X (the weights), prop_expl_var$X and vip(). Generated
#             with mixOmics 6.26.0 (Bioconductor 3.18, the newest for R
#             4.3), as available where these fixtures were made.
#
# Each case holds:
#   name, counts, groups      The table (samples in rows) and its groups.
#   norm, pseudocount, scale  plsda()'s arguments; the references get the
#                             table after the same normalization, in base R.
#   k                         Components requested.
#   x                         The normalized table given to the references.
#   expected$sklearn, expected$sklearn_default, expected$mixomics
#                             Reference results (NULL where not run).

library(vegan)

data(dune, dune.env)

dir.create(skdir <- file.path('data-raw', 'expected', 'plsda'), showWarnings = FALSE)

dune <- as.matrix(dune)

set.seed(20261118)
n   <- 110
big <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.4), n, 40)
dimnames(big) <- list(sprintf('S%03d', seq_len(n)), sprintf('F%02d', 1:40))
site <- factor(sample(c('a', 'b', 'c'), n, TRUE))
big[site == 'b', 2:6] <- big[site == 'b', 2:6] + 2
big[site == 'c', 7:9] <- big[site == 'c', 7:9] + 3


# The normalizations plsda() offers, in base R.
transform <- function (m, norm, pseudocount = NULL) {
  m <- m + 0
  switch(
    EXPR      = norm,
    none      = m,
    hellinger = sqrt(m / rowSums(m)),
    clr       = { l <- log(m + pseudocount); l - rowMeans(l) } )
}

cs <- function (name, counts, groups, k, norm = 'none', scale = TRUE,
                pseudocount = NULL, mixomics = TRUE)
  list(name = name, counts = counts, groups = groups, norm = norm,
       pseudocount = pseudocount, scale = scale, k = k,
       x = transform(counts, norm, pseudocount), run_mixomics = mixomics)

two <- factor(ifelse(dune.env$Management == 'NM', 'natural', 'managed'))

cases <- list(
  cs('dune_mgmt',      dune, dune.env$Management, 3),
  cs('dune_two',       dune, two, 2),
  cs('dune_const',     cbind(dune, Const = 3, Zero = 0), dune.env$Management, 3, mixomics = FALSE),
  cs('dune_hellinger', dune, dune.env$Management, 3, norm = 'hellinger', scale = FALSE),
  cs('dune_clr',       dune, dune.env$Moisture, 4, norm = 'clr', pseudocount = 1),
  cs('big_site',       big,  site, 5) )


# Inputs for plsda.py.
manifest <- data.frame(
  name  = vapply(cases, `[[`, '', 'name'),
  scale = vapply(cases, `[[`, TRUE, 'scale'),
  k     = vapply(cases, `[[`, 0, 'k') )
write.csv(manifest, file.path(skdir, 'cases.csv'), row.names = FALSE)
for (c in cases) {
  write.csv(c$x, file.path(skdir, paste0(c$name, '_x.csv')))
  write.csv(data.frame(group = as.character(c$groups), row.names = rownames(c$x)),
            file.path(skdir, paste0(c$name, '_groups.csv')))
}


# scikit-learn's results, if plsda.py has run.
rd <- function (f) {
  f <- file.path(skdir, f)
  if (!file.exists(f)) return (NULL)
  unname(as.matrix(read.csv(f, header = FALSE)))
}
sk <- function (name, suffix = '') {
  scores <- rd(paste0(name, '_sk', suffix, '_x_scores.csv'))
  if (is.null(scores)) return (NULL)
  list(
    scores     = scores,
    weights    = rd(paste0(name, '_sk', suffix, '_x_weights.csv')),
    loadings   = rd(paste0(name, '_sk', suffix, '_x_loadings.csv')),
    y_loadings = rd(paste0(name, '_sk', suffix, '_y_loadings.csv')) )
}


# mixOmics' results.
mix <- function (c) {
  if (!c$run_mixomics) return (NULL)
  suppressMessages(library(mixOmics))
  m <- mixOmics::plsda(c$x, c$groups, ncomp = c$k, scale = c$scale, tol = 1e-15, max.iter = 1e5)
  list(
    scores  = unname(as.matrix(m$variates$X)),
    weights = unname(as.matrix(m$loadings$X)),
    prop    = unname(m$prop_expl_var$X),
    vip     = unname(as.matrix(mixOmics::vip(m))) )
}


for (i in seq_along(cases)) {
  c <- cases[[i]]
  c$expected <- list(
    sklearn         = sk(c$name),
    sklearn_default = sk(c$name, '_default'),
    mixomics        = mix(c) )
  c$run_mixomics <- NULL
  cases[[i]] <- c
}

if (is.null(cases[[1]]$expected$sklearn))
  message('No scikit-learn results yet: run data-raw/expected/plsda.py, then this script again.')

saveRDS(
  list(cases = cases),
  file.path('tests', 'tinytests', 'expected', 'plsda.rds'),
  version = 2 )
