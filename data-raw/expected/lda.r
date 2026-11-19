# Writes tests/tinytests/expected/lda.rds: inputs and reference results for
# test-lda.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/lda.r      # writes the inputs
#     python data-raw/expected/lda.py      # scikit-learn's results
#     Rscript data-raw/expected/lda.r      # folds them in
#
# References:
#
#   sklearn   scikit-learn 1.9.1's LinearDiscriminantAnalysis(solver = "svd")
#             (BSD-3), which lda() ports, fitted to the normalized table and
#             the groups (with the same priors), from the CSV files that
#             data-raw/expected/lda.py writes to data-raw/expected/lda/:
#             transform() scores, scalings_, explained_variance_ratio_,
#             means_ and xbar_. Its scores and scalings_ are larger than
#             MASS's (and lda()'s) by sqrt(n / (n - g)).
#   mass      MASS::lda() (GPL; outputs only, never its source), fitted to
#             the normalized table: predict()$x, scaling, svd, means and
#             prior; and MASS::lda(CV = TRUE) on the first m columns, for
#             every m: the leave-one-out classes, and the posterior for all
#             columns (lda_loo()'s reference). MASS 7.3-60.0.1, as shipped
#             with R 4.3.3.
#
# Each case holds:
#   name, counts, groups      The table (samples in rows) and its groups.
#   norm, pseudocount, prior  lda()'s arguments; the references get the table
#                             after the same normalization, in base R.
#   x                         The normalized table given to the references.
#   expected$sklearn, expected$mass
#                             Reference results.

library(vegan)
library(MASS)

data(dune, dune.env)
data(iris)

dir.create(skdir <- file.path('data-raw', 'expected', 'lda'), showWarnings = FALSE)

dune <- as.matrix(dune)

iris_x <- as.matrix(iris[, 1:4])
rownames(iris_x) <- sprintf('i%03d', seq_len(nrow(iris_x)))

# Sparse counts: 4 groups of unequal sizes (a singleton among them), with
# group effects on a few features.
set.seed(20261119)
n   <- 151
big <- matrix(rpois(n * 12, 4) * rbinom(n * 12, 1, 0.6), n, 12)
dimnames(big) <- list(sprintf('S%03d', seq_len(n)), sprintf('F%02d', 1:12))
site <- factor(c(sample(c('a', 'b', 'c'), n - 1, TRUE, prob = c(.5, .3, .2)), 'd'))
big[site == 'b', 2:4] <- big[site == 'b', 2:4] + 2
big[site == 'c', 5:6] <- big[site == 'c', 5:6] + 3
big[site == 'd', 7]   <- big[site == 'd', 7] + 5

# Principal coordinates of Bray-Curtis distances (the input of capda()).
pco <- cmdscale(vegdist(dune), k = 6, eig = TRUE)$points
dimnames(pco) <- list(rownames(dune), paste0('PCoA', 1:6))

# Eight common species (every site has some).
sp8 <- dune[, order(-colSums(dune > 0))[1:8]]
stopifnot(all(rowSums(sp8) > 0))


# The normalizations lda() offers, in base R.
transform <- function (m, norm, pseudocount = NULL) {
  m <- m + 0
  switch(
    EXPR      = norm,
    none      = m,
    percent   = m / rowSums(m) * 100,
    hellinger = sqrt(m / rowSums(m)) )
}

cs <- function (name, counts, groups, norm = 'none', pseudocount = NULL, prior = NULL)
  list(name = name, counts = counts, groups = factor(groups), norm = norm,
       pseudocount = pseudocount, prior = prior,
       x = transform(counts, norm, pseudocount))

cases <- list(
  cs('iris',           iris_x, iris$Species),
  cs('iris_prior',     iris_x, iris$Species, prior = c(setosa = 0.2, versicolor = 0.3, virginica = 0.5)),
  cs('iris_two',       iris_x[51:150, ], droplevels(iris$Species[51:150])),
  cs('dune_pcoa',      pco,    dune.env$Management),
  cs('dune_hellinger', sp8,    dune.env$Moisture, norm = 'hellinger'),
  cs('big_site',       big,    site) )


# Inputs for lda.py.
manifest <- data.frame(name = vapply(cases, `[[`, '', 'name'))
write.csv(manifest, file.path(skdir, 'cases.csv'), row.names = FALSE)
for (c in cases) {
  write.csv(c$x, file.path(skdir, paste0(c$name, '_x.csv')))
  write.csv(data.frame(group = as.character(c$groups), row.names = rownames(c$x)),
            file.path(skdir, paste0(c$name, '_groups.csv')))
  pr <- if (is.null(c$prior)) as.vector(table(c$groups)) / length(c$groups) else c$prior[levels(c$groups)]
  write.csv(data.frame(prior = unname(pr)), file.path(skdir, paste0(c$name, '_prior.csv')), row.names = FALSE)
}


# scikit-learn's results, if lda.py has run.
rd <- function (f) {
  f <- file.path(skdir, f)
  if (!file.exists(f)) return (NULL)
  unname(as.matrix(read.csv(f, header = FALSE)))
}
sk <- function (name) {
  scores <- rd(paste0(name, '_sk_scores.csv'))
  if (is.null(scores)) return (NULL)
  list(
    scores   = scores,
    scalings = rd(paste0(name, '_sk_scalings.csv')),
    evr      = as.vector(rd(paste0(name, '_sk_evr.csv'))),
    means    = rd(paste0(name, '_sk_means.csv')),
    xbar     = as.vector(rd(paste0(name, '_sk_xbar.csv'))) )
}


# MASS's results.
mass <- function (c) {
  pr  <- if (is.null(c$prior)) as.vector(table(c$groups)) / length(c$groups) else unname(c$prior[levels(c$groups)])
  f   <- MASS::lda(c$x, c$groups, prior = pr)
  loo <- vapply(seq_len(ncol(c$x)), function (m)
    as.integer(MASS::lda(c$x[, seq_len(m), drop = FALSE], c$groups, prior = pr, CV = TRUE)$class),
    integer(nrow(c$x)) )
  cv  <- MASS::lda(c$x, c$groups, prior = pr, CV = TRUE)
  list(
    scores    = unname(predict(f)$x),
    scaling   = unname(f$scaling),
    svd       = f$svd,
    means     = unname(f$means),
    prior     = unname(f$prior),
    loo_class = loo,
    loo_post  = unname(cv$posterior) )
}


for (i in seq_along(cases)) {
  c <- cases[[i]]
  c$expected <- list(sklearn = sk(c$name), mass = mass(c))
  cases[[i]] <- c
}

if (is.null(cases[[1]]$expected$sklearn))
  message('No scikit-learn results yet: run data-raw/expected/lda.py, then this script again.')

saveRDS(
  list(cases = cases),
  file.path('tests', 'tinytests', 'expected', 'lda.rds'),
  version = 2 )
