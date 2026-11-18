# Writes tests/tinytests/expected/ca.rds: inputs and reference results for
# test-ca.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/ca.r      # writes scikit-bio's inputs
#     python data-raw/expected/ca.py      # scikit-bio's outputs
#     Rscript data-raw/expected/ca.r      # folds them in
#
# vegan and ade4 are GPL: only their outputs are used, never their source.
# Generated with vegan 2.7.8 and ade4 1.7-24.
#
# Each case holds:
#   name, counts, k   A label, a samples x features matrix, and the axes
#                     requested.
#   warn              The warning ca() gives when it returns fewer than
#                     min(k, n - 1) axes, or NULL.
#   expected          vegan::ca()'s results (the same as vegan::cca(counts)),
#                     named after ecodive's attributes: sites
#                     (scores(display = "sites"), scaling 2), loadings
#                     (scores(display = "species"), scaling 2, with 0 for
#                     features with no observations, which vegan leaves
#                     out), eig, prop ("Proportion Explained") and trace
#                     (tot.chi).
#   ade4              ade4::dudi.coa()'s eig (every eigenvalue), l1 (row
#                     scores of unit weighted variance) and co (column
#                     principal coordinates), on the features with
#                     observations, with the same number of axes.
#
# And `skbio`: scikit-bio 0.7.4's ca() (BSD-3) in scaling 2, of the dune and
# varespec tables, from the CSV files that data-raw/expected/ca.py writes to
# data-raw/expected/ca/: every axis's sites and species scores, and its
# eigenvalues.

library(vegan)

data(dune, varespec)

dune <- as.matrix(dune)
vare <- as.matrix(varespec)

set.seed(20261121)
n   <- 110
big <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.4), n, 40)
big[, 1] <- big[, 1] + 1
dimnames(big) <- list(sprintf('S%03d', seq_len(n)), sprintf('F%02d', 1:40))

# A table of rank 3 (4 features), and one with an all-zero feature.
small <- dune[1:10, c('Achimill', 'Agrostol', 'Bromhord', 'Poatriv')]
small <- small + 1
dune_zero <- cbind(dune, Zero = 0)


cs <- function (name, counts, k = 3) list(name = name, counts = counts, k = k)

cases <- list(
  cs('dune',                 dune),
  cs('dune, every axis',     dune, k = 19),
  cs('dune, zero feature',   dune_zero, k = 4),
  cs('dune, k > n - 1',      dune[1:8, ], k = 10),
  cs('dune, n = 3',          dune[c(1, 7, 15), ], k = 3),
  cs('varespec',             vare, k = 6),
  cs('varespec, every axis', vare, k = 23),
  cs('n = 110',              big, k = 10),
  cs('rank 3',               small, k = 8) )


for (i in seq_along(cases)) {

  cs <- cases[[i]]
  Y  <- cs$counts

  fit <- vegan::ca(Y)
  ev  <- eigenvals(fit)
  ev  <- ev[ev > 0]
  k   <- min(cs$k, length(ev))
  sc  <- scores(fit, choices = seq_len(k), display = c('sites', 'species'),
                scaling = 2, droplist = FALSE)

  k1      <- min(cs$k, nrow(Y) - 1)
  cs$warn <- if (k < k1) sprintf('only %d of the first %d', k, k1)

  species <- matrix(0, ncol(Y), k, dimnames = list(colnames(Y), colnames(sc$species)))
  species[rownames(sc$species), ] <- sc$species

  cs$expected <- list(
    sites    = unclass(sc$sites),
    loadings = species,
    eig      = unname(ev[seq_len(k)]),
    prop     = unname(summary(eigenvals(fit))['Proportion Explained', names(ev)[seq_len(k)]]),
    trace    = fit$tot.chi )

  Yo <- Y[, colSums(Y) > 0, drop = FALSE]
  a  <- ade4::dudi.coa(as.data.frame(Yo), scannf = FALSE, nf = k)
  cs$ade4 <- list(
    eig = a$eig,
    l1  = unname(as.matrix(a$l1)),
    co  = unname(as.matrix(a$co)) )

  cases[[i]] <- cs
}


# scikit-bio: write its inputs, and read its outputs if ca.py has run.
skdir <- file.path('data-raw', 'expected', 'ca')
dir.create(skdir, showWarnings = FALSE, recursive = TRUE)
write.csv(dune, file.path(skdir, 'dune.csv'))
write.csv(vare, file.path(skdir, 'varespec.csv'))

rd    <- function (f) as.matrix(read.csv(file.path(skdir, f), row.names = 1, check.names = FALSE))
skbio <- NULL
if (file.exists(file.path(skdir, 'skbio_dune_samples.csv'))) {
  skbio <- list()
  for (nm in c('dune', 'varespec')) {
    Y <- if (nm == 'dune') dune else vare
    m <- nrow(Y) - 1 # axes with non-zero eigenvalues
    skbio[[nm]] <- list(
      counts  = Y,
      sites   = rd(sprintf('skbio_%s_samples.csv',  nm))[, seq_len(m)],
      species = rd(sprintf('skbio_%s_features.csv', nm))[, seq_len(m)],
      eigvals = rd(sprintf('skbio_%s_eigvals.csv',  nm))[seq_len(m), 1] )
  }
}


out <- 'tests/tinytests/expected/ca.rds'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  object = list(
    cases    = cases,
    skbio    = skbio,
    versions = list(
      vegan = format(packageVersion('vegan')),
      ade4  = format(packageVersion('ade4')),
      R     = R.version.string ) ),
  file     = out,
  version  = 2,
  compress = 'xz' )

cat(out, ':', file.size(out), 'bytes\n')
