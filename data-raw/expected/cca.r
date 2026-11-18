# Writes tests/tinytests/expected/cca.rds: inputs and vegan::cca() results
# for test-cca.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/cca.r
#
# vegan is GPL: only its outputs are used, never its source. Generated with
# vegan 2.7.8 and permute 0.10.0.
#
# Each case holds:
#   counts, data      A samples x features matrix and its sample metadata.
#   constraints       Right-hand side of the constraints, e.g. "Management + A1",
#                     or NULL (conditions only, or no model).
#   conditions        Right-hand side of Condition(), or NULL.
#   k                 Axes requested.
#   warn              The warning cca() gives when it returns fewer than
#                     min(k, n - 1) axes, or NULL.
#   permutations      A shuffleSet() matrix, for the adjusted R-squared.
#   expected          vegan's results, named after ecodive's attributes:
#                     sites (scores(display = "sites"), scaling 2), eig, prop
#                     ("Proportion Explained"), trace, conditional,
#                     constrained, unconstrained, constrained_eig,
#                     r_squared and adj_r_squared (RsquareAdj() with
#                     `permutations`; numeric(0) without constraints),
#                     biplot, centroids, and loadings (scores(display =
#                     "species"), scaling 2, with 0 for the features with no
#                     observations, which vegan leaves out).
#
# And `skbio`: scikit-bio 0.7.4's cca() (BSD-3) of the dune table on
# `~ Management + A1`, in scaling 2, from the CSV files that
# data-raw/expected/cca.py writes to data-raw/expected/cca/ (this script
# writes its inputs there first): every axis's sites, species and biplot
# scores, and its eigenvalues.

library(vegan)

data(dune, dune.env, varespec, varechem)

set.seed(20261120)
n   <- 110
big <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.4), n, 40)
big[, 1] <- big[, 1] + 1
dimnames(big) <- list(sprintf('S%03d', seq_len(n)), sprintf('F%02d', 1:40))
big_env <- data.frame(
  row.names = rownames(big),
  group     = factor(sample(c('ctrl', 'trt'), n, TRUE)),
  site      = factor(sample(c('a', 'b', 'c'), n, TRUE)),
  ph        = round(rnorm(n, 7, 0.5), 2) )
big_env$ph <- big_env$ph + 0.3 * (big_env$group == 'trt')
big[big_env$group == 'trt', 2:6] <- big[big_env$group == 'trt', 2:6] + 2

dune <- as.matrix(dune)
vare <- as.matrix(varespec)

# An all-zero feature.
dune_zero <- cbind(dune, Zero = 0)


cs <- function (name, counts, data, constraints, conditions = NULL, k = 3)
  list(name = name, counts = counts, data = data, constraints = constraints,
       conditions = conditions, k = k)

cases <- list(
  cs('dune, factor + numeric',             dune, dune.env, 'Management + A1'),
  cs('dune, factor',                       dune, dune.env, 'Management'),
  cs('dune, numeric',                      dune, dune.env, 'A1', k = 4),
  cs('dune, partial',                      dune, dune.env, 'Management', 'A1'),
  cs('dune, partial, factor cond',         dune, dune.env, 'Moisture', 'Use', k = 5),
  cs('dune, ordered, partial',             dune, dune.env, 'Manure + A1', 'Moisture', k = 4),
  cs('dune, many terms',                   dune, dune.env, 'Moisture + Management + Use + A1', k = 6),
  cs('dune, interaction',                  dune, dune.env, 'Management * A1', k = 4),
  # k beyond the constrained axes, and beyond n - 1.
  cs('dune, k > constrained axes',         dune, dune.env, 'Use', k = 6),
  cs('dune, k > n - 1',                    dune[1:8, ], dune.env[1:8, ], 'Management', k = 10),
  # Conditions only, and no model at all (CA).
  cs('dune, conditions only',              dune, dune.env, NULL, 'A1'),
  cs('dune, conditions only, factor',      dune, dune.env, NULL, 'Management + A1', k = 4),
  cs('dune, unconstrained',                dune, dune.env, NULL, k = 5),
  # Aliased columns, and no residual degrees of freedom.
  cs('dune, aliased',                      dune, cbind(dune.env, A1x2 = 2 * dune.env$A1), 'Management + A1 + A1x2'),
  cs('dune, aliased, partial',             dune, dune.env, 'Moisture + Management', 'Use + Manure'),
  cs('dune, saturated',                    dune, cbind(dune.env, id = factor(1:20)), 'id'),
  # An all-zero feature.
  cs('dune, zero feature',                 dune_zero, dune.env, 'Management + A1'),
  # More features than samples; numeric constraints.
  cs('varespec, 5 numeric',                vare, varechem, 'N + P + K + Al + pH', k = 8),
  cs('varespec, 12 numeric',               vare, varechem, 'N + P + K + Ca + Mg + S + Al + Fe + Mn + Zn + Mo + Baresoil', k = 14),
  cs('varespec, partial',                  vare, varechem, 'N + P + K', 'Baresoil + pH', k = 5),
  cs('varespec, unconstrained',            vare, varechem, NULL, k = 6),
  # n = 110, fewer features than samples.
  cs('n = 110, 2 groups',                  big, big_env, 'group', k = 4),
  cs('n = 110, partial',                   big, big_env, 'group + ph', 'site', k = 5),
  cs('n = 110, k > rank',                  big, big_env, 'group', k = 45) )


# The right-hand side of a case's formula.
case_rhs <- function (cs) {
  rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
  if (length(rhs) == 0) '1' else paste(rhs, collapse = ' + ')
}


for (i in seq_along(cases)) {

  cs <- cases[[i]]
  Y  <- cs$counts

  fit <- if (is.null(cs$constraints) && is.null(cs$conditions)) {
    vegan::cca(Y)
  } else {
    vegan::cca(as.formula(paste('Y ~', case_rhs(cs))), data = cs$data)
  }

  # vegan 2.7.8 needs shuffleSet()'s object, not a plain matrix.
  set.seed(i)
  ss <- shuffleSet(nrow(Y), 99)
  cs$permutations <- matrix(as.integer(ss), nrow(ss))

  r2 <- if (is.null(cs$constraints)) list(r.squared = numeric(0), adj.r.squared = numeric(0))
        else RsquareAdj(fit, permutations = ss)

  ev  <- eigenvals(fit)
  ev  <- ev[ev > 0]
  k   <- min(cs$k, length(ev))
  sc  <- scores(fit, choices = seq_len(k), display = c('sites', 'species', 'bp', 'cn'),
                scaling = 2, droplist = FALSE)

  k1      <- min(cs$k, nrow(Y) - 1)
  cs$warn <- if (k < k1) sprintf('only %d of the first %d', k, k1)

  # Features with no observations: vegan leaves them out; cca() gives 0.
  species <- matrix(0, ncol(Y), k, dimnames = list(colnames(Y), colnames(sc$species)))
  species[rownames(sc$species), ] <- sc$species

  cs$expected <- list(
    sites           = unclass(sc$sites),
    loadings        = species,
    eig             = unname(ev[seq_len(k)]),
    prop            = unname(summary(eigenvals(fit))['Proportion Explained', names(ev)[seq_len(k)]]),
    trace           = fit$tot.chi,
    conditional     = if (is.null(fit$pCCA)) 0 else fit$pCCA$tot.chi,
    constrained     = if (is.null(fit$CCA)) 0 else fit$CCA$tot.chi,
    unconstrained   = if (is.null(fit$CA$tot.chi)) 0 else fit$CA$tot.chi,
    constrained_eig = if (is.null(fit$CCA)) numeric(0) else unname(fit$CCA$eig),
    r_squared       = r2$r.squared,
    adj_r_squared   = r2$adj.r.squared,
    biplot          = sc$biplot,
    centroids       = sc$centroids )

  cases[[i]] <- cs
}


# scikit-bio: write its inputs, and read its outputs if cca.py has run.
skdir <- file.path('data-raw', 'expected', 'cca')
dir.create(skdir, showWarnings = FALSE, recursive = TRUE)
write.csv(dune, file.path(skdir, 'dune.csv'))
skx <- model.matrix(~ Management + A1, dune.env)[, -1]
write.csv(skx, file.path(skdir, 'dune_env.csv'))

rd    <- function (f) as.matrix(read.csv(file.path(skdir, f), row.names = 1, check.names = FALSE))
skbio <- NULL
if (file.exists(file.path(skdir, 'skbio_samples.csv'))) {
  m <- min(19, nrow(dune) - 1) # axes with non-zero eigenvalues
  skbio <- list(
    counts   = dune,
    x        = as.data.frame(skx),
    sites    = rd('skbio_samples.csv')[, seq_len(m)],
    species  = rd('skbio_features.csv')[, seq_len(m)],
    biplot   = rd('skbio_biplot.csv'),
    eigvals  = rd('skbio_eigvals.csv')[seq_len(m), 1] )
}


out <- 'tests/tinytests/expected/cca.rds'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  object = list(
    cases    = cases,
    skbio    = skbio,
    versions = list(
      vegan   = format(packageVersion('vegan')),
      permute = format(packageVersion('permute')),
      R       = R.version.string ) ),
  file     = out,
  version  = 2,
  compress = 'xz' )

cat(out, ':', file.size(out), 'bytes\n')
