# Writes tests/tinytests/expected/rda.rds: inputs and vegan::rda() results
# for test-rda.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/rda.r
#
# vegan is GPL: only its outputs are used, never its source. Generated with
# vegan 2.7.8 and permute 0.10.0. (vegan 2.6-x permutes partial models
# differently, so its anova.cca() p-values for `perm$partial` differ from
# permanova()'s; 2.7-x matches them.)
#
# Each case holds:
#   counts, data      A samples x features matrix and its sample metadata.
#   norm, pseudocount, scale
#                     rda()'s arguments. vegan::rda() is given the table
#                     after the same normalization, done here in base R.
#   constraints       Right-hand side of the constraints, e.g. "Management + A1",
#                     or NULL (conditions only, or no model).
#   conditions        Right-hand side of Condition(), or NULL.
#   k                 Axes requested.
#   warn              The warning rda() gives when it returns fewer than
#                     min(k, n - 1) axes, or NULL.
#   expected          vegan's results, named after ecodive's attributes:
#                     sites (scores(display = "sites"), scaling 2), eig, prop
#                     ("Proportion Explained"), trace, conditional,
#                     constrained, unconstrained, constrained_eig,
#                     r_squared and adj_r_squared (RsquareAdj(); numeric(0)
#                     without constraints), biplot, centroids, and
#                     loadings (scores(display = "species"), scaling 2).
#
# And `perm`: one model with anova.cca()'s pseudo-F and p-value for a fixed
# permutation matrix, which permanova() of Euclidean distances must match.
#
# And `skbio`: scikit-bio 0.7.4's rda() (BSD-3) of the Hellinger-transformed
# dune table on `~ Management + A1`, in scaling 2, from the CSV files that
# data-raw/expected/rda.py writes to data-raw/expected/rda/ (this script
# writes its inputs there first): every axis's sites, species and biplot
# scores, and its "eigvals", which are singular values.

library(vegan)

data(dune, dune.env, varespec, varechem)

set.seed(20261118)
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

# A constant feature and an all-zero feature.
dune_const <- cbind(dune, Const = 3, Zero = 0)


# The normalizations rda() offers, in base R.
transform <- function (m, norm, pseudocount = NULL) {
  m <- m + 0
  switch(
    EXPR      = norm,
    none      = m,
    percent   = m / rowSums(m),
    chord     = m / sqrt(rowSums(m^2)),
    hellinger = sqrt(m / rowSums(m)),
    binary    = (m > 0) + 0,
    clr       = { l <- log(m + pseudocount); l - rowMeans(l) } )
}


cs <- function (name, counts, data, constraints, conditions = NULL, k = 3,
                norm = 'none', scale = FALSE, pseudocount = NULL)
  list(name = name, counts = counts, data = data, norm = norm,
       pseudocount = pseudocount, scale = scale, constraints = constraints,
       conditions = conditions, k = k)

cases <- list(
  cs('dune, hellinger, factor + numeric',  dune, dune.env, 'Management + A1', norm = 'hellinger'),
  cs('dune, chord, factor',                dune, dune.env, 'Management', norm = 'chord'),
  cs('dune, none, factor + numeric',       dune, dune.env, 'Management + A1'),
  cs('dune, scale',                        dune, dune.env, 'Management + A1', scale = TRUE),
  cs('dune, hellinger, scale',             dune, dune.env, 'Moisture', norm = 'hellinger', scale = TRUE, k = 4),
  cs('dune, hellinger, partial',           dune, dune.env, 'Management', 'A1', norm = 'hellinger'),
  cs('dune, none, partial, factor cond',   dune, dune.env, 'Moisture', 'Use', k = 5),
  cs('dune, hellinger, ordered, partial',  dune, dune.env, 'Manure + A1', 'Moisture', norm = 'hellinger', k = 4),
  cs('dune, hellinger, many terms',        dune, dune.env, 'Moisture + Management + Use + A1', norm = 'hellinger', k = 6),
  cs('dune, percent',                      dune, dune.env, 'Management', norm = 'percent'),
  cs('dune, binary',                       dune, dune.env, 'Management + A1', norm = 'binary'),
  cs('dune, clr',                          dune, dune.env, 'Management + A1', norm = 'clr', pseudocount = 1),
  cs('dune, interaction',                  dune, dune.env, 'Management * A1', norm = 'hellinger', k = 4),
  # k beyond the constrained axes, and beyond n - 1.
  cs('dune, k > constrained axes',         dune, dune.env, 'Use', norm = 'hellinger', k = 6),
  cs('dune, k > n - 1',                    dune[1:8, ], dune.env[1:8, ], 'Management', norm = 'hellinger', k = 10),
  # Conditions only, and no model at all (PCA).
  cs('dune, conditions only',              dune, dune.env, NULL, 'A1', norm = 'hellinger'),
  cs('dune, conditions only, factor',      dune, dune.env, NULL, 'Management + A1', k = 4),
  cs('dune, unconstrained',                dune, dune.env, NULL, norm = 'hellinger', k = 4),
  cs('dune, unconstrained, none',          dune, dune.env, NULL, k = 3),
  # Aliased columns, and no residual degrees of freedom.
  cs('dune, aliased',                      dune, cbind(dune.env, A1x2 = 2 * dune.env$A1), 'Management + A1 + A1x2', norm = 'hellinger'),
  cs('dune, aliased, partial',             dune, dune.env, 'Moisture + Management', 'Use + Manure', norm = 'hellinger'),
  cs('dune, saturated',                    dune, cbind(dune.env, id = factor(1:20)), 'id', norm = 'hellinger'),
  # A constant and an all-zero feature.
  cs('dune, constant features',            dune_const, dune.env, 'Management + A1', norm = 'none'),
  # More features than samples, 12 numeric constraints.
  cs('varespec, hellinger, 12 numeric',    vare, varechem, 'N + P + K + Ca + Mg + S + Al + Fe + Mn + Zn + Mo + Baresoil', norm = 'hellinger', k = 14),
  cs('varespec, scale, partial',           vare, varechem, 'N + P + K', 'Baresoil + pH', scale = TRUE, k = 5),
  # n = 110, fewer features than samples.
  cs('n = 110, hellinger, 2 groups',       big, big_env, 'group', norm = 'hellinger', k = 4),
  cs('n = 110, hellinger, partial',        big, big_env, 'group + ph', 'site', norm = 'hellinger', k = 3),
  cs('n = 110, none, partial',             big, big_env, 'group + ph', 'site', k = 5),
  cs('n = 110, chord, k > rank',           big, big_env, 'group', norm = 'chord', k = 45) )


# The right-hand side of a case's formula.
case_rhs <- function (cs) {
  rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
  if (length(rhs) == 0) '1' else paste(rhs, collapse = ' + ')
}


for (i in seq_along(cases)) {

  cs <- cases[[i]]
  Y  <- transform(cs$counts, cs$norm, cs$pseudocount)

  # vegan stops on constant features with scale = TRUE; rda() leaves them
  # out, so give vegan the table without them.
  if (cs$scale) Y <- Y[, apply(Y, 2, function (x) diff(range(x)) > 0), drop = FALSE]

  fit <- if (is.null(cs$constraints) && is.null(cs$conditions)) {
    vegan::rda(Y, scale = cs$scale)
  } else {
    vegan::rda(as.formula(paste('Y ~', case_rhs(cs))), data = cs$data, scale = cs$scale)
  }

  r2  <- RsquareAdj(fit)
  if (is.null(cs$constraints)) r2 <- list(r.squared = numeric(0), adj.r.squared = numeric(0))

  ev  <- eigenvals(fit)
  ev  <- ev[ev > 0]
  k   <- min(cs$k, length(ev))
  sc  <- scores(fit, choices = seq_len(k), display = c('sites', 'species', 'bp', 'cn'),
                scaling = 2, droplist = FALSE)

  k1      <- min(cs$k, nrow(Y) - 1)
  cs$warn <- if (k < k1) sprintf('only %d of the first %d', k, k1)

  species <- sc$species
  if (cs$scale) { # Features left out have loadings of 0.
    full <- matrix(0, ncol(cs$counts), k, dimnames = list(colnames(cs$counts), colnames(species)))
    full[rownames(species), ] <- species
    species <- full
  }

  cs$expected <- list(
    sites           = unclass(sc$sites),
    loadings        = unclass(species),
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


# The permutation test: anova.cca() on a tb-RDA, with a fixed permutation
# matrix, and the same for a partial model (permutations of the residuals
# under the reduced model, which permanova() also does).
set.seed(20261119)
perm_mtx <- shuffleSet(nrow(dune), 199)
h        <- transform(dune, 'hellinger')
fit1     <- vegan::rda(h ~ Management + A1, data = dune.env)
fit2     <- vegan::rda(h ~ Management + Condition(A1), data = dune.env)
a1       <- anova(fit1, permutations = perm_mtx)
a2       <- anova(fit2, permutations = perm_mtx)

perm <- list(
  counts       = dune,
  data         = dune.env,
  permutations = unclass(perm_mtx)[, , drop = FALSE],
  full         = list(formula = '~ Management + A1', f = a1$F[1], p = a1$`Pr(>F)`[1]),
  partial      = list(formula = '~ Management + Condition(A1)', f = a2$F[1], p = a2$`Pr(>F)`[1]) )


# scikit-bio: write its inputs, and read its outputs if rda.py has run.
skdir <- file.path('data-raw', 'expected', 'rda')
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
    biplot   = rd('skbio_biplot.csv')[, seq_len(m)],
    eigvals  = rd('skbio_eigvals.csv')[seq_len(m), 1] )
}


out <- 'tests/tinytests/expected/rda.rds'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  object = list(
    cases    = cases,
    perm     = perm,
    skbio    = skbio,
    versions = list(
      vegan   = format(packageVersion('vegan')),
      permute = format(packageVersion('permute')),
      R       = R.version.string ) ),
  file     = out,
  version  = 2,
  compress = 'xz' )

cat(out, ':', file.size(out), 'bytes\n')
