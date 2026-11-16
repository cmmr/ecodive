# Writes tests/tinytests/expected/dbrda.rds: inputs and vegan::dbrda() results
# for test-dbrda.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/dbrda.r
#
# Each case holds:
#   dm, data          A dist object and a data.frame of sample metadata.
#   constraints       Right-hand side of the constraints, e.g. "Management + A1".
#   conditions        Right-hand side of Condition(), or NULL.
#   k                 Axes requested.
#   x, z              The numeric design (model.matrix() without intercept)
#                     for the constraints and conditions, as vegan builds it.
#   expected          vegan's results, named after ecodive's attributes:
#                     sites (scores(display = "sites"), scaling 2), eig, prop
#                     ("Proportion Explained"), trace, conditional,
#                     constrained, unconstrained, constrained_eig, biplot,
#                     centroids.

library(vegan)

data(dune, dune.env)

set.seed(20261116)
n   <- 110
big <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.4), n, 40)
big[, 1] <- big[, 1] + 1
rownames(big) <- sprintf('S%03d', seq_len(n))
big_env <- data.frame(
  row.names = rownames(big),
  group     = factor(sample(c('ctrl', 'trt'), n, TRUE)),
  site      = factor(sample(c('a', 'b', 'c'), n, TRUE)),
  ph        = round(rnorm(n, 7, 0.5), 2) )
big_env$ph <- big_env$ph + 0.3 * (big_env$group == 'trt')


cases <- list(
  list(name = 'bray, factor + continuous',  dm = vegdist(dune), data = dune.env, constraints = 'Management + A1', conditions = NULL,         k = 3),
  list(name = 'bray, factor only',          dm = vegdist(dune), data = dune.env, constraints = 'Management',      conditions = NULL,         k = 3),
  list(name = 'bray, 2 constrained axes',   dm = vegdist(dune), data = dune.env, constraints = 'Use',             conditions = NULL,         k = 3),
  list(name = 'bray, continuous only',      dm = vegdist(dune), data = dune.env, constraints = 'A1',              conditions = NULL,         k = 4),
  list(name = 'bray, partial',              dm = vegdist(dune), data = dune.env, constraints = 'Management',      conditions = 'A1',         k = 3),
  list(name = 'euclidean (rescaled), partial', dm = dist(dune), data = dune.env, constraints = 'Moisture',        conditions = 'Use',        k = 5),
  list(name = 'jaccard, many terms',        dm = vegdist(dune, 'jaccard'), data = dune.env, constraints = 'Moisture + Management + Use + A1', conditions = NULL, k = 6),
  list(name = 'bray, n = 110, 2 groups',    dm = vegdist(big),  data = big_env,  constraints = 'group',           conditions = NULL,         k = 4),
  list(name = 'bray, n = 110, partial',     dm = vegdist(big),  data = big_env,  constraints = 'group + ph',      conditions = 'site',       k = 3) )


for (i in seq_along(cases)) {

  cs  <- cases[[i]]
  dm  <- cs$dm
  rhs <- cs$constraints
  if (!is.null(cs$conditions)) rhs <- paste0(rhs, ' + Condition(', cs$conditions, ')')
  fit <- dbrda(as.formula(paste('dm ~', rhs)), data = cs$data)

  ev  <- eigenvals(fit)
  ev  <- ev[ev > 0]
  k   <- min(cs$k, length(ev))
  sc  <- scores(fit, choices = seq_len(k), display = c('sites', 'bp', 'cn'), droplist = FALSE)

  cs$x <- model.matrix(as.formula(paste('~', cs$constraints)), cs$data)[, -1, drop = FALSE]
  cs$z <- if (!is.null(cs$conditions))
    model.matrix(as.formula(paste('~', cs$conditions)), cs$data)[, -1, drop = FALSE]
  attr(cs$x, 'assign') <- attr(cs$x, 'contrasts') <- NULL
  if (!is.null(cs$z)) attr(cs$z, 'assign') <- attr(cs$z, 'contrasts') <- NULL

  cs$expected <- list(
    sites           = unclass(sc$sites),
    eig             = unname(ev[seq_len(k)]),
    prop            = unname(summary(eigenvals(fit))['Proportion Explained', names(ev)[seq_len(k)]]),
    trace           = fit$tot.chi,
    conditional     = if (is.null(fit$pCCA)) 0 else fit$pCCA$tot.chi,
    constrained     = fit$CCA$tot.chi,
    unconstrained   = fit$CA$tot.chi,
    constrained_eig = unname(fit$CCA$eig),
    biplot          = sc$biplot,
    centroids       = sc$centroids )

  cases[[i]] <- cs
}


out <- 'tests/tinytests/expected/dbrda.rds'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  object = list(
    cases    = cases,
    versions = list(
      vegan = format(packageVersion('vegan')),
      R     = R.version.string ) ),
  file     = out,
  version  = 2,
  compress = 'xz' )

cat(out, ':', file.size(out), 'bytes\n')
