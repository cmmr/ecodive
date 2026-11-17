# Writes tests/tinytests/expected/capscale.rds: inputs and vegan::capscale()
# results for test-capscale.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/capscale.r
#
# Each case holds:
#   dm, data          A dist object and a data.frame of sample metadata.
#   constraints       Right-hand side of the constraints, e.g. "Management + A1",
#                     or NULL (conditions only).
#   conditions        Right-hand side of Condition(), or NULL.
#   k                 Axes requested.
#   warn              The warning capscale() gives when it returns fewer than
#                     min(k, n - 1) axes, or NULL.
#   euclidean         Whether G has no eigenvalue below vegan's tolerance (no
#                     imaginary inertia).
#   expected          vegan's results, named after ecodive's attributes:
#                     sites (scores(display = "sites"), scaling 2), eig, prop
#                     ("Proportion Explained"), trace (tot.chi), conditional,
#                     constrained, unconstrained, imaginary,
#                     constrained_eig, r_squared and adj_r_squared
#                     (RsquareAdj(); numeric(0) without constraints), biplot,
#                     centroids.

library(vegan)

data(dune, dune.env, varespec, varechem)

set.seed(20261117)
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

# Hellinger distances: Euclidean, all below 4.
hell <- dist(sqrt(dune / rowSums(dune)))

# A few negative eigenvalues among many positive ones: Gower matrix
# X X' - Y Y', with X of full rank. capscale() then finds the eigenvectors
# of the few, not of the many.
set.seed(7)
nf  <- 60
X   <- scale(matrix(rnorm(nf * 150), nf), scale = FALSE)
Y   <- scale(matrix(rnorm(nf * 3), nf), scale = FALSE)
G   <- tcrossprod(X) - 3 * tcrossprod(Y)
few <- as.dist(sqrt(outer(diag(G), diag(G), '+') - 2 * G))
attr(few, 'Labels') <- sprintf('F%02d', seq_len(nf))
few_env <- data.frame(
  row.names = attr(few, 'Labels'),
  g = factor(rep(c('a', 'b', 'c'), length.out = nf)),
  x = round(rnorm(nf), 3),
  z = round(rnorm(nf), 3) )

dune_alias <- cbind(dune.env, A1x2 = 2 * dune.env$A1)
dune_sat   <- cbind(dune.env, id = factor(1:20))


cases <- list(
  # Bray-Curtis: negative eigenvalues.
  list(name = 'bray, factor + continuous',  dm = vegdist(dune), data = dune.env, constraints = 'Management + A1', conditions = NULL,         k = 3),
  list(name = 'bray, factor only',          dm = vegdist(dune), data = dune.env, constraints = 'Management',      conditions = NULL,         k = 3),
  list(name = 'bray, continuous only',      dm = vegdist(dune), data = dune.env, constraints = 'A1',              conditions = NULL,         k = 4),
  list(name = 'bray, partial',              dm = vegdist(dune), data = dune.env, constraints = 'Management',      conditions = 'A1',         k = 3),
  list(name = 'bray, ordered, partial',     dm = vegdist(dune), data = dune.env, constraints = 'Manure + A1',     conditions = 'Moisture',   k = 4),
  list(name = 'bray, interaction',          dm = vegdist(dune), data = dune.env, constraints = 'Use * A1',        conditions = NULL,         k = 3),
  list(name = 'bray, n = 110, partial',     dm = vegdist(big),  data = big_env,  constraints = 'group + ph',      conditions = 'site',       k = 3),
  list(name = 'bray, 12 numeric constraints', dm = vegdist(varespec), data = varechem, constraints = 'N + P + K + Ca + Mg + S + Al + Fe + Mn + Zn + Mo + Baresoil', conditions = NULL, k = 14),
  list(name = 'jaccard, many terms',        dm = vegdist(dune, 'jaccard'), data = dune.env, constraints = 'Moisture + Management + Use + A1', conditions = NULL, k = 6),
  # Manhattan distances of 4 or more (rescaled), with negative eigenvalues.
  list(name = 'manhattan (rescaled), partial', dm = vegdist(varespec, 'manhattan'), data = varechem, constraints = 'N + P', conditions = 'K', k = 5),
  # Euclidean: no imaginary inertia.
  list(name = 'euclidean (rescaled), partial', dm = dist(dune), data = dune.env, constraints = 'Moisture', conditions = 'Use',  k = 5),
  list(name = 'hellinger, factor',          dm = hell,          data = dune.env, constraints = 'Management',      conditions = NULL,         k = 4),
  list(name = 'euclidean (rescaled), n = 110', dm = dist(big),  data = big_env,  constraints = 'group + ph',      conditions = 'site',       k = 4),
  # A few negative eigenvalues.
  list(name = 'few negative, partial',      dm = few,           data = few_env,  constraints = 'g + x',           conditions = 'z',          k = 6),
  # Aliased (redundant) columns.
  list(name = 'bray, aliased',              dm = vegdist(dune), data = dune_alias, constraints = 'Management + A1 + A1x2', conditions = NULL, k = 3),
  list(name = 'bray, aliased, partial',     dm = vegdist(dune), data = dune.env, constraints = 'Moisture + Management', conditions = 'Use + Manure', k = 3),
  # More axes than positive eigenvalues; and a saturated model.
  list(name = 'bray, k > positive eigenvalues', dm = vegdist(dune), data = dune.env, constraints = 'Management', conditions = NULL, k = 19),
  list(name = 'bray, saturated',            dm = vegdist(dune), data = dune_sat, constraints = 'id',              conditions = NULL,         k = 3),
  # Conditions only.
  list(name = 'bray, conditions only',      dm = vegdist(dune), data = dune.env, constraints = NULL, conditions = 'Management + A1', k = 4) )

# The right-hand side of a case's formula.
case_rhs <- function (cs) {
  rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
  paste(rhs, collapse = ' + ')
}


for (i in seq_along(cases)) {

  cs  <- cases[[i]]
  dm  <- cs$dm
  fit <- capscale(as.formula(paste('dm ~', case_rhs(cs))), data = cs$data)

  # numeric(0) for both without constraints (conditions only).
  r2  <- RsquareAdj(fit)

  ev  <- eigenvals(fit)
  ev  <- ev[ev > 0]
  k   <- min(cs$k, length(ev))
  sc  <- scores(fit, choices = seq_len(k), display = c('sites', 'bp', 'cn'), droplist = FALSE)

  # capscale() warns when fewer than min(k, n - 1) axes have positive
  # eigenvalues.
  k1      <- min(cs$k, attr(dm, 'Size') - 1)
  cs$warn <- if (k < k1) sprintf('only %d of the first %d', k, k1)

  imaginary    <- if (is.null(fit$CA$imaginary.chi)) 0 else fit$CA$imaginary.chi
  cs$euclidean <- imaginary == 0

  cs$expected <- list(
    sites           = unclass(sc$sites),
    eig             = unname(ev[seq_len(k)]),
    prop            = unname(summary(eigenvals(fit))['Proportion Explained', names(ev)[seq_len(k)]]),
    trace           = fit$tot.chi,
    conditional     = if (is.null(fit$pCCA)) 0 else fit$pCCA$tot.chi,
    constrained     = if (is.null(fit$CCA)) 0 else fit$CCA$tot.chi,
    unconstrained   = if (is.null(fit$CA$tot.chi)) 0 else fit$CA$tot.chi,
    imaginary       = imaginary,
    constrained_eig = if (is.null(fit$CCA)) numeric(0) else unname(fit$CCA$eig),
    r_squared       = r2$r.squared,
    adj_r_squared   = r2$adj.r.squared,
    biplot          = sc$biplot,
    centroids       = sc$centroids )

  cases[[i]] <- cs
}


out <- 'tests/tinytests/expected/capscale.rds'
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
