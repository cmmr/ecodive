# Writes tests/tinytests/expected/varpart.rds: inputs and vegan::varpart()
# results for test-varpart.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/varpart.r
#
# vegan is GPL: only its outputs are used, never its source. Generated with
# vegan 2.7.8 and permute 0.10.0.
#
# `data` holds the tables: vegan's mite, mite.env and mite.pcnm, varespec
# and varechem, and `noise`, random variables unrelated to anything.
#
# Each case holds:
#   x          The community: list(table = name in `data`, kind), kind being
#              'hellinger' (RDA of the Hellinger-transformed table; vegan's
#              transfo = "hel"), 'bray' (db-RDA of vegdist()'s Bray-Curtis
#              distances, ecodive's bray()), 'euclidean' (db-RDA of the
#              Euclidean distances of the Hellinger-transformed table), or
#              'chisquare' (CCA; vegan's chisquare = TRUE).
#   sets       Two to four sets, each list(formula = "~ a + b") over
#              `env`, or list(table = name in `data`, cols = column
#              indices), given to vegan as a data frame.
#   env        The name in `data` of the metadata for the formulas.
#   permutations  For 'chisquare': a shuffleSet() matrix.
#   expected   vegan's results: indfract (the Venn fractions' adjusted
#              R-squared), fract (the sets and unions: Df, R.square,
#              Adj.R.square, named as vegan names them), contr (the
#              adjusted R-squared of contr1, and of contr2 for four sets,
#              named as vegan names them), SS.Y, and warning (vegan's
#              warning, if any).

library(vegan)

data(mite, mite.env, mite.pcnm, varespec, varechem)

set.seed(20261119)
noise <- data.frame(row.names = rownames(varechem), z1 = rnorm(24), z2 = rnorm(24))
noise_m <- data.frame(row.names = rownames(mite), z1 = rnorm(70), z2 = rnorm(70), z3 = rnorm(70))

data <- list(
  mite      = as.matrix(mite),
  mite.env  = mite.env,
  mite.pcnm = mite.pcnm,
  varespec  = as.matrix(varespec),
  varechem  = varechem,
  noise     = noise,
  noise_m   = noise_m )

f  <- function (s) list(formula = s)
tb <- function (t, cols = NULL) list(table = t, cols = cols)


cases <- list(

  # mite: two, three and four sets; Hellinger table and Bray-Curtis.
  list(name = 'mite2_hel',  x = list(table = 'mite', kind = 'hellinger'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'))),
  list(name = 'mite3_hel',  x = list(table = 'mite', kind = 'hellinger'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'), tb('mite.pcnm'))),
  list(name = 'mite4_hel',  x = list(table = 'mite', kind = 'hellinger'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'),
                   tb('mite.pcnm', 1:5), tb('mite.pcnm', 6:10))),
  list(name = 'mite2_bray', x = list(table = 'mite', kind = 'bray'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'))),
  list(name = 'mite3_bray', x = list(table = 'mite', kind = 'bray'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'), tb('mite.pcnm'))),
  list(name = 'mite4_bray', x = list(table = 'mite', kind = 'bray'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'),
                   tb('mite.pcnm', 1:5), tb('mite.pcnm', 6:10))),
  list(name = 'mite3_eucl', x = list(table = 'mite', kind = 'euclidean'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'), tb('mite.pcnm'))),

  # A data frame with factors (dummy-coded), and a set that explains nothing.
  list(name = 'mite3_df',   x = list(table = 'mite', kind = 'hellinger'), env = 'mite.env',
       sets = list(tb('mite.env', 1:3), tb('mite.pcnm', 1:8), tb('noise_m'))),
  list(name = 'mite3_df_bray', x = list(table = 'mite', kind = 'bray'), env = 'mite.env',
       sets = list(tb('mite.env', 1:3), tb('mite.pcnm', 1:8), tb('noise_m'))),

  # varespec: two, three and four sets.
  list(name = 'vare2_hel',  x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N + P + K'), f('~ Al + Fe + Mn'))),
  list(name = 'vare3_hel',  x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N + P + K'), f('~ Al + Fe + Mn'), f('~ pH + Humdepth'))),
  list(name = 'vare4_hel',  x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N + P'), f('~ K + Ca'), f('~ Al + Fe'), f('~ Humdepth + pH'))),
  list(name = 'vare2_bray', x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       sets = list(f('~ N + P + K'), f('~ Al + Fe + Mn'))),
  list(name = 'vare3_bray', x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       sets = list(f('~ N + P + K'), f('~ Al + Fe + Mn'), f('~ pH + Humdepth'))),
  list(name = 'vare4_bray', x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       sets = list(f('~ N + P'), f('~ K + Ca'), f('~ Al + Fe'), f('~ Humdepth + pH'))),

  # One variable per set.
  list(name = 'vare2_one',  x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N'), f('~ pH'))),
  list(name = 'vare2_one_bray', x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       sets = list(f('~ N'), f('~ pH'))),

  # Collinear sets: a variable in two sets, and a set within another.
  list(name = 'vare3_coll', x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N + P'), f('~ P + K'), f('~ Ca'))),
  list(name = 'vare2_nested', x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       sets = list(f('~ N + P + K'), f('~ N'))),

  # Sets that explain nothing: negative fractions.
  list(name = 'vare3_noise', x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       sets = list(f('~ N + P'), tb('noise', 1), tb('noise', 2))),

  # CCA (chi-square), with permutations.
  list(name = 'mite2_chi',  x = list(table = 'mite', kind = 'chisquare'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'))),
  list(name = 'mite3_chi',  x = list(table = 'mite', kind = 'chisquare'), env = 'mite.env',
       sets = list(f('~ SubsDens + WatrCont'), f('~ Substrate + Shrub + Topo'), tb('mite.pcnm', 1:6))),
  list(name = 'vare3_chi',  x = list(table = 'varespec', kind = 'chisquare'), env = 'varechem',
       sets = list(f('~ N + P'), f('~ K + Ca'), f('~ Al + Fe'))) )

# Not four sets with chisquare = TRUE: there vegan 2.7.8's results disagree
# with its own cca() (its SS.Y for mite is 13.8, where cca()'s total inertia,
# and its SS.Y for two or three sets, is 1.70), so test-varpart.r checks
# four sets against ecodive's cca() instead.


vegan_set <- function (s, env) {
  if (!is.null(s$formula)) return (stats::as.formula(s$formula))
  t <- data[[s$table]]
  if (!is.null(s$cols)) t <- t[, s$cols, drop = FALSE]
  as.data.frame(t)
}

for (i in seq_along(cases)) {

  cs  <- cases[[i]]
  Y   <- data[[cs$x$table]]
  env <- data[[cs$env]]
  X   <- lapply(cs$sets, vegan_set)
  n   <- nrow(Y)

  args <- switch(
    cs$x$kind,
    hellinger = list(Y = Y, transfo = 'hel'),
    bray      = list(Y = vegdist(Y)),
    euclidean = list(Y = dist(decostand(Y, 'hellinger'))),
    chisquare = {
      set.seed(i)
      cases[[i]]$permutations <- shuffleSet(n, 99)
      list(Y = Y, chisquare = TRUE, permutations = cases[[i]]$permutations)
    } )

  warn <- NULL
  v <- withCallingHandlers(
    do.call(varpart, c(args[1], X, list(data = env), args[-1])),
    warning = function (w) { warn <<- c(warn, conditionMessage(w)); invokeRestart('muffleWarning') } )

  p   <- v$part
  adj <- function (d) if (is.null(d$Adj.R.square)) d$Adj.R.squared else d$Adj.R.square
  r2  <- function (d) if (is.null(d$R.square)) d$R.squared else d$R.square

  contr <- c(
    stats::setNames(adj(p$contr1), rownames(p$contr1)),
    if (!is.null(p$contr2)) stats::setNames(adj(p$contr2), rownames(p$contr2)) )

  cases[[i]]$expected <- list(
    indfract = stats::setNames(adj(p$indfract), rownames(p$indfract)),
    fract    = data.frame(
      Df = p$fract$Df, R.square = r2(p$fract), Adj.R.square = adj(p$fract),
      row.names = rownames(p$fract) ),
    contr    = contr,
    SS.Y     = p$SS.Y,
    warning  = warn )
}

ref <- list(data = data, cases = cases)
saveRDS(ref, 'tests/tinytests/expected/varpart.rds', version = 2)
cat('Wrote', length(cases), 'cases;',
    file.size('tests/tinytests/expected/varpart.rds'), 'bytes\n')
