# Writes tests/tinytests/expected/forward_select.rds: inputs and
# vegan::ordiR2step() and adespatial::forward.sel() results for
# test-forward_select.r, so that the tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/forward_select.r
#
# vegan and adespatial are GPL: only their outputs are used, never their
# source. Generated with vegan 2.7.8, permute 0.10.0 and adespatial 0.3-30.
#
# `data` holds the tables: vegan's mite, mite.env and mite.pcnm (its first
# six columns, joined to mite.env as `mite.vars`), and varespec and
# varechem; and the permutation matrices, perm_mite and perm_vare.
#
# Each case holds:
#   x             list(table = name in `data`, kind): 'hellinger' (RDA of
#                 the Hellinger-transformed table; decostand("hel")) or
#                 'bray' (db-RDA of vegdist()'s Bray-Curtis distances,
#                 ecodive's bray()).
#   formula       The scope, as forward_select()'s `formula` over `env`.
#   env           The name in `data` of the metadata.
#   permutations  The name in `data` of a shuffleSet() matrix, given to
#                 ordiR2step().
#   expected      ordiR2step()'s anova table: `variable`, `adj` (R2.adj),
#                 `df`, `F`, and `p` for each step, and `global` (the
#                 "<All variables>" R2.adj).
#
# `forward_sel`: adespatial::forward.sel() of the Hellinger-transformed
# mite table on the numeric columns of mite.vars, with adjR2thresh set to
# the global model's adjusted R-squared: its table (variables, R2, R2Cum,
# AdjR2Cum, F), whose p-values come from its own random permutations.

library(vegan)
library(adespatial)

data(mite, mite.env, mite.pcnm, varespec, varechem)

mite.vars <- cbind(mite.env, mite.pcnm[, 1:6])

data <- list(
  mite      = as.matrix(mite),
  mite.vars = mite.vars,
  varespec  = as.matrix(varespec),
  varechem  = varechem )

set.seed(20261119)
data$perm_mite <- shuffleSet(70, 199)
data$perm_vare <- shuffleSet(24, 499)

cases <- list(
  list(name = 'mite_hel',   x = list(table = 'mite', kind = 'hellinger'), env = 'mite.vars',
       formula = '~ .', permutations = 'perm_mite'),
  list(name = 'mite_bray',  x = list(table = 'mite', kind = 'bray'), env = 'mite.vars',
       formula = '~ .', permutations = 'perm_mite'),
  list(name = 'mite_cond',  x = list(table = 'mite', kind = 'hellinger'), env = 'mite.vars',
       formula = '~ SubsDens + WatrCont + Substrate + Topo + V2 + V3 + V4 + V5 + V6 + Condition(V1 + Shrub)',
       permutations = 'perm_mite'),
  list(name = 'mite_cond_bray', x = list(table = 'mite', kind = 'bray'), env = 'mite.vars',
       formula = '~ SubsDens + WatrCont + Substrate + V2 + V3 + V4 + Condition(Topo)',
       permutations = 'perm_mite'),
  list(name = 'mite_inter', x = list(table = 'mite', kind = 'hellinger'), env = 'mite.vars',
       formula = '~ SubsDens * Topo + WatrCont + Shrub', permutations = 'perm_mite'),
  list(name = 'vare_hel',   x = list(table = 'varespec', kind = 'hellinger'), env = 'varechem',
       formula = '~ .', permutations = 'perm_vare'),
  list(name = 'vare_bray',  x = list(table = 'varespec', kind = 'bray'), env = 'varechem',
       formula = '~ .', permutations = 'perm_vare') )


# vegan's model of the case, with the scope's terms (or none) as constraints.
fit <- function (cs, rhs) {
  y   <- data[[cs$x$table]]
  env <- data[[cs$env]]
  f   <- switch(
    cs$x$kind,
    hellinger = { Y <- decostand(y, 'hel'); stats::as.formula(paste('Y', rhs)) },
    bray      = { Y <- vegdist(y);          stats::as.formula(paste('Y', rhs)) } )
  environment(f) <- environment()
  if (cs$x$kind == 'bray') dbrda(f, data = env) else rda(f, data = env)
}

for (i in seq_along(cases)) {
  cs   <- cases[[i]]
  full <- cs$formula
  if (full == '~ .') full <- paste('~', paste(names(data[[cs$env]]), collapse = ' + '))
  cond <- regmatches(full, regexpr('Condition\\(.*\\)', full))
  base <- if (length(cond)) paste('~', cond) else '~ 1'

  r <- ordiR2step(fit(cs, base), fit(cs, full), permutations = data[[cs$permutations]], trace = FALSE)
  a <- r$anova
  s <- grepl('^\\+ ', rownames(a))
  cases[[i]]$expected <- list(
    variable = sub('^\\+ ', '', rownames(a)[s]),
    adj      = a$R2.adj[s],
    df       = a$Df[s],
    F        = a$F[s],
    p        = a$`Pr(>F)`[s],
    global   = a$R2.adj[rownames(a) == '<All variables>'] )
}


# adespatial::forward.sel(), as a second check of the order.
num <- mite.vars[, vapply(mite.vars, is.numeric, TRUE)]
hel <- decostand(mite, 'hel')
set.seed(1)
fw  <- forward.sel(
  hel, num, nperm = 999, verbose = FALSE,
  adjR2thresh = RsquareAdj(rda(hel ~ ., num))$adj.r.squared )
forward_sel <- data.frame(
  variable  = fw$variables,
  r2        = fw$R2,
  r2_cum    = fw$R2Cum,
  adj_cum   = fw$AdjR2Cum,
  F         = fw$F,
  stringsAsFactors = FALSE )


saveRDS(
  list(data = data, cases = cases, forward_sel = forward_sel,
       vegan = format(packageVersion('vegan')), adespatial = format(packageVersion('adespatial'))),
  file.path('tests', 'tinytests', 'expected', 'forward_select.rds'),
  version = 2 )
