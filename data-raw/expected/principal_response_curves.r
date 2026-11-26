# Writes tests/tinytests/expected/principal_response_curves.rds: inputs and
# vegan::prc() results for test-principal_response_curves.r, so that the
# tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/principal_response_curves.r
#
# vegan is GPL: only its outputs are used, never its source. Generated with
# vegan 2.7.8 and permute 0.10-1.
#
# `data` holds vegan's pyrifos table (log(10 x + 1) abundances of 178 taxa in
# 132 samples: 12 ditches, 11 weeks, 5 doses; as a matrix), its `week`,
# `dose` and `ditch` factors (built as in ?pyrifos), and three shuffleSet()
# matrices of 99 permutations each: `perm_plots` (whole ditches, Plots(type =
# "free"), Within(type = "none")), `perm_free` (free) and `perm_week`
# (within weeks, how(blocks = week)).
#
# Each case holds:
#   rows          The samples used (indices into pyrifos), or NULL for all.
#   transform     'none', 'hellinger' (decostand(x, "hellinger") given to
#                 prc(); ecodive's norm = 'hellinger') or 'counts' (the
#                 table given to ecodive is round((exp(x) - 1) / 10), and
#                 prc() is given its clr with pseudocount 1:
#                 log(x + 1) - rowMeans(log(x + 1)); ecodive's norm = 'clr',
#                 pseudocount = 1).
#   treatment     'dose', or 'dose2' (two doses, 0 and 44).
#   control       The control level (prc() is given the treatment with it
#                 as the first level).
#   permutations  The name in `data` of the permutation matrix (rows
#                 restricted to `rows` are not used: cases with `rows` run
#                 no tests).
#   expected      For axes 1 and 2: summary(prc, axis, scaling)'s
#                 coefficients (treatments x weeks) and sp, for scaling
#                 1, 2 and 3 ("sites", "species", "symmetric"):
#                 coef[[scaling]][[axis]], sp[[scaling]][[axis]]; eig, the
#                 first constrained eigenvalues; and, with permutations,
#                 test_first and test_all: anova(prc, first = TRUE) and
#                 anova(prc) as list(F, p, df, residual, F_perm (from
#                 permutest())).

library(vegan)
library(permute)

data(pyrifos)
week  <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose  <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
ditch <- gl(12, 1, length = 132)

data <- list(
  pyrifos = as.matrix(pyrifos),
  week    = week,
  dose    = dose,
  ditch   = ditch )

# vegan needs the shuffleSet() objects; the fixture stores plain matrices.
set.seed(20261125)
sets <- list(
  perm_plots = shuffleSet(132, control = how(
    plots = Plots(strata = ditch, type = 'free'), within = Within(type = 'none'), nperm = 99 )),
  perm_free  = shuffleSet(132, 99),
  perm_week  = shuffleSet(132, control = how(blocks = week, nperm = 99)) )
for (nm in names(sets)) {
  p <- sets[[nm]][, , drop = FALSE]
  attributes(p) <- list(dim = dim(p))
  data[[nm]] <- p
}

# Unbalanced: 20 samples dropped, keeping control samples at every week.
set.seed(3)
repeat {
  drop <- sort(sample(132, 20))
  keep <- setdiff(seq_len(132), drop)
  if (all(table(week[keep][dose[keep] == '0']) > 0)) break
}

cases <- list(
  list(name = 'plots',     rows = NULL, transform = 'none', treatment = 'dose', control = '0', permutations = 'perm_plots'),
  list(name = 'free',      rows = NULL, transform = 'none', treatment = 'dose', control = '0', permutations = 'perm_free'),
  list(name = 'week',      rows = NULL, transform = 'none', treatment = 'dose', control = '0', permutations = 'perm_week'),
  list(name = 'hellinger', rows = NULL, transform = 'hellinger', treatment = 'dose', control = '0', permutations = 'perm_plots'),
  list(name = 'clr',       rows = NULL, transform = 'counts', treatment = 'dose', control = '0', permutations = 'perm_free'),
  list(name = 'control44', rows = NULL, transform = 'none', treatment = 'dose', control = '44', permutations = 'perm_plots'),
  list(name = 'two_doses', rows = which(dose %in% c('0', '44')), transform = 'none', treatment = 'dose', control = '0', permutations = NULL),
  list(name = 'empty_cell', rows = which(!(dose == '44' & week == '4')), transform = 'none', treatment = 'dose', control = '0', permutations = NULL),
  list(name = 'unbalanced', rows = keep, transform = 'none', treatment = 'dose', control = '0', permutations = NULL) )

for (i in seq_along(cases)) {

  cs   <- cases[[i]]
  rows <- if (is.null(cs$rows)) seq_len(132) else cs$rows
  y    <- data$pyrifos[rows, ]
  y    <- switch(
    cs$transform,
    none      = y,
    hellinger = as.matrix(decostand(y, 'hellinger')),
    counts    = { l <- log(round((exp(y) - 1) / 10) + 1); l - rowMeans(l) } )
  trt <- droplevels(relevel(dose[rows], cs$control))
  tim <- droplevels(week[rows])

  m <- prc(y, trt, tim)

  ex <- list(coef = list(), sp = list(), eig = unname(m$CCA$eig[1:2]))
  for (sc in c('sites', 'species', 'symmetric')) {
    ex$coef[[sc]] <- lapply(1:2, function (a) unclass(summary(m, axis = a, scaling = sc)$coefficients))
    ex$sp[[sc]]   <- lapply(1:2, function (a) summary(m, axis = a, scaling = sc)$sp)
  }

  if (!is.null(cs$permutations)) {
    P  <- sets[[cs$permutations]]
    a1 <- anova(m, permutations = P, first = TRUE)
    a0 <- anova(m, permutations = P)
    p1 <- permutest(m, permutations = P, first = TRUE)
    p0 <- permutest(m, permutations = P)
    ex$test_first <- list(F = a1$F[[1]], p = a1$`Pr(>F)`[[1]], df = a1$Df, residual = a1$Variance[[2]], F_perm = as.vector(p1$F.perm))
    ex$test_all   <- list(F = a0$F[[1]], p = a0$`Pr(>F)`[[1]], df = a0$Df, residual = a0$Variance[[2]], F_perm = as.vector(p0$F.perm))
  }

  cases[[i]]$expected <- ex
}

saveRDS(
  list(data = data, cases = cases),
  'tests/tinytests/expected/principal_response_curves.rds',
  version = 2, compress = 'xz' )
