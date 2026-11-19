# Writes tests/tinytests/expected/bioenv.rds: inputs and reference results for
# test-bioenv.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/bioenv.r      # writes scikit-bio's inputs
#     python data-raw/expected/bioenv.py      # scikit-bio's outputs
#     Rscript data-raw/expected/bioenv.r      # folds them in
#
# vegan is GPL: only its outputs are used, never its source.
# Generated with vegan 2.7.8 and permute 0.10.0.
#
# Inputs: varespec's Bray-Curtis dissimilarities (vegdist()'s default, as
# ecodive's bray()), as `dm`; varechem's 14 variables, as `env`; and
# Euclidean distances on its pH, as `given` (vegan's `partial`).
#
# Each case holds:
#   method, metric, max_vars, partial, env   The arguments: vegan's method,
#                     metric, upto, partial (TRUE: `given`) and the
#                     variables' table (`env`, unless given here).
#   expected          summary() of vegan::bioenv(): size, variables (joined
#                     by ", ") and correlation, as a data frame.
#
# vegan refuses method = "kendall" with `partial`, so those cases are left
# out; the tests check them against a direct computation instead.
#
# And `skbio`: scikit-bio 0.7.4's bioenv() (BSD-3) on the default options
# (Spearman, Euclidean), of the same dm and env, from the CSV file that
# data-raw/expected/bioenv.py writes to data-raw/expected/bioenv/.

library(vegan)

data(varespec, varechem)

dm    <- vegdist(varespec)
env   <- varechem
given <- dist(varechem$pH)

dir <- file.path('data-raw', 'expected', 'bioenv')
dir.create(dir, showWarnings = FALSE)

write.csv(as.matrix(dm), file.path(dir, 'dm.csv'))
write.csv(env, file.path(dir, 'env.csv'))

vegan_case <- function (method, metric, max_vars = ncol(env), partial = FALSE, env_ = NULL) {
  e <- if (is.null(env_)) env else env_
  b <- suppressMessages(bioenv(
    dm, e, method = method, metric = metric, upto = max_vars,
    partial = if (partial) given ))
  s <- summary(b)
  list(
    method = method, metric = metric, max_vars = max_vars, partial = partial,
    env = env_,
    expected = data.frame(
      size        = s$size,
      variables   = gsub(' ', ', ', s$variables, fixed = TRUE),
      correlation = s$correlation,
      stringsAsFactors = FALSE ))
}

cases <- list()
for (method in c('spearman', 'pearson', 'kendall'))
  for (metric in c('euclidean', 'mahalanobis', 'manhattan', 'gower'))
    for (partial in c(FALSE, TRUE)) {
      if (partial && method == 'kendall') next
      cases[[length(cases) + 1]] <- vegan_case(method, metric, partial = partial)
      cases[[length(cases) + 1]] <- vegan_case(method, metric, 3L, partial)
    }

# Ties: a copy of a variable ties with it in every subset it is in, so the
# first (in the order of the columns) is kept.
tied <- env[, c('N', 'P', 'K', 'Ca', 'pH')]
tied$Kcopy <- tied$K
tied <- tied[, c('Kcopy', 'N', 'P', 'K', 'Ca', 'pH')]
cases[[length(cases) + 1]] <- vegan_case('spearman', 'euclidean', env_ = tied)
cases[[length(cases) + 1]] <- vegan_case('pearson',  'manhattan', env_ = tied)


skbio <- NULL
f <- file.path(dir, 'skbio.csv')
if (file.exists(f)) {
  s <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  skbio <- data.frame(
    size        = as.integer(s$size),
    variables   = s$vars,
    correlation = as.double(s$correlation),
    stringsAsFactors = FALSE )
}
if (is.null(skbio)) message('Run data-raw/expected/bioenv.py, then this script again.')

saveRDS(
  list(dm = dm, env = env, given = given, cases = cases, skbio = skbio),
  file.path('tests', 'tinytests', 'expected', 'bioenv.rds'), version = 2 )
