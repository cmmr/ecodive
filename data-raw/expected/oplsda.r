# Writes tests/tinytests/expected/oplsda.rds: inputs and reference results
# for test-oplsda.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/oplsda.r      # writes the inputs
#     python data-raw/expected/oplsda.py      # pyopls's results
#     Rscript data-raw/expected/oplsda.r      # folds them in
#
# References:
#
#   pyopls    pyopls 20.3.post1's OPLS (MIT), which oplsda() ports, fitted
#             to the normalized table and the 0/1 indicator of the second
#             group, then scikit-learn 1.9.1's PLSRegression(1, scale =
#             False) fitted to the filtered table, from the CSV files that
#             data-raw/expected/oplsda.py writes to
#             data-raw/expected/oplsda/.
#   ropls     ropls::opls(predI = 1, orthoI = k) (CeCILL; outputs only,
#             never its source), fitted to the normalized table:
#             scores, loadings, weights, VIP (vipVn), R2X
#             per component and R2Y(cum), which ropls rounds to 3
#             significant digits; and, unrounded, R2X from its scores,
#             loadings and scaled table, and R2Y from its fitted values. Generated with ropls 1.34.0 (Bioconductor 3.18, the
#             newest for R 4.3), as available where these fixtures were
#             made, not the plan's 1.44.
#
# Each case holds:
#   name, counts, groups      The table (samples in rows) and its groups.
#   norm, pseudocount, scale  oplsda()'s arguments; the references get the
#                             table after the same normalization, in base R.
#   k                         Orthogonal components requested.
#   x                         The normalized table given to the references.
#   expected$pyopls, expected$ropls
#                             Reference results (NULL where not run).

library(vegan)

data(dune, dune.env)

dir.create(pydir <- file.path('data-raw', 'expected', 'oplsda'), showWarnings = FALSE)

dune <- as.matrix(dune)

set.seed(20261119)
n   <- 110
big <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.4), n, 40)
dimnames(big) <- list(sprintf('S%03d', seq_len(n)), sprintf('F%02d', 1:40))
site <- factor(sample(c('a', 'b'), n, TRUE))
big[site == 'b', 2:6] <- big[site == 'b', 2:6] + 2
big[, 10:14] <- big[, 10:14] + rpois(n, 4)   # variation unrelated to site


# The normalizations oplsda() offers, in base R.
transform <- function (m, norm, pseudocount = NULL) {
  m <- m + 0
  switch(
    EXPR      = norm,
    none      = m,
    hellinger = sqrt(m / rowSums(m)),
    clr       = { l <- log(m + pseudocount); l - rowMeans(l) } )
}

cs <- function (name, counts, groups, k, norm = 'none', scale = TRUE,
                pseudocount = NULL, ropls = TRUE)
  list(name = name, counts = counts, groups = groups, norm = norm,
       pseudocount = pseudocount, scale = scale, k = k,
       x = transform(counts, norm, pseudocount), run_ropls = ropls)

two   <- factor(ifelse(dune.env$Management == 'NM', 'natural', 'managed'))
hay   <- factor(ifelse(dune.env$Use == 'Hayfield', 'hay', 'other'))
moist <- factor(ifelse(as.integer(dune.env$Moisture) > 2, 'wet', 'dry'))

cases <- list(
  cs('dune_two',       dune, two,   2),
  cs('dune_one',       dune, two,   1),
  cs('dune_const',     cbind(dune, Const = 3, Zero = 0), two, 2, ropls = FALSE),
  cs('dune_hellinger', dune, hay,   2, norm = 'hellinger', scale = FALSE),
  cs('dune_clr',       dune, moist, 3, norm = 'clr', pseudocount = 1),
  cs('big_site',       big,  site,  4) )


# Inputs for oplsda.py.
manifest <- data.frame(
  name  = vapply(cases, `[[`, '', 'name'),
  scale = vapply(cases, `[[`, TRUE, 'scale'),
  k     = vapply(cases, `[[`, 0, 'k') )
write.csv(manifest, file.path(pydir, 'cases.csv'), row.names = FALSE)
for (c in cases) {
  write.csv(c$x, file.path(pydir, paste0(c$name, '_x.csv')))
  write.csv(data.frame(group = as.integer(c$groups) - 1L, row.names = rownames(c$x)),
            file.path(pydir, paste0(c$name, '_groups.csv')))
}


# pyopls's results, if oplsda.py has run.
rd <- function (f) {
  f <- file.path(pydir, f)
  if (!file.exists(f)) return (NULL)
  unname(as.matrix(read.csv(f, header = FALSE)))
}
py <- function (name) {
  t_ortho <- rd(paste0(name, '_py_t_ortho.csv'))
  if (is.null(t_ortho)) return (NULL)
  list(
    t_ortho   = t_ortho,
    w_ortho   = rd(paste0(name, '_py_w_ortho.csv')),
    p_ortho   = rd(paste0(name, '_py_p_ortho.csv')),
    t_pred    = rd(paste0(name, '_py_t_pred.csv')),
    w_pred    = rd(paste0(name, '_py_w_pred.csv')),
    p_pred    = rd(paste0(name, '_py_p_pred.csv')),
    q_pred    = drop(rd(paste0(name, '_py_q_pred.csv'))),
    r2y       = drop(rd(paste0(name, '_py_r2y.csv'))) )
}


# ropls's results.
rop <- function (c) {
  if (!c$run_ropls) return (NULL)
  suppressMessages(library(ropls))
  m <- ropls::opls(
    c$x, c$groups, predI = 1, orthoI = c$k,
    scaleC = if (c$scale) 'standard' else 'center',
    crossvalI = min(7, nrow(c$x)), permI = 0,
    fig.pdfC = 'none', info.txtC = 'none' )
  scores   <- unname(cbind(ropls::getScoreMN(m), ropls::getScoreMN(m, orthoL = TRUE)))
  loadings <- unname(cbind(ropls::getLoadingMN(m), ropls::getLoadingMN(m, orthoL = TRUE)))
  # Unrounded R2X and R2Y from ropls' own outputs: each component's
  # |t|^2 |p|^2 over the sum of squares of its scaled table (xModelMN); and
  # the R^2 of its fitted values (yPreMN) for y, unscaled from yModelMN.
  sup <- m@suppLs
  ym  <- sup$yModelMN
  y   <- ym * attr(ym, 'scaled:scale') + attr(ym, 'scaled:center')
  list(
    version       = as.character(packageVersion('ropls')),
    scores        = scores,
    loadings      = loadings,
    weights       = unname(cbind(ropls::getWeightMN(m), ropls::getWeightMN(m, orthoL = TRUE))),
    vip           = unname(ropls::getVipVn(m)),
    r2x           = unname(m@modelDF[, 'R2X']),
    r2x_cum       = unname(m@summaryDF[, 'R2X(cum)']),
    r2y_cum       = unname(m@summaryDF[, 'R2Y(cum)']),
    r2x_exact     = colSums(scores^2) * colSums(loadings^2) / sum(sup$xModelMN^2),
    r2y_exact     = 1 - sum((y - sup$yPreMN)^2) / sum((y - mean(y))^2) )
}


for (i in seq_along(cases)) {
  c <- cases[[i]]
  c$expected <- list(pyopls = py(c$name), ropls = rop(c))
  c$run_ropls <- NULL
  cases[[i]] <- c
}

if (is.null(cases[[1]]$expected$pyopls))
  message('No pyopls results yet: run data-raw/expected/oplsda.py, then this script again.')

saveRDS(
  list(cases = cases),
  file.path('tests', 'tinytests', 'expected', 'oplsda.rds'),
  version = 2 )
