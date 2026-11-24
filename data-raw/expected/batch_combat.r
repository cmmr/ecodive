# Writes tests/tinytests/expected/batch_combat.rds: inputs and reference
# results for test-batch_correct.r, so that the tests need no reference
# package.
#
# Run from the package root, with sva installed:
#
#     Rscript data-raw/expected/batch_combat.r
#
# The reference is sva's ComBat() (Artistic License 2.0). The committed
# fixture comes from sva 3.60.0 (Bioconductor 3.23), installed from its
# source tarball into R 4.3.3 with its dependencies from Bioconductor 3.18.
#
# Tables (samples x features): log2(counts + 1) of a simulated count table
# with planted batch shifts and a covariate, two features constant within a
# batch (left unchanged by ComBat), and a batch of one sample (mean only);
# and the clr of a count table, computed here in base R, for
# batch_combat(norm = 'clr').

suppressPackageStartupMessages(library(sva))
cat('sva', format(packageVersion('sva')), '\n')

set.seed(20261124)
n <- 45; D <- 25
meta <- data.frame(
  batch = factor(rep(c('b1', 'b2', 'b3'), times = c(15, 18, 12))),
  cond  = rep(c('ctl', 'case'), length.out = n),
  age   = round(runif(n, 20, 70)),
  row.names = sprintf('S%02d', seq_len(n)),
  stringsAsFactors = FALSE )
mu <- matrix(rexp(D, 1 / 60), n, D, byrow = TRUE)
mu[meta$batch == 'b2', 1:10]  <- mu[meta$batch == 'b2', 1:10] * 3
mu[meta$batch == 'b3', 6:14]  <- mu[meta$batch == 'b3', 6:14] / 2
mu[meta$cond == 'case', 16:20] <- mu[meta$cond == 'case', 16:20] * 2
counts <- matrix(rnbinom(n * D, mu = mu, size = 3), n, D,
                 dimnames = list(rownames(meta), sprintf('F%02d', seq_len(D))))
counts[meta$batch == 'b3', 24] <- 0           # constant within b3
logx <- log2(counts + 1)
logx[meta$batch == 'b1', 25] <- 1.5           # constant within b1

# A batch of one sample.
meta1 <- meta[1:31, ]
meta1$batch <- factor(c(as.character(meta$batch[1:30]), 'b4'))
logx1 <- logx[1:31, ]

# clr with pseudocount 1, as normalize(norm = 'clr').
lc  <- log(counts + 1)
clr <- lc - rowMeans(lc)

combat <- function (x, md, cov, par, mo = FALSE, ref = NULL) {
  mod <- if (is.null(cov)) NULL else model.matrix(stats::reformulate(cov), md)
  out <- suppressMessages(utils::capture.output(
    res <- ComBat(t(x), md$batch, mod = mod, par.prior = par, mean.only = mo, ref.batch = ref) ))
  t(res)
}

cases <- list(
  par            = list(table = 'logx', cov = NULL,            par = TRUE),
  par_cov        = list(table = 'logx', cov = c('cond', 'age'), par = TRUE),
  nonpar         = list(table = 'logx', cov = NULL,            par = FALSE),
  nonpar_cov     = list(table = 'logx', cov = 'cond',          par = FALSE),
  par_mean       = list(table = 'logx', cov = 'cond',          par = TRUE,  mo = TRUE),
  nonpar_mean    = list(table = 'logx', cov = NULL,            par = FALSE, mo = TRUE),
  par_ref        = list(table = 'logx', cov = 'cond',          par = TRUE,  ref = 'b2'),
  nonpar_ref     = list(table = 'logx', cov = NULL,            par = FALSE, ref = 'b1'),
  par_ref_mean   = list(table = 'logx', cov = NULL,            par = TRUE,  mo = TRUE, ref = 'b3'),
  single         = list(table = 'logx1', cov = NULL,           par = TRUE),
  clr            = list(table = 'clr',  cov = 'cond',          par = TRUE) )

tables <- list(logx = logx, logx1 = logx1, clr = clr)
metas  <- list(logx = meta, logx1 = meta1, clr = meta)

results <- lapply(cases, function (cs) {
  combat(
    tables[[cs$table]], metas[[cs$table]], cs$cov, cs$par,
    isTRUE(cs$mo), cs$ref )
})

saveRDS(
  list(
    version = format(packageVersion('sva')),
    tables  = tables,
    counts  = counts,
    metas   = metas,
    cases   = cases,
    results = results ),
  file.path('tests', 'tinytests', 'expected', 'batch_combat.rds'),
  version = 2, compress = 'xz' )
