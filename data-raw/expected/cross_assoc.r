# Writes tests/tinytests/expected/cross_assoc.rds: inputs and reference
# results for test-cross_assoc.r, so that the tests need no reference
# package.
#
# Run from the package root, with ecodive installed (for normalize() only):
#
#     Rscript data-raw/expected/cross_assoc.r    # writes SciPy's inputs
#     python data-raw/expected/cross_assoc.py    # SciPy's outputs
#     Rscript data-raw/expected/cross_assoc.r    # folds them in
#
# Hmisc, psych, ppcor and microeco are GPL: only their outputs are used,
# never their source. Generated with R 4.3.3, Hmisc 5.3-1, psych 2.6.9,
# ppcor 1.1 and microeco (version in `versions`), and SciPy 1.18.1.
#
# Inputs (simulated): `counts`, 40 samples x 25 taxa (negative-binomial
# counts with many zeros; taxon 25 is all zero), and `env`, 6 variables:
#   ph     follows taxon 1 (a strong association);
#   temp   rounded, so with ties;
#   salt   3 NAs;
#   const  constant (NA correlations);
#   depth  integers 1 to 5 (many ties), 2 NAs (a different pattern);
#   od     unrelated, no NAs.
# and `z`, two covariates (age, numeric; sex, 0/1) for partial correlations.
#
# The relative abundances that cross_assoc() correlates (norm = 'percent')
# are ecodive::normalize(counts, 'percent'), written as hexadecimal floats
# for SciPy.
#
# Fields (each a 25 x 6 matrix, taxa x variables, unless said):
#   scipy     spearman, pearson, kendall and their *_p: scipy.stats.spearmanr(),
#             pearsonr(), kendalltau(method = 'asymptotic') on each pair's
#             complete samples; n.
#   cortest   The same from stats::cor() and cor.test(exact = FALSE)
#             (Kendall: stats::cor.test()'s own tie-corrected variance).
#   hmisc     Hmisc::rcorr(cbind(rel, env), type = 'spearman'): r, n, P (the
#             25 x 6 block), with the constant column left out (rcorr()
#             refuses nothing, but gives NaN).
#   psych     psych::corr.test(rel, env, method = 'spearman', adjust = 'BH'):
#             r, n, p, p.adj.
#   ppcor     ppcor::pcor.test(rel[, i], env[, j], z, method = m) for
#             m = 'spearman', 'pearson', on the complete samples of env's
#             columns without NAs (ph, temp, od): estimate, p.value,
#             statistic (3 x 25 x 2 values each, as 25 x 3 matrices).
#   microeco  trans_env$cal_cor(use_data = 'Genus', method = 'spearman',
#             p_adjust_type = 'All', p_adjust_method = 'fdr') on ph, temp
#             and od (its relative abundances of taxa 1-24): its res_cor
#             table, as a data frame. It calls cor.test() with its default
#             exact = NULL, so its p-values are exact (AS 89) for the taxa
#             without ties (no zeros).

suppressPackageStartupMessages(library(ecodive))

dir <- file.path('data-raw', 'expected', 'cross_assoc')
dir.create(dir, showWarnings = FALSE)

set.seed(20261124)
n <- 40; p <- 25
mu     <- exp(rnorm(p, 2, 1))
counts <- matrix(
  rnbinom(n * p, mu = rep(mu, each = n), size = 0.6), n, p,
  dimnames = list(sprintf('S%02d', seq_len(n)), sprintf('Taxon%02d', seq_len(p))) )
counts[, 1] <- counts[, 1] + rpois(n, 5)
counts[, 25] <- 0L
stopifnot(all(rowSums(counts) > 0))

rel <- normalize(counts, 'percent')

env <- data.frame(
  ph    = as.vector(scale(log1p(rel[, 1] * 100))) + rnorm(n, sd = 0.5),
  temp  = round(rnorm(n, 15, 2)),
  salt  = rnorm(n, 30, 3),
  const = 7,
  depth = sample(1:5, n, TRUE),
  od    = rnorm(n),
  row.names = rownames(counts) )
env$salt[c(4, 17, 33)] <- NA
env$depth[c(9, 17)]    <- NA

z <- data.frame(age = round(rnorm(n, 50, 10)), sex = rbinom(n, 1, 0.5),
                row.names = rownames(counts))

hexmat <- function (m) {
  out <- matrix(sprintf('%a', m), nrow(m), ncol(m), dimnames = list(NULL, colnames(m)))
  out[is.na(m)] <- 'NA'
  out
}
write.csv(hexmat(rel), file.path(dir, 'x.csv'), row.names = FALSE, quote = FALSE)
write.csv(hexmat(as.matrix(env)), file.path(dir, 'y.csv'), row.names = FALSE, quote = FALSE)

scipy_file <- file.path(dir, 'scipy.csv')
if (!file.exists(scipy_file)) {
  message('Now run: python data-raw/expected/cross_assoc.py; then this script again.')
  quit(save = 'no')
}


#--------------------------------------------------
# SciPy.
#--------------------------------------------------
sc  <- read.csv(scipy_file, colClasses = 'character')
q   <- ncol(env)
mat <- function (v) {
  v <- as.numeric(ifelse(v == 'nan', NA, v))
  m <- matrix(NA_real_, p, q, dimnames = list(colnames(rel), names(env)))
  m[cbind(as.integer(sc$i), as.integer(sc$j))] <- v
  m
}
scipy <- lapply(
  c(spearman = 'spearman', spearman_p = 'spearman_p', pearson = 'pearson',
    pearson_p = 'pearson_p', kendall = 'kendall', kendall_p = 'kendall_p', n = 'n'),
  function (f) mat(sc[[f]]) )


#--------------------------------------------------
# stats::cor.test(exact = FALSE).
#--------------------------------------------------
cortest <- list()
for (m in c('spearman', 'pearson', 'kendall')) {
  r <- pv <- matrix(NA_real_, p, q, dimnames = list(colnames(rel), names(env)))
  for (i in seq_len(p)) for (j in seq_len(q)) {
    ok <- !is.na(env[[j]])
    a  <- rel[ok, i]; b <- env[[j]][ok]
    if (length(unique(a)) < 2 || length(unique(b)) < 2) next
    ct <- suppressWarnings(cor.test(a, b, method = m, exact = FALSE))
    r[i, j]  <- ct$estimate
    pv[i, j] <- ct$p.value
  }
  cortest[[m]] <- r
  cortest[[paste0(m, '_p')]] <- pv
}


#--------------------------------------------------
# Hmisc, psych, ppcor.
#--------------------------------------------------
versions <- c(R = as.character(getRversion()))

if (requireNamespace('Hmisc', quietly = TRUE)) {
  versions[['Hmisc']] <- as.character(packageVersion('Hmisc'))
  e  <- as.matrix(env[, names(env) != 'const'])
  rc <- Hmisc::rcorr(cbind(rel[, -25], e), type = 'spearman')
  blk <- function (m) m[colnames(rel)[-25], colnames(e)]
  hmisc <- list(r = blk(rc$r), n = blk(rc$n), P = blk(rc$P))
}

if (requireNamespace('psych', quietly = TRUE)) {
  versions[['psych']] <- as.character(packageVersion('psych'))
  ct <- suppressWarnings(psych::corr.test(
    rel[, -25], env[, names(env) != 'const'], method = 'spearman', adjust = 'BH', ci = FALSE ))
  psych <- list(r = ct$r, n = ct$n, p = ct$p, p.adj = ct$p.adj)
}

if (requireNamespace('ppcor', quietly = TRUE)) {
  versions[['ppcor']] <- as.character(packageVersion('ppcor'))
  vars  <- c('ph', 'temp', 'od')
  ppcor <- list()
  for (m in c('spearman', 'pearson')) {
    est <- pv <- st <- matrix(NA_real_, p - 1, length(vars),
                              dimnames = list(colnames(rel)[-25], vars))
    for (i in seq_len(p - 1)) for (j in vars) {
      pt_ <- ppcor::pcor.test(rel[, i], env[[j]], z, method = m)
      est[i, j] <- pt_$estimate; pv[i, j] <- pt_$p.value; st[i, j] <- pt_$statistic
    }
    ppcor[[m]] <- list(estimate = est, p.value = pv, statistic = st)
  }
}


#--------------------------------------------------
# microeco: Spearman with p-values adjusted over all pairs.
#--------------------------------------------------
microeco <- NULL
if (requireNamespace('microeco', quietly = TRUE)) {
  versions[['microeco']] <- as.character(packageVersion('microeco'))
  otu <- as.data.frame(t(counts[, -25]))
  tax <- data.frame(
    Kingdom = 'k__Bacteria', Phylum = 'p__P', Class = 'c__C', Order = 'o__O',
    Family = 'f__F', Genus = paste0('g__', rownames(otu)), row.names = rownames(otu) )
  sam <- data.frame(SampleID = colnames(otu), row.names = colnames(otu))
  mt  <- microeco::microtable$new(otu_table = otu, tax_table = tax, sample_table = sam)
  mt$tidy_dataset()
  mt$cal_abund()
  e  <- env[, c('ph', 'temp', 'od')]
  te <- microeco::trans_env$new(dataset = mt, add_data = e)
  te$cal_cor(use_data = 'Genus', method = 'spearman', taxa_name_full = FALSE,
             p_adjust_method = 'fdr', p_adjust_type = 'All')
  microeco <- as.data.frame(te$res_cor)
}


saveRDS(
  list(counts = counts, env = env, z = z, scipy = scipy, cortest = cortest,
       hmisc = hmisc, psych = psych, ppcor = ppcor, microeco = microeco,
       versions = versions),
  file.path('tests', 'tinytests', 'expected', 'cross_assoc.rds'),
  version = 2 )
