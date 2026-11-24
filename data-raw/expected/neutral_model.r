# Writes tests/tinytests/expected/neutral_model.rds: inputs and reference
# results for test-neutral_model.r, so that the tests need no reference
# package.
#
# Run from the package root, after data-raw/expected/neutral_model.py:
#
#     NEUFIT_DIR=/path/to/neufit python data-raw/expected/neutral_model.py
#     TYRA_DIR=/path/to/tyRa Rscript data-raw/expected/neutral_model.r
#
# References, all run as black boxes (outputs only; their code is not read
# or ported):
#
#   * iCAMP 1.9.1 snm() (GPL; its documentation says it holds Burns et
#     al.'s (2016) original code, with metacommunity settings added), with
#     minpack.lm 1.2.4 and Hmisc 5.3.1 from CRAN.
#   * tyRa 0.1.0 fit_sncm() (github.com/DanielSprockett/tyRa, R/modeling.R
#     with the MD5 below; its documentation says the function's author is
#     Adam Burns and that it was published in Burns et al. 2016; its
#     licence file is missing, so it is run, never read).
#   * Hmisc 5.3.1 binconf(method = 'wilson') at the edges.
#
# The Python script's outputs (neufit's, and the 40-digit mpmath
# reference) are read from data-raw/expected/neutral_model/; every real
# number there is a hexadecimal float, read back exactly.

stopifnot(
  packageVersion('iCAMP')      == '1.9.1',
  packageVersion('minpack.lm') == '1.2.4',
  packageVersion('Hmisc')      == '5.3.1' )
suppressPackageStartupMessages({
  library(minpack.lm)
  library(Hmisc)
  library(stats4)
})

tyra_file <- file.path(Sys.getenv('TYRA_DIR'), 'R', 'modeling.R')
stopifnot(unname(tools::md5sum(tyra_file)) == 'bdd275a1b8e5b8817fdfa1132fb53963')
tyra <- new.env()
sys.source(tyra_file, envir = tyra)

dir <- file.path('data-raw', 'expected', 'neutral_model')
hex <- function (x) as.numeric(x)
csv <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

# Tables, samples x features, integer counts.
read_table <- function (nm) {
  x <- as.matrix(read.delim(
    file.path(dir, sprintf('table-%s.tsv', nm)), row.names = 1, check.names = FALSE ))
  x <- t(x)
  storage.mode(x) <- 'integer'
  x
}
tables <- lapply(
  setNames(nm = c('sim', 'dm', 'uneven', 'uneven_rarefied', 'sat', 'pool')),
  read_table )


# neufit's outputs.
nstats <- csv('neufit-stats.csv')
neufit <- lapply(setNames(nm = nstats$table), function (nm) {
  s <- nstats[nstats$table == nm, ]
  d <- csv(sprintf('neufit-%s.csv', nm))
  list(
    n_reads   = as.integer(s$n_reads),
    n_samples = as.integer(s$n_samples),
    m         = hex(s$m),
    m_stderr  = if (s$m_stderr == 'NA') NA_real_ else hex(s$m_stderr),
    r_squared = hex(s$r_squared),
    m_tight   = hex(s$m_tight),
    table     = data.frame(
      feature        = d$feature,
      mean_abundance = hex(d$mean_abundance),
      occurrence     = hex(d$occurrence),
      predicted      = hex(d$predicted),
      lower          = hex(d$lower),
      upper          = hex(d$upper),
      stringsAsFactors = FALSE ))
})

hp <- csv('hp-stats.csv')
hp <- data.frame(
  table   = hp$table,
  K       = as.integer(hp$K),
  m       = hex(hp$m),
  m_lower = hex(hp$m_lower),
  m_upper = hex(hp$m_upper),
  cutoff  = hex(hp$cutoff),
  stringsAsFactors = FALSE )

w <- csv('wilson.csv')
wilson_sm <- data.frame(
  n          = as.integer(w$n),
  conf_level = as.numeric(w$conf_level),
  count      = hex(w$count),
  lower      = hex(w$lower),
  upper      = hex(w$upper) )


# iCAMP's snm() and tyRa's fit_sncm() on each case.
cases <- list(
  sim             = list(counts = 'sim'),
  dm              = list(counts = 'dm'),
  uneven          = list(counts = 'uneven'),
  uneven_rarefied = list(counts = 'uneven_rarefied'),
  sat             = list(counts = 'sat'),
  dm_pool         = list(counts = 'dm', pool = 'pool') )

burns <- lapply(cases, function (cs) {
  x    <- tables[[cs$counts]]
  pool <- if (!is.null(cs$pool)) tables[[cs$pool]]
  ic   <- suppressWarnings(iCAMP::snm(x, meta.com = pool))
  ty   <- suppressWarnings(tyra$fit_sncm(x, pool = pool))
  list(
    icamp = list(stats = ic$stats, detail = data.frame(
      feature = rownames(ic$detail), ic$detail, stringsAsFactors = FALSE )),
    tyra  = list(stats = ty[[1]], detail = data.frame(
      feature = rownames(ty[[2]]), ty[[2]], stringsAsFactors = FALSE )) )
})


# Hmisc's Wilson interval at and near the edges.
wilson_hmisc <- do.call(rbind, lapply(c(10, 50), function (n) {
  do.call(rbind, lapply(c(0.95, 0.9), function (cl) {
    x <- c(0, 1e-12, 0.5, n / 3, n - 0.5, n - 1e-12, n)
    b <- Hmisc::binconf(x, n, alpha = 1 - cl, method = 'wilson')
    data.frame(n = n, conf_level = cl, count = x, lower = b[, 'Lower'], upper = b[, 'Upper'])
  }))
}))


saveRDS(
  object  = list(
    tables       = tables,
    neufit       = neufit,
    hp           = hp,
    burns        = burns,
    wilson_sm    = wilson_sm,
    wilson_hmisc = wilson_hmisc ),
  file    = file.path('tests', 'tinytests', 'expected', 'neutral_model.rds'),
  version = 2 )
