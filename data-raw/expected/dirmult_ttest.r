# Writes tests/tinytests/expected/dirmult_ttest.rds: inputs and reference
# results for test-dirmult_ttest.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     python data-raw/expected/dirmult_ttest.py   # scikit-bio's outputs
#     Rscript data-raw/expected/dirmult_ttest.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): dirmult_ttest() on given Gamma draws (NumPy's,
# written out), with every p_adjust; each draw's Welch statistics; its run
# with its own draws on a deep table; and SciPy's stdtr() and stdtrit() on
# a grid. Every real number in the files that the Python script writes to
# data-raw/expected/dirmult_ttest/ is a hexadecimal float, read back
# exactly.

dir <- file.path('data-raw', 'expected', 'dirmult_ttest')

hex   <- function (x) as.numeric(x)
words <- function (x) strsplit(x, ' ', fixed = TRUE)[[1]]
csv   <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

read_table <- function (nm) {
  m <- as.matrix(read.csv(file.path(dir, sprintf('table-%s.csv', nm)), header = FALSE))
  dimnames(m) <- NULL
  storage.mode(m) <- 'integer'
  m
}

cases   <- csv('cases.csv')
results <- csv('results.csv')
perdraw <- csv('perdraw.csv')
sampled <- csv('sampled.csv')
tdist   <- csv('tdist.csv')

tables <- lapply(setNames(nm = unique(c(cases$table, sampled$table))), read_table)

cols <- c('t', 'log2_fc', 'ci_lower', 'ci_upper', 'p_value', 'p_adjusted')

as_result <- function (r) {
  out <- lapply(setNames(nm = cols), function (k) hex(words(r[[k]])))
  out$reject <- words(r$reject) == '1'
  out
}

ref <- lapply(setNames(nm = cases$case), function (cs) {

  cc <- cases[cases$case == cs, ]
  x  <- tables[[cc$table]]
  nd <- as.integer(cc$n_draws)

  # The Gamma draws, n x D x n_draws.
  g <- csv(sprintf('gamma-%s.csv', cs))
  G <- array(NA_real_, c(nrow(x), ncol(x), nd))
  for (k in seq_len(nrow(g)))
    G[as.integer(g$sample[k]) + 1L, , as.integer(g$draw[k]) + 1L] <- hex(unlist(g[k, -(1:2)]))
  stopifnot(!anyNA(G))

  r <- results[results$case == cs, ]
  p <- perdraw[perdraw$case == cs, ]
  draws <- lapply(
    setNames(nm = c('diff', 'se', 'dof', 't', 'p_tail', 'lower', 'upper')),
    function (k) {
      m <- matrix(NA_real_, nd, ncol(x))
      m[cbind(as.integer(p$draw) + 1L, as.integer(p$feature) + 1L)] <- hex(p[[k]])
      m
    })

  list(
    table     = cc$table,
    groups    = words(cc$groups),
    contrast  = if (nzchar(cc$treatment))
                  c(cc$treatment, if (nzchar(cc$reference)) cc$reference else NA),
    n_draws   = nd,
    gamma     = G,
    results   = lapply(setNames(nm = r$p_adjust), function (m) as_result(r[r$p_adjust == m, ])),
    draws     = draws )
})

samp <- list(
  table    = sampled$table,
  groups   = words(sampled$groups),
  contrast = c(sampled$treatment, sampled$reference),
  result   = as_result(sampled[1, ]) )

saveRDS(
  object   = list(
    tables  = tables,
    cases   = ref,
    sampled = samp,
    tdist   = data.frame(
      what  = tdist$what,
      df    = hex(tdist$df),
      x     = hex(tdist$x),
      value = hex(tdist$value),
      stringsAsFactors = FALSE ) ),
  file     = file.path('tests', 'tinytests', 'expected', 'dirmult_ttest.rds'),
  version  = 2 )
