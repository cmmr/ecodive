# Writes tests/tinytests/expected/ancom.rds: inputs and reference results
# for test-ancom.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     python data-raw/expected/ancom.py   # scikit-bio's outputs
#     Rscript data-raw/expected/ancom.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): ancom() on each table plus its default
# pseudocount, with SciPy's f_oneway() and kruskal() (sig_test 'anova' and
# 'kruskal') and every p_adjust; each pair's statistic and p-value; each
# feature's adjusted p-values (tables of up to 15 features); its cutoff
# rule on crafted W; and SciPy's F and chi-squared tail probabilities on a
# grid. Every number in the files that the Python script writes to
# data-raw/expected/ancom/ is a hexadecimal float, read back exactly.

dir <- file.path('data-raw', 'expected', 'ancom')

hex <- function (x) as.numeric(x)

read_table <- function (nm) {
  m <- as.matrix(read.csv(
    file.path(dir, sprintf('table-%s.csv', nm)),
    header = FALSE, colClasses = 'character' ))
  m <- matrix(hex(m), nrow(m), ncol(m))
  if (all(m == round(m))) storage.mode(m) <- 'integer'
  m
}

csv <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

tabs     <- csv('tables.csv')
results  <- csv('results.csv')
pairs    <- csv('pairs.csv')
adjusted <- csv('adjusted.csv')
dist     <- csv('dist.csv')
cutoff   <- csv('cutoff.csv')

words <- function (x) strsplit(x, ' ', fixed = TRUE)[[1]]

ref <- lapply(setNames(tabs$table, tabs$table), function (nm) {

  r <- list(
    counts = read_table(nm),
    pc     = hex(tabs$pc[tabs$table == nm]),
    groups = words(tabs$groups[tabs$table == nm]) )

  for (st in c('anova', 'kruskal')) {
    p <- pairs[pairs$table == nm & pairs$sig_test == st, ]
    p <- p[order(as.integer(p$k)), ]
    r[[st]] <- list(statistic = hex(p$statistic), p_value = hex(p$p_value))

    a <- adjusted[adjusted$table == nm & adjusted$sig_test == st, ]
    if (nrow(a)) {
      D <- ncol(r$counts)
      r[[st]]$adjusted <- lapply(split(a, a$p_adjust), function (b) {
        m <- matrix(NA_real_, D, D)
        m[cbind(as.integer(b$i), as.integer(b$j))] <- hex(b$value)
        m
      })
    }
  }
  r
})

runs <- data.frame(
  table    = results$table,
  sig_test = results$sig_test,
  p_adjust = results$p_adjust,
  alpha    = hex(results$alpha),
  tau      = hex(results$tau),
  theta    = hex(results$theta),
  stringsAsFactors = FALSE )
runs$W      <- lapply(results$W,      function (x) as.integer(words(x)))
runs$signif <- lapply(results$signif, function (x) as.logical(as.integer(words(x))))

dists <- data.frame(
  dist = dist$dist,
  df1  = as.integer(dist$df1),
  df2  = as.integer(dist$df2),
  x    = hex(dist$x),
  p    = hex(dist$p),
  stringsAsFactors = FALSE )

cuts <- data.frame(
  D     = as.integer(cutoff$D),
  tau   = hex(cutoff$tau),
  theta = hex(cutoff$theta),
  quirk = cutoff$quirk == '1',
  stringsAsFactors = FALSE )
cuts$W      <- lapply(cutoff$W,      function (x) as.integer(words(x)))
cuts$signif <- lapply(cutoff$signif, function (x) as.logical(as.integer(words(x))))

saveRDS(
  list(tables = ref, runs = runs, dist = dists, cutoff = cuts),
  file.path('tests', 'tinytests', 'expected', 'ancom.rds'),
  version = 2 )
