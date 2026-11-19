# Writes tests/tinytests/expected/log_ratio.rds: inputs and reference results
# for test-log_ratio.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     python data-raw/expected/log_ratio.py   # scikit-bio's and Qurro's outputs
#     Rscript data-raw/expected/log_ratio.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): ilr() with sbp_basis() of a sequential binary
# partition whose first row is the balance, and the closed form in numpy.
# Qurro 0.9.0 (BSD-3): qarcoal()'s log ratio of summed counts. All from the
# files that the Python script writes to data-raw/expected/log_ratio/, where
# every number is a hexadecimal float, read back exactly.
#
# compositions (GPL; output only, never its source), version 2.0.9:
# balance(acomp(x + pc), ~ (F1*F2)/(F3*F4)) of the same cases, the same
# balance with the same sign, as a cross-check.

dir <- file.path('data-raw', 'expected', 'log_ratio')

hex <- function (x) suppressWarnings(as.numeric(x))  # 'NA' -> NA

read_table <- function (nm) {
  m <- as.matrix(read.csv(
    file.path(dir, sprintf('table-%s.csv', nm)),
    header = FALSE, colClasses = 'character' ))
  m <- matrix(hex(m), nrow(m), ncol(m))
  if (all(m == round(m))) storage.mode(m) <- 'integer'
  m
}

cases  <- read.csv(file.path(dir, 'cases.csv'),  colClasses = 'character')
values <- read.csv(file.path(dir, 'values.csv'), colClasses = 'character')

tables <- list()
for (nm in unique(cases$table))
  tables[[nm]] <- read_table(nm)

ints <- function (s) as.integer(strsplit(s, ' ')[[1]])

field <- function (id, f) {
  v <- values[values$id == id & values$field == f, ]
  out <- rep(NA_real_, max(as.integer(v$sample)))
  out[as.integer(v$sample)] <- hex(v$value)
  out
}

ref <- lapply(seq_len(nrow(cases)), function (i) {
  cs <- cases[i, ]
  list(
    table  = cs$table,
    num    = ints(cs$num),
    den    = ints(cs$den),
    pc     = hex(cs$pc),
    skbio  = field(cs$id, 'skbio'),
    closed = field(cs$id, 'closed'),
    qurro  = field(cs$id, 'qurro'),
    lr_pc1 = field(cs$id, 'lr_pc1') )
})


# compositions, output only: balance(acomp(x + pc), ~ num / den), with the
# parts named and grouped by `*` in the formula.
compositions <- NULL
if (requireNamespace('compositions', quietly = TRUE)) {
  compositions <- list(version = as.character(packageVersion('compositions')))
  compositions$values <- lapply(ref, function (cs) {
    x <- tables[[cs$table]] + cs$pc
    colnames(x) <- paste0('F', seq_len(ncol(x)))
    f <- stats::as.formula(sprintf(
      '~ (%s) / (%s)',
      paste0('F', cs$num, collapse = '*'),
      paste0('F', cs$den, collapse = '*') ))
    as.vector(unclass(compositions::balance(compositions::acomp(x), f)))
  })
}


saveRDS(
  list(tables = tables, cases = ref, compositions = compositions),
  file.path('tests', 'tinytests', 'expected', 'log_ratio.rds'),
  version = 2 )
