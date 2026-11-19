# Writes tests/tinytests/expected/proportionality.rds: inputs and reference
# results for test-proportionality.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     python data-raw/expected/proportionality.py   # scikit-bio's outputs
#     Rscript data-raw/expected/proportionality.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): pairwise_vlr() and vlr() of each table plus its
# default pseudocount, from the files that the Python script writes to
# data-raw/expected/proportionality/, where every number is a hexadecimal
# float, read back exactly.
#
# rho_p() and phi_s() are tested against base-R versions of their published
# formulas in the test file itself. propr (GPL; output only) could not be
# installed when the fixtures were made, so no propr outputs are stored.

dir <- file.path('data-raw', 'expected', 'proportionality')

hex <- function (x) as.numeric(x)

read_table <- function (nm) {
  m <- as.matrix(read.csv(
    file.path(dir, sprintf('table-%s.csv', nm)),
    header = FALSE, colClasses = 'character' ))
  m <- matrix(hex(m), nrow(m), ncol(m))
  if (all(m == round(m))) storage.mode(m) <- 'integer'
  m
}

tabs   <- read.csv(file.path(dir, 'tables.csv'), colClasses = 'character')
values <- read.csv(file.path(dir, 'values.csv'), colClasses = 'character')

field <- function (nm, f) {
  v <- values[values$table == nm & values$field == f, ]
  hex(v$value)[order(as.integer(v$k))]
}

ref <- lapply(setNames(tabs$table, tabs$table), function (nm) list(
  counts   = read_table(nm),
  pc       = hex(tabs$pc[tabs$table == nm]),
  pairwise = field(nm, 'pairwise'),
  vlr      = field(nm, 'vlr') ))

saveRDS(
  ref, file.path('tests', 'tinytests', 'expected', 'proportionality.rds'),
  version = 2 )
