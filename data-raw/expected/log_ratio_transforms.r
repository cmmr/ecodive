# Writes tests/tinytests/expected/log_ratio_transforms.rds: inputs and
# reference results for test-log_ratio_transforms.r, so that the tests need
# no reference package.
#
# Run from the package root:
#
#     python data-raw/expected/log_ratio_transforms.py   # scikit-bio's outputs
#     Rscript data-raw/expected/log_ratio_transforms.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): alr(), ilr() (default basis, an SBP converted by
# sbp_basis(), a random orthonormal basis), _gram_schmidt_basis() and
# sbp_basis(), from the files that the Python script writes to
# data-raw/expected/log_ratio_transforms/. Every number is a hexadecimal
# float there, so it is read back exactly.
#
# compositions (GPL; output only, never its source), version 2.0.9: alr()
# and ilr() of the same tables plus their pseudocounts, as a cross-check.
# compositions' alr() divides by the last part; its default ilr basis
# (ilrBase()) is scikit-bio's with the opposite sign.

dir <- file.path('data-raw', 'expected', 'log_ratio_transforms')

hex <- function (x) as.numeric(x)

read_hex <- function (file) {
  m <- as.matrix(read.csv(file.path(dir, file), header = FALSE, colClasses = 'character'))
  matrix(hex(m), nrow(m), ncol(m))
}

read_int <- function (file)
  unname(as.matrix(read.csv(file.path(dir, file), header = FALSE)))

cases  <- read.csv(file.path(dir, 'cases.csv'),  colClasses = 'character')
values <- read.csv(file.path(dir, 'values.csv'), colClasses = 'character')

tables <- list()
for (nm in unique(cases$table))
  tables[[nm]] <- read_int(sprintf('counts-%s.csv', nm))

field <- function (id, f) {
  v <- values[values$id == id & values$field == f, ]
  if (nrow(v) == 0) return (NULL)
  r <- as.integer(v$row); k <- as.integer(v$col)
  m <- matrix(NA_real_, max(r), max(k))
  m[cbind(r, k)] <- hex(v$value)
  m
}

skbio <- lapply(seq_len(nrow(cases)), function (i) {
  cs <- cases[i, ]
  list(
    table   = cs$table,
    fn      = cs$fn,
    ref_idx = if (nzchar(cs$ref_idx)) as.integer(cs$ref_idx),
    kind    = if (nzchar(cs$kind)) cs$kind,
    pc      = hex(cs$pc),
    basis   = field(cs$id, 'basis'),
    result  = field(cs$id, 'result') )
})

helmert <- list()
for (D in c(2, 3, 4, 7, 20, 60))
  helmert[[as.character(D)]] <- read_hex(sprintf('helmert-%d.csv', D))

sbp <- list()
for (nm in c('doc', 'random20'))
  sbp[[nm]] <- list(
    sbp   = read_int(sprintf('sbp-%s.csv', nm)),
    basis = read_hex(sprintf('sbp-%s-basis.csv', nm)) )


# compositions, output only.
compositions <- NULL
if (requireNamespace('compositions', quietly = TRUE)) {
  pcs <- vapply(skbio, function (s) s$pc, 0)[!duplicated(vapply(skbio, function (s) s$table, ''))]
  names(pcs) <- unique(vapply(skbio, function (s) s$table, ''))
  compositions <- list(version = as.character(packageVersion('compositions')))
  for (nm in names(tables)) {
    x <- tables[[nm]] + pcs[[nm]]
    compositions[[nm]] <- list(
      alr = unclass(compositions::alr(compositions::acomp(x))),
      ilr = unclass(compositions::ilr(compositions::acomp(x))) )
    attributes(compositions[[nm]]$alr) <- list(dim = dim(compositions[[nm]]$alr))
    attributes(compositions[[nm]]$ilr) <- list(dim = dim(compositions[[nm]]$ilr))
  }
}

saveRDS(
  list(tables = tables, skbio = skbio, helmert = helmert, sbp = sbp,
       compositions = compositions),
  file.path('tests', 'tinytests', 'expected', 'log_ratio_transforms.rds'),
  version = 2 )
