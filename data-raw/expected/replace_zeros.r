# Writes tests/tinytests/expected/replace_zeros.rds: inputs and reference
# results for test-replace_zeros.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     python data-raw/expected/replace_zeros.py   # scikit-bio's outputs
#     Rscript data-raw/expected/replace_zeros.r   # folds them in
#
# scikit-bio 0.7.4 (BSD-3): multi_replace() of every case that the Python
# script writes to data-raw/expected/replace_zeros/. Every number is a
# hexadecimal float there, so it is read back exactly.
#
# zCompositions (GPL; output only, never its source), version 1.6.2:
# multRepl() of the same tables, as a cross-check. Given closed proportions,
# frac = 1 and a detection limit of delta for every part, it
# replaces each zero with delta and scales the non-zero parts as
# multi_replace() does. A per-sample delta is a matrix `dl`, one row per
# sample. Its z.warning / z.delete options are off (`z.warning = 1`,
# `z.delete = FALSE`), as ecodive never drops a feature.

dir <- file.path('data-raw', 'expected', 'replace_zeros')

read_hex <- function (file) {
  m <- as.matrix(read.csv(file.path(dir, file), header = FALSE, colClasses = 'character'))
  matrix(as.numeric(m), nrow(m), ncol(m))
}

cases <- read.csv(file.path(dir, 'cases.csv'), colClasses = 'character')

skbio <- lapply(seq_len(nrow(cases)), function (i) {
  cs <- cases[i, ]
  list(
    id     = cs$id,
    x      = read_hex(paste0(cs$id, '-x.csv')),
    delta  = if (cs$delta == 'per_sample') as.vector(read_hex(paste0(cs$id, '-d.csv')))
             else if (nzchar(cs$delta)) as.numeric(cs$delta),
    result = read_hex(paste0(cs$id, '-r.csv')) )
})
names(skbio) <- cases$id


# zCompositions::multRepl(), output only. ====

zcomp <- NULL
if (requireNamespace('zCompositions', quietly = TRUE)) {

  stopifnot(as.character(packageVersion('zCompositions')) == '1.6.2')

  zcomp <- lapply(skbio, function (cs) {
    x <- cs$x
    if (!any(x == 0) || any(colSums(x) == 0)) return (NULL) # multRepl() needs zeros, and no empty part
    p <- x / rowSums(x)
    D <- ncol(x)
    d <- if (is.null(cs$delta)) (1 / D)^2 else cs$delta
    dl <- if (length(d) > 1) matrix(d, nrow(x), D) else rep(d, D)
    r <- suppressMessages(zCompositions::multRepl(
      p, label = 0, dl = dl, frac = 1, z.warning = 1, z.delete = FALSE ))
    unname(as.matrix(r))
  })
  zcomp <- zcomp[!vapply(zcomp, is.null, logical(1))]
}


saveRDS(
  list(skbio = skbio, zcomp = zcomp),
  file.path('tests', 'tinytests', 'expected', 'replace_zeros.rds'),
  version = 2 )
