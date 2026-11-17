# Expected results for tests/tinytests/test-nmds.r, from scikit-learn's
# non-metric SMACOF. Run from the package root:
#
#     Rscript data-raw/expected/nmds.r
#
# Writes tests/tinytests/expected/nmds.rds.
#
# The reference values themselves come from data-raw/expected/nmds.py (Python
# with scikit-learn 1.9.1), whose text outputs are committed under
# data-raw/expected/nmds/: rerun it only to change the cases. This script
# only converts those files, so it needs neither Python nor ecodive.
#
# nmds.rds is a list of:
#   dists   Named list of `dist` objects (bray30, bray60, dup12, line40,
#           ties50), also used by data-raw/expected/mds.r.
#   cases   Data frame, one row per case: id, case, dist, k, n_iter (the
#           iterations scikit-learn ran), skl_random and skl_cmds (the
#           stress-1 of scikit-learn's own best configurations from 8
#           random starts and from classical MDS).
#   values  List, one element per case id, of:
#             init   n x k starting configuration;
#             trace  data frame (iter, stress): the raw stress after every
#                    10th iteration and after the last;
#             x      n x k configuration after the last iteration.

smacof_fixtures <- function (dir, prefix) {

  cases  <- read.csv(file.path(dir, paste0(prefix, '-cases.csv')),  stringsAsFactors = FALSE)
  vals   <- read.csv(file.path(dir, paste0(prefix, '-values.csv')), stringsAsFactors = FALSE)
  dists  <- nmds_dists()

  field <- function (id, f) {
    v <- vals[vals$id == id & vals$field == f, ]
    v[order(v$index), c('index', 'value')]
  }

  values <- lapply(seq_len(nrow(cases)), function (i) {
    cs <- cases[i, ]
    n  <- attr(dists[[cs$dist]], 'Size')
    tr <- field(cs$id, 's')
    list(
      init  = matrix(field(cs$id, 'i')$value, n, cs$k),
      trace = data.frame(iter = as.integer(tr$index), stress = tr$value),
      x     = matrix(field(cs$id, 'x')$value, n, cs$k) )
  })

  list(dists = dists, cases = cases, values = values)
}


# Dissimilarities, one per line in dist order, as nmds.py wrote them.
nmds_dists <- function (dir = file.path('data-raw', 'expected', 'nmds')) {
  names <- c('bray30', 'bray60', 'dup12', 'line40', 'ties50')
  sapply(names, simplify = FALSE, function (name) {
    v <- as.double(readLines(file.path(dir, paste0('nmds-dist-', name, '.txt'))))
    n <- as.integer(round((1 + sqrt(1 + 8 * length(v))) / 2))
    structure(v, class = 'dist', Size = n, Diag = FALSE, Upper = FALSE)
  })
}


if (sys.nframe() == 0L) {
  saveRDS(
    smacof_fixtures(file.path('data-raw', 'expected', 'nmds'), 'nmds'),
    file.path('tests', 'tinytests', 'expected', 'nmds.rds'),
    version = 2 )
}
