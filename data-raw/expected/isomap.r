# Expected results for tests/tinytests/test-isomap.r, from scikit-learn's
# Isomap. Run from the package root:
#
#     Rscript data-raw/expected/isomap.r
#
# Converts the text files in data-raw/expected/isomap/, written by
# data-raw/expected/isomap.py (scikit-learn 1.9.1), into
# tests/tinytests/expected/isomap.rds: a list with one element per dataset.
# There is no Python where the package is developed, so the committed text
# files are the reference; this script only reads them.
#
# Each element holds the input distances `dm` (a dist object, computed here
# with stats::dist() from the coordinates when the file gives coordinates,
# with the same arithmetic as the generator), `k`, `n_neighbors` and
# `radius` (NA when not used), and scikit-learn's results: `components`,
# `membership`, `edges` (undirected, i < j, 1-based, sorted, joining edges
# included), `eig`, `embedding`, and `geo` (the packed geodesic distances,
# or NULL), plus `sklearn_same`, whether Isomap itself, end to end, agreed
# with the reference (see isomap.py).

read_fixture <- function (path) {

  lines <- readLines(path)
  lines <- lines[!startsWith(lines, '#')]
  parts <- strsplit(lines, ' ', fixed = TRUE)
  fx    <- lapply(parts, function (p) p[-1])
  names(fx) <- vapply(parts, `[[`, '', 1)

  text <- c('name', 'metric')
  fx[!names(fx) %in% text] <- lapply(fx[!names(fx) %in% text], function (x)
    if (identical(x, 'NA')) NA else as.numeric(x) )

  n  <- as.integer(fx$n)
  dm <- if (is.null(fx$dist)) {
    dist(matrix(fx$coords, nrow = n, ncol = fx$ncol), method = fx$metric)
  } else {
    structure(fx$dist, Size = n, Diag = FALSE, Upper = FALSE, class = 'dist')
  }

  list(
    name         = fx$name,
    dm           = dm,
    k            = as.integer(fx$k),
    n_neighbors  = as.integer(fx$n_neighbors),
    radius       = as.double(fx$radius),
    sklearn_same = as.logical(fx$sklearn_same),
    components   = as.integer(fx$components),
    membership   = as.integer(fx$membership),
    edges        = matrix(as.integer(fx$edges), ncol = 2, byrow = TRUE),
    eig          = fx$eig,
    embedding    = matrix(fx$embedding, nrow = n),
    geo          = fx$geo )
}

names <- c(
  'swiss', 'swiss150', 'radius', 'radius_split', 'ties', 'clusters',
  'cluster_ties', 'duplicates', 'gradient', 'zero_join' )

expected <- lapply(names, function (nm)
  read_fixture(file.path('data-raw', 'expected', 'isomap', paste0(nm, '.txt'))) )
names(expected) <- names

stopifnot(identical(unname(vapply(expected, `[[`, '', 'name')), names))

saveRDS(expected, 'tests/tinytests/expected/isomap.rds', version = 2)
