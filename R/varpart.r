# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Variation Partitioning
#'
#' Splits the variation in a community among two to four sets of explanatory
#' variables (such as environment, space, and time): the fractions that each
#' set explains on its own, the fractions they explain jointly, and what
#' none of them explains. Each fraction is an adjusted \eqn{R^2} from
#' (partial) redundancy analysis of a count table ([rda()]), distance-based
#' redundancy analysis of distances ([dbrda()]), or canonical correspondence
#' analysis ([cca()]), as in `vegan::varpart()`, with the same results.
#'
#' @inheritParams documentation
#' @inheritParams rda
#' @inheritParams cca
#' @inherit pseudocount_section
#' @inheritSection documentation Input Types
#'
#' @param x  The community: a `dist` object of distances among the samples,
#'        such as from [beta_div()] or [bray()], for db-RDA; or a count
#'        table (samples \eqn{\times} features, or any of the table types
#'        listed under "Input Types"), for RDA, or with `chisquare = TRUE`
#'        for CCA. A plain matrix is always a table; give distances as a
#'        `dist` object (`as.dist(m)`).
#'
#' @param ...  Two to four sets of explanatory variables, each one of:
#'   * a one-sided formula over the columns of `data`, such as
#'     `~ ph + temp`;
#'   * a data frame, whose factor and character columns are coded as dummy
#'     variables, as `model.matrix()` codes them;
#'   * a numeric matrix, used as it is (such as spatial eigenvectors); or
#'   * a vector or factor (one variable).
#'
#'   The kinds may be mixed. Rows of a data frame, matrix, or named vector
#'   are matched to the samples of `x` by name when both have names (other
#'   rows are ignored); without names, they must be in the samples' order.
#'   Name the sets to label them, as in `varpart(x, env = ~ ph + temp, space
#'   = pcnm)`; otherwise they are `X1`, `X2`, and so on.
#'
#' @param data  Sample metadata for the sets given as formulas: a data frame
#'        with one row per sample, matched to the samples of `x` by row
#'        name (or, without row names, in the same order). Default: `NULL`.
#'
#' @param norm  Normalize the count table first, as in [rda()]: `'none'`,
#'        `'percent'`, `'chord'`, `'hellinger'`, `'binary'`, `'clr'`, or
#'        `'rclr'`. `vegan::varpart(transfo = 'hel')` is `norm =
#'        'hellinger'`. Not used for distances, and only `'none'` with
#'        `chisquare = TRUE`. Default: `'none'`.
#'
#' @param chisquare  Partition the chi-square inertia of the count table by
#'        CCA, rather than its variance by RDA, as `vegan::varpart(chisquare
#'        = TRUE)` does. Adjusted \eqn{R^2} values are then estimated by
#'        permutation, as in [cca()]. Default: `FALSE`.
#'
#' @param permutations  With `chisquare = TRUE`: the number of permutations
#'        for each adjusted \eqn{R^2} (at least 1), or a matrix with one
#'        permutation of the samples per row (such as from
#'        `permute::shuffleSet()`), so that results can be compared exactly
#'        with vegan's. The same permutations are used for every set and
#'        union of sets. Not used otherwise. Default: `1000`.
#'
#' @param seed  With `chisquare = TRUE` and a number of `permutations`: the
#'        random seed for the permutations. Not used otherwise. Default: `0`.
#'
#' @return A data frame of class `c('ecodive_varpart', 'data.frame')`, with
#'   one row per fraction of the Venn diagram, named as vegan names them:
#'   `[a]` to `[d]` for two sets, `[a]` to `[h]` for three, and `[a]` to
#'   `[p]` for four, the last being the residuals. The first rows are the
#'   fractions unique to each set (`[a]` for the first set, `[b]` for the
#'   second, ...). Columns:
#'
#'   * `sets`: Which sets the fraction belongs to, such as `"env"` (unique
#'     to `env`), `"env & space"` (shared by those two only), or
#'     `"residuals"`.
#'   * `df`: For a unique fraction, the number of model columns that its
#'     set adds to all the others (its rank, less any it shares with them).
#'     `NA` for shared fractions and the residuals.
#'   * `r_squared`: The fraction of the inertia, unadjusted, which is biased
#'     upwards (vegan shows `NA`).
#'   * `adj_r_squared`: The adjusted fraction (Peres-Neto et al. 2006). A
#'     shared fraction may be negative.
#'   * `testable`: Whether a (partial) constrained ordination can test the
#'     fraction: only the unique fractions, each by the set conditioned on
#'     all the others, as in `rda(x, formula = ~ set + Condition(others))`.
#'
#'   Attributes:
#'
#'   * `unions`: A data frame of every set and union of sets (vegan's
#'     "Partition table"), with columns `fractions` (the Venn fractions it
#'     comprises, such as `"[a+d+f+g]"`), `df` (its rank), `r_squared`,
#'     `adj_r_squared`, and `testable` (all `TRUE`). Rows are named by the
#'     sets, joined by `+`.
#'   * `conditional`: For three or four sets, a data frame of each set
#'     conditioned on some, but not all, of the others (vegan's
#'     "Controlling 1 table" and "Controlling 2 tables"), with rows named
#'     like `"env | space"` and the same columns as `unions`. `NULL` for
#'     two sets.
#'   * `method`: `'rda'`, `'dbrda'`, or `'cca'`.
#'   * `norm`: The normalization, for RDA.
#'   * `total`: The total inertia: for RDA, the sum of squares of the
#'     centred (normalized) table (vegan's "Total variation (SS)"); for
#'     db-RDA, \eqn{\sum d_{ij}^2 / n} over pairs of samples; for CCA, the
#'     total chi-square inertia.
#'   * `n`: The number of samples.
#'   * `set_labels`: The sets as given, named by their names.
#'   * `n_perm`: For CCA, the number of permutations per adjusted
#'     \eqn{R^2}.
#'
#' @details
#'
#' **Method**
#'
#' For each set and each union of sets, the constrained ordination of `x`
#' on that set's model columns gives \eqn{R^2}, the fraction of the total
#' inertia it explains, and its adjusted \eqn{R^2}. For RDA and db-RDA that
#' is Ezekiel's adjustment, \eqn{1 - (1 - R^2)(n - 1)/(n - m - 1)} for
#' \eqn{n} samples and \eqn{m} model columns (the rank); for CCA it is
#' estimated by permutation, as [cca()] does (Peres-Neto et al. 2006). The
#' fractions of the Venn diagram are then sums and differences of those
#' values: with two sets X1 and X2, \eqn{[a] = R^2_{12} - R^2_2},
#' \eqn{[b] = R^2_{12} - R^2_1}, \eqn{[c] = R^2_1 + R^2_2 - R^2_{12}}, and
#' the residuals \eqn{[d] = 1 - R^2_{12}}; with more sets, by inclusion and
#' exclusion. Model columns that are linear combinations of others (within
#' or between sets) are dropped from each union, as `qr()` decides with a
#' tolerance of `1e-7`, so shared columns move into the shared fractions.
#'
#' Only the unique fractions can be tested, each by a partial constrained
#' ordination of the set given all the others, as with [permanova()] on
#' `x`'s distances and `formula = ~ set + Condition(others)`.
#'
#'
#' **Speed**
#'
#' For RDA and db-RDA, one pass over the table's non-zeros (or the
#' distances) gives the cross products of all the sets' columns, from
#' which every union's \eqn{R^2} follows without another pass, so the cost
#' is about that of a single [rda()] or [dbrda()] of all the sets together.
#' No \eqn{n \times n}{n x n} matrix is formed, and sparse tables stay
#' sparse. For CCA, each of the \eqn{2^k - 1} unions of \eqn{k} sets is a
#' [cca()] with its permutations. Results are identical for any number of
#' `cpus`.
#'
#'
#' **Differences from vegan::varpart()**
#'
#' * The result is one data frame of the Venn fractions, with vegan's
#'   other tables as attributes; there is no plot.
#' * `r_squared` is given for the Venn fractions too.
#' * `df` is `NA`, not `0`, for shared fractions and the residuals.
#' * Normalizations are ecodive's (`norm`), not `decostand()`'s
#'   `transfo`; there is no `scale`, and no `add` for distances.
#' * Samples are matched by name between `x`, `data`, and the sets.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::varpart()` has the same name. When both packages are attached,
#' the one attached last masks the other; call `ecodive::varpart()` to be
#' sure.
#'
#' @references
#' Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the
#' spatial component of ecological variation. *Ecology*, 73, 1045-1055.
#' \doi{10.2307/1940179}
#'
#' Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006).
#' Variation partitioning of species data matrices: estimation and
#' comparison of fractions. *Ecology*, 87, 2614-2625.
#' \doi{10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2}
#'
#' @seealso [rda()], [dbrda()], [cca()], [permanova()]
#'
#' @export
#' @examples
#'     # 40 samples: two environmental variables, and coordinates.
#'     set.seed(1)
#'     meta <- data.frame(
#'       ph   = runif(40, 5, 8),
#'       temp = runif(40, 10, 20),
#'       x    = runif(40),
#'       y    = runif(40) )
#'     counts <- matrix(rpois(40 * 15, 4), 40, 15)
#'     counts[, 1:5] <- counts[, 1:5] + round(3 * meta$ph + 4 * meta$x)
#'
#'     # tb-RDA on the Hellinger-transformed table.
#'     vp <- varpart(
#'       counts, env = ~ ph + temp, space = ~ x + y,
#'       data = meta, norm = 'hellinger' )
#'     vp
#'     plot(vp)
#'
#'     # db-RDA on Bray-Curtis distances; a set may be a matrix.
#'     varpart(bray(counts), ~ ph + temp, as.matrix(meta[, c('x', 'y')]), data = meta)
#'
#'     # Test the fraction unique to `env`.
#'     permanova(
#'       euclidean(counts, norm = 'hellinger'),
#'       data = meta, formula = ~ ph + temp + Condition(x + y) )
#'
varpart <- function (
    x, ..., data = NULL, norm = 'none', pseudocount = NULL, margin = 1L,
    chisquare = FALSE, permutations = 1000L, seed = 0L, cpus = n_cpus() ) {

  env <- parent.frame()

  validate_cpus()

  if (!isTRUE(chisquare) && !isFALSE(chisquare))
    stop('`chisquare` must be TRUE or FALSE.', call. = FALSE)


  #--------------------------------------------------------
  # The response.
  #--------------------------------------------------------
  y <- varpart_response(x, norm, pseudocount, margin, chisquare)
  n <- y$n

  if (n < 3)
    stop('`x` must have at least 3 samples.', call. = FALSE)


  #--------------------------------------------------------
  # The sets, with their labels.
  #--------------------------------------------------------
  sets  <- list(...)
  exprs <- as.list(substitute(list(...)))[-1]
  k     <- length(sets)

  if (k < 2 || k > 4)
    stop('Give two to four sets of explanatory variables in `...`; got ', k, '.', call. = FALSE)

  nms <- names(sets)
  if (is.null(nms)) nms <- rep('', k)
  nms[is.na(nms) | nms == ''] <- paste0('X', seq_len(k))[is.na(nms) | nms == '']
  if (anyDuplicated(nms))
    stop('The sets must have different names: ', name_list(nms[duplicated(nms)]), call. = FALSE)

  set_labels <- vapply(exprs, function (e) paste(deparse(e, width.cutoff = 500L), collapse = ' '), '')
  names(set_labels) <- nms


  #--------------------------------------------------------
  # Each set's model columns, aligned with the samples.
  #--------------------------------------------------------
  if (!is.null(data) && !is.data.frame(data) && !is.matrix(data))
    stop('`data` must be a data.frame.', call. = FALSE)

  cols <- vector('list', k)
  for (i in seq_len(k))
    cols[[i]] <- varpart_set(sets[[i]], nms[[i]], data, n, y$labels, env)

  X     <- do.call(cbind, cols)
  owner <- rep(seq_len(k), vapply(cols, ncol, 0L))
  storage.mode(X) <- 'double'


  #--------------------------------------------------------
  # R-squared, adjusted R-squared, and rank of every union
  # of sets, indexed by bit mask (1 .. 2^k - 1).
  #--------------------------------------------------------
  masks <- seq_len(2^k - 1)
  in_mask <- function (mask) bitwAnd(mask, 2^(seq_len(k) - 1)) > 0
  mask_cols <- function (mask) which(owner %in% which(in_mask(mask)))

  if (y$kind == 'cca') {
    fits <- varpart_cca(y, X, masks, mask_cols, permutations, seed, cpus)
  }
  else {
    model <- response_cross(y, X, cpus)
    if (!(model$total > 0))
      stop('`x` has no variation to partition.', call. = FALSE)
    fits <- lapply(masks, function (mask) {
      fit <- subset_inertia(model, mask_cols(mask))
      r2  <- fit$inertia / model$total
      list(r2 = r2, adj = adj_r_squared(r2, n, fit$rank), rank = fit$rank)
    })
  }

  r2   <- vapply(fits, `[[`, 0, 'r2')
  adj  <- vapply(fits, `[[`, 0, 'adj')
  rank <- vapply(fits, `[[`, 0L, 'rank')

  single <- 2^(seq_len(k) - 1)
  if (any(empty <- rank[single] == 0))
    stop('Set ', name_list(nms[empty]), ' has no variation among the samples.', call. = FALSE)

  if (sum(rank[single]) > rank[[2^k - 1]])
    warning(
      'Some sets are collinear (their columns are linear combinations of ',
      'each other\'s): ', sum(rank[single]) - rank[[2^k - 1]], ' of ',
      sum(rank[single]), ' model columns are shared, and are counted in ',
      'the shared fractions.', call. = FALSE )

  out <- varpart_fractions(k, nms, r2, adj, rank)

  attr(out, 'method')     <- y$kind
  attr(out, 'norm')       <- if (y$kind == 'rda') y$norm_name
  attr(out, 'total')      <- if (y$kind == 'cca') fits[[1]]$total else model$total
  attr(out, 'n')          <- n
  attr(out, 'set_labels') <- set_labels
  if (y$kind == 'cca') attr(out, 'n_perm') <- fits[[1]]$n_perm
  class(out) <- c('ecodive_varpart', 'data.frame')
  out
}



#------------------------------------------------------------------------------
# The community as varpart() and forward_select() read it: distances (kind
# 'dbrda'), a count table ('rda'), or, with `chisquare`, a count table for
# CCA ('cca'). Returns a list with `kind`, `n`, `labels`, and what
# response_cross() or C_cca() needs: `dm` (packed lower triangle), or
# `counts`, `margin`, `norm` (code), `norm_name`, and `pseudocount`.
#------------------------------------------------------------------------------
varpart_response <- function (x, norm, pseudocount, margin, chisquare) {

  if (inherits(x, 'formula'))
    stop('`x` must be a dist object or a count table, not a formula. ',
         'Give the sets of explanatory variables after it.', call. = FALSE)

  if (inherits(x, 'dist')) {

    if (chisquare)
      stop('`chisquare = TRUE` needs a count table, not distances.', call. = FALSE)
    if (!identical(norm, 'none'))
      stop('`norm` applies to a count table, not distances.', call. = FALSE)

    dv <- dist_vector(x)
    return (list(kind = 'dbrda', n = dv$n, labels = dv$labels, dm = dv$dm))
  }

  counts <- x
  validate_counts()
  norm_name <- if (is.character(norm) && length(norm) == 1) tolower(norm) else norm
  validate_norm()
  validate_pseudocount()
  validate_margin()

  if (chisquare && norm != NORM_NONE)
    stop('With `chisquare = TRUE`, `norm` must be \'none\': CCA analyzes ',
         'the counts themselves.', call. = FALSE)

  dims <- cca_dims(counts, margin)

  list(
    kind        = if (chisquare) 'cca' else 'rda',
    n           = dims$n,
    labels      = dims$dn[[1]],
    counts      = counts,
    margin      = margin,
    norm        = norm,
    norm_name   = match.arg(norm_name, c('none', 'percent', 'chord', 'hellinger', 'binary', 'clr', 'rclr')),
    pseudocount = pseudocount )
}



#------------------------------------------------------------------------------
# One set of explanatory variables as an n x m numeric matrix, aligned with
# the samples: a one-sided formula over `data`, a data frame (dummy-coded),
# a numeric matrix, or a vector or factor.
#------------------------------------------------------------------------------
varpart_set <- function (s, name, data, n, labels, env) {

  what <- paste0('Set ', name_list(name))

  if (inherits(s, 'formula')) {

    if (length(s) == 3)
      stop(what, ' must be a one-sided formula, such as `~ ph + temp`.', call. = FALSE)
    if (is.null(data))
      stop(what, ' is a formula, so `data` is required.', call. = FALSE)

    d <- model_design(data, s, labels = labels, n = n, env = env)
    if (!is.null(d$z))
      stop(what, ' has Condition() terms. varpart() conditions each set on ',
           'the others itself; give the variables as a set of their own.', call. = FALSE)
    return (d$x)
  }

  if (is.atomic(s) && is.null(dim(s))) {
    if (!is.numeric(s) && !is.factor(s) && !is.character(s) && !is.logical(s))
      stop(what, ' must be a formula, data frame, matrix, or vector.', call. = FALSE)
    nms <- names(s)
    if (!is.null(nms) && (anyNA(nms) || anyDuplicated(nms))) nms <- NULL
    s <- data.frame(s, stringsAsFactors = FALSE)
    names(s) <- name
    if (!is.null(nms)) rownames(s) <- nms
  }

  if (is.matrix(s)) {

    if (!is.numeric(s) && !is.logical(s))
      stop(what, ' is a matrix, so it must be numeric.', call. = FALSE)

    idx <- match_rows(labels, n, rownames(s), nrow(s), ref = '`x`', what = what)
    m   <- s[idx, , drop = FALSE]
    storage.mode(m) <- 'double'
    if (is.null(colnames(m))) colnames(m) <- paste0(name, '.', seq_len(ncol(m)))
    rownames(m) <- NULL

    if (anyNA(m))
      stop(what, ' has missing values (NA).', call. = FALSE)
    if (!all(is.finite(m)))
      stop(what, ' has infinite values.', call. = FALSE)
    return (m)
  }

  if (is.data.frame(s)) {

    if (ncol(s) == 0)
      stop(what, ' has no columns.', call. = FALSE)

    rn  <- if (.row_names_info(s) > 0) rownames(s) else NULL
    idx <- match_rows(labels, n, rn, nrow(s), ref = '`x`', what = what)
    s   <- droplevels(s[idx, , drop = FALSE])
    rownames(s) <- NULL

    na <- names(s)[vapply(s, anyNA, logical(1))]
    if (length(na))
      stop(what, ' has missing values (NA) in ', name_list(na), '.', call. = FALSE)
    inf <- names(s)[vapply(s, function (v) is.numeric(v) && !all(is.finite(v)), logical(1))]
    if (length(inf))
      stop(what, ' has infinite values in ', name_list(inf), '.', call. = FALSE)

    for (j in seq_along(s)) if (is.character(s[[j]])) s[[j]] <- factor(s[[j]])

    tt <- stats::terms(~ ., data = s)
    return (centred_model_matrix(tt, s)$x)
  }

  stop(what, ' must be a formula, data frame, matrix, or vector.', call. = FALSE)
}



#------------------------------------------------------------------------------
# Shared with forward_select(): the cross products of the design's centred
# columns with the community, Xc' B Xc, and the total inertia trace(B), from
# one pass over the data (C_cord_cross()). B is Gower's matrix of the
# distances (y$kind 'dbrda') or A A' for the centred, normalized table A
# ('rda'), both unscaled. Returns list(cross, total, xc), xc being the
# centred design.
#------------------------------------------------------------------------------
response_cross <- function (y, X, cpus) {

  xc <- X - rep(colMeans(X), each = nrow(X))
  attributes(xc) <- list(dim = dim(X), dimnames = list(NULL, colnames(X)))

  res <- if (y$kind == 'dbrda')
    .Call(C_cord_cross, 0L, y$dm, y$n, NULL, 0L, 0, xc, cpus)
  else
    .Call(C_cord_cross, 1L, y$counts, y$n, y$margin, y$norm, y$pseudocount, xc, cpus)

  list(cross = res$cross, total = res$total, xc = xc)
}



#------------------------------------------------------------------------------
# Shared with forward_select(): the inertia that columns `cols` of the design
# explain, trace(Q' B Q) for an orthonormal basis Q of those (centred)
# columns, from response_cross()'s `model`. The rank is qr()'s, with
# tolerance 1e-7, as vegan decides it; with Xc_S P = Q R, Q = Xc_S P R^-1,
# so trace(Q' B Q) = trace(R^-T (P' Xc_S' B Xc_S P) R^-1). Returns
# list(inertia, rank, kept), `kept` being the columns of the basis.
#------------------------------------------------------------------------------
subset_inertia <- function (model, cols, tol = 1e-7) {

  if (length(cols) == 0)
    return (list(inertia = 0, rank = 0L, kept = integer(0)))

  q <- qr(model$xc[, cols, drop = FALSE], tol = tol)
  r <- q$rank
  if (r == 0)
    return (list(inertia = 0, rank = 0L, kept = integer(0)))

  kept <- cols[q$pivot[seq_len(r)]]
  ri   <- backsolve(qr.R(q)[seq_len(r), seq_len(r), drop = FALSE], diag(r))
  list(
    inertia = sum(ri * (model$cross[kept, kept, drop = FALSE] %*% ri)),
    rank    = r,
    kept    = kept )
}



# Shared with forward_select(): Ezekiel's adjusted R-squared for n samples
# and m model columns, as vegan::RsquareAdj() gives it for rda() and dbrda();
# NA without residual degrees of freedom.
adj_r_squared <- function (r2, n, m) {
  if (n - m - 1 <= 0) return (NA_real_)
  1 - (1 - r2) * (n - 1) / (n - m - 1)
}



#------------------------------------------------------------------------------
# varpart(chisquare = TRUE): each union's R-squared and permutation adjusted
# R-squared from C_cca(), with the same permutations for every union, as
# cca() computes them for an unconditioned model.
#------------------------------------------------------------------------------
varpart_cca <- function (y, X, masks, mask_cols, permutations, seed, cpus) {

  validate_seed()
  n  <- y$n
  pm <- perm_matrix(permutations, n)
  if (pm$n_perm < 1)
    stop('With `chisquare = TRUE`, `permutations` must be at least 1: ',
         'the adjusted R-squared is estimated by permutation.', call. = FALSE)

  lapply(masks, function (mask) {

    d   <- cord_design(X[, mask_cols(mask), drop = FALSE], NULL, n)
    res <- .Call(
      C_cca, y$counts, y$margin, d$design, 0L, 1L,
      pm$perms, pm$n_perm, seed, TRUE, cpus )

    if (!is.null(res$empty))
      stop(cca_empty_msg(res$empty, y$labels), call. = FALSE)

    r2 <- cca_r_squared(res, pm$n_perm, d, n)
    list(
      r2 = r2[[1]], adj = r2[[2]], rank = as.integer(sum(res$kept)),
      total = res$total, n_perm = pm$n_perm )
  })
}



# The sets of each lettered fraction of k sets, in vegan's order (the
# residuals, last, are not included).
varpart_regions <- function (k) {
  switch(
    k - 1L,
    list(1, 2, 1:2),
    list(1, 2, 3, 1:2, 2:3, c(1, 3), 1:3),
    list(1, 2, 3, 4, 1:2, 2:3, c(1, 3), c(1, 4), c(2, 4), 3:4,
         c(1, 2, 4), 1:3, 2:4, c(1, 3, 4), 1:4) )
}


#------------------------------------------------------------------------------
# The Venn fractions, from every union's R-squared, adjusted R-squared, and
# rank (indexed by bit mask). Fraction W (the region shared by exactly the
# sets in W) is sum over S in W of (-1)^(|W| - |S| + 1) V(all \ S), with
# V(empty) = 0. Letters follow vegan's.
#------------------------------------------------------------------------------
varpart_fractions <- function (k, nms, r2, adj, rank) {

  regions <- varpart_regions(k)

  all_mask <- 2^k - 1
  mask_of  <- function (s) sum(2^(s - 1))
  value    <- function (v, mask) if (mask == 0) 0 else v[[mask]]
  subsets  <- function (w) {
    # Indices, as combn(2, 1) would be combn(1:2, 1).
    unlist(lapply(0:length(w), function (j) {
      lapply(utils::combn(length(w), j, simplify = FALSE), function (i) w[i])
    }), recursive = FALSE)
  }
  fraction <- function (v, w) {
    sum(vapply(subsets(w), function (s) {
      (-1)^(length(w) - length(s) + 1) * value(v, all_mask - mask_of(s))
    }, 0))
  }

  letters_ <- letters[seq_len(length(regions) + 1L)]
  # The fractions (letters) of the regions for which `hit` is TRUE.
  join <- function (hit) paste0('[', paste(letters_[which(hit)], collapse = '+'), ']')
  covers <- function (sets) join(vapply(regions, function (r) any(r %in% sets), TRUE))

  unique_ <- vapply(regions, length, 0L) == 1L
  df      <- vapply(seq_along(regions), function (i) {
    if (!unique_[[i]]) return (NA_integer_)
    as.integer(rank[[all_mask]] - value(rank, all_mask - mask_of(regions[[i]])))
  }, 0L)

  out <- data.frame(
    sets          = c(vapply(regions, function (r) paste(nms[r], collapse = ' & '), ''), 'residuals'),
    df            = c(df, NA_integer_),
    r_squared     = c(vapply(regions, function (w) fraction(r2,  w), 0), 1 - r2[[all_mask]]),
    adj_r_squared = c(vapply(regions, function (w) fraction(adj, w), 0), 1 - adj[[all_mask]]),
    testable      = c(unique_, FALSE),
    row.names     = paste0('[', letters_, ']'),
    stringsAsFactors = FALSE )


  # Every union (vegan's "Partition table"), by size, then in order.
  unions <- unlist(lapply(seq_len(k), function (j) utils::combn(k, j, simplify = FALSE)), recursive = FALSE)
  attr(out, 'unions') <- data.frame(
    fractions     = vapply(unions, covers, ''),
    df            = vapply(unions, function (s) as.integer(rank[[mask_of(s)]]), 0L),
    r_squared     = vapply(unions, function (s) r2[[mask_of(s)]], 0),
    adj_r_squared = vapply(unions, function (s) adj[[mask_of(s)]], 0),
    testable      = TRUE,
    row.names     = vapply(unions, function (s) paste(nms[s], collapse = '+'), ''),
    stringsAsFactors = FALSE )


  # Each set given some, but not all, of the others (vegan's "Controlling
  # 1 table" and "Controlling 2 tables").
  cond <- NULL
  if (k > 2) {
    pairs <- list()
    for (j in seq_len(k - 2L))
      for (i in seq_len(k))
        for (z in utils::combn(k - 1L, j, simplify = FALSE))
          pairs[[length(pairs) + 1L]] <- list(i = i, z = setdiff(seq_len(k), i)[z])

    cond <- data.frame(
      fractions     = vapply(pairs, function (p) {
        join(vapply(regions, function (r) p$i %in% r && !any(p$z %in% r), TRUE))
      }, ''),
      df            = vapply(pairs, function (p) as.integer(rank[[mask_of(c(p$i, p$z))]] - rank[[mask_of(p$z)]]), 0L),
      r_squared     = vapply(pairs, function (p) r2[[mask_of(c(p$i, p$z))]] - r2[[mask_of(p$z)]], 0),
      adj_r_squared = vapply(pairs, function (p) adj[[mask_of(c(p$i, p$z))]] - adj[[mask_of(p$z)]], 0),
      testable      = TRUE,
      row.names     = vapply(pairs, function (p) paste(nms[p$i], '|', paste(nms[p$z], collapse = '+')), ''),
      stringsAsFactors = FALSE )
  }
  attr(out, 'conditional') <- cond

  out
}



#' @export
print.ecodive_varpart <- function (x, digits = 4L, ...) {

  method <- attr(x, 'method')
  n      <- attr(x, 'n')
  labs   <- attr(x, 'set_labels')
  nms    <- names(labs)
  k      <- length(labs)
  unions <- attr(x, 'unions')
  total  <- attr(x, 'total')

  how <- switch(
    method,
    rda   = sprintf('RDA%s', if (!identical(attr(x, 'norm'), 'none')) sprintf(' (%s)', attr(x, 'norm')) else ''),
    dbrda = 'db-RDA',
    cca   = sprintf('CCA (adjusted R-squared from %d permutations)', attr(x, 'n_perm')) )

  cat(sprintf('Variation partitioning by %s: %d samples, %d sets\n', how, n, k))
  cat(sprintf(
    'Total inertia: %s%s\n\n', format(signif(total, 6)),
    if (method == 'rda') sprintf(' (sum of squares; variance %s)', format(signif(total / (n - 1), 6))) else '' ))

  fmt  <- function (v) if (is.na(v)) 'NA' else formatC(v, format = 'f', digits = digits)
  w    <- max(nchar(c(nms, 'All')))
  wl   <- max(nchar(labs))
  all_ <- unions[nrow(unions), ]
  for (i in seq_len(k + 1L)) {
    s <- if (i > k) all_ else unions[nms[[i]], ]
    cat(sprintf('  %-*s  %-*s  %3d df, adj. R-squared %s\n',
                w, if (i > k) 'All' else nms[[i]], wl, if (i > k) '' else labs[[i]],
                s$df, fmt(s$adj_r_squared)))
  }
  cat('\n')

  venn <- varpart_venn(x, nms)
  if (length(venn)) cat(venn, sep = '\n')

  df <- as.data.frame(unclass(x)[c('sets', 'df', 'r_squared', 'adj_r_squared', 'testable')],
                      row.names = rownames(x), stringsAsFactors = FALSE)
  df$r_squared     <- round(df$r_squared, digits)
  df$adj_r_squared <- round(df$adj_r_squared, digits)
  print(df, na.print = '')

  invisible(x)
}



# A text Venn diagram of the adjusted fractions, for two or three sets
# (character lines; none for four).
varpart_venn <- function (x, nms) {

  k   <- length(nms)
  if (k > 3) return (character(0))

  val <- function (i) {
    v <- x$adj_r_squared[[i]]
    if (is.na(v)) '     NA' else sprintf('%7.4f', v)
  }
  short <- function (s) if (nchar(s) > 10) paste0(substr(s, 1, 9), '~') else s

  if (k == 2) {
    h <- 8; w <- 46
    boxes  <- list(c(0, 5, 0, 28), c(1, 6, 16, 44))
    labels <- list(
      list(2, 3, '[a]'), list(3, 3, val(1)),
      list(2, 18, '[c]'), list(3, 18, val(3)),
      list(3, 31, '[b]'), list(4, 31, val(2)) )
    names_ <- list(list(0, 2, short(nms[[1]])), list(6, 18, short(nms[[2]])))
    resid  <- 4
  }
  else {
    h <- 14; w <- 48
    boxes  <- list(c(0, 9, 0, 30), c(2, 8, 18, 46), c(5, 12, 6, 38))
    labels <- list(
      list(2, 3, '[a]'),   list(3, 3, val(1)),
      list(3, 20, '[d]'),  list(4, 20, val(4)),
      list(3, 33, '[b]'),  list(4, 33, val(2)),
      list(6, 20, '[g]'),  list(7, 20, val(7)),
      list(6, 31, '[e]'),  list(7, 31, val(5)),
      list(6, 9, '[f]'),   list(7, 9, val(6)),
      list(10, 12, paste('[c]', val(3))) )
    names_ <- list(
      list(0, 2, short(nms[[1]])), list(2, 34, short(nms[[2]])),
      list(12, 9, short(nms[[3]])) )
    resid  <- 8
  }

  g <- matrix(' ', h, w)
  for (b in boxes) {
    r <- b[1:2] + 1; cc <- b[3:4] + 1
    for (j in cc[1]:cc[2]) for (i in r) g[i, j] <- if (g[i, j] == '|') '+' else '-'
    for (i in r[1]:r[2]) for (j in cc) g[i, j] <- if (g[i, j] == '-') '+' else '|'
    for (i in r) for (j in cc) g[i, j] <- '+'
  }
  put <- function (row, col, s) {
    ch <- strsplit(s, '')[[1]]
    g[row + 1, col + seq_along(ch)] <<- ch
  }
  for (l in names_) put(l[[1]], l[[2]], paste0(' ', l[[3]], ' '))
  for (l in labels) put(l[[1]], l[[2]], l[[3]])

  lines <- sub(' +$', '', apply(g, 1L, paste, collapse = ''))
  while (length(lines) && lines[[length(lines)]] == '') lines <- lines[-length(lines)]
  c(paste0('  ', lines),
    sprintf('  Residuals [%s] %s', letters[nrow(x)], trimws(val(nrow(x)))), '')
}
