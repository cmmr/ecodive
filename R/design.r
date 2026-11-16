# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Helpers shared by the methods that test or constrain a distance matrix with
# sample metadata: permanova() (and dbrda()).
#
#   model_design()   metadata + formula -> centred design matrices, with the
#                    rows in the same order as the dist's samples.
#   design_basis()   orthonormal basis (or, for partial models, the kept
#                    columns) of a design, with R's qr() rank decisions.
#   condition_basis() orthonormal basis of the Condition() columns.
#   perm_matrix()    checks a `permutations` argument.
#   perm_ss()        per-column sums of squares, trace(Q' G Q), for the
#                    observed design and every permutation (C_permanova).



#------------------------------------------------------------------------------
# Aligns `data` with the samples of `dm` and builds the centred model matrices
# of `formula`'s terms (x) and of its Condition() terms (z).
#
#   data     A data.frame (or matrix) of sample metadata, or a single vector
#            or factor (taken as one variable named `group`).
#   formula  A one-sided formula, or NULL for every column of `data`. A
#            left-hand side is ignored.
#   labels   The dist's sample names, or NULL.
#   n        The number of samples.
#   strata   NULL, the name of a column of `data`, or a vector with one
#            value per row of `data`, in the same order.
#
# Rows of `data` are matched to `labels` by row name; rows for other samples
# are left out (so their NAs don't matter). Unnamed rows (or an unlabelled
# dist) are matched by position instead, and must be one per sample.
#
# Returns a list:
#   x            n x p centred model matrix, without the intercept.
#   assign       Each column of x's term, an index into term_labels.
#   term_labels  The terms of x.
#   terms        The terms object of x (for drop.scope()).
#   z            n x pz centred matrix of the Condition() variables, or NULL.
#   strata       Integer codes 0 .. n_strata - 1, or NULL.
#   data         `data`, with rows aligned to the samples.
#------------------------------------------------------------------------------
model_design <- function (data, formula = NULL, labels = NULL, n, strata = NULL) {

  if (missing(data) || is.null(data))
    stop('`data` is required: a data.frame of sample metadata.', call. = FALSE)


  # A vector or factor is one variable, `group`.
  if (is.null(dim(data)) && is.atomic(data)) {
    nms  <- names(data)
    data <- data.frame(group = data, stringsAsFactors = FALSE)
    if (!is.null(nms) && !anyDuplicated(nms)) rownames(data) <- nms
  }

  if (!is.data.frame(data)) {
    if (!is.matrix(data))
      stop('`data` must be a data.frame, matrix, vector, or factor.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }


  # `strata` as a column name, or a vector parallel to the rows of `data`.
  if (!is.null(strata)) {
    if (is.character(strata) && length(strata) == 1 && strata %in% names(data)) {
      strata <- data[[strata]]
    }
    else if (length(strata) != nrow(data) || !is.null(dim(strata))) {
      stop(
        '`strata` must be the name of a column of `data`, or have one value ',
        'per row of `data` (', nrow(data), ').', call. = FALSE )
    }
  }


  #--------------------------------------------------------
  # Align the rows of `data` with the dist's samples.
  #--------------------------------------------------------
  # Row names that are only numbers, as from subsetting a data frame without
  # row names, count as names only when they match the labels.
  named <- !is.null(labels) && .row_names_info(data) > 0
  if (named && all(grepl("^[0-9]+$", rownames(data))))
    named <- all(labels %in% rownames(data))

  if (named) {

    if (anyDuplicated(labels))
      stop('`dm` has duplicate sample names, so `data` cannot be matched to them.', call. = FALSE)
    if (anyDuplicated(rownames(data)))
      stop('`data` has duplicate row names.', call. = FALSE)

    idx <- match(labels, rownames(data))

    if (anyNA(idx)) {
      miss <- labels[is.na(idx)]
      stop(
        sprintf('%d sample%s in `dm` %s no row in `data`: %s',
          length(miss), if (length(miss) == 1) '' else 's',
          if (length(miss) == 1) 'has' else 'have', name_list(miss) ),
        call. = FALSE )
    }

    # Rows for samples not in `dm` (e.g. dropped by rarefying) are left out.
  }
  else {
    if (nrow(data) != n)
      stop(
        sprintf('`data` has %d rows, but `dm` has %d samples.', nrow(data), n),
        if (!is.null(labels)) '\nGive `data` row names to match its rows to the samples by name.',
        call. = FALSE )
    idx <- seq_len(n)
  }

  data <- data[idx, , drop = FALSE]
  if (!is.null(labels)) rownames(data) <- labels

  if (!is.null(strata)) {
    strata <- strata[idx]
    if (anyNA(strata)) stop('`strata` must not contain NAs.', call. = FALSE)
    strata <- as.integer(factor(strata)) - 1L
  }


  #--------------------------------------------------------
  # Parse the formula; split off Condition() terms.
  #--------------------------------------------------------
  if (is.null(formula)) formula <- ~ .
  if (is.character(formula)) formula <- stats::as.formula(formula, env = parent.frame(2))
  if (!inherits(formula, 'formula'))
    stop('`formula` must be a formula, such as `~ treatment + ph`.', call. = FALSE)

  env <- environment(formula)
  if (length(formula) == 3) formula <- formula[-2]
  environment(formula) <- if (is.null(env)) globalenv() else env

  tt <- tryCatch(
    stats::terms(formula, data = data, specials = 'Condition'),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  if (attr(tt, 'response') > 0)
    tt <- stats::delete.response(tt)

  # Every variable must be a column of `data` (or a constant), so that it
  # lines up with the samples.
  vars    <- all.vars(stats::delete.response(tt))
  missing <- vars[!vars %in% names(data)]
  missing <- missing[!vapply(missing, function (v) {
    obj <- get0(v, envir = environment(formula), inherits = TRUE)
    !is.null(obj) && !is.function(obj) && length(obj) == 1 }, logical(1))]
  if (length(missing))
    stop('Not found in `data`: ', name_list(missing), call. = FALSE)

  used <- intersect(vars, names(data))
  na   <- used[vapply(data[used], anyNA, logical(1))]
  if (length(na))
    stop('Missing values (NA) in `data` column', if (length(na) > 1) 's', ': ',
         name_list(na), call. = FALSE)

  term_labels <- attr(tt, 'term.labels')
  cond_vars   <- attr(tt, 'specials')$Condition  # rows of attr(tt, 'factors')
  cond_terms  <- integer(0)
  z           <- NULL

  if (length(cond_vars)) {

    fac        <- attr(tt, 'factors')
    cond_terms <- which(colSums(fac[cond_vars, , drop = FALSE] != 0) > 0)

    if (any(colSums(fac[, cond_terms, drop = FALSE] != 0) > 1))
      stop('Condition() terms cannot be part of an interaction.', call. = FALSE)

    # Condition(a + b, c) -> ~ a + b + c
    vlist <- as.list(attr(tt, 'variables'))[-1]
    args  <- unlist(lapply(vlist[cond_vars], function (v) as.list(v)[-1]), recursive = FALSE)
    rhs   <- Reduce(function (a, b) call('+', a, b), args)
    zf    <- stats::as.formula(call('~', rhs), env = environment(formula))
    z     <- centred_model_matrix(stats::terms(zf, data = data), data)$x

    if (ncol(z) == 0) z <- NULL
  }

  if (length(cond_terms) == length(term_labels))
    stop('`formula` has no terms to test.', call. = FALSE)

  xt <- if (length(cond_terms)) stats::drop.terms(tt, cond_terms) else tt
  attr(xt, 'specials') <- NULL

  mm <- centred_model_matrix(xt, data)

  list(
    x           = mm$x,
    assign      = mm$assign,
    term_labels = attr(xt, 'term.labels'),
    terms       = xt,
    z           = z,
    strata      = strata,
    data        = data )
}


# Model matrix without its intercept, with each column centred.
centred_model_matrix <- function (tt, data) {

  mf <- stats::model.frame(tt, data, na.action = stats::na.fail)
  mm <- stats::model.matrix(tt, mf)

  assign <- attr(mm, 'assign')
  keep   <- assign > 0
  x      <- mm[, keep, drop = FALSE]
  assign <- assign[keep]

  x <- x - rep(colMeans(x), each = nrow(x))
  attributes(x) <- list(dim = dim(x), dimnames = list(rownames(data), colnames(x)))

  list(x = x, assign = assign)
}


# 'a', 'b', 'c', ... for error messages.
name_list <- function (x, max = 5) {
  x <- paste0("'", x, "'")
  if (length(x) > max) x <- c(x[seq_len(max)], '...')
  paste(x, collapse = ', ')
}



#------------------------------------------------------------------------------
# Orthonormal, centred basis of the Condition() columns z, using qr()'s rank
# decision. Returns an n x rank matrix, or NULL.
#------------------------------------------------------------------------------
condition_basis <- function (z, tol = 1e-7) {
  if (is.null(z) || ncol(z) == 0) return(NULL)
  qz   <- qr(z, tol = tol)
  rank <- qz$rank
  if (rank == 0) return(NULL)
  q <- qr.Q(qz)[, seq_len(rank), drop = FALSE]
  attr(q, 'z') <- z[, qz$pivot[seq_len(rank)], drop = FALSE]
  q
}



#------------------------------------------------------------------------------
# The columns of design x that span its model, with qr()'s rank decision and
# column order (as vegan uses).
#
# Without conditions (qz = NULL), `cols` is an orthonormal basis, in the order
# of x's kept columns. Then each column's quadratic form is that column's
# sequential sum of squares.
#
# With conditions, `cols` is x's kept columns themselves: C_permanova()
# permutes them, makes them orthogonal to qz, and orthonormalizes them.
#
# Returns list(cols, keep (indices into x's columns), rank).
#------------------------------------------------------------------------------
design_basis <- function (x, qz = NULL, tol = 1e-7) {

  n <- nrow(x)

  if (ncol(x) == 0)
    return(list(cols = matrix(0, n, 0), keep = integer(0), rank = 0L))

  if (is.null(qz)) {
    qx   <- qr(x, tol = tol)
    rank <- qx$rank
    keep <- qx$pivot[seq_len(rank)]
    cols <- qr.Q(qx)[, seq_len(rank), drop = FALSE]
  }
  else {
    zk   <- attr(qz, 'z')
    pz   <- ncol(zk)
    qx   <- qr(cbind(zk, x), tol = tol)
    rank <- qx$rank - pz
    keep <- qx$pivot[pz + seq_len(rank)] - pz
    cols <- x[, keep, drop = FALSE]
  }

  colnames(cols) <- colnames(x)[keep]
  list(cols = cols, keep = keep, rank = rank)
}



#------------------------------------------------------------------------------
# Checks `permutations`: a count, or a matrix with one permutation of
# 1 .. n per row. Returns list(n_perm, perms), where perms is NULL (random
# permutations) or an n x n_perm integer matrix of 0-based indices.
#------------------------------------------------------------------------------
perm_matrix <- function (permutations, n) {

  if (is.matrix(permutations)) {

    p <- permutations
    if (!is.numeric(p) || anyNA(p) || ncol(p) != n || any(p != round(p)))
      stop('A `permutations` matrix must have ', n, ' columns (one per sample) ',
           'of whole numbers, with one permutation of 1 to ', n, ' per row.', call. = FALSE)

    storage.mode(p) <- 'integer'
    ok <- apply(p, 1L, function (r) all(sort.int(r, method = 'radix') == seq_len(n)))
    if (!all(ok))
      stop('Row ', which(!ok)[[1]], ' of `permutations` is not a permutation of 1 to ', n, '.', call. = FALSE)

    return(list(n_perm = nrow(p), perms = t(p) - 1L))
  }

  if (!is.numeric(permutations) || length(permutations) != 1 || is.na(permutations) ||
      permutations < 0 || permutations %% 1 != 0 || permutations > .Machine$integer.max - 1)
    stop('`permutations` must be a whole number of at least 0, or a matrix ',
         'with one permutation per row.', call. = FALSE)

  list(n_perm = as.integer(permutations), perms = NULL)
}



#------------------------------------------------------------------------------
# Sums of squares, trace(Q' G Q) column by column, for the identity (row 1)
# and each permutation (rows 2 .. n_perm + 1).
#
#   dm       Packed lower triangle (double), from dist_vector().
#   n        Samples.
#   cols     n x m design columns (see design_basis()), or NULL with `grp`.
#   sets     Each column's set: columns are orthonormalized within a set.
#   qz       Orthonormal basis of the conditions, or NULL.
#   grp      One factor instead of cols: integer group codes 0 .. g - 1.
#   pm       From perm_matrix().
#   strata   Integer codes 0 .. n_strata - 1, or NULL.
#
# Returns list(ss = (n_perm + 1) x m matrix (x 1 for grp: the model's sum
# of squares), total = trace(G)).
#------------------------------------------------------------------------------
perm_ss <- function (
    dm, n, cols = NULL, sets = NULL, qz = NULL, grp = NULL,
    pm, strata = NULL, seed = 0L, cpus = 1L ) {

  if (!is.null(cols)) {
    storage.mode(cols) <- 'double'
    if (is.null(sets)) sets <- rep(1L, ncol(cols))
    sets <- as.integer(sets)
  }
  if (!is.null(qz)) {
    qz <- unname(qz[, , drop = FALSE])
    attributes(qz) <- list(dim = dim(qz))
    storage.mode(qz) <- 'double'
  }
  if (!is.null(grp)) grp <- as.integer(grp)
  if (!is.null(strata)) strata <- as.integer(strata)

  res <- .Call(
    C_permanova, dm, as.integer(n), cols, sets, qz, grp, pm$perms, strata,
    as.integer(c(pm$n_perm, seed, cpus)) )

  list(ss = t(res$ss), total = res$total)
}
