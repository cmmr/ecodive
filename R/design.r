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
#   match_rows()     which row of a table holds each sample (also envfit()).
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
#            value per row of `data`, in the same order. A column named
#            here is left out of the default formula (`~ .`).
#   allow_empty  Whether a formula with no terms to test (only Condition()
#            terms) is allowed: dbrda() then gives the unconstrained
#            ordination of the conditioned distances.
#   env      Where a character `formula`'s variables are looked up.
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
#   idx          The row of `data` (as given) for each sample.
#------------------------------------------------------------------------------
model_design <- function (
    data, formula = NULL, labels = NULL, n, strata = NULL,
    allow_empty = FALSE, env = parent.frame(2) ) {

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
  strata_col <- NULL
  if (!is.null(strata)) {
    if (is.character(strata) && length(strata) == 1 && strata %in% names(data)) {
      strata_col <- strata
      strata     <- data[[strata]]
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
  # Rows for samples not in `dm` (e.g. dropped by rarefying) are left out.
  dn  <- if (.row_names_info(data) > 0) rownames(data) else NULL
  idx <- match_rows(labels, n, dn, nrow(data))

  # Factor levels with no samples left (e.g. only in rows left out) would
  # give all-zero, aliased model matrix columns.
  data <- droplevels(data[idx, , drop = FALSE])
  if (!is.null(labels) && !anyNA(labels) && !anyDuplicated(labels))
    rownames(data) <- labels

  if (!is.null(strata)) {
    strata <- strata[idx]
    if (anyNA(strata)) stop('`strata` must not contain NAs.', call. = FALSE)
    strata <- as.integer(factor(strata, exclude = NULL)) - 1L
  }


  #--------------------------------------------------------
  # Parse the formula; split off Condition() terms.
  #--------------------------------------------------------
  if (is.null(formula)) {
    formula <- ~ .
    # A strata column is a blocking factor, not a term to test.
    if (!is.null(strata_col))
      formula <- stats::as.formula(call('~', call('-', quote(.), as.name(strata_col))))
  }
  if (is.character(formula)) {
    # Parse without evaluating: R 3.6's as.formula('trt') evaluates `trt`.
    f <- tryCatch(str2lang(formula), error = function (e) NULL)
    if (!is.call(f) || !identical(f[[1]], as.name('~')))
      stop('`formula` must be a formula, such as `~ treatment + ph`.', call. = FALSE)
    formula <- eval(f)
    environment(formula) <- env
  }
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

  inf <- used[vapply(data[used], function (v) is.numeric(v) && !all(is.finite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in `data` column', if (length(inf) > 1) 's', ': ',
         name_list(inf), call. = FALSE)

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

  if (length(cond_terms) == length(term_labels)) {

    if (!allow_empty)
      stop('`formula` has no terms to test.', call. = FALSE)

    # Only conditions: nothing to constrain.
    return (list(
      x           = matrix(0, nrow(data), 0, dimnames = list(rownames(data), NULL)),
      assign      = integer(0),
      term_labels = character(0),
      terms       = NULL,
      z           = z,
      strata      = strata,
      data        = data,
      idx         = idx ))
  }

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
    data        = data,
    idx         = idx )
}



#------------------------------------------------------------------------------
# Which row of a table (`what`, with row names `names` or NULL, and `n_rows`
# rows) holds each of the n samples of `ref` (named `labels`, or NULL).
# Returns n indices; rows for other samples are left out.
#
# Matched by name when both have names. Names that are only numbers, as from
# subsetting a data frame without row names, count as names only when they
# match every label. Otherwise matched by position, so the counts must agree.
# `unit` is what the table's rows are called in messages.
#------------------------------------------------------------------------------
match_rows <- function (
    labels, n, names, n_rows, ref = '`dm`', what = '`data`', unit = 'row' ) {

  named <- !is.null(labels) && !is.null(names)
  if (named && all(grepl("^[0-9]+$", names)))
    named <- all(labels %in% names)

  if (named) {

    if (anyDuplicated(labels))
      stop(ref, ' has duplicate sample names, so ', what, ' cannot be matched to them.', call. = FALSE)
    if (anyDuplicated(names))
      stop(what, ' has duplicate ', unit, ' names.', call. = FALSE)

    idx <- match(labels, names)

    if (anyNA(idx)) {
      miss <- labels[is.na(idx)]
      stop(
        sprintf('%d sample%s in %s %s no %s in %s: %s',
          length(miss), if (length(miss) == 1) '' else 's', ref,
          if (length(miss) == 1) 'has' else 'have', unit, what, name_list(miss) ),
        call. = FALSE )
    }

    return (idx)
  }

  if (n_rows != n)
    stop(
      sprintf('%s has %d %s%s, but %s has %d sample%s.', what, n_rows, unit,
              if (n_rows == 1) '' else 's', ref, n, if (n == 1) '' else 's'),
      if (!is.null(labels) && unit == 'value')
        sprintf('\nName the values of %s to match them to the samples by name.', what)
      else if (!is.null(labels))
        sprintf('\nGive %s %s names to match its %ss to the samples by name.', what, unit, unit),
      call. = FALSE )

  seq_len(n)
}



#------------------------------------------------------------------------------
# One rule for every per-sample argument (`groups`, `y`, `given`, `strata`):
#
#   * With `data`: the name of one of its columns, or a vector with one value
#     per row of `data`, in the order of its rows.
#   * Without `data`: a vector with one value per sample, matched to the
#     samples by name when both are named, else in the samples' order.
#
#   v       The argument's value: NULL, a column name, or an atomic vector.
#   arg     Its name, for messages.
#   data    NULL, or the metadata as a data frame, with its rows as given.
#   idx     With `data`: the row of `data` for each sample (match_rows()).
#   labels  The samples' names, or NULL.
#   n       The number of samples.
#   ref     What the samples belong to, for messages.
#
# Returns the values in the samples' order (NULL for NULL).
#------------------------------------------------------------------------------
per_sample <- function (v, arg, data, idx, labels, n, ref = '`dm`') {

  if (is.null(v)) return (NULL)

  if (!is.atomic(v) || !is.null(dim(v)))
    stop(
      '`', arg, '` must be a vector with one value per sample',
      if (!is.null(data)) ', or the name of a column of `data`', '.', call. = FALSE )

  one_name <- is.character(v) && length(v) == 1

  if (!is.null(data)) {

    if (one_name && !is.na(v) && v %in% names(data))
      return (data[[v]][idx])

    if (one_name && nrow(data) != 1)
      stop(
        '`', arg, '` must be the name of a column of `data`, but ',
        name_list(v), ' is not one.', call. = FALSE )

    if (length(v) != nrow(data))
      stop(
        'With `data`, `', arg, '` must be the name of a column of `data`, ',
        'or have one value per row of `data` (', nrow(data), ').', call. = FALSE )

    return (v[idx])
  }

  if (one_name && n != 1)
    stop('`', arg, '` names a column, but `data` is not given.', call. = FALSE)

  nms <- names(v)
  if (!is.null(nms) && (anyNA(nms) || anyDuplicated(nms))) nms <- NULL

  v[match_rows(labels, n, nms, length(v), ref = ref, what = paste0('`', arg, '`'), unit = 'value')]
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
    if (!is.numeric(p) || anyNA(p) || ncol(p) != n || any(p != round(p)) ||
        any(p < 1 | p > n))
      stop('A `permutations` matrix must have ', n, ' columns (one per sample) ',
           'of whole numbers, with one permutation of 1 to ', n, ' per row.', call. = FALSE)

    storage.mode(p) <- 'integer'
    ok <- apply(p, 1L, function (r) all(sort.int(r, method = 'radix') == seq_len(n)))
    if (!all(ok))
      stop('Row ', which(!ok)[[1]], ' of `permutations` is not a permutation of 1 to ', n, '.', call. = FALSE)

    return(list(n_perm = nrow(p), perms = t(p) - 1L))
  }

  if (!is.numeric(permutations) || length(permutations) != 1 || !is.finite(permutations) ||
      permutations < 0 || permutations %% 1 != 0)
    stop('`permutations` must be a whole number of at least 0, or a matrix ',
         'with one permutation per row.', call. = FALSE)

  # Each permutation's statistics are kept, and memory for them is taken
  # before the first check for an interrupt.
  if (permutations > 1e8)
    stop('`permutations` must be at most 1e8 (100 million).', call. = FALSE)

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
