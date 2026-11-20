# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#------------------------------------------------------------------------------
# The design of the differential-abundance model methods (ancombc(),
# ancombc2(), da_lm()): `data` aligned with the samples, `formula` turned
# into a model matrix with an intercept (unless removed) and treatment
# contrasts for every factor, as patsy's in scikit-bio.
#
# `groups` follows the convention of da_test(): with `formula`, the name of
# the term of interest (a column of `data` and a term of `formula`); alone,
# the same as `formula = ~ groups`: a column name of `data`, or one value
# per sample, which becomes the `data` column `group`.
#
#   formula, data, groups  As the methods' arguments.
#   labels, n              The samples' names (or NULL) and number.
#   what                   The method's name, for messages.
#   standardize            TRUE: centre and scale every numeric variable of
#                          `formula` (scale()) before the model matrix is
#                          made, as LinDA does (linda()).
#
# Returns a list:
#   X          n x P model matrix (assign and contrasts attributes kept).
#   formula    The formula used.
#   groups     The term of interest, or NULL.
#   intercept  0-based column of the intercept, or -1.
#   factors    One integer vector per factor with a main-effect term, for
#              C_ancombc: c(standard, rebase, L, level_col, codes) (see
#              src/ancombc.c).
#------------------------------------------------------------------------------
da_design <- function (formula, data, groups, labels, n, what, standardize = FALSE) {

  if (is.null(formula) && is.null(groups))
    stop(
      '`formula` or `groups` is required, such as `formula = ~ status + age` ',
      'or `groups = "status"`.', call. = FALSE )

  if (!is.null(data) && is.null(dim(data)) && is.atomic(data)) {
    nms  <- names(data)
    data <- data.frame(group = data, stringsAsFactors = FALSE)
    if (!is.null(nms) && !anyDuplicated(nms)) rownames(data) <- nms
  }
  if (!is.null(data) && !is.data.frame(data)) {
    if (!is.matrix(data))
      stop('`data` must be a data.frame, matrix, vector, or factor.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  idx <- NULL
  if (!is.null(data)) {
    dn  <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx <- match_rows(labels, n, dn, nrow(data), ref = '`counts`')
  }


  #--------------------------------------------------------
  # `groups`: a column name, or values (which become the
  # column `group`).
  #--------------------------------------------------------
  group_vals <- NULL
  if (!is.null(groups)) {

    if (!is.atomic(groups) || !is.null(dim(groups)))
      stop(
        '`groups` must be the name of a column of `data`, or a vector with ',
        'one value per sample.', call. = FALSE )

    is_name <- is.character(groups) && length(groups) == 1L && !is.na(groups) &&
      !is.null(data) && groups %in% names(data)

    if (!is_name) {
      if (!is.null(formula))
        stop(
          'With `formula`, `groups` must name its term of interest, a ',
          'column of `data` in `formula`.', call. = FALSE )
      if (!is.null(data) && 'group' %in% names(data))
        stop(
          '`data` already has a column named "group": give `groups` as the ',
          'name of a column of `data`.', call. = FALSE )
      group_vals <- per_sample(groups, 'groups', data, idx, labels, n, ref = '`counts`')
      groups     <- 'group'
    }

    if (is.null(formula))
      formula <- stats::as.formula(call('~', as.name(groups)), env = globalenv())
  }

  if (is.null(data)) {
    if (is.null(group_vals))
      stop('`data` is required: a data.frame of sample metadata.', call. = FALSE)
    data <- data.frame(group = group_vals, stringsAsFactors = FALSE)
  }
  else {
    data <- data[idx, , drop = FALSE]
    if (!is.null(group_vals)) data[['group']] <- group_vals
  }


  #--------------------------------------------------------
  # The formula.
  #--------------------------------------------------------
  env <- globalenv()
  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    txt <- if (grepl('~', formula, fixed = TRUE)) formula else paste('~', formula)
    f   <- tryCatch(str2lang(txt), error = function (e) NULL)
    if (!is.call(f) || !identical(f[[1]], as.name('~')))
      stop('`formula` must be a formula, such as `~ status + age`.', call. = FALSE)
    formula <- eval(f)
    environment(formula) <- env
  }
  if (!inherits(formula, 'formula'))
    stop('`formula` must be a formula, such as `~ status + age`.', call. = FALSE)
  if (length(formula) == 3L) {
    env     <- environment(formula)
    formula <- formula[-2L]
    environment(formula) <- if (is.null(env)) globalenv() else env
  }

  tt <- tryCatch(
    stats::terms(formula, data = data),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  vars    <- all.vars(tt)
  missing <- vars[!vars %in% names(data)]
  if (length(missing))
    stop('Not found in `data`: ', name_list(missing), call. = FALSE)

  na <- vars[vapply(data[vars], anyNA, logical(1))]
  if (length(na))
    stop('Missing values (NA) in `data` column', if (length(na) > 1) 's', ': ',
         name_list(na), call. = FALSE)

  inf <- vars[vapply(data[vars], function (v) is.numeric(v) && !all(is.finite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in `data` column', if (length(inf) > 1) 's', ': ',
         name_list(inf), call. = FALSE)


  if (standardize) {
    const <- vars[vapply(data[vars], function (v) is.numeric(v) && !isTRUE(stats::sd(v) > 0), logical(1))]
    if (length(const))
      stop(
        'A numeric variable in `formula` has a single value, so it can\'t be ',
        'standardized: ', name_list(const), call. = FALSE )
    for (v in vars) if (is.numeric(data[[v]])) data[[v]] <- as.vector(scale(data[[v]]))
  }


  #--------------------------------------------------------
  # Categorical variables: factors with treatment contrasts,
  # levels in their own order (factors) or sorted.
  #--------------------------------------------------------
  as_levels <- function (x) {
    if (is.factor(x))    return (factor(as.character(x), levels = levels(x)[levels(x) %in% as.character(x)]))
    if (is.logical(x))   return (factor(x, levels = c(FALSE, TRUE)[c(FALSE, TRUE) %in% x]))
    if (is.character(x)) return (factor(x, levels = sort(unique(x), method = 'radix')))
    x
  }
  for (v in vars) data[[v]] <- as_levels(data[[v]])

  mf <- if (standardize) suppressWarnings(stats::model.frame(tt, data, na.action = stats::na.pass))
        else stats::model.frame(tt, data, na.action = stats::na.pass)
  for (i in seq_along(mf)) mf[[i]] <- as_levels(mf[[i]])

  is_fac <- vapply(mf, is.factor, logical(1))
  one    <- names(mf)[is_fac][vapply(mf[is_fac], nlevels, integer(1)) < 2L]
  if (length(one))
    stop('A categorical variable in `formula` has a single value: ', name_list(one), call. = FALSE)

  ctr <- stats::setNames(rep(list('contr.treatment'), sum(is_fac)), names(mf)[is_fac])
  X   <- tryCatch(
    stats::model.matrix(tt, mf, contrasts.arg = if (length(ctr)) ctr),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )
  if (standardize && !all(is.finite(X)))
    stop(
      'The design has non-finite values. With `standardize = TRUE`, ',
      'functions in `formula` (such as log(x)) apply to the standardized ',
      'variables: use `standardize = FALSE`.', call. = FALSE )

  P <- ncol(X)
  if (P == 0L)
    stop('`formula` has no terms.', call. = FALSE)
  if (n <= P)
    stop(
      what, ' needs more samples than design columns: ', n, ' sample',
      if (n != 1L) 's', ', ', P, ' column', if (P != 1L) 's', '.', call. = FALSE )
  if (qr(X)$rank < P)
    warning(
      'The design is rank deficient: some of its columns are combinations of ',
      'others, so their coefficients are not unique.', call. = FALSE )


  #--------------------------------------------------------
  # Factors with a main-effect term (_r_fit_info()).
  #--------------------------------------------------------
  asg       <- attr(X, 'assign')
  fac       <- attr(tt, 'factors')
  ord       <- attr(tt, 'order')
  intercept <- if (attr(tt, 'intercept') == 1L) 0L else -1L
  factors   <- list()

  for (v in names(mf)[is_fac]) {
    if (!v %in% rownames(fac)) next
    main <- which(fac[v, ] > 0 & colSums(fac > 0) == 1L)
    if (length(main) != 1L) next
    cols <- which(asg == main) - 1L
    L    <- nlevels(mf[[v]])
    if      (length(cols) == L - 1L) { standard <- 1L; lcol <- c(-1L, cols) }
    else if (length(cols) == L)      { standard <- 0L; lcol <- cols }
    else next
    in_ia  <- any(fac[v, ] > 0 & ord > 1L)
    rebase <- as.integer(standard == 1L && !in_ia && intercept >= 0L)
    factors[[length(factors) + 1L]] <- as.integer(c(
      standard, rebase, L, lcol, as.integer(mf[[v]]) - 1L ))
  }


  #--------------------------------------------------------
  # The term of interest. The post-hoc tests check that it
  # is categorical, with at least three levels.
  #--------------------------------------------------------
  if (!is.null(groups)) {
    term <- match(groups, attr(tt, 'term.labels'))
    if (is.na(term))
      stop('`groups` must be a term in `formula`: ', name_list(groups), ' is not.', call. = FALSE)
  }

  # For the post-hoc tests: the terms, and the levels of each factor.
  attr(X, 'term.labels') <- attr(tt, 'term.labels')
  attr(X, 'xlevels')     <- lapply(mf[is_fac], levels)

  list(
    X         = X,
    formula   = formula,
    groups    = groups,
    intercept = intercept,
    factors   = factors )
}
