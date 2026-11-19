# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' How a result was made
#'
#' ecodive's results record each step that made them in one attribute per
#' step, so that the record travels with the result through later steps:
#'
#' * `dist`: how a distance matrix was computed, set by [beta_div()], its
#'   metric functions (such as [bray()]), [vlr()], [rho_p()] and [phi_s()].
#' * `ord`: how an ordination was computed, set by [pcoa()], [pca()],
#'   [dbrda()], [procrustes()] and every other ordination.
#'
#' Each record is a small list of class `ecodive_step`, which `str()` and
#' `print()` show on one line. Its first element, `name`, is the step's
#' display name (such as `"Bray-Curtis Dissimilarity"` or `"Principal
#' Coordinates Analysis"`), so that any package can label it without
#' ecodive's tables; the rest are the settings it was made with. Read them
#' as a list: `attr(x, "dist")$metric`.
#'
#' The records accumulate as results pass along: an ordination of a
#' distance matrix keeps the matrix's `dist` record beside its own `ord`
#' record, and the permutation tests ([permanova()], ...) keep their
#' distance matrix's `dist`, and [envfit()] its ordination's `dist` and
#' `ord`. `plot()` titles its plots with them, such as "Principal
#' Coordinates Analysis of Bray-Curtis Dissimilarity".
#'
#' @section The `dist` record:
#' * `name`: The metric's display name, as [list_metrics()] gives it.
#' * `metric`: Its ID, such as `"bray"`.
#' * `norm`: The normalization applied before computing it (by name, such
#'   as `"percent"`), whether chosen with `norm` or fixed by the metric
#'   (such as [hellinger()]'s). Absent for the UniFrac metrics.
#' * `pseudocount`: The pseudocount used (`0` for none), for the metrics
#'   that take one. `denominator`, for [rho_p()] and [phi_s()].
#' * `power`, `alpha`: For [minkowski()] and [generalized_unifrac()].
#' * `clust`, `order`: With the `clust` argument: the clustering method,
#'   and each sample's position in the input.
#'
#' @section The `ord` record:
#' * `name`: The method's display name, as [list_dist_ords()] and
#'   [list_count_ords()] give it.
#' * `method`: Its ID, such as `"pcoa"` (also the ordination's `method`
#'   attribute).
#' * The other arguments it was called with that are single values, such
#'   as `n_axes`, `correction`, `seed`, `norm` or `formula` (as text). Data
#'   arguments (`dm`, `counts`, `data`, ...) and `cpus`, which never changes
#'   the result, are left out.
#'
#' @name ecodive_steps
#' @examples
#'     ord <- pcoa(bray(iris[, 1:4], norm = 'percent'))
#'     str(ord)
#'     attr(ord, 'dist')
#'     attr(ord, 'dist')$name
#'     attr(ord, 'ord')$n_axes
NULL



# A step record: `name`, then the settings in `...` and in the list
# `settings`, without the NULL ones.
ecodive_step <- function (name, ..., settings = list()) {
  v <- c(list(...), settings)
  v <- v[!vapply(v, is.null, logical(1))]
  structure(c(list(name = name), v), class = 'ecodive_step')
}


#' @export
format.ecodive_step <- function (x, ...) {
  v  <- unclass(x)
  nm <- v$name
  v$name <- NULL
  if (length(v) == 0) return (paste0('', nm))
  val <- vapply(v, function (z) {
    if (inherits(z, 'ecodive_step')) return (paste0('<', z$name, '>'))
    if (is.atomic(z) && length(z) <= 3) return (paste(format(z, trim = TRUE), collapse = ' '))
    sprintf('%s [1:%d]', if (is.atomic(z)) typeof(z) else class(z)[[1]], length(z))
  }, '')
  sprintf('%s (%s)', nm, paste(names(v), val, sep = ' = ', collapse = ', '))
}


#' @export
print.ecodive_step <- function (x, ...) {
  cat(format(x), '\n', sep = '')
  invisible(x)
}


# One line in str(), for records in the attributes of a result.
#' @export
str.ecodive_step <- function (object, ...) {
  cat(' ', format(object), '\n', sep = '')
  invisible()
}



# The display name of a step record, or NULL.
step_name <- function (step) {
  if (!inherits(step, 'ecodive_step')) return (NULL)
  nm <- step$name
  if (is.character(nm) && length(nm) == 1 && !is.na(nm)) nm
}


# A title for a result: `what`, then " of " and the name of the step it was
# made from (`from`, a record), if any.
step_title <- function (what, from) {
  nm <- step_name(from)
  if (is.null(nm) || is.null(what)) return (what)
  paste(what, 'of', nm)
}



# The `ord` record of an ordination by method `method` (its ID, the name of
# its function): its display name, and the single-valued arguments that the
# function was called with, read from its frame, which is found up the call
# stack. `...` adds settings.
ord_step <- function (method, ...) {

  f <- if (is.character(method) && length(method) == 1 &&
           exists(method, envir = ENV, mode = 'function', inherits = FALSE))
    get(method, envir = ENV, mode = 'function', inherits = FALSE)

  settings <- list()
  if (!is.null(f)) {
    for (i in rev(seq_len(sys.nframe()))) {
      if (identical(sys.function(i), f)) {
        settings <- step_settings(sys.frame(i), names(formals(f)))
        break
      }
    }
  }
  settings$method <- NULL

  ecodive_step(ord_method_name(method), method = method, ..., settings = settings)
}


# The arguments of a step worth recording, from its function's frame `env`:
# the single values among its arguments `nms`, formulas as text,
# normalizations by name, and whole numbers as integers. Data arguments and
# `cpus` are left out.
STEP_SKIP <- c('counts', 'dm', 'data', 'ord', 'target', 'tree', 'cpus', '...', 'x', 'y', 'given')

step_settings <- function (env, nms) {
  out <- list()
  for (nm in setdiff(nms, STEP_SKIP)) {
    if (!exists(nm, envir = env, inherits = FALSE)) next
    v <- tryCatch(get(nm, envir = env, inherits = FALSE), error = function (e) NULL)
    if (inherits(v, 'formula')) v <- paste(deparse(v, width.cutoff = 500L), collapse = ' ')
    if (identical(nm, 'norm') && is.numeric(v)) v <- norm_name(v)
    if (is.atomic(v) && length(v) == 1 && is.null(dim(v))) {
      if (is.factor(v)) v <- as.character(v)
      # Whole numbers as integers, so that n_axes = 3 and 3L record alike.
      if (is.double(v) && is.finite(v) && v == round(v) && abs(v) < .Machine$integer.max)
        v <- as.integer(v)
      out[[nm]] <- v
    }
  }
  out
}
