# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Titles of results
#'
#' ecodive's distance matrices, ordinations and tests carry one attribute
#' for ecodive's own `plot()` and `print()`: `display`, a list holding
#' `title`, the title that `plot()` gives the result unless `main` is
#' given, and for a distance matrix made with `clust`, the clustering
#' method `clust`, which `plot()` keeps and names in its subtitle.
#'
#' A title is the method's name, then " of " and the title of the result it
#' was computed from, when that has one:
#'
#' * **Distance matrices**: the metric's name, as [list_metrics()] gives
#'   it, such as `"Bray-Curtis Dissimilarity"`.
#' * **Ordinations**: the method's name, as [list_dist_ords()] and
#'   [list_count_ords()] give it, and for an ordination of distances, the
#'   distance matrix's title: `"Principal Coordinates Analysis of
#'   Bray-Curtis Dissimilarity"`.
#' * **Permutation tests**: the test's name and the distance matrix's
#'   title; for [mantel()], both matrices' titles when both have one.
#'   `print()` shows it as the heading.
#' * [envfit()]: "Fitted variables on" and the ordination's title.
#'
#' Results don't record how they were made: neither the arguments they were
#' made with, nor earlier steps, such as normalizing, rarefying or
#' transforming the counts, which return plain tables. A title names only
#' what ecodive knows to be true of the result.
#'
#' Arithmetic on a distance matrix, such as `sqrt(d)` or `1 - d`, gives new
#' values, so the result loses its title and its `method` attribute (the
#' metric's ID), and plots as "Distances"; it keeps `clust` and `order`.
#' To title a result yourself, set `attr(x, "display")$title`, or give
#' `plot()` a `main`.
#'
#' A result's other attributes are results in their own right, such as an
#' ordination's `eig`, or tell other functions how to read it, such as a
#' distance matrix's `method` and, with `clust`, `order`. Each function's
#' page describes them.
#'
#' @name ecodive_display
#' @examples
#'     d   <- bray(iris[, 1:4])
#'     ord <- pcoa(d)
#'     ord
#'     attr(ord, 'display')$title
#'
#'     # Arithmetic drops the title.
#'     attr(sqrt(d), 'display')$title
#'
#'     # Your own title.
#'     attr(ord, 'display')$title <- 'Iris flowers'
#'     plot(ord, 'Species', iris)
NULL



# `x` with the values in `...` set in its `display` list (a NULL value
# removes one). The list is kept last among the attributes, and dropped
# when empty.
set_display <- function (x, ...) {
  disp <- attr(x, 'display', exact = TRUE)
  if (!is.list(disp)) disp <- list()
  new <- list(...)
  for (nm in names(new)) {
    if (is.null(new[[nm]])) disp[[nm]] <- NULL
    else                    disp[[nm]] <- new[[nm]]
  }
  attr(x, 'display') <- NULL
  if (length(disp)) attr(x, 'display') <- disp
  x
}


# One value of the `display` list of `x`, or NULL.
get_display <- function (x, name) {
  disp <- attr(x, 'display', exact = TRUE)
  if (is.list(disp)) disp[[name]]
}


# The title of `x`: a single non-empty string, or NULL.
display_title <- function (x) {
  title <- get_display(x, 'title')
  if (is.character(title) && length(title) == 1 && !is.na(title) && nzchar(title)) title
}


# The title for plot(): the title of `x`, or else `fallback`.
plot_title <- function (x, fallback = NULL) {
  title <- display_title(x)
  if (is.null(title)) fallback else title
}


# The title of a result made by `what` (a display name) from an input
# titled `from` (a string, or NULL): "<what> of <from>", or `what` alone.
title_of <- function (what, from = NULL) {
  if (is.character(from) && length(from) == 1 && !is.na(from) && nzchar(from) && !is.null(what))
    paste(what, 'of', from)
  else what
}



# Arithmetic on a distance matrix gives new values, which the metric's title
# and ID no longer describe; `clust` and `order` still hold (see
# ?ecodive_display).
#' @export
Ops.ecodive_dist <- function (e1, e2) {
  dist_unlabel(NextMethod())
}

#' @export
Math.ecodive_dist <- function (x, ...) {
  dist_unlabel(NextMethod())
}

dist_unlabel <- function (d) {
  attr(d, 'method') <- NULL
  set_display(d, title = NULL)
}
