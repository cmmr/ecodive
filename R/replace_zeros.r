# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# replace_zeros_mult() follows multi_replace() and closure() from
# scikit-bio 0.7.4 (skbio/stats/composition/_base.py), Copyright (c) 2013--,
# scikit-bio development team, under the BSD-3 license given in
# src/replace_zeros.c and inst/COPYRIGHTS.
#
# replace_zeros_czm(), _gbm(), _sq() and _bl() are a clean-room
# implementation of zCompositions' cmultRepl() methods (GPL): their formulas
# come from its documentation and its outputs, never its source (see
# plans/zero_replacement.md and data-raw/expected/zero_replacement_probe.r).



# Zero replacement methods. The function for row `id` is
# replace_zeros_<id>(counts, margin = 1L, ...): each takes the table first
# and `margin` second, checks it with zero_replacement_input() and builds
# its result with zero_replacement_result(), so that a new method is a new
# row and a new function, with no change to the wrapper.
ZERO_REPLACEMENTS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                    id,   alt_ids
      Multiplicative Replacement,              mult, multiplicative multreplacement multirepl multireplace
      Count Zero Multiplicative,               czm,  countzeromultiplicative
      Geometric Bayesian-Multiplicative,       gbm,  geometricbayesianmultiplicative cmultrepl
      Square-Root Bayesian-Multiplicative,     sq,   squareroot sqrtbm
      Bayes-Laplace Bayesian-Multiplicative,   bl,   bayeslaplace laplace
  " )

  return (df)
})



#' Zero Replacement
#'
#' Replaces the zeros of a count table with small positive values, so that
#' log ratios can be taken without a pseudocount. Runs any of ecodive's zero
#' replacement methods by name, as [count_ord()] does for ordinations.
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts or proportions.
#'
#' @param method   The name of a zero replacement method: `"mult"`, `"czm"`,
#'        `"gbm"`, `"sq"` or `"bl"`. Flexible matching is supported (see
#'        below). Programmatic access via `list_zero_replacements()`.
#'        Required: there is no default, because the choice of method affects
#'        every downstream log ratio. `"gbm"` is the common choice for count
#'        tables (it is zCompositions' default).
#'
#' @param ...   Further arguments for the method's function, such as
#'        `margin`, [replace_zeros_mult()]'s `delta`, or `frac`, `threshold`
#'        and `adjust` (see [replace_zeros_gbm()]). Each method keeps its own
#'        defaults.
#'
#' @return The method's result: a dense numeric `matrix` of proportions with
#'   no zeros, in the orientation of `counts` (samples \eqn{\times} features
#'   for `margin = 1`, features \eqn{\times} samples for `margin = 2`), with
#'   its names. Each sample sums to 1. Attributes `method` (the method's id),
#'   then the method's results, ending with `n_replaced`, the number of
#'   zeros replaced. `replace_zeros(counts, 'gbm', ...)` is the same as
#'   `replace_zeros_gbm(counts, ...)`.
#'
#' @details
#'
#' Logs of zero are undefined, so every log-ratio method, such as
#' `norm = 'clr'`, [aitchison()], [alr()] and [ilr()], needs the zeros dealt
#' with. A `pseudocount` added to every value changes the ratios between the
#' non-zero parts, and acts as an arbitrary weighting. A zero replacement
#' instead gives each zero a small value and scales the sample's non-zero
#' parts down to make room for it, keeping their ratios.
#'
#' The result has no zeros, so it can be given to any of those functions
#' without a pseudocount (and without the warning that `pseudocount = NULL`
#' gives for zeros). The result is in proportions: every log-ratio method is
#' scale invariant, so nothing is lost.
#'
#' `mult` gives every zero the same proportion, whatever the sample's depth.
#' The other methods are for counts: `czm` imputes a fraction of one read in
#' each sample, and the Bayesian-multiplicative methods (`gbm`, `sq`, `bl`)
#' also use how abundant the feature is in the other samples. `gbm` is the
#' common default for count tables.
#'
#' Conclusions that depend on how zeros were treated are fragile. Repeating
#' an analysis with each method of `list_zero_replacements('id')` is a
#' simple check.
#'
#' **List of Zero Replacement Methods**
#'
#' | Option / Function Name | Method Name                                |
#' | :--------------------- | :----------------------------------------- |
#' | `mult`                 | Multiplicative replacement, with a fixed value \eqn{\delta} ([replace_zeros_mult()]) |
#' | `czm`                  | Count zero multiplicative ([replace_zeros_czm()]) |
#' | `gbm`                  | Geometric Bayesian-multiplicative ([replace_zeros_gbm()]) |
#' | `sq`                   | Square-root Bayesian-multiplicative ([replace_zeros_sq()]) |
#' | `bl`                   | Bayes-Laplace Bayesian-multiplicative ([replace_zeros_bl()]) |
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters. E.g.
#' `method = 'MULT'`, `'multiplicative'` or `'multiplicative replacement'`
#' selects `mult`, and `'GBM'` or `'geometric Bayesian-multiplicative'`
#' selects `gbm`. Common other names also work: `'multRepl'` (from
#' zCompositions) and `'multi_replace'` (from scikit-bio) for `mult`;
#' `'cmultRepl'` (zCompositions' function, whose default method is GBM) for
#' `gbm`; `'square root'` for `sq`; and `'Laplace'` for `bl`.
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#'
#' **Same names in other packages**
#'
#' None known. zCompositions' `multRepl()` and `cmultRepl()`, and
#' scikit-bio's `multi_replace()`, do the same job under other names.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V.
#' (2003). Dealing with zeros and missing values in compositional data sets
#' using nonparametric imputation. *Mathematical Geology*, 35, 253–278.
#' \doi{10.1023/A:1023866030544}
#'
#' Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
#' Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of count
#' zeros in compositional data sets. *Statistical Modelling*, 15, 134–158.
#' \doi{10.1177/1471082X14535524}
#'
#' Lubbe, S., Filzmoser, P., & Templ, M. (2021). Comparison of zero
#' replacement strategies for compositional data with large numbers of
#' zeros. *Chemometrics and Intelligent Laboratory Systems*, 210, 104248.
#' \doi{10.1016/j.chemolab.2021.104248}
#'
#' @seealso [replace_zeros_mult()], [replace_zeros_gbm()],
#'   [list_zero_replacements()]; [aitchison()], [alr()], [ilr()] and
#'   [normalize()] (`'clr'`), which take the result.
#'
#' @export
#' @examples
#'     x <- replace_zeros(ex_counts, 'gbm')
#'     x[1:3, 1:4]
#'     attr(x, 'n_replaced')
#'
#'     # No zeros left, so no pseudocount is needed.
#'     aitchison(x)
#'
#'     # The same analysis with every method.
#'     lapply(list_zero_replacements('id'), function (m) {
#'       aitchison(replace_zeros(ex_counts, m))
#'     })
#'
replace_zeros <- function (counts, method, ...) {

  if (missing(method))
    stop(
      '`method` is required: the zero replacement method, one of ',
      paste(collapse = ', ', sort(ZERO_REPLACEMENTS$id)), '. There is no ',
      'default, because the choice affects every downstream log ratio.',
      call. = FALSE )

  method <- match_zero_replacement(method)

  # quote(counts), so that a large table isn't copied into the call.
  do.call(method$func, c(list(quote(counts)), list(...)))
}



#' Find and Browse Zero Replacement Methods
#'
#' Programmatic access to the list of zero replacement methods that
#' [replace_zeros()] runs, and their associated functions.
#'
#' @param method   The name of a zero replacement method to search for.
#'        Supports partial matching. All non-alpha characters are ignored.
#'
#' @param val   Sets the return value for this function call. See "Value"
#'        section below. Default: `"data.frame"`
#'
#' @param nm   What value to use for the names of the returned object.
#'        Default is `"id"` when `val` is `"list"` or `"func"`, otherwise the
#'        default is `NA` (no name).
#'
#' @return
#'
#' **`match_zero_replacement()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Method name, e.g. `"Multiplicative Replacement"`
#' * `id` : Method ID, e.g. `"mult"`
#' * `func` : The function for this method, e.g. `ecodive::replace_zeros_mult`
#' * `params` : Formal args for `func`, e.g. `c("counts", "margin", "delta")`
#'
#'
#' **`list_zero_replacements()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_zero_replacement()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of method IDs.
#' * `val = "name"` : A character vector of method names.
#'
#' If `nm` is set, then the names of the vector or list will be the method ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname zero_replacements
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available methods.
#'     list_zero_replacements()
#'
#'     # Find the method for zCompositions' 'multRepl'.
#'     m <- match_zero_replacement('multRepl')
#'     m$id
#'
#'     # zCompositions' 'cmultRepl' defaults to GBM.
#'     match_zero_replacement('cmultRepl')$id
#'
list_zero_replacements <- function (
    val = c('data.frame', 'list', 'func', 'id', 'name'),
    nm  = c(NA, 'id', 'name') ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }

  df <- ZERO_REPLACEMENTS

  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) zero_replacement_info(df[df$id == id,]))
    if (!is.na(nm)) names(result) <- df[[nm]]
    return (result)
  }

  if (val == 'func') {
    result <- mget(paste0('replace_zeros_', df$id), ENV)
    names(result) <- if (is.na(nm)) NULL else df[[nm]]
    return (result)
  }

  result <- df[[val]]
  if (!is.na(nm)) names(result) <- df[[nm]]
  return (result)
}


#' @rdname zero_replacements
#' @export
#'
match_zero_replacement <- function (method) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')
  if (is.na(method))         stop('`method` must not be NA')

  df <- ZERO_REPLACEMENTS

  # Exact match
  if (method %in% df$id)
    return (zero_replacement_info(df[df$id == method,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  match <- match_name(needle, 'ZERO_REPLACEMENTS', df$id)

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  zero_replacement_info(df[df$id == match,])
}


# One row of ZERO_REPLACEMENTS as a list, with the method's function
# (replace_zeros_<id>) and its args.
zero_replacement_info <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(paste0('replace_zeros_', info$id), ENV)
  info$params  <- names(formals(info$func))
  info
}



#' Multiplicative Zero Replacement
#'
#' Replaces every zero with a small proportion \eqn{\delta}, and scales each
#' sample's non-zero proportions down by the same factor so that the sample
#' still sums to 1 (Martín-Fernández et al. 2003). The ratios between
#' non-zero parts are unchanged.
#'
#' @inheritParams replace_zeros
#' @inheritParams documentation
#'
#' @param delta  The proportion that every zero becomes: a single value, or
#'        one value per sample (matched to the samples by name when both are
#'        named, else in order), each between 0 and 1. `NULL` (the default)
#'        uses \eqn{1 / D^2}, with \eqn{D} the number of features, as
#'        scikit-bio does. A sample's \eqn{\delta} times its number of zeros
#'        must be below 1.
#'
#' @return A dense numeric `matrix` of proportions with no zeros, in the
#'   orientation of `counts` (samples \eqn{\times} features for `margin =
#'   1`, features \eqn{\times} samples for `margin = 2`, and for a single
#'   sample given as a vector), with its names. Each sample sums to 1.
#'   Attributes:
#'
#'   * `method`: `'mult'`.
#'   * `delta`: the \eqn{\delta} used: one value, or one per sample (named
#'     after the samples).
#'   * `n_replaced`: the number of zeros replaced.
#'
#' @details
#' For a sample with proportions \eqn{x_j = c_j / \sum_k c_k} and
#' \eqn{Z} zeros,
#' \deqn{r_j = \delta \text{ if } x_j = 0, \qquad
#'       r_j = (1 - Z\delta)\, x_j \text{ otherwise.}}
#'
#' The \eqn{\delta} is the same for every sample (unless given per sample),
#' whatever its depth: a zero in a sample of 1,000 reads and one in a sample
#' of 100,000 become the same proportion. Choose \eqn{\delta} below the
#' smallest proportion that the samples could detect, e.g. a fraction of
#' `1 / max(depths)`, or give each sample its own, e.g.
#' `delta = 0.5 / rowSums(counts)` (half a read).
#'
#' A feature that is zero in every sample is replaced like any other zero.
#' Samples with no non-zero values are an error: there is nothing to scale.
#'
#' The result equals scikit-bio's `multi_replace()` (version 0.7.4) to the
#' last bit: each sample's total is summed in numpy's (pairwise) order.
#' Unlike scikit-bio, a \eqn{\delta} that would leave a sample's non-zero
#' proportions at exactly 0 (\eqn{Z\delta = 1}) is an error too, so that
#' the result has no zeros.
#'
#' The result feeds straight into [aitchison()], [alr()], [ilr()] and any
#' function with `norm = 'clr'`, with no pseudocount.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V.
#' (2003). Dealing with zeros and missing values in compositional data sets
#' using nonparametric imputation. *Mathematical Geology*, 35, 253–278.
#' \doi{10.1023/A:1023866030544}
#'
#' @seealso [replace_zeros()], [list_zero_replacements()]
#'
#' @export
#' @examples
#'     x <- replace_zeros_mult(ex_counts)
#'     x[1:3, 1:4]
#'     attr(x, 'delta')
#'     rowSums(x)
#'
#'     # The ratios between non-zero parts are kept.
#'     nz <- ex_counts[1, ] > 0
#'     all.equal(x[1, nz] / sum(x[1, nz]), ex_counts[1, nz] / sum(ex_counts[1, nz]))
#'
#'     # Half a read in each sample.
#'     x <- replace_zeros_mult(ex_counts, delta = 0.5 / rowSums(ex_counts))
#'     attr(x, 'delta')
#'
replace_zeros_mult <- function (counts, margin = 1L, delta = NULL) {

  info <- zero_replacement_input()

  if (is.null(delta)) {
    delta <- .Call(C_mult_default_delta, info$D)
  }
  else {
    if (!is.numeric(delta) || !is.null(dim(delta)) || length(delta) == 0 ||
        anyNA(delta) || !all(is.finite(delta)) || !all(delta > 0 & delta < 1))
      stop(
        '`delta` must be a proportion between 0 and 1 (exclusive): one value, ',
        'or one per sample.', call. = FALSE )
    storage.mode(delta) <- 'double'  # keeping its names
    if (length(delta) > 1L) {
      delta <- per_sample(
        delta, 'delta', NULL, NULL, info$samples, info$n, ref = '`counts`' )
      names(delta) <- info$samples
    }
    else {
      delta <- unname(delta)
    }
  }

  # A sample's non-zero proportions are scaled by 1 - (its zeros) * delta.
  bad <- which(1 - info$sample_zeros * delta <= 0)
  if (length(bad))
    stop(
      sprintf(
        paste0(
          'Multiplicative replacement created negative proportions (or zeros) ',
          'in %d sample%s (%s): `delta` times a sample\'s number of zeros must ',
          'be below 1. Consider using a smaller `delta`.' ),
        length(bad), if (length(bad) == 1) '' else 's',
        name_list(if (is.null(info$samples)) bad else info$samples[bad]) ),
      call. = FALSE )

  res <- .Call(C_replace_zeros_mult, counts, margin, delta)

  zero_replacement_result(
    res, info, margin, 'mult', list(delta = delta) )
}



#' Bayesian-Multiplicative and Count Zero Replacement
#'
#' Replaces the zeros of a count table with small proportions that depend on
#' each sample's depth, and scales each sample's non-zero proportions down by
#' the same factor, keeping their ratios. The methods of zCompositions'
#' `cmultRepl()` (Martín-Fernández et al. 2015): geometric
#' Bayesian-multiplicative (`gbm`, the usual choice), square-root
#' Bayesian-multiplicative (`sq`), Bayes-Laplace Bayesian-multiplicative
#' (`bl`), and count zero multiplicative (`czm`).
#'
#' @inheritParams replace_zeros
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of non-negative counts (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Each sample's total is its number of reads, so give counts,
#'        not proportions (other values give a warning).
#'
#' @param frac  For `czm`, the fraction of the upper limit (`threshold`
#'        reads) that a zero becomes; and, with `adjust`, the fraction of a
#'        feature's smallest observed proportion that an imputed value above
#'        it becomes. Between 0 and 1 (exclusive). Default: `0.65`.
#'
#' @param threshold  For `czm`, the upper limit of a zero, in reads: a zero
#'        in a sample of \eqn{n} reads is below `threshold / n`. Above 0 and
#'        at most 1. Default: `0.5`.
#'
#' @param adjust  If `TRUE` (the default), an imputed proportion above the
#'        smallest proportion observed for that feature in any sample is
#'        replaced by `frac` times that proportion, so that no zero is
#'        imputed above what was actually seen.
#'
#' @return A dense numeric `matrix` of proportions with no zeros, in the
#'   orientation of `counts` (samples \eqn{\times} features for `margin =
#'   1`, features \eqn{\times} samples for `margin = 2`, and for a single
#'   sample given as a vector), with its names. Each sample sums to 1.
#'   Attributes:
#'
#'   * `method`: `'gbm'`, `'sq'`, `'bl'` or `'czm'`.
#'   * `n_adjusted`: the number of imputed values that `adjust` lowered.
#'   * `n_replaced`: the number of zeros replaced.
#'
#' @details
#' For sample \eqn{i} with counts \eqn{c_{ij}}, depth
#' \eqn{n_i = \sum_j c_{ij}}, proportions \eqn{x_{ij} = c_{ij} / n_i} and
#' \eqn{D} features, each zero becomes a proportion \eqn{r_{ij}}, and each
#' non-zero proportion is scaled so that the sample sums to 1:
#' \deqn{r_{ij} = x_{ij} \Bigl(1 - \sum_{k:\, c_{ik} = 0} r_{ik}\Bigr)
#'       \quad \text{if } c_{ij} > 0.}
#'
#' **Bayesian-multiplicative methods** (`gbm`, `sq`, `bl`). A zero is the
#' posterior estimate of a multinomial probability under a Dirichlet prior
#' with expectation \eqn{t_{ij}} and strength \eqn{s_i}:
#' \deqn{r_{ij} = \frac{s_i\, t_{ij}}{n_i + s_i} \quad \text{if } c_{ij} = 0.}
#' The prior expectation is the feature's pooled proportion in the other
#' samples,
#' \deqn{t_{ij} = \frac{\sum_{k \ne i} c_{kj}}{\sum_{k \ne i} n_k},}
#' so a feature that is abundant elsewhere gets a larger value. The methods
#' differ in the strength:
#'
#' | Method | \eqn{s_i} |
#' | :----- | :-------- |
#' | `gbm`  | \eqn{1 / g_i}, with \eqn{g_i} the geometric mean of \eqn{t_{i1}, \ldots, t_{iD}} |
#' | `sq`   | \eqn{\sqrt{n_i}} |
#' | `bl`   | \eqn{D} |
#'
#' **Count zero multiplicative** (`czm`). Every zero of sample \eqn{i}
#' becomes \eqn{r_{ij} = frac \cdot threshold / n_i}: by default 65% of
#' half a read. The other samples are used only by `adjust`.
#'
#' **Adjustment.** With `adjust = TRUE`, an \eqn{r_{ij}} above
#' \eqn{m_j = \min_{k:\, c_{kj} > 0} x_{kj}}, the smallest proportion of
#' feature \eqn{j} observed in any sample, becomes \eqn{frac \cdot m_j}
#' (before the non-zero parts are scaled). This happens mostly in shallow
#' samples.
#'
#' **Requirements.** No sample may be all zeros. The Bayesian methods need
#' at least two samples, and no feature may be zero in every sample (its
#' prior would be 0). `gbm` also needs every feature to be non-zero in at
#' least two samples: a feature seen in a single sample makes that sample's
#' \eqn{t_{ij}} zero, its geometric mean zero and its prior strength
#' infinite. Filter such features out first (e.g.
#' `counts[, colSums(counts > 0) >= 2]`), or use `sq` or `bl`, which handle
#' them. `czm` works on a single sample, and replaces the zeros of a
#' feature that is zero everywhere like any other.
#'
#' The proportion of zeros is not checked: features that are zeros in most
#' samples, as are most features of a microbiome table, are imputed mostly
#' from their prior. zCompositions warns about (or drops) rows and columns
#' that are at least 80% zeros; ecodive never drops any.
#'
#' **Compared with zCompositions.** The results equal `cmultRepl(counts,
#' method = ..., output = "prop", z.delete = FALSE)` of zCompositions 1.6.2,
#' to within about \eqn{10^{-15}} (relative). Unlike `cmultRepl()`, a table
#' without zeros gives its proportions instead of an error, `czm` accepts a
#' single sample and features that are zero everywhere, and `czm` refuses
#' zeros whose imputed values would leave the sample's non-zero parts at
#' zero or below.
#'
#' The result feeds straight into [aitchison()], [alr()], [ilr()] and any
#' function with `norm = 'clr'`, with no pseudocount.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
#' Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of count
#' zeros in compositional data sets. *Statistical Modelling*, 15, 134–158.
#' \doi{10.1177/1471082X14535524}
#'
#' Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015). zCompositions
#' — R package for multivariate imputation of left-censored data under a
#' compositional approach. *Chemometrics and Intelligent Laboratory
#' Systems*, 143, 85–96. \doi{10.1016/j.chemolab.2015.02.019}
#'
#' @seealso [replace_zeros()], [replace_zeros_mult()],
#'   [list_zero_replacements()]
#'
#' @rdname replace_zeros_cmult
#' @export
#' @examples
#'     # Geometric Bayesian-multiplicative replacement.
#'     x <- replace_zeros_gbm(ex_counts)
#'     x
#'
#'     # The ratios between non-zero parts are kept.
#'     nz <- ex_counts[1, ] > 0
#'     all.equal(x[1, nz] / sum(x[1, nz]), ex_counts[1, nz] / sum(ex_counts[1, nz]))
#'
#'     # The other methods.
#'     replace_zeros_sq(ex_counts)
#'     replace_zeros_bl(ex_counts)
#'     replace_zeros_czm(ex_counts)
#'
#'     # Without the adjustment.
#'     attr(replace_zeros_czm(ex_counts), 'n_adjusted')
#'     attr(replace_zeros_czm(ex_counts, adjust = FALSE), 'n_adjusted')
#'
replace_zeros_gbm <- function (counts, margin = 1L, frac = 0.65, adjust = TRUE) {
  cmult_replace('gbm', frac, adjust = adjust)
}


#' @rdname replace_zeros_cmult
#' @export
replace_zeros_sq <- function (counts, margin = 1L, frac = 0.65, adjust = TRUE) {
  cmult_replace('sq', frac, adjust = adjust)
}


#' @rdname replace_zeros_cmult
#' @export
replace_zeros_bl <- function (counts, margin = 1L, frac = 0.65, adjust = TRUE) {
  cmult_replace('bl', frac, adjust = adjust)
}


#' @rdname replace_zeros_cmult
#' @export
replace_zeros_czm <- function (
    counts, margin = 1L, frac = 0.65, threshold = 0.5, adjust = TRUE ) {
  cmult_replace('czm', frac, threshold, adjust)
}


# The shared body of replace_zeros_czm(), _gbm(), _sq() and _bl(), with the
# caller's `counts` and `margin`. See src/replace_zeros.c for the formulas.
cmult_replace <- function (
    method, frac, threshold = 0.5, adjust = TRUE, env = parent.frame() ) {

  if (!is.numeric(frac) || length(frac) != 1 || is.na(frac) || !(frac > 0 && frac < 1))
    stop('`frac` must be a single number between 0 and 1 (exclusive).', call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
      !(threshold > 0 && threshold <= 1))
    stop('`threshold` must be a single number above 0 and at most 1.', call. = FALSE)
  if (!is.logical(adjust) || length(adjust) != 1 || is.na(adjust))
    stop('`adjust` must be TRUE or FALSE.', call. = FALSE)

  bayes <- method != 'czm'
  info  <- zero_replacement_input(
    empty_features_ok = !bayes, whole = TRUE,
    min_samples = if (bayes) 2L else 1L, env = env )

  if (method == 'gbm') {
    single <- which(info$feature_nonzeros == 1L)
    if (length(single))
      stop(
        sprintf(
          paste0(
            '`counts` has %d feature%s that %s non-zero in only one sample ',
            '(%s). GBM cannot be computed for that sample: the feature\'s ',
            'prior there is zero, so the prior strength is infinite. ',
            'Remove %s first (e.g. keep the features found in at least 2 ',
            'samples), or use `sq` or `bl`.' ),
          length(single), if (length(single) == 1) '' else 's',
          if (length(single) == 1) 'is' else 'are',
          name_list(if (is.null(info$features)) single else info$features[single]),
          if (length(single) == 1) 'it' else 'them' ),
        call. = FALSE )
  }

  code <- match(method, c('czm', 'bl', 'sq', 'gbm')) - 1L
  res  <- .Call(
    C_replace_zeros_cmult, env$counts, env$margin, code,
    as.double(frac), as.double(threshold), adjust )

  # Only czm can impute more than a sample holds: a Bayesian imputation is
  # below the sample's share n / (n + s) of its zero parts' priors.
  bad <- which(!(res$kept > 0))
  if (length(bad))
    stop(
      sprintf(
        paste0(
          'Count zero multiplicative replacement created negative proportions ',
          '(or zeros) in %d sample%s (%s): the values imputed for %s zeros ',
          'add up to 1 or more. Use a smaller `frac` or `threshold`.' ),
        length(bad), if (length(bad) == 1) '' else 's',
        name_list(if (is.null(info$samples)) bad else info$samples[bad]),
        if (length(bad) == 1) 'its' else 'their' ),
      call. = FALSE )

  zero_replacement_result(
    res$result, info, env$margin, method, list(n_adjusted = res$n_adjusted) )
}



# Checks the table of a zero replacement in the calling function, and sets
# its `counts` (a table the C code reads) and `margin`. Returns the table's
# sizes and names, and its zeros:
#
#   D, n               features and samples.
#   features, samples  their names, or NULL.
#   sample_zeros       the number of zeros in each sample.
#   feature_nonzeros   the number of samples where each feature is non-zero.
#
# Every method refuses negative, non-finite and NA values, and samples with
# no non-zero value. The options are for methods that need more:
#
#   empty_features_ok  FALSE: refuse features that are zero in every sample
#                      (nothing can be estimated for them).
#   whole              TRUE: warn if the values aren't whole numbers (the
#                      method uses each sample's depth).
#   mostly_zero        TRUE: warn when any sample or feature is at least 80%
#                      zeros.
#   min_samples        2: refuse a single sample (the method estimates each
#                      sample's prior from the others).
zero_replacement_input <- function (
    empty_features_ok = TRUE, whole = FALSE, mostly_zero = FALSE,
    min_samples = 1L, env = parent.frame() ) {

  validate_counts(env)
  validate_margin(env)

  counts <- env$counts
  margin <- env$margin

  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  pkg  <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  if (whole && any(vals != round(vals)))
    warning(
      '`counts` has values that are not whole numbers. This method uses each ',
      "sample's total as its number of reads, so give it counts, not ",
      'proportions.', call. = FALSE )
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])
  samples  <- if (length(dn[[margin]]))      as.character(dn[[margin]])

  if (n < min_samples)
    stop(
      'This method needs at least ', min_samples, ' samples: the prior of ',
      'each sample is estimated from the others. Use `czm` for a single ',
      'sample.', call. = FALSE )

  zc <- .Call(C_zero_counts, counts, margin)

  which_names <- function (i, nms) name_list(if (is.null(nms)) i else nms[i])

  empty <- which(zc$sample_zeros == D)
  if (length(empty))
    stop(
      sprintf(
        '`counts` has %d sample%s with no non-zero values (%s): there is nothing to scale. Remove %s first.',
        length(empty), if (length(empty) == 1) '' else 's',
        which_names(empty, samples), if (length(empty) == 1) 'it' else 'them' ),
      call. = FALSE )

  if (!empty_features_ok) {
    empty <- which(zc$feature_nonzeros == 0L)
    if (length(empty))
      stop(
        sprintf(
          '`counts` has %d feature%s that %s zero in every sample (%s): nothing can be estimated for %s. Remove %s first.',
          length(empty), if (length(empty) == 1) '' else 's',
          if (length(empty) == 1) 'is' else 'are',
          which_names(empty, features),
          if (length(empty) == 1) 'it' else 'them',
          if (length(empty) == 1) 'it' else 'them' ),
        call. = FALSE )
  }

  if (mostly_zero) {
    s <- which(zc$sample_zeros >= 0.8 * D)
    f <- which(n - zc$feature_nonzeros >= 0.8 * n)
    if (length(s) || length(f))
      warning(
        paste0(
          if (length(s))
            sprintf('%d sample%s (%s)', length(s), if (length(s) == 1) '' else 's',
                    which_names(s, samples)),
          if (length(s) && length(f)) ' and ',
          if (length(f))
            sprintf('%d feature%s (%s)', length(f), if (length(f) == 1) '' else 's',
                    which_names(f, features)),
          if (length(s) + length(f) == 1) ' is' else ' are',
          ' at least 80% zeros, so most of their values will be imputed. ',
          'Consider filtering them out first.' ),
        call. = FALSE )
  }

  if (as.double(n) * D > 1e8)
    warning(
      sprintf('The result is a dense %s x %s matrix (%.1f GB).',
        n, D, as.double(n) * D * 8 / 1e9), call. = FALSE )

  list(
    D                = D,
    n                = n,
    features         = features,
    samples          = samples,
    sample_zeros     = zc$sample_zeros,
    feature_nonzeros = zc$feature_nonzeros )
}



# The result of a zero replacement: `res` (the dense matrix from C, in the
# input's orientation) with the input's names, and the attributes `method`,
# then the method's parameters (`params`, a named list), then `n_replaced`.
zero_replacement_result <- function (res, info, margin, method, params) {

  dn <- if (margin == 1L) list(info$samples, info$features)
        else              list(info$features, info$samples)
  if (!all(vapply(dn, is.null, logical(1)))) dimnames(res) <- dn

  attr(res, 'method') <- method
  for (p in names(params)) attr(res, p) <- params[[p]]
  attr(res, 'n_replaced') <- sum(as.double(info$sample_zeros))
  res
}
