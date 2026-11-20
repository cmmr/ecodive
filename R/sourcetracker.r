# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# sourcetracker() follows SourceTracker2 2.0.1-dev (BSD-3; Copyright (c)
# 2016, Biota Technology): gibbs_helper() in sourcetracker/_gibbs.py, and
# collapse_source_data(), gibbs() and collate_gibbs_results() in
# sourcetracker/_sourcetracker.py. Its sampler is ported to C in
# src/sourcetracker.c, whose header holds that license in full. See
# inst/COPYRIGHTS.



# The source-tracking methods: their IDs (the results' `method`) and display
# names, for plot titles and print() headers. feast() joins it later.
SOURCE_TRACKERS <- data.frame(
  id   = 'sourcetracker',
  name = 'SourceTracker',
  stringsAsFactors = FALSE )

source_tracker_name <- function (method) {
  i <- match(method, SOURCE_TRACKERS$id)
  if (length(i) != 1 || is.na(i)) NULL else SOURCE_TRACKERS$name[[i]]
}



#' Microbial source tracking (SourceTracker)
#'
#' Estimates what share of each **sink** sample's community came from each
#' of several **source** environments, and from an **unknown** source, with
#' SourceTracker's Bayesian Gibbs sampler (Knights et al. 2011), as
#' SourceTracker2 implements it.
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of whole-number counts (samples
#'        \eqn{\times} features), sources and sinks together, or any table
#'        type listed under "Input Types".
#'
#' @param sources  Each sample's source environment: a vector with one
#'        value per sample (matched by name when both are named), or the
#'        name of a column of `data`. **`NA` marks a sink.** Samples with the
#'        same label are pooled into one environment. To leave a sample out
#'        altogether, remove it from `counts`.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample (matched to the samples by row name when both have names),
#'        for `sources` to name a column of.
#'
#' @param depth  The rarefaction depth of the pooled source environments, or
#'        `NULL` not to rarefy them. Default: `1000`.
#'
#' @param sink_depth  The rarefaction depth of the sinks, or `NULL` not to
#'        rarefy them. Default: `depth`.
#'
#' @param alpha1  The prior count added to every feature of every known
#'        source environment. Higher values trust the sources less and smooth
#'        their feature distributions. Default: `0.001`.
#'
#' @param alpha2  The prior count of every feature in the unknown source, as
#'        a fraction of the sink's depth. Higher values make the unknown
#'        source smoother. Must be positive. Default: `0.1`.
#'
#' @param beta  The prior count of the sink's reads in each source (the
#'        unknown included). Must be positive. Default: `10`.
#'
#' @param n_restarts  The number of independent Gibbs chains for each sink.
#'        Default: `10`.
#'
#' @param n_draws  The number of draws taken from each chain. Default: `1`.
#'
#' @param n_burnin  The number of passes over the sink's reads before the
#'        first draw. SourceTracker calls this `burnin`. Default: `100`.
#'
#' @param thin  The number of passes between draws. SourceTracker calls this
#'        `delay`. Default: `1`, as in SourceTracker2 (SourceTracker 1's
#'        default is 10).
#'
#' @param loo  Leave one out: instead of the sinks, track each pooled and
#'        rarefied source environment in turn against the others. Samples
#'        with `NA` in `sources` are then ignored. Default: `FALSE`.
#'
#' @param features  Also return each sink's reads by source and feature
#'        (attribute `features`). Default: `FALSE`.
#'
#' @param seed  Random seed for the rarefaction and the Gibbs chains: any
#'        integer. The same `seed` gives the same result for any `cpus`.
#'        Default: `0`.
#'
#' @return A numeric matrix of class `ecodive_source`: one row per sink, and
#'   one column per source environment (in the order of `levels(sources)`
#'   for a factor, otherwise sorted by name in the C locale), then
#'   `"Unknown"`. Each row holds the sink's mean source proportions, and
#'   sums to 1. With `loo = TRUE`, there is one row per environment, and 0
#'   in each environment's own column.
#'
#'   Attributes:
#'   \describe{
#'     \item{`sd`}{The standard deviations of the proportions over the
#'       draws, in the same shape. With the defaults they rest on only 10
#'       draws (`n_restarts` \eqn{\times} `n_draws`).}
#'     \item{`n_draws`}{The number of draws in total for each sink:
#'       `n_restarts * n_draws`.}
#'     \item{`features`}{With `features = TRUE`: a list with one matrix per
#'       sink, of the sink's reads assigned to each source (rows, as the
#'       columns of the result) from each feature (columns), averaged over
#'       the draws. Each matrix sums to the sink's (rarefied) depth.}
#'     \item{`method`}{`"sourcetracker"`.}
#'   }
#'
#'   It is still a matrix: `x[, "Unknown"]`, `rowSums(x)` and
#'   `barplot(t(x))` work as for any other. [plot()][plot.ecodive_source]
#'   draws the usual stacked bar chart.
#'
#' @details
#' **Sinks, sources and the unknown source.** A *sink* is a community whose
#' origins are in question (say an infant's gut, or a hospital surface);
#' the *sources* are samples of the environments it may have come from
#' (the mother's skin, vagina and stool; patients' skin; tap water). Each
#' environment's samples are pooled, and each read of the sink is assigned,
#' by the Gibbs sampler, to one of the environments or to an *unknown*
#' source, whose composition is learned from the reads assigned to it. The
#' proportions are the shares of the sink's reads assigned to each, averaged
#' over the draws.
#'
#' The model assumes that a sink is a mixture of its sources (and the
#' unknown). It is not when communities change after they arrive (growth,
#' competition, a different habitat), which can make sources look more or
#' less important than they were (Shenhav et al. 2019). In their simulations,
#' SourceTracker also tended to underestimate the unknown source.
#'
#' **The pipeline**, as SourceTracker2's:
#' 1. Each environment's samples are pooled by their **mean**, truncated to
#'    whole numbers (as SourceTracker2 does, although its README says
#'    "rounded").
#' 2. The pooled environments are rarefied to `depth`, and the sinks to
#'    `sink_depth`, without replacement. A source environment or sink with
#'    fewer reads than its depth is an error. With `loo = TRUE`, the
#'    held-out environment is already rarefied, and `sink_depth` is unused.
#' 3. For each sink, `n_restarts` chains each start from a uniform random
#'    assignment of the reads, run `n_burnin + (n_draws - 1) * thin + 1`
#'    passes over the reads (each in a new random order), and take a draw
#'    (the reads' current assignment) after pass `n_burnin + 1` and every
#'    `thin` passes after that.
#'
#' A read of feature \eqn{t} is assigned to known source \eqn{v} with
#' weight \eqn{(m_{vt} + \alpha_1)/(m_v + \tau\alpha_1)\,(c_v + \beta)},
#' and to the unknown source with weight \eqn{(u_t + \alpha_2 n)/(U +
#' \alpha_2 n\tau)\,(c_U + \beta)}, where \eqn{m} are the rarefied source
#' counts, \eqn{c} the sink's reads currently assigned to each source,
#' \eqn{u} (sum \eqn{U}) those assigned to the unknown by feature, \eqn{n}
#' the sink's depth, and \eqn{\tau} the number of features: **every column
#' of `counts`**, including features that are zero everywhere, as in
#' SourceTracker2. Removing empty features therefore changes the result a
#' little.
#'
#' **Rarefaction depth matters.** Rarefying the sources to `depth` sets how
#' much evidence each environment's feature distribution carries; rarefying
#' the sinks to `sink_depth` sets how many reads are assigned, and so how
#' smooth the unknown source is (`alpha2` is per read of the sink). Results
#' at different depths aren't directly comparable, and a deeper `depth`
#' makes rare features of the sources count.
#'
#' **Reproducibility.** Rarefying sink \eqn{i} uses the same random stream
#' as [rarefy()] would for sample \eqn{i} of `counts`, so that each sink
#' is rarefied exactly as `rarefy(counts, sink_depth, seed = seed)`
#' rarefies it. Each chain has its own stream, chosen by the sink's
#' position in `counts` and the restart. The result is therefore the same
#' for any `cpus`, and a sink's result doesn't change when sinks after it
#' in `counts` are added or removed.
#'
#' **The stacked bar chart.** `plot(x)` draws one bar per sink, divided by
#' source, with the unknown in grey; `barplot(t(x), legend.text = TRUE)`
#' draws much the same with base graphics alone.
#'
#' @inheritSection documentation Input Types
#'
#' @section Differences from SourceTracker2:
#' * Each sink and restart has its own random stream of `seed`, so results
#'   are reproducible and identical for any `cpus`. SourceTracker2's
#'   `--jobs` workers either share one random state or start unseeded, so
#'   its results depend on `jobs` and the operating system. Its random
#'   numbers (NumPy's, and scikit-bio's subsampling) aren't reproduced, so
#'   results agree with SourceTracker2's in distribution, not draw for draw.
#' * One `sources` vector, with `NA` for sinks, replaces SourceTracker2's
#'   `SourceSink` and `Env` metadata columns.
#' * `sd` is the standard deviation of the per-draw proportions. SourceTracker2
#'   divides each draw's counts by the total over *all* draws before taking
#'   the standard deviation, so its reported values are these divided by
#'   the number of draws.
#' * `features` holds the mean over the draws. SourceTracker2's feature
#'   tables (its docstring says averages) are sums over the draws: multiply
#'   by `attr(x, "n_draws")` to get them.
#' * Counts must be whole numbers; SourceTracker2 truncates others.
#' * `loo = TRUE` works. SourceTracker2's leave-one-out (2.0.1-dev, with
#'   current pandas) fails: it keeps the held-out environment among the
#'   sources. ecodive does what it intends, holding it out.
#'
#' @section Differences from SourceTracker 1:
#' The original R script's `delay` (here `thin`) defaults to 10, not 1.
#' Its leave-one-out holds out individual source samples; SourceTracker2's
#' (and so `loo = TRUE` here) holds out whole pooled environments.
#' SourceTracker 1's `tune.st()` (choosing `alpha1` and `alpha2` by
#' leave-one-out) isn't provided.
#'
#' @section Same names in other packages:
#' SourceTracker 1's R script defines `sourcetracker()`, the same method
#' (it isn't a package on CRAN). FEAST (Shenhav et al. 2019; package
#' 'FEAST' on GitHub) is the other common source-tracking method.
#'
#' @references
#' Knights, D., Kuczynski, J., Charlson, E. S., Zaneveld, J., Mozer, M. C.,
#' Collman, R. G., Bushman, F. D., Knight, R., & Kelley, S. T. (2011).
#' Bayesian community-wide culture-independent microbial source tracking.
#' *Nature Methods*, 8, 761-763. \doi{10.1038/nmeth.1650}
#'
#' Shenhav, L., Thompson, M., Joseph, T. A., Briscoe, L., Furman, O.,
#' Bogumil, D., Mizrahi, I., Pe'er, I., & Halperin, E. (2019). FEAST: fast
#' expectation-maximization for microbial source tracking. *Nature
#' Methods*, 16, 627-632. \doi{10.1038/s41592-019-0431-x}
#'
#' @export
#' @examples
#'     # Two source environments, three samples of each, and three sinks:
#'     # a mix of soil and water, mostly soil, and mostly something else.
#'     soil  <- c(300, 200, 100, 50, 0, 0, 0, 0)
#'     water <- c(0, 0, 50, 100, 200, 300, 0, 0)
#'     other <- c(0, 0, 0, 0, 0, 0, 300, 350)
#'     counts <- rbind(
#'       soil1  = soil + c(10, 0, 5, 0, 0, 0, 0, 0),
#'       soil2  = soil + c(0, 20, 0, 5, 0, 0, 0, 0),
#'       soil3  = soil,
#'       water1 = water + c(0, 0, 0, 0, 10, 0, 0, 0),
#'       water2 = water + c(0, 0, 0, 10, 0, 20, 0, 0),
#'       water3 = water,
#'       mix    = soil + water,
#'       soily  = 3 * soil + water,
#'       novel  = soil + 5 * other )
#'     sources <- c(rep('soil', 3), rep('water', 3), NA, NA, NA)
#'
#'     st <- sourcetracker(counts, sources, depth = 500, n_burnin = 20)
#'     st
#'
#'     plot(st)
#'
#'     # Leave one out: how distinct are the environments?
#'     sourcetracker(counts, sources, depth = 500, n_burnin = 20, loo = TRUE)
sourcetracker <- function (
    counts, sources, data = NULL, margin = 1L,
    depth = 1000L, sink_depth = depth,
    alpha1 = 1e-3, alpha2 = 0.1, beta = 10,
    n_restarts = 10L, n_draws = 1L, n_burnin = 100L, thin = 1L,
    loo = FALSE, features = FALSE, seed = 0L, cpus = n_cpus() ) {

  force(sink_depth)
  if (missing(sources)) sources <- NULL
  sourcetracker_fit(environment())$result
}



# sourcetracker() for the tests: returns list(result, collapsed, sources,
# sinks, draws, passes): the pooled environments before and after
# rarefaction (E x F), the rarefied sinks (sinks x F), the per-draw read
# counts (V x draws x sinks) and each pass's draw flag.
sourcetracker_details <- function (
    counts, sources, data = NULL, margin = 1L,
    depth = 1000L, sink_depth = depth,
    alpha1 = 1e-3, alpha2 = 0.1, beta = 10,
    n_restarts = 10L, n_draws = 1L, n_burnin = 100L, thin = 1L,
    loo = FALSE, features = FALSE, seed = 0L, cpus = n_cpus() ) {

  force(sink_depth)
  if (missing(sources)) sources <- NULL
  sourcetracker_fit(environment(), details = TRUE)
}



# Which passes of a chain are followed by a draw, as SourceTracker2's
# gibbs_sampler(): n_burnin + (n_draws - 1) * thin + 1 passes, and a draw
# after pass p (1-based) when p > n_burnin and (p - n_burnin - 1) %% thin is 0.
st_draw_passes <- function (n_burnin, n_draws, thin) {
  p <- seq_len(n_burnin + (n_draws - 1L) * thin + 1L)
  as.integer(p > n_burnin & (p - n_burnin - 1L) %% thin == 0L)
}



sourcetracker_fit <- function (env, details = FALSE) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)
  validate_seed(env)

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
  if (!all(vals %% 1 == 0))
    stop('`counts` must be whole numbers (integers).', call. = FALSE)
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  n        <- as.integer(dims[[margin]])
  F        <- as.integer(dims[[3L - margin]])
  samples  <- if (length(dn[[margin]]))      as.character(dn[[margin]])
  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])


  #--------------------------------------------------------
  # Settings.
  #--------------------------------------------------------
  one_num <- function (arg, ok, what) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) || !ok(v))
      stop('`', arg, '` must be ', what, '.', call. = FALSE)
    as.double(v)
  }
  one_int <- function (arg, min, max = .Machine$integer.max) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v %% 1 != 0 || v < min || v > max)
      stop(
        '`', arg, '` must be a single integer ',
        if (min == 0) 'of at least 0' else 'greater than 0',
        if (max < .Machine$integer.max) paste(' and at most', max), '.', call. = FALSE )
    as.integer(v)
  }
  one_depth <- function (arg) {
    v <- env[[arg]]
    if (is.null(v)) return (0L)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v %% 1 != 0 ||
        v < 1 || v > .Machine$integer.max)
      stop('`', arg, '` must be a single positive integer, or NULL for no rarefaction.', call. = FALSE)
    as.integer(v)
  }
  one_flag <- function (arg) {
    v <- env[[arg]]
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)
    v
  }

  alpha1     <- one_num('alpha1', function (v) v >= 0, 'a single non-negative number')
  alpha2     <- one_num('alpha2', function (v) v > 0,  'a single positive number')
  beta       <- one_num('beta',   function (v) v > 0,  'a single positive number')
  n_restarts <- one_int('n_restarts', 1, 2^24 - 1)
  n_draws    <- one_int('n_draws', 1)
  n_burnin   <- one_int('n_burnin', 0)
  thin       <- one_int('thin', 1)
  depth      <- one_depth('depth')
  sink_depth <- one_depth('sink_depth')
  loo        <- one_flag('loo')
  features_  <- one_flag('features')
  seed       <- env$seed
  cpus       <- env$cpus

  n_passes <- n_burnin + (n_draws - 1) * thin + 1
  if (n_passes > .Machine$integer.max || (as.double(n_restarts) * n_draws) > .Machine$integer.max)
    stop('`n_burnin`, `n_draws`, `thin` and `n_restarts` give too many passes or draws.', call. = FALSE)
  passes <- st_draw_passes(n_burnin, n_draws, thin)


  #--------------------------------------------------------
  # Sources and sinks.
  #--------------------------------------------------------
  data <- env$data
  idx  <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn_data <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(samples, n, dn_data, nrow(data), ref = '`counts`')
  }
  if (is.null(env$sources))
    stop('`sources` is required: one value per sample (NA for sinks), or a column name of `data`.', call. = FALSE)
  lab <- per_sample(env$sources, 'sources', data, idx, samples, n, ref = '`counts`')

  if (is.factor(lab)) {
    envs <- levels(droplevels(lab[!is.na(lab)]))
    lab  <- as.character(lab)
  }
  else {
    lab  <- as.character(lab)
    envs <- sort(unique(lab[!is.na(lab)]), method = 'radix')
  }
  E <- length(envs)
  if (E == 0L)
    stop('`sources` has no source samples: give each source sample its environment, and NA for sinks.', call. = FALSE)
  if ('Unknown' %in% envs)
    stop('"Unknown" can\'t be a source environment: it names the unknown source.', call. = FALSE)

  env_idx <- match(lab, envs) - 1L
  env_idx[is.na(env_idx)] <- -1L

  sink_rows <- which(is.na(lab))
  if (!loo && length(sink_rows) == 0L)
    stop('`sources` has no sinks (NA values). Use `loo = TRUE` to track each environment against the others.', call. = FALSE)

  sample_names <- if (is.null(samples)) as.character(seq_len(n)) else samples

  n_sinks <- if (loo) E else length(sink_rows)
  if (as.double(n_sinks) * n_restarts * n_draws * (E + 1) > .Machine$integer.max)
    stop('Too many sinks and draws: reduce `n_restarts` or `n_draws`.', call. = FALSE)


  #--------------------------------------------------------
  # Pool the environments, and check the depths.
  #--------------------------------------------------------
  collapsed <- .Call(C_st_collapse, counts, margin, env_idx, E)
  env_sums  <- rowSums(collapsed)

  if (any(env_sums > .Machine$integer.max))
    stop('A source environment has more than 2^31 - 1 reads after pooling.', call. = FALSE)
  if (depth > 0L && any(env_sums < depth)) {
    low <- envs[env_sums < depth]
    stop(sprintf(
      '%d source environment%s %s fewer than `depth` (%d) reads after pooling (the shallowest has %s): %s',
      length(low), if (length(low) == 1) '' else 's', if (length(low) == 1) 'has' else 'have',
      depth, format(min(env_sums)), name_list(low) ), call. = FALSE)
  }
  if (any(env_sums == 0) && (loo || alpha1 == 0))
    stop(
      'Source environment', if (sum(env_sums == 0) == 1) ' ' else 's ',
      name_list(envs[env_sums == 0]), ' ha', if (sum(env_sums == 0) == 1) 's' else 've',
      ' no reads after pooling', if (!loo) ', so `alpha1` must be positive', '.', call. = FALSE )

  if (!loo) {
    sums <- sample_sums(counts, margin)[sink_rows]
    if (any(sums > .Machine$integer.max))
      stop('A sink has more than 2^31 - 1 reads.', call. = FALSE)
    low <- if (sink_depth > 0L) sums < sink_depth else sums == 0
    if (any(low))
      stop(sprintf(
        '%d sink%s %s %s: %s',
        sum(low), if (sum(low) == 1) '' else 's', if (sum(low) == 1) 'has' else 'have',
        if (sink_depth > 0L) sprintf('fewer than `sink_depth` (%d) reads', sink_depth) else 'no reads',
        name_list(sample_names[sink_rows][low]) ), call. = FALSE)
  }


  #--------------------------------------------------------
  # Sample.
  #--------------------------------------------------------
  res <- .Call(
    C_sourcetracker, counts, margin, sink_rows - 1L, collapsed, depth,
    sink_depth, alpha1, alpha2, beta, n_restarts, passes, loo, features_,
    seed, details, cpus )

  draws <- res$draws                 # V x D x sinks
  D     <- dim(draws)[[2]]
  K     <- dim(draws)[[3]]
  V     <- dim(draws)[[1]]
  cols  <- c(envs, 'Unknown')
  rows  <- if (loo) envs else sample_names[sink_rows]

  prop <- matrix(0, K, E + 1L, dimnames = list(rows, cols))
  sds  <- prop
  for (j in seq_len(K)) {
    d   <- matrix(as.double(draws[, , j]), V, D)
    tot <- colSums(d)
    p   <- d / rep(tot, each = V)
    m   <- rowSums(d) / sum(tot)
    s   <- sqrt(rowSums((p - m)^2) / D)
    at  <- if (loo) c(seq_len(E)[-j], E + 1L) else seq_len(E + 1L)
    prop[j, at] <- m
    sds[j, at]  <- s
  }

  out <- prop
  attr(out, 'sd')      <- sds
  attr(out, 'n_draws') <- D

  if (features_) {
    fts <- lapply(seq_len(K), function (j) {
      f  <- res$features[[j]] / D
      at <- if (loo) c(seq_len(E)[-j], E + 1L) else seq_len(E + 1L)
      m  <- matrix(0, E + 1L, F, dimnames = list(cols, features))
      m[at, ] <- f
      m
    })
    names(fts) <- rows
    attr(out, 'features') <- fts
  }

  attr(out, 'method') <- 'sourcetracker'
  class(out) <- c('ecodive_source', 'matrix', 'array')

  if (!details) return (list(result = out))

  dimnames(collapsed)   <- list(envs, features)
  dimnames(res$sources) <- list(envs, features)
  if (!loo) dimnames(res$sinks) <- list(rows, features)
  list(
    result    = out,
    collapsed = collapsed,
    sources   = res$sources,
    sinks     = if (loo) res$sources else res$sinks,
    draws     = draws,
    passes    = passes )
}



#------------------------------------------------------------------------------
# Printing and plotting: shared by every source-tracking method (method
# `sourcetracker`, and later `feast`).
#------------------------------------------------------------------------------

#' @rdname print_ecodive
#' @export
print.ecodive_source <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  att <- attributes(x)
  sym <- print_sym()
  k   <- ncol(x)
  nm  <- source_tracker_name(att$method)

  head <- c(
    sprintf('A matrix: %d sink%s %s %d source%s', nrow(x), if (nrow(x) == 1) '' else 's',
            sym$times, k, if (k == 1) '' else 's'),
    if (!is.null(nm)) sprintf('Source contributions (%s)', nm),
    if (is.numeric(att$n_draws)) sprintf('Mean of %d draw%s', att$n_draws, if (att$n_draws == 1) '' else 's') )

  m <- unclass(x)
  attributes(m) <- att[c('dim', 'dimnames')]
  df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
  if (is.null(colnames(m))) names(df) <- sprintf('[,%d]', seq_len(k))

  cat(
    print_table_lines(
      df, head, n = n, width = width, digits = digits, types = FALSE,
      rows = 'sinks', attrs = print_attrs(x, c('dim', 'dimnames', 'class')) ),
    sep = '\n' )
  invisible(x)
}


#' Plot source contributions
#'
#' The usual stacked bar chart of a source-tracking result: one bar per
#' sink, divided into its estimated proportions from each source, with the
#' unknown source last, in grey.
#'
#' @param x  The result of [sourcetracker()].
#'
#' @param horiz  Draw the bars horizontally. Default: `FALSE`.
#'
#' @param legend  Where to put the legend of sources: a position keyword for
#'        [graphics::legend()] (such as `"topright"`), or `NULL` for none.
#'        Default: `"right"`, outside the plot.
#'
#' @param ...  Passed on to [graphics::barplot()], such as `main`, `col`,
#'        `names.arg`, `las` or `border`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     counts <- rbind(
#'       a1 = c(90, 10, 0, 0), a2 = c(80, 20, 0, 0),
#'       b1 = c(0, 0, 70, 30), b2 = c(0, 5, 60, 35),
#'       s1 = c(45, 10, 30, 15), s2 = c(80, 15, 5, 0) )
#'     st <- sourcetracker(
#'       counts, c('a', 'a', 'b', 'b', NA, NA), depth = NULL, n_burnin = 10 )
#'     plot(st)
plot.ecodive_source <- function (x, horiz = FALSE, legend = 'right', ...) {

  args <- list(...)
  pos  <- plot_legend_pos(legend)
  m    <- unclass(x)
  attributes(m) <- attributes(x)[c('dim', 'dimnames')]
  k    <- ncol(m)
  srcs <- colnames(m)
  if (is.null(srcs)) srcs <- sprintf('Source %d', seq_len(k))
  unk  <- identical(srcs[[k]], 'Unknown')
  col  <- if (unk) c(plot_palette(k - 1L), 'grey80') else plot_palette(k)
  col  <- plot_arg(args, 'col', col)

  nm   <- source_tracker_name(attr(x, 'method', exact = TRUE))
  main <- if (is.null(nm)) 'Source contributions' else sprintf('Source contributions (%s)', nm)

  # Room for the legend outside the plot, on the right.
  outside <- identical(pos, 'right')
  if (outside) {
    w   <- max(graphics::strwidth(srcs, units = 'inches')) + 0.6
    old <- graphics::par(mai = graphics::par('mai') + c(0, 0, 0, w))
    on.exit(graphics::par(old), add = TRUE)
  }

  plot_do(
    graphics::barplot, c(list(t(m)), plot_drop(args, 'col')),
    list(
      col = col, horiz = horiz, las = 1, border = NA,
      main = main,
      xlab = if (horiz) 'Proportion' else '',
      ylab = if (horiz) '' else 'Proportion',
      xlim = if (horiz) c(0, 1),
      ylim = if (!horiz) c(0, 1) ))

  if (!is.null(pos)) {
    if (outside) {
      usr <- graphics::par('usr')
      graphics::legend(
        x = usr[[2]] + 0.02 * (usr[[2]] - usr[[1]]), y = usr[[4]],
        legend = rev(srcs), fill = rev(rep_len(col, k)), border = NA,
        bty = 'n', xpd = TRUE )
    }
    else
      graphics::legend(pos, legend = rev(srcs), fill = rev(rep_len(col, k)), border = NA, bty = 'n')
  }

  invisible(x)
}
