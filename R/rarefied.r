# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The averaging options (mean, median, medoid) and the medoid's definition
# (the iteration with the smallest sum of Euclidean distances between
# condensed distance matrices) follow q2-boots (Raspet et al. 2025; BSD-3,
# (c) 2024 Caporaso Lab). No code is ported. See plans/rarefied_diversity.md.



#' Diversity Averaged over Repeated Rarefactions
#'
#' Rarefy a count table many times, compute alpha or beta diversity from
#' each rarefied table, and average the results. A single rarefaction throws
#' reads away at random, so its diversity values carry sampling noise;
#' averaging many keeps rarefaction's control for uneven sequencing depth
#' and removes most of that noise. This is mothur's `dist.shared(subsample =
#' T)`, vegan's `avgdist()`, and QIIME 2's q2-boots.
#'
#' * `beta_div_rarefied()`: a distance matrix, as from [beta_div()].
#' * `alpha_div_rarefied()`: one value per sample, as from [alpha_div()].
#'
#' Iteration \eqn{t} (\eqn{t = 1, \ldots,} `times`) rarefies `counts` with
#' [rarefy()], with seed `seed + t - 1`, and computes the metric from it,
#' exactly as `beta_div(rarefy(counts, depth, seed = seed + t - 1), metric)`
#' would. So `times = 1` gives the same result as that call.
#'
#' @inheritParams beta_div
#' @inheritParams alpha_div
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features) of whole
#'        numbers; see **Input Types**. A `phyloseq`, `rbiom` or
#'        `TreeSummarizedExperiment` object's tree is used for the
#'        phylogenetic metrics.
#'
#' @param metric  The name of a diversity metric: for `beta_div_rarefied()`,
#'        one of [beta_div()]'s metrics, such as `'bray'`; for
#'        `alpha_div_rarefied()`, one of [alpha_div()]'s, such as
#'        `'shannon'`. Matched as those functions match it (case-insensitive,
#'        with prefixes and aliases). See [list_metrics()].
#'
#' @param depth  The number of observations to rarefy every sample to: a
#'        single whole number of at least 1. Samples with fewer observations
#'        are dropped (with a warning), and are listed in the result's
#'        `dropped` attribute. `NULL` (the default) chooses the depth as
#'        [rarefy()] does (see "Choosing the depth").
#'
#' @param times  The number of rarefactions to average, a whole number of
#'        at least 1. Default: `100`, as `vegan::avgdist()`.
#'
#' @param average  How the rarefactions' values are combined: `"mean"`
#'        (the default), `"median"`, or for `beta_div_rarefied()` only,
#'        `"medoid"`. Case-insensitive, and a prefix is enough. See "Choosing
#'        how to average".
#'
#' @param sd  Also return the standard deviation of each value over the
#'        `times` rarefactions, as attribute `sd`. Default: `FALSE`.
#'
#' @param pairs  Which pairs of samples to compute distances for, as in
#'        [beta_div()]: indexes into the distance matrix of the samples that
#'        are **kept** (not dropped for having fewer than `depth`
#'        observations), as with `beta_div(rarefy(counts, depth), metric,
#'        pairs = pairs)`. Default: `NULL` (all pairs).
#'
#' @param seed  An integer seed for the rarefactions: iteration \eqn{t} uses
#'        `seed + t - 1`, as `rarefy(times = )` does. The same seed gives the
#'        same result, for any `cpus`. Default: `0`.
#'
#' @param margin  The margin containing samples: `1` if samples are rows,
#'        `2` if samples are columns. Ignored for `phyloseq`, `rbiom` and
#'        `SummarizedExperiment` objects. Default: `1`.
#'
#' @return
#' `beta_div_rarefied()`: a distance matrix of class `c("ecodive_dist",
#' "dist")` over the kept samples, as [beta_div()] returns, with attribute
#' `method` (the metric's ID). Its title (in its `display` attribute; see
#' [ecodive_display]) names the metric and the rarefactions, such as
#' `"Bray-Curtis Dissimilarity (rarefied 100×)"`, so ordinations and tests
#' of it are titled after it. Further attributes:
#' * `depth`: the depth every sample was rarefied to.
#' * `dropped`: the names (or, without names, the positions) of the samples
#'   dropped for having fewer than `depth` observations; empty if none.
#' * `sd` (with `sd = TRUE`): the standard deviations, a distance matrix
#'   over the same samples, in the same order.
#' * `medoid` (with `average = "medoid"`): the iteration \eqn{t} whose
#'   distance matrix was returned.
#' * `order` (with `clust`): see [beta_div()].
#'
#' `alpha_div_rarefied()`: a named numeric vector over the kept samples, as
#' [alpha_div()] returns, with attributes `depth`, `dropped`, and with
#' `sd = TRUE`, `sd` (a named numeric vector).
#'
#' @section When to use it:
#' Rarefying controls for uneven sequencing depth: every sample is compared
#' at the same number of reads, so that differences in depth don't pass for
#' differences in diversity. McMurdie & Holmes (2014) called rarefying
#' "inadmissible" because a single rarefaction discards data and adds
#' noise; averaging many rarefactions answers the second objection. Schloss
#' (2024) found rarefaction better than the alternatives at controlling for
#' uneven sequencing effort in amplicon data, and Schloss (2026) found it
#' better than robust Aitchison PCA ([rpca()]) and other compositional
#' methods. To compare both views on your own data, ordinate both and test
#' their agreement (see the examples).
#'
#' @section Choosing the depth:
#' Every sample with fewer than `depth` observations is dropped, the same
#' samples from every iteration. A higher depth keeps more of each sample's
#' information; a lower one keeps more samples. With `depth = NULL`, the
#' depth is the smallest library size at which the observations kept (that
#' depth times the number of samples with at least that many) make up at
#' least 10% of all the observations in the table, as in [rarefy()]. Check
#' `attr(x, "dropped")`, and choose `depth` yourself when the samples
#' dropped matter.
#'
#' @section Choosing how to average:
#' * `"mean"` (the default, as in `vegan::avgdist()` and mothur): the
#'   average of each value over the rarefactions. It needs one pass over
#'   them, and memory for a few distance matrices.
#' * `"median"`: the median of each value; less affected by an unusual
#'   rarefaction.
#' * `"medoid"` (beta only): the distance matrix of the one rarefaction
#'   closest to all the others (smallest sum of Euclidean distances between
#'   the iterations' distances). Unlike the mean and the median, the result
#'   is the distance matrix of one real rarefied table, so it keeps
#'   whatever properties the metric has (such as being Euclidean), which a
#'   cell-by-cell average can lose.
#'
#' q2-boots' defaults are `average = "medoid"` for beta diversity and
#' `average = "median"` for alpha diversity; set them to get the same.
#'
#' **Memory and time.** The mean rarefies `times` times. The median and the
#' medoid need every rarefaction's value of each pair of samples at once,
#' so for many samples they work through the pairs in blocks that fit in
#' 256 MB (`8 * times` bytes per pair), rarefying the table `times` times for
#' each block, with the same seeds and so the same draws. With
#' `times = 100`, 500 samples (124,750 pairs) fit in one block, 1,000
#' samples need two, and 10,000 samples about 150. The medoid also needs a
#' `times` \eqn{\times} `times` matrix, and rarefies once more at the end,
#' for the medoid's matrix.
#'
#' @section Is an ordination driven by depth?:
#' [envfit()] fits library size onto any ordination, with a permutation
#' p-value: `envfit(ord, data.frame(log_depth = log(rowSums(counts))))`.
#' A large \eqn{r^2} means that the ordination still follows sequencing
#' depth. `cor(ord, log(rowSums(counts)))` gives the correlation with each
#' axis. See the examples.
#'
#' @section Warnings:
#' A warning raised by the metric (such as the Aitchison distance's choice
#' of pseudocount) is given once, after all the rarefactions, however many
#' of them raised it.
#'
#' @section Same names in other packages:
#' None, but `vegan::avgdist()` computes the same as `beta_div_rarefied()`
#' with `average = "mean"` (or `"median"`), from R's random number
#' generator, so their results agree only up to Monte Carlo error.
#'
#' @section Differences from q2-boots:
#' * No bootstrap mode: [rarefy()] draws only without replacement, so there
#'   is no equivalent of q2-boots' `replacement = True`.
#' * The default `average` is the mean, for alpha and beta diversity alike
#'   (q2-boots: the median for alpha, the medoid for beta).
#' * The rarefactions are seeded with ecodive's own PCG generator, so the
#'   draws differ from q2-boots' (and from vegan's) for the same seed.
#'
#' @references
#' Raspet, I., Gehret, E., Herman, C., Meilander, J., Manley, A., Simard, A.,
#' Bolyen, E., & Caporaso, J. G. (2025). Facilitating bootstrapped and
#' rarefaction-based microbiome diversity analysis with q2-boots.
#' *F1000Research*, 14, 87. \doi{10.12688/f1000research.156295.1}
#'
#' Schloss, P. D. (2024). Rarefaction is currently the best approach to
#' control for uneven sequencing effort in amplicon sequence analyses.
#' *mSphere*, 9, e00354-23. \doi{10.1128/msphere.00354-23}
#'
#' Schloss, P. D. (2026). Rarefaction is better than robust Aitchison PCA and
#' other compositional data analysis methods at controlling for uneven
#' sequencing effort. *mSystems*, 11, e00930-26.
#' \doi{10.1128/msystems.00930-26}
#'
#' McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why rarefying
#' microbiome data is inadmissible. *PLoS Computational Biology*, 10,
#' e1003531. \doi{10.1371/journal.pcbi.1003531}
#'
#' @seealso [rarefy()], [beta_div()], [alpha_div()], [rpca()], and
#'   [rarefaction_curve()] to choose a depth.
#'
#' @name rarefied_diversity
#' @export
#' @examples
#'     # Bray-Curtis dissimilarity, averaged over 100 rarefactions to the
#'     # automatically chosen depth.
#'     dm <- beta_div_rarefied(ex_counts, 'bray')
#'     dm
#'     attr(dm, 'depth')
#'
#'     # With the spread over the rarefactions.
#'     dm <- beta_div_rarefied(ex_counts, 'bray', depth = 100, sd = TRUE)
#'     attr(dm, 'sd')
#'
#'     # q2-boots' defaults: the medoid for beta, the median for alpha.
#'     beta_div_rarefied(ex_counts, 'bray', depth = 100, average = 'medoid')
#'     alpha_div_rarefied(ex_counts, 'shannon', depth = 100, average = 'median')
#'
#'     # Ordinations and tests are titled after the rarefied matrix.
#'     ord <- pcoa(dm, n_axes = 2)
#'     ord
#'
#'     # Does the ordination still follow sequencing depth?
#'     kept <- attr(dm, 'Labels')
#'     envfit(ord, data.frame(log_depth = log(rowSums(ex_counts)[kept])))
#'
beta_div_rarefied <- function (
    counts,
    metric,
    depth       = NULL,
    times       = 100L,
    average     = 'mean',
    sd          = FALSE,
    margin      = 1L,
    norm        = 'none',
    pseudocount = NULL,
    power       = 1.5,
    alpha       = 0.5,
    tree        = NULL,
    pairs       = NULL,
    clust       = NULL,
    seed        = 0L,
    cpus        = n_cpus() ) {

  rarefied_beta_fit(
    counts = counts, metric = metric, depth = depth, times = times,
    average = average, sd = sd, margin = margin, norm = norm,
    pseudocount = pseudocount, power = power, alpha = alpha, tree = tree,
    pairs = pairs, clust = clust, seed = seed, cpus = cpus )
}



#' @rdname rarefied_diversity
#' @export
alpha_div_rarefied <- function (
    counts,
    metric,
    depth   = NULL,
    times   = 100L,
    average = 'mean',
    sd      = FALSE,
    cutoff  = 10L,
    digits  = 3L,
    tree    = NULL,
    margin  = 1L,
    q       = 2,
    bias_corrected = TRUE,
    seed    = 0L,
    cpus    = n_cpus() ) {

  s <- rarefied_setup('alpha')

  validate_cutoff()
  validate_digits()
  validate_q()
  validate_bias_corrected()

  args <- list(
    cutoff = cutoff, digits = digits, tree = s$tree, margin = s$margin,
    q = q, bias_corrected = bias_corrected, cpus = s$cpus )
  args <- args[intersect(names(args), s$metric$params)]

  times <- s$times
  warns <- new.env(parent = emptyenv())
  M     <- NULL

  rarefied_quietly(warns, {
    for (t in seq_len(times)) {
      a <- do.call(s$metric$func, c(list(counts = rarefied_draw(s, t)), args))
      if (t == 1L) {
        nms <- names(a)
        M   <- matrix(NA_real_, length(a), times)
        acc <- as.vector(a)
      }
      else {
        acc <- acc + as.vector(a)
      }
      M[, t] <- a
    }
  })

  val <- switch(
    EXPR   = s$average,
    mean   = acc / times,
    median = .Call(C_row_medians, M, s$cpus) )
  names(val) <- nms

  attr(val, 'depth')   <- s$depth
  attr(val, 'dropped') <- s$dropped
  if (s$sd) {
    sds <- apply(M, 1L, stats::sd)
    names(sds) <- nms
    attr(val, 'sd') <- sds
  }

  rarefied_rewarn(warns)
  val
}



# beta_div_rarefied(), with `budget`: the most memory, in bytes, that the
# median's and medoid's blocks of pairs (8 * times bytes each) may take.
# Tests set a tiny budget to force many blocks.
rarefied_beta_fit <- function (
    counts, metric, depth, times, average, sd, margin, norm, pseudocount,
    power, alpha, tree, pairs, clust, seed, cpus, budget = 256 * 2^20 ) {

  validate_clust()
  s <- rarefied_setup('beta')

  times  <- s$times
  n_keep <- sum(s$keep)
  n_dist <- n_keep * (n_keep - 1) / 2

  # `pairs` indexes the distance matrix of the kept samples.
  if (!is.null(pairs)) {
    env <- new.env()
    env$counts <- matrix(0, n_keep, 1L)
    env$margin <- 1L
    env$pairs  <- pairs
    validate_pairs(env)
    pairs <- env$pairs
  }

  # The arguments the metric takes, besides counts and pairs.
  args <- list(
    margin = s$margin, norm = norm, pseudocount = pseudocount,
    power = power, alpha = alpha, tree = s$tree, clust = NULL,
    cpus = s$cpus )
  args <- args[intersect(names(args), s$metric$params)]

  # The distance matrix of iteration t, for `prs` (NULL: all pairs).
  dist_t <- function (t, prs) {
    a <- args
    a['pairs']  <- list(prs)
    a$counts    <- rarefied_draw(s, t)
    do.call(s$metric$func, a[intersect(s$metric$params, names(a))])
  }

  warns  <- new.env(parent = emptyenv())
  medoid <- NULL

  rarefied_quietly(warns, {

    if (s$average == 'mean') {

      # One pass: a running sum (and Welford's mean and M2, for sd).
      for (t in seq_len(times)) {
        d <- dist_t(t, pairs)
        v <- as.vector(d)
        if (t == 1L) {
          template <- d
          acc <- v
          if (s$sd) { w_mean <- v; w_m2 <- numeric(length(v)) }
        }
        else {
          acc <- acc + v
          if (s$sd) {
            delta  <- v - w_mean
            w_mean <- w_mean + delta / t
            w_m2   <- w_m2 + delta * (v - w_mean)
          }
        }
      }
      val <- acc / times
      if (s$sd)
        sds <- if (times > 1L) sqrt(w_m2 / (times - 1L)) else rep(NA_real_, length(val))
    }

    else {

      # Median and medoid: every iteration's value of a pair at once, in
      # blocks of pairs, rerunning the (identical) rarefactions per block.
      todo <- if (is.null(pairs)) seq_len(n_dist) else sort(unique(pairs))
      size <- max(1, floor(budget / (8 * times)))
      blocks <- if (length(todo)) split(todo, ceiling(seq_along(todo) / size)) else list(integer(0))

      val <- rep(NA_real_, n_dist)
      sds <- rep(NA_real_, n_dist)
      ssq <- matrix(0, times, times)

      for (b in seq_along(blocks)) {

        blk <- blocks[[b]]
        M   <- matrix(NA_real_, length(blk), times)
        for (t in seq_len(times)) {
          d <- dist_t(t, blk)
          if (b == 1L && t == 1L) template <- d
          M[, t] <- d[blk]
        }

        if (s$average == 'median') val[blk] <- .Call(C_row_medians, M, s$cpus)
        else                       ssq      <- .Call(C_medoid_accum, M, ssq, s$cpus)

        if (s$sd && times > 1L) {
          mu <- rowMeans(M)
          s2 <- numeric(length(blk))
          for (t in seq_len(times)) s2 <- s2 + (M[, t] - mu)^2
          sds[blk] <- sqrt(s2 / (times - 1L))
        }
        remove('M')
      }

      if (s$average == 'medoid') {
        # The iteration with the smallest sum of Euclidean distances to the
        # others (the first, if tied), and its whole distance matrix.
        medoid   <- which.min(rowSums(sqrt(ssq)))
        template <- dist_t(medoid, pairs)
        val      <- as.vector(template)
      }
    }
  })


  # The result: the template's labels and attributes, with the new values.
  out <- val
  attributes(out) <- attributes(template)
  title <- paste0(
    METRICS$name[match(s$metric$id, METRICS$id)],
    ' (rarefied ', format(times, scientific = FALSE), '\u00d7)' )

  out <- dist_clust(out, clust, s$cpus)

  attr(out, 'depth')   <- s$depth
  attr(out, 'dropped') <- s$dropped

  if (s$sd) {
    sdd <- sds
    attributes(sdd) <- attributes(template)
    attr(sdd, 'method') <- NULL
    ord <- attr(out, 'order', exact = TRUE)
    if (!is.null(ord)) {
      sdd <- dist_reorder(sdd, ord)
      attr(sdd, 'order') <- ord
    }
    sdd <- set_display(sdd, title = title_of('SD', title), clust = get_display(out, 'clust'))
    attr(out, 'sd') <- sdd
  }

  if (!is.null(medoid)) attr(out, 'medoid') <- medoid

  out <- set_display(out, title = title)

  rarefied_rewarn(warns)
  out
}



# Checks the arguments of beta_div_rarefied() or alpha_div_rarefied() (`div`
# is 'beta' or 'alpha') in the calling function's frame, and returns what
# the iterations need: the counts (as doubles, as rarefy() makes them),
# margin, tree, depth, times, average, sd, seeds, cpus, the matched metric,
# which samples are kept, and the names (or positions) of those dropped.
rarefied_setup <- function (div, env = parent.frame()) {

  # The metric first, so that a misspelt metric is the error reported.
  metric <- match_metric(env$metric, div = div)

  counts <- rarefied_counts(env)
  margin <- env$margin
  times  <- rarefied_times(env$times)


  # depth: one for every sample.
  depth <- env$depth
  if (!is.null(depth)) {
    if (is.numeric(depth) && length(depth) > 1)
      stop(
        '`depth` must be a single whole number: one depth for every sample. ',
        'To rarefy each sample to a depth of its own, use rarefy().', call. = FALSE )
    if (!is.numeric(depth) || !is.null(dim(depth)) || length(depth) != 1 ||
        is.na(depth) || depth < 1 || depth %% 1 != 0 || depth > .Machine$integer.max)
      stop('`depth` must be NULL or a single whole number of at least 1.', call. = FALSE)
    depth <- as.integer(depth)
  }


  average <- rarefied_average(env$average, div)


  # sd
  sd <- env$sd
  if (!(identical(sd, TRUE) || identical(sd, FALSE)))
    stop('`sd` must be TRUE or FALSE.', call. = FALSE)


  # The samples kept: the same in every iteration, as rarefy() drops them.
  sums <- sample_sums(counts, margin)
  if (is.null(depth)) depth <- auto_depth(sums)
  keep <- sums >= depth

  info    <- count_samples(counts, margin)
  dropped <- if (is.null(info$labels)) which(!keep) else info$labels[!keep]
  n_keep  <- sum(keep)
  min_n   <- if (div == 'beta') 2L else 1L

  if (n_keep < min_n)
    stop(sprintf(
      'Only %d of the %d samples %s at least %d observations; %s needed.',
      n_keep, length(keep), if (n_keep == 1) 'has' else 'have', depth,
      if (div == 'beta') 'at least two are' else 'at least one is' ), call. = FALSE)

  if (n_keep < length(keep))
    warning(sprintf(
      '%d of the %d samples %s fewer than %d observations and %s dropped.',
      sum(!keep), length(keep), if (sum(!keep) == 1) 'has' else 'have', depth,
      if (sum(!keep) == 1) 'was' else 'were' ), call. = FALSE)


  list(
    counts  = counts,
    margin  = margin,
    tree    = env$tree,
    depth   = depth,
    times   = times,
    average = average,
    sd      = sd,
    seeds   = rarefied_seeds(env$seed, times),
    cpus    = env$cpus,
    metric  = metric,
    keep    = keep,
    dropped = dropped )
}



# The checks of `counts`, `margin`, `seed` and `cpus` in `env` that the
# rarefied functions and rarefaction_curve() share. Returns the counts as
# doubles, as rarefy() makes them, so that each draw is exactly rarefy()'s.
rarefied_counts <- function (env) {

  validate_counts(env)
  validate_margin(env)
  assert_integer_counts(env)
  validate_seed(env)
  validate_cpus(env)

  counts <- env$counts

  if (isTRUE(any(stored_values(counts) < 0)))
    stop('`counts` must not contain negative values.', call. = FALSE)

  if (is.matrix(counts)) {
    if (!is.double(counts)) counts[] <- as.double(counts)
  }
  else if (inherits(counts, 'simple_triplet_matrix')) {
    if (!is.double(counts$v)) counts$v <- as.double(counts$v)
  }

  counts
}


# `times`, checked, as an integer.
rarefied_times <- function (times) {
  if (!is.numeric(times) || length(times) != 1 || is.na(times) ||
      times < 1 || times %% 1 != 0 || times > .Machine$integer.max)
    stop('`times` must be a whole number of at least 1.', call. = FALSE)
  as.integer(times)
}


# `average`, matched (exactly, else by a unique prefix, ignoring case) to
# the choices for `div` ('beta' or 'alpha').
rarefied_average <- function (average, div) {
  choices <- if (div == 'beta') c('mean', 'median', 'medoid') else c('mean', 'median')
  if (!is.character(average) || length(average) != 1 || is.na(average) || !nzchar(average))
    stop('`average` must be one of ', paste0('"', choices, '"', collapse = ', '), '.', call. = FALSE)
  needle  <- tolower(average)
  hit     <- if (needle %in% choices) needle else choices[startsWith(choices, needle)]
  if (length(hit) != 1) {
    if (div == 'alpha' && nchar(needle) >= 4 && startsWith('medoid', needle))
      stop(
        '`average = "medoid"` is only for beta diversity; ',
        'use "mean" or "median" for alpha diversity.', call. = FALSE )
    stop(
      '`average = "', average, '"` must be one of ',
      paste0('"', choices, '"', collapse = ', '),
      if (length(hit) > 1) ' (it matches several)', '.', call. = FALSE )
  }
  hit
}


# The seeds of `times` rarefactions, as rarefy(times =) makes them: `seed`,
# `seed + 1`, ..., wrapped to R's integer range.
rarefied_seeds <- function (seed, times) {
  ((seed + 2**31 - 1 + seq_len(times)) %% 2**32) - 2**31
}



# The rarefied table of iteration t: exactly rarefy(counts, depth, seed =
# seeds[t], margin = margin, cpus = cpus, warn = FALSE).
rarefied_draw <- function (s, t) {
  m <- .Call(C_rarefy, s$counts, s$depth, s$seeds[[t]], s$margin, s$cpus)
  if (s$margin == 1L) m[s$keep, , drop = FALSE]
  else                m[, s$keep, drop = FALSE]
}



# Evaluates `expr`, holding back its warnings: each distinct message is
# recorded once in `store`, for rarefied_rewarn() to give after the loop.
rarefied_quietly <- function (store, expr) {
  withCallingHandlers(
    expr,
    warning = function (w) {
      store$msgs <- unique(c(store$msgs, conditionMessage(w)))
      invokeRestart('muffleWarning')
    })
}

rarefied_rewarn <- function (store) {
  for (msg in store$msgs) warning(msg, call. = FALSE)
}
