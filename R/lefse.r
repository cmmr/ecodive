# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# lefse() is a port of LEfSe 1.1.2 (lefse/lefse.py, lefse_run.py,
# lefse_format_input.py and, for the plot, lefse_plot_res.py), under this
# license:
#
#   Copyright (c) 2012, Nicola Segata and Curtis Huttenhower
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
#   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
#   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# The statistics LEfSe gets from R's kruskal.test(), coin's wilcox_test()
# and MASS's lda() (GPL; their code was not read) are reproduced from their
# documentation and outputs: see src/lefse.c.



#' LEfSe: LDA Effect Size Biomarker Discovery
#'
#' Finds the features that differ between groups of samples with LEfSe
#' (Segata et al. 2011): a Kruskal-Wallis test of each feature, pairwise
#' Wilcoxon tests between subgroups for those that pass, and a bootstrapped
#' linear discriminant analysis whose effect size, the "LDA score", ranks
#' them. A port of LEfSe 1.1.2, the program papers cite.
#'
#' @inheritParams da_wilcoxon
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**.
#'
#' @param groups  The classes to compare: the name of a column of `data`,
#'        or a vector with one value per sample. Two or more classes. A
#'        factor's levels give their order; other values are sorted.
#'
#' @param subgroups  Subclasses within the classes, for the Wilcoxon step:
#'        the name of a column of `data`, or one value per sample. `NULL`
#'        (the default) makes each class one subgroup.
#'
#' @param total  Each sample is rescaled to sum to `total` before testing:
#'        the default, `1e6`, gives counts per million, as LEfSe's Galaxy
#'        instructions and lefser do (LEfSe's `-o 1000000`), and is the scale
#'        the `lda_threshold` of 2 was set for. `NULL` uses the values as
#'        given (LEfSe's command-line default).
#'
#' @param alpha  The significance level of the Kruskal-Wallis step: a
#'        feature passes when its adjusted p-value is below `alpha` (LEfSe's
#'        `-a`). Default `0.05`.
#'
#' @param wilcoxon  Whether to run the Wilcoxon step (LEfSe's `--wilc`).
#'        Default `TRUE`.
#'
#' @param wilcoxon_alpha  The significance level of the Wilcoxon step
#'        (LEfSe's `-w`): each comparison's two-sided p-value must be below
#'        `2 * wilcoxon_alpha`, as in LEfSe. Default `0.05`.
#'
#' @param lda_threshold  Features are reported (`reject`) when their LDA
#'        score is above this (LEfSe's `-l`). Default `2`.
#'
#' @param n_boots  The number of bootstrap subsamples the LDA is fitted to
#'        (LEfSe's `-b`). Default `30`.
#'
#' @param boot_fraction  The size of each bootstrap subsample, as a fraction
#'        of the samples (LEfSe's `-f`), above 0 and at most 1. Default
#'        `0.67`.
#'
#' @param min_subgroup  Subgroups with fewer samples than this are compared
#'        by their medians alone, without a test (LEfSe's `--min_c`).
#'        Default `10`.
#'
#' @param one_against_one  With three or more classes: `FALSE` (the
#'        default) passes a feature when one class differs from every other
#'        one ("one against all"); `TRUE` needs every pair of classes to
#'        differ (LEfSe's `-y 1`, stricter).
#'
#' @param same_subgroups  With `subgroups`: whether to compare only the
#'        subgroups of the same name in different classes (LEfSe's `-e`).
#'        Default `FALSE`: every subgroup of one class against every
#'        subgroup of the other.
#'
#' @param p_adjust  How to adjust the Kruskal-Wallis p-values for testing
#'        every feature: any method of [stats::p.adjust()]: `'holm'` (the
#'        default, as ecodive's other DA methods), `'BH'` (recommended, as
#'        the benchmarks), or `'none'` (LEfSe, which has no adjustment).
#'
#' @param seed  An integer seed for the jitter and the bootstrap subsamples.
#'        The same `seed` gives the same scores, for any `cpus`.
#'
#' @return A data frame with a row per feature, in the order of `counts`, of
#'   class `c('ecodive_lefse', 'data.frame')`, with columns:
#'
#'   * `feature`: the feature's name, or its index as text.
#'   * `enriched`: the class where its mean is largest (LEfSe's output
#'     names it for the reported features).
#'   * `log_max_mean`: \eqn{\log_{10}} of that mean, at least 0 (LEfSe's
#'     second column).
#'   * `kw_stat`, `p_value`: the Kruskal-Wallis statistic and p-value.
#'   * `p_adjusted`: the p-value adjusted with `p_adjust`.
#'   * `wilcoxon`: whether the feature passed the Wilcoxon step; `NA` when
#'     it failed the Kruskal-Wallis step or `wilcoxon = FALSE`.
#'   * `lda_score`: the LDA score, for the features that passed both steps;
#'     `NaN` for the others.
#'   * `reject`: `TRUE` for the features LEfSe reports: passed both steps,
#'     with `lda_score > lda_threshold`.
#'
#'   Missing values are `NaN` (all of a feature's values tied gives `NaN`
#'   Kruskal-Wallis statistics and p-values, and the feature fails).
#'   Attributes: `method` (`"lefse"`), `groups` (the classes, in order),
#'   `subgroups` (the subgroups' names, or `NULL`), `alpha`,
#'   `wilcoxon_alpha`, `lda_threshold`, `n_boots`, `boot_fraction`, `total`
#'   (`NA` for none), `p_adjust`, `log_base` (`NA`: the score is a
#'   \eqn{\log_{10}} of a difference, not a log fold change), and `means`
#'   (the class \eqn{\times} feature matrix of class means of the rescaled
#'   table).
#'
#' @details
#' **The steps.** Each sample is first rescaled to sum to `total`. Then:
#'
#' 1. **Kruskal-Wallis.** Each feature is tested across the classes, as
#'    [stats::kruskal.test()] does (average ranks, the tie correction). It
#'    passes when its adjusted p-value is below `alpha`.
#' 2. **Wilcoxon** (when `wilcoxon = TRUE`), for the features that passed.
#'    For each pair of classes, every subgroup of one is compared with every
#'    subgroup of the other (with `same_subgroups`, only those of the same
#'    name). When both subgroups have at least `min_subgroup` samples, the
#'    comparison is a two-sided Wilcoxon rank-sum test (the normal
#'    approximation with the variance corrected for ties and no continuity
#'    correction, as `coin::wilcox_test()`'s asymptotic test), which must
#'    give a p-value below `2 * wilcoxon_alpha`; otherwise only the medians
#'    are compared. Every comparison must find the medians in the same
#'    direction, and differ. The pair of classes differs when all its
#'    comparisons do. With two classes the feature passes when they differ;
#'    with more, when one class differs from all others (or, with
#'    `one_against_one`, every pair does).
#' 3. **LDA effect size**, for the features that passed. In each class where
#'    a feature has few distinct values (at most half the class's samples,
#'    or 4), each of its values \eqn{v} is jittered to \eqn{|v + e|},
#'    \eqn{e \sim N(0, \max(0.05 v, 0.01)^2)}. Then, `n_boots` times, a
#'    subsample of `boot_fraction` of the samples is drawn with replacement
#'    (drawn again, up to 1000 times, while a class has too few samples or
#'    a feature too few distinct values in a class), and a linear
#'    discriminant analysis of the classes on all these features is fitted
#'    to it, as `MASS::lda(tol = 1e-10)` would. With \eqn{u} its first
#'    discriminant's coefficients scaled to unit length, each pair of
#'    classes \eqn{(a, b)} gives feature \eqn{j}
#'    \deqn{\frac{|\bar{x}_{aj} - \bar{x}_{bj}| + |u_j| \, |(\bar{x}_a - \bar{x}_b)^\top u|}{2},}
#'    with the subsample's class means. These are averaged over the
#'    subsamples; \eqn{m_j} is the largest over the pairs of classes, and
#'    the LDA score is \eqn{\log_{10}(1 + m_j)}.
#'
#' A feature is reported (`reject`) when it passed both tests and its score
#' is above `lda_threshold`.
#'
#' The LDA score is LEfSe's effect size, not one of [lda()]'s discriminant
#' axes or eigenvalues: half the difference between class means, in the
#' units of `total`, plus a term of the same size from the discriminant,
#' on a log scale. It grows with a feature's abundance.
#'
#' **Use with care.** LEfSe is the most used differential-abundance method
#' in microbiome papers, and ecodive offers it because readers and reviewers
#' expect it, not because it is recommended. In Nearing et al. (2022, 38
#' datasets), it was among the methods calling the most features
#' significant, with among the highest false discovery rates when samples
#' were split at random, and its results were among the least consistent
#' with the other methods'; the authors advised against using it without
#' a correction of the p-values (`p_adjust = 'BH'`). LEfSe itself adjusts
#' nothing (`p_adjust = 'none'` reproduces it); the Wilcoxon step's level
#' is twice `wilcoxon_alpha`; with fewer than `min_subgroup` samples a
#' class (or subgroup) is not tested at all in that step, only its median
#' compared; the score is larger for abundant features; and it compares
#' relative abundances (each sample rescaled to `total`), not absolute ones.
#' The scores, though not the two tests, depend on `seed`. The plain
#' per-feature tests, [da_wilcoxon()], [da_ttest()] and [da_lm()], did
#' better in the benchmarks; [da_sensitivity()] shows whether the features
#' found survive other ways of handling the zeros.
#'
#' **Differences from LEfSe.**
#'
#' * Random numbers come from ecodive's own generator (PCG), not Python's,
#'   so the jitter and the subsamples, and with them the scores, differ as
#'   they would with another seed. Given the same jitter and subsamples, the
#'   scores agree to rounding.
#' * Tables are flat: LEfSe sums the features of a hierarchical table into
#'   their parent clades, and sums only the top level for the sample totals;
#'   here every column is a feature, and a sample's total is the sum of all
#'   of them. To test clades, add their columns to the table.
#' * `p_adjust` is new, and its default (`'holm'`) adjusts the
#'   Kruskal-Wallis p-values; LEfSe adjusts nothing (`'none'`).
#' * Not offered: LEfSe's `-s` (which *widens* the Wilcoxon level to
#'   \eqn{1 - (1 - \alpha)^k} or \eqn{k\alpha} for \eqn{k} comparisons), `-c`
#'   (a beta-version Wilcoxon rule), `-r svm` (disabled in LEfSe 1.1.2) and
#'   subjects (`-u`, used only to sort samples).
#' * When a feature is constant within every class of a subsample, MASS's
#'   `lda()` stops, and LEfSe with it; here the feature is left out of that
#'   subsample's discriminant (its coefficient is 0).
#' * With `same_subgroups`, subgroups are matched by name. LEfSe matches
#'   the ends of their names after dropping as many characters as the
#'   class's name has, which also matches some subgroups of different names
#'   when those names occur in one class only.
#'
#' **Differences from lefser** (lefser 1.12.1, Bioconductor 3.18): lefser
#' allows two classes only; its Wilcoxon step runs only with subclasses
#' (`blockCol`), tests every pair of subclasses whatever their size (no
#' `min_subgroup` median rule), and asks the z statistics' signs to agree
#' instead of the medians'; it fits `MASS::lda()` with its default
#' tolerance (\eqn{10^{-4}}, LEfSe's is \eqn{10^{-10}}); its jitter's
#' standard deviation is one value per class (from the class's largest
#' value), not one per value; its subsamples are \eqn{2/3} of the samples;
#' its scores are signed, from the difference of the second class's mean
#' and the first's (so a feature richer in the first class gets a smaller
#' magnitude than LEfSe gives it); it keeps scores *at least* the threshold
#' (LEfSe: above it); its `relativeAb()` gives proportions, not counts per
#' million; and it uses R's random numbers.
#'
#' **Computation.** The tests run in parallel C over features, from a
#' sparse table without a dense copy (only the features that reach the LDA
#' are made dense); the subsamples' discriminant analyses run in parallel,
#' each subsample with its own random stream, so results are identical for
#' any `cpus`. The discriminant takes a singular value decomposition of the
#' scaled within-class residuals, keeping the directions above the
#' tolerance, so more features than samples are fine (MASS warns that
#' "variables are collinear" there, and so does nothing here); with such
#' collinear features the direction found depends on rounding, in MASS too.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Segata, N., Izard, J., Waldron, L., Gevers, D., Miropolsky, L.,
#' Garrett, W. S., & Huttenhower, C. (2011). Metagenomic biomarker discovery
#' and explanation. *Genome Biology*, 12, R60.
#' \doi{10.1186/gb-2011-12-6-r60}
#'
#' Khleborodova, A., Gamboa-Tuz, S. D., Ramos, M., Segata, N., Waldron, L.,
#' & Oh, S. (2024). lefser: implementation of metagenomic biomarker
#' discovery tool, LEfSe, in R. *Bioinformatics*, 40(12), btae707.
#' \doi{10.1093/bioinformatics/btae707}
#'
#' Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
#' differential abundance methods produce different results across 38
#' datasets. *Nature Communications*, 13, 342.
#' \doi{10.1038/s41467-022-28034-z}
#'
#' Fisher, R. A. (1936). The use of multiple measurements in taxonomic
#' problems. *Annals of Eugenics*, 7, 179-188.
#' \doi{10.1111/j.1469-1809.1936.tb02137.x}
#'
#' @seealso [da_test()], [da_wilcoxon()], [da_lm()], [da_sensitivity()],
#'   [lda()] (a discriminant analysis of the samples, not LEfSe's score),
#'   [plot.ecodive_lefse()].
#'
#' @export
#' @examples
#'     # Features 1 to 4 are richer in group 'b'.
#'     set.seed(1)
#'     x <- matrix(rpois(40 * 12, 50), 40, 12, dimnames = list(NULL, paste0('F', 1:12)))
#'     x[21:40, 1:4] <- x[21:40, 1:4] * 4
#'     grp <- rep(c('a', 'b'), each = 20)
#'
#'     res <- lefse(x, grp, cpus = 1)
#'     res
#'     plot(res)
#'
#'     # As LEfSe reports it, without adjusting the p-values.
#'     lefse(x, grp, p_adjust = 'none', cpus = 1)
#'
lefse <- function (
    counts, groups, data = NULL, subgroups = NULL, margin = 1L, total = 1e6,
    alpha = 0.05, wilcoxon = TRUE, wilcoxon_alpha = 0.05, lda_threshold = 2,
    n_boots = 30L, boot_fraction = 0.67, min_subgroup = 10L,
    one_against_one = FALSE, same_subgroups = FALSE, p_adjust = 'holm',
    seed = 0L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  lefse_run(environment())$result
}



# lefse(), and the details behind it, for tests: `jitter` (n x p, the
# jittered values of the features that reach the LDA, in table order) and
# `boots` (m x n_boots, 1-based samples of each bootstrap) replace the
# random draws. Returns list(result, features (the 1-based features scored),
# x (their jittered values), boots (the samples drawn, 1-based), lda (p x
# pairs of classes: each pair's averaged effect), pairs, wilcoxon (each
# comparison's p-value, n_kw x comparisons, NA where only medians were
# compared)).
lefse_details <- function (
    counts, groups, data = NULL, subgroups = NULL, margin = 1L, total = 1e6,
    alpha = 0.05, wilcoxon = TRUE, wilcoxon_alpha = 0.05, lda_threshold = 2,
    n_boots = 30L, boot_fraction = 0.67, min_subgroup = 10L,
    one_against_one = FALSE, same_subgroups = FALSE, p_adjust = 'holm',
    seed = 0L, cpus = n_cpus(), jitter = NULL, boots = NULL ) {

  if (missing(groups)) groups <- NULL
  lefse_run(environment())
}



#------------------------------------------------------------------------------
# The checks, the three steps, and the result.
#------------------------------------------------------------------------------
lefse_run <- function (env) {

  info <- lefse_input(env)
  n    <- info$n
  D    <- info$D
  cpus <- env$cpus

  #--------------------------------------------------------
  # Classes and subgroups.
  #--------------------------------------------------------
  des <- lefse_groups(env, n, info$samples)
  cls <- des$cls
  lev <- des$levels
  k   <- length(lev)
  sub <- des$sub

  m <- as.integer(n * env$boot_fraction)
  if (m <= k)
    stop(
      '`boot_fraction` leaves ', m, ' samples per bootstrap subsample: ',
      'it must leave more than the ', k, ' groups.', call. = FALSE )

  total <- if (is.null(env$total)) NA_real_ else env$total

  #--------------------------------------------------------
  # Kruskal-Wallis.
  #--------------------------------------------------------
  kw <- .Call(C_lefse_kw, env$counts, env$margin, total, cls - 1L, k, cpus)
  H  <- kw$stat
  H[is.na(H)] <- NaN
  p  <- stats::pchisq(H, k - 1L, lower.tail = FALSE)
  p[is.na(p)] <- NaN

  padj <- p
  ok   <- !is.na(p)
  padj[ok] <- stats::p.adjust(p[ok], env$p_adjust)
  kw_pass <- !is.na(padj) & padj < env$alpha

  means <- kw$means
  dimnames(means) <- list(lev, info$features)
  big      <- apply(means, 2L, max)
  enriched <- lev[max.col(t(means), ties.method = 'first')]

  #--------------------------------------------------------
  # Wilcoxon.
  #--------------------------------------------------------
  wil  <- rep(NA, D)
  wdet <- NULL
  if (env$wilcoxon && any(kw_pass)) {
    fi   <- which(kw_pass)
    w    <- lefse_wilcoxon(env, info, des, fi, total)
    wil[fi] <- w$pass
    wdet <- w$p
  }
  pass <- kw_pass & (!env$wilcoxon | wil %in% TRUE)

  #--------------------------------------------------------
  # LDA effect size.
  #--------------------------------------------------------
  score <- rep(NaN, D)
  lda   <- NULL
  fl    <- which(pass)
  pairs <- utils::combn(k, 2L)
  if (length(fl)) {

    x <- .Call(C_lefse_columns, env$counts, env$margin, total, fl - 1L, cpus)
    jittered <- FALSE
    if (!is.null(env$jitter)) {
      x <- env$jitter
      if (!is.numeric(x) || !identical(dim(x), c(n, length(fl))) || any(!is.finite(x)))
        stop('`jitter` must be a finite ', n, ' x ', length(fl), ' matrix.', call. = FALSE)
      storage.mode(x) <- 'double'
      dimnames(x) <- NULL
      jittered <- TRUE
    }

    bts <- env$boots
    if (!is.null(bts)) {
      if (!is.numeric(bts) || !identical(dim(bts), c(m, as.integer(env$n_boots))) ||
          anyNA(bts) || any(bts < 1 | bts > n | bts %% 1 != 0))
        stop('`boots` must be a ', m, ' x ', env$n_boots, ' matrix of sample indices.', call. = FALSE)
      bts <- matrix(as.integer(bts) - 1L, m)
    }

    min_cl <- max(as.integer(min(tabulate(cls, k)) * env$boot_fraction * env$boot_fraction * 0.5), 1L)

    lda <- .Call(
      C_lefse_lda, x, fl - 1L, jittered, cls - 1L, k,
      pairs[1L, ] - 1L, pairs[2L, ] - 1L, env$n_boots, m, min_cl, bts,
      env$seed, cpus )

    mx <- apply(lda$means, 1L, max)
    score[fl] <- log10(1 + mx)
  }

  reject <- pass & !is.na(score) & score > env$lda_threshold

  out <- data.frame(
    feature      = info$features,
    enriched     = enriched,
    log_max_mean = log10(pmax(big, 1)),
    kw_stat      = H,
    p_value      = p,
    p_adjusted   = padj,
    wilcoxon     = wil,
    lda_score    = score,
    reject       = reject,
    stringsAsFactors = FALSE )

  attr(out, 'method')         <- 'lefse'
  attr(out, 'groups')         <- lev
  attr(out, 'subgroups')      <- des$sub_names
  attr(out, 'alpha')          <- env$alpha
  attr(out, 'wilcoxon_alpha') <- env$wilcoxon_alpha
  attr(out, 'lda_threshold')  <- env$lda_threshold
  attr(out, 'n_boots')        <- env$n_boots
  attr(out, 'boot_fraction')  <- env$boot_fraction
  attr(out, 'total')          <- total
  attr(out, 'p_adjust')       <- env$p_adjust
  attr(out, 'log_base')       <- NA_real_
  attr(out, 'means')          <- means
  class(out) <- c('ecodive_lefse', 'data.frame')

  list(
    result   = out,
    features = fl,
    x        = if (length(fl)) lda$x,
    boots    = if (length(fl)) lda$boots + 1L,
    lda      = if (length(fl)) lda$means,
    pairs    = pairs,
    wilcoxon = wdet )
}



#------------------------------------------------------------------------------
# The arguments' checks. Returns list(D, n, features, samples).
#------------------------------------------------------------------------------
lefse_input <- function (env) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)
  validate_seed(env)

  counts <- env$counts
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
  remove('vals')

  num <- function (arg, lo, hi, lo_ok = FALSE, hi_ok = FALSE, what) {
    v <- env[[arg]]
    ok <- is.numeric(v) && length(v) == 1L && !is.na(v) && is.finite(v) &&
      (if (lo_ok) v >= lo else v > lo) && (if (hi_ok) v <= hi else v < hi)
    if (!ok) stop('`', arg, '` must be ', what, '.', call. = FALSE)
    env[[arg]] <- as.double(v)
  }
  num('alpha',          0, 1, what = 'a single number strictly between 0 and 1')
  num('wilcoxon_alpha', 0, 1, what = 'a single number strictly between 0 and 1')
  num('boot_fraction',  0, 1, hi_ok = TRUE, what = 'a single number above 0 and at most 1')
  num('lda_threshold', -Inf, Inf, what = 'a single finite number')

  if (!is.null(env$total))
    num('total', 0, Inf, what = 'NULL or a single positive number')

  int <- function (arg, lo) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) || v %% 1 != 0 || v < lo || v > .Machine$integer.max)
      stop('`', arg, '` must be a whole number of at least ', lo, '.', call. = FALSE)
    env[[arg]] <- as.integer(v)
  }
  int('n_boots', 1L)
  int('min_subgroup', 0L)

  for (arg in c('wilcoxon', 'one_against_one', 'same_subgroups'))
    if (!isTRUE(env[[arg]]) && !isFALSE(env[[arg]]))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)

  env$p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)

  margin <- env$margin
  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  if (D < 1L) stop('`counts` must have at least one feature.', call. = FALSE)

  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])
  list(
    D        = D,
    n        = n,
    features = if (is.null(features)) as.character(seq_len(D)) else features,
    samples  = if (length(dn[[margin]])) as.character(dn[[margin]]) )
}



#------------------------------------------------------------------------------
# Classes and subgroups. Returns list(cls (1-based class of each sample),
# levels, sub (1-based subgroup of each sample), sub_class (each subgroup's
# class), sub_name (each subgroup's name, NA without subgroups), sub_names
# (NULL or the distinct names)).
#------------------------------------------------------------------------------
lefse_groups <- function (env, n, samples) {

  ref <- '`counts`'
  val <- test_group_values(env$groups, env$data, n, samples, ref = ref)$groups
  lev <- if (is.factor(val)) levels(val) else sort(unique(val), method = 'radix')
  cls <- match(val, lev)
  lev <- as.character(lev)
  if (length(lev) < 2L)
    stop('`groups` must have at least two groups.', call. = FALSE)

  if (is.null(env$subgroups)) {
    return (list(
      cls = cls, levels = lev, sub = cls, sub_class = seq_along(lev),
      sub_name = rep(NA_character_, length(lev)), sub_names = NULL ))
  }

  data <- env$data
  idx  <- NULL
  if (!is.null(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(samples, n, dn, nrow(data), ref = ref)
  }
  sv <- per_sample(env$subgroups, 'subgroups', data, idx, samples, n, ref = ref)
  if (anyNA(sv))
    stop('`subgroups` must not contain NAs.', call. = FALSE)

  slev <- if (is.factor(sv)) levels(droplevels(sv)) else sort(unique(sv), method = 'radix')
  scode <- match(sv, slev)
  slev  <- as.character(slev)

  # Subgroups: each class's subgroup values, in order.
  key  <- (cls - 1L) * length(slev) + scode
  keys <- sort(unique(key))
  list(
    cls       = cls,
    levels    = lev,
    sub       = match(key, keys),
    sub_class = (keys - 1L) %/% length(slev) + 1L,
    sub_name  = slev[(keys - 1L) %% length(slev) + 1L],
    sub_names = slev )
}



#------------------------------------------------------------------------------
# test_rep_wilcoxon_r(), for the features `fi` (1-based). Returns
# list(pass (logical, per feature), p (the comparisons' p-values, NA where
# only medians were compared)).
#------------------------------------------------------------------------------
lefse_wilcoxon <- function (env, info, des, fi, total) {

  k   <- length(des$levels)
  sc  <- des$sub_class
  ns  <- tabulate(des$sub, length(sc))
  nf  <- length(fi)
  has <- !is.null(des$sub_names)

  # Every comparison: pair of classes, subgroups, and whether it is made.
  cmp <- NULL
  for (a in seq_len(k - 1L)) for (b in (a + 1L):k)
    for (i in which(sc == a)) for (j in which(sc == b))
      cmp <- rbind(cmp, c(a, b, i, j, !(has && env$same_subgroups && des$sub_name[[i]] != des$sub_name[[j]])))

  made <- cmp[, 5L] == 1L
  st <- .Call(
    C_lefse_wilcoxon, env$counts, env$margin, total, fi - 1L, des$sub - 1L,
    length(sc), cmp[made, 3L] - 1L, cmp[made, 4L] - 1L, env$cpus )

  # The comparisons' p-values: coin::wilcox_test()'s asymptotic test.
  n1 <- ns[cmp[made, 3L]]
  n2 <- ns[cmp[made, 4L]]
  N  <- n1 + n2
  z  <- t((t(st$w) - n1 * n2 / 2) / sqrt(n1 * n2 / 12 * ((N + 1) - t(st$ties) / (N * (N - 1)))))
  pv <- 2 * stats::pnorm(-abs(z))
  med <- n1 < env$min_subgroup | n2 < env$min_subgroup
  pv[, med] <- NA

  # The decisions, comparison by comparison, as LEfSe makes them.
  npair <- k * (k - 1L) / 2L
  diff  <- matrix(TRUE, nf, npair)
  pidx  <- 0L
  col   <- 0L
  for (a in seq_len(k - 1L)) for (b in (a + 1L):k) {
    pidx  <- pidx + 1L
    rows  <- which(cmp[, 1L] == a & cmp[, 2L] == b)
    brk   <- rep(FALSE, nf)
    first <- rep(TRUE, nf)
    dirn  <- rep(NA, nf)
    for (r in rows) {
      if (!made[[r]]) next
      col <- col + 1L
      sx  <- st$med1[, col]
      sy  <- st$med2[, col]
      lt  <- sx < sy
      eq  <- sx == sy
      live <- !brk
      # Both subgroups hold one and the same value: fails.
      same <- st$same[, col]
      brk[live & same] <- TRUE
      live <- live & !same
      tres <- if (med[[col]]) rep(TRUE, nf) else pv[, col] < 2 * env$wilcoxon_alpha
      tres[is.na(tres)] <- FALSE
      f1 <- live & first
      if (any(f1)) {
        dirn[f1 & tres] <- lt[f1 & tres]
        brk[f1 & !tres] <- TRUE
        first[f1] <- FALSE
      }
      later <- live & !f1
      brk[later & (!tres | lt != dirn | eq)] <- TRUE
    }
    diff[, pidx] <- !brk
  }

  pass <- if (k == 2L) diff[, 1L]
  else if (env$one_against_one) apply(diff, 1L, all)
  else {
    pr <- utils::combn(k, 2L)
    apply(diff, 1L, function (d) any(vapply(seq_len(k), function (cl) {
      all(d[pr[1L, ] == cl | pr[2L, ] == cl]) }, logical(1))))
  }

  colnames(pv) <- NULL
  list(pass = pass, p = pv)
}



#' Plot a LEfSe result
#'
#' LEfSe's bar chart of LDA scores: one horizontal bar per reported feature
#' (`reject`), its length the feature's LDA score, coloured by the class
#' it is enriched in.
#'
#' @param x  The result of [lefse()].
#'
#' @param n_features  How many reported features to show, those with the
#'        largest scores. Default `Inf`: all of them.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `xlab`, `sub` or `xlim`.
#'
#' @details
#' As LEfSe draws it: with two classes, the first class's bars point left
#' (negative) and the second's right, each sorted by score; with more, the
#' bars are grouped by class, in the classes' order, and sorted by score
#' within each. The dashed lines are at `lda_threshold`. The features that
#' passed both tests but scored at or below the threshold are not shown.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     x <- matrix(rpois(40 * 12, 50), 40, 12, dimnames = list(NULL, paste0('F', 1:12)))
#'     x[21:40, 1:4] <- x[21:40, 1:4] * 4
#'     x[1:20, 5:6]  <- x[1:20, 5:6] * 3
#'     plot(lefse(x, rep(c('a', 'b'), each = 20), cpus = 1))
plot.ecodive_lefse <- function (x, n_features = Inf, ...) {

  n   <- plot_count(n_features, 'n_features')
  grp <- attr(x, 'groups', exact = TRUE)
  th  <- attr(x, 'lda_threshold', exact = TRUE)
  k   <- length(grp)
  two <- k == 2L

  rows <- which(x$reject %in% TRUE & is.finite(x$lda_score))
  n_rep <- length(rows)
  rows <- utils::head(rows[order(-x$lda_score[rows])], n)

  cls <- match(x$enriched[rows], grp)
  val <- x$lda_score[rows]
  if (two) val[cls == 1L] <- -val[cls == 1L]

  # Bottom to top, as LEfSe's lefse_plot_res.py orders them.
  o    <- if (two) order(val) else order(cls, val)
  rows <- rows[o]; cls <- cls[o]; val <- val[o]
  kk   <- length(rows)
  col  <- plot_palette(k)

  sub <- sprintf(
    '%s: %d of %d features reported', paste(grp, collapse = ' vs '), n_rep, nrow(x) )
  if (kk < n_rep) sub <- sprintf('%s; top %d shown', sub, kk)

  opar <- graphics::par('mai')
  on.exit(graphics::par(mai = opar))
  linch <- max(graphics::strwidth(x$feature[rows], 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)

  lim  <- max(abs(val), th, 0, na.rm = TRUE)
  xlim <- if (two) c(-lim, lim) else c(0, lim)
  y    <- seq_len(kk)

  plot_do(
    graphics::plot.default, c(list(0, 0), list(...)),
    list(
      type = 'n', yaxt = 'n', xlim = xlim, ylim = c(0.5, max(kk, 1) + 0.5),
      main = paste('LDA scores of', da_method_name(attr(x, 'method', exact = TRUE))),
      xlab = 'LDA score (log 10)', ylab = '', sub = sub ))

  if (kk) {
    graphics::rect(0, y - 0.4, val, y + 0.4, col = col[cls], border = NA)
    graphics::axis(2, at = y, labels = x$feature[rows], tick = FALSE, las = 1, cex.axis = 0.8)
  }
  else {
    graphics::text(mean(xlim), 1, 'No features reported', col = PLOT_MUTED)
  }
  graphics::abline(v = 0, col = 'grey70')
  if (is.finite(th) && th > 0)
    graphics::abline(v = if (two) c(-th, th) else th, lty = 2, col = 'grey50')
  graphics::legend(
    'bottomright', legend = grp, fill = col, border = NA, bty = 'n', cex = 0.8 )

  invisible(x)
}
