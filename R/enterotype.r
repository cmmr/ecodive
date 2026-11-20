# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Written from the published algorithms (clean room): PAM from Kaufman &
# Rousseeuw (1990), silhouettes from Rousseeuw (1987), and the
# Calinski-Harabasz index from Calinski & Harabasz (1974). cluster::pam(),
# cluster::silhouette() and clusterSim::index.G1() (all GPL) were only run, to
# compare outputs. See src/pam.c.


#' Enterotypes
#'
#' Clusters samples into community types as Arumugam et al. (2011) did for
#' the human gut: partitioning around medoids (PAM) on Jensen-Shannon
#' distances, with the number of types chosen by the Calinski-Harabasz (CH)
#' index, and each sample's fit given by its silhouette.
#'
#' @param x  A count table (or proportions), from which [jensen()]
#'        distances are computed; or a `dist` object, such as one from
#'        [beta_div()], used as it is, so that other distances can be
#'        tried. Samples are rows unless `margin = 2`.
#'
#' @param n_types  The numbers of types to try (each from 2 to the number
#'        of samples minus 1); a single value fixes it. Default `2:10`, cut
#'        to `2:(n - 1)` for fewer than 11 samples.
#'
#' @param ch  Which form of the Calinski-Harabasz index chooses the number
#'        of types: `"medoids"` (the default), as
#'        `clusterSim::index.G1(centrotypes = "medoids")` computes it in
#'        Arumugam et al.'s published script, or `"centroids"`, the pseudo-F
#'        of the distances. See Details.
#'
#' @param counts  Only with distances as `x`: an optional count table of
#'        the same samples (matched by name, or else by position), from which
#'        each type's mean composition (`profiles`) and `top_feature` are
#'        computed. Without it, a `dist` gives no profiles.
#'
#' @param margin  If your samples are in the table's rows, set to `1L`
#'        (the default). If they are in columns, set to `2L`.
#'
#' @param cpus  How many parallel processing threads to use. Results are
#'        the same for any number. Default: all available cores.
#'
#' @details
#' **PAM** (Kaufman & Rousseeuw 1990, ch. 2) chooses \eqn{k} samples as
#' medoids, and assigns every sample to its nearest medoid. BUILD picks the
#' medoids one at a time: first the sample with the smallest total distance
#' to all others, then each time the sample that most reduces the total
#' distance to the nearest medoid. SWAP then makes the best exchange of a
#' medoid with a non-medoid, as long as one lowers the total. Ties go to the
#' Types are numbered in the order of their medoids'
#' positions in the input. Where samples tie, PAM breaks the tie as
#' `cluster::pam()` does in nearly every case tested: in BUILD, the highest
#' sample index; in SWAP, the first swap found; and a sample equally near
#' two medoids joins the lower-indexed one.
#'
#' **The Calinski-Harabasz index** compares the spread between types with
#' the spread within them; the number of types with the largest index is
#' chosen (the smallest, on ties). Both forms are
#' \deqn{CH = \frac{(T - W)/(k - 1)}{W/(n - k)}}
#' for \eqn{n} samples in \eqn{k} types, with \eqn{T} the total spread
#' and \eqn{W} the spread within types.
#'
#' * `ch = "medoids"` follows `clusterSim::index.G1(x, cl, d, centrotypes
#'   = "medoids")`, which Arumugam et al.'s script calls with the
#'   Jensen-Shannon distances as `d` and the samples' relative abundances
#'   as `x`. The medoid \eqn{m_c} of each type, and \eqn{M} of all
#'   samples, is the member with the smallest total distance to the others.
#'   The spreads are then measured by squared Euclidean distances
#'   \eqn{e} between the samples' proportions:
#'   \eqn{T = \sum_i e(i, M)^2} and \eqn{W = \sum_i e(i, m_{c(i)})^2}.
#'   For distances given without `counts`, \eqn{e} is the distance itself;
#'   this equals `index.G1()` when the distances are Euclidean distances
#'   of its `x`.
#' * `ch = "centroids"` is the pseudo-F of the distances (the F statistic
#'   of [permanova()]):
#'   \eqn{T = \frac1n\sum_{i<j}d_{ij}^2} and
#'   \eqn{W = \sum_c \frac{1}{n_c}\sum_{i<j \in c} d_{ij}^2}, for types
#'   of \eqn{n_c} samples. For Euclidean distances, this is exactly
#'   Calinski & Harabasz's (1974) index. It uses only the distances.
#'
#' **The silhouette** of sample \eqn{i} (Rousseeuw 1987) is
#' \eqn{s_i = (b_i - a_i) / \max(a_i, b_i)}, where \eqn{a_i} is its mean
#' distance to the rest of its type and \eqn{b_i} its smallest mean
#' distance to another type; 0 for a sample alone in its type. Values near
#' 1 mean a sample sits well within its type; values near 0 or below, that
#' it lies between types. Arumugam et al. checked their clusters with the
#' mean silhouette, which `selection` reports for every number of types.
#'
#' **Distances.** The paper used the Jensen-Shannon *distance*, the square
#' root of the divergence: [jensen()], not [jsd()]. The base of the
#' logarithm doesn't matter, since PAM, the CH index and silhouettes are
#' unchanged when every distance is scaled. Zeros need no pseudocount
#' (\eqn{0 \log 0 = 0}); some published scripts add a small one (such as
#' \eqn{10^{-6}}), which changes the distances, and so the results,
#' slightly. Arumugam et al. used genus abundances, so aggregate features
#' to genus first to follow them.
#'
#' **A caution.** Costea et al. (2018), with the authors of the original
#' paper, describe enterotypes as a coarse stratification of what is a
#' continuous landscape of community composition, rather than discrete
#' types. Check the silhouettes, and consider ordinations (e.g. [pcoa()])
#' alongside them.
#'
#' @section Differences from cluster and clusterSim:
#' * Medoids, assignments, objectives and silhouettes match
#'   `cluster::pam()` (with its default `pamonce = FALSE`) and
#'   `cluster::silhouette()` to rounding, and the CH index matches
#'   `clusterSim::index.G1()`, as described above.
#' * `cluster::pam()` numbers its clusters by their first appearance in the
#'   data; `enterotype()` numbers types by their medoids' positions. The
#'   partitions are the same.
#' * Exact ties can still be broken differently: `cluster::pam()` breaks a
#'   few ties in its first BUILD step otherwise, which can lead to another,
#'   equally good, set of medoids.
#'
#' @return A data frame of class
#'   `c("ecodive_enterotype", "ecodive_types", "data.frame")`, one row per
#'   sample (see [ecodive_types]), with columns `sample`, `type` (a factor)
#'   and `silhouette`. Its attributes:
#'   * `types`: a data frame, one row per type, with `type`, `n_samples`,
#'     `top_feature` (the feature of largest mean proportion), `medoid` (the
#'     medoid sample's name, or its position) and `silhouette` (the mean of
#'     the type's samples).
#'   * `profiles`: a types \eqn{\times} features matrix of mean proportions,
#'     or `NULL` for distances without `counts`.
#'   * `selection`: a data frame, one row per number of types tried, with
#'     `n_types`, `ch`, `silhouette` (the mean), `objective` (PAM's mean
#'     distance to the nearest medoid) and `selected`.
#'   * `objective`: PAM's mean distance from each sample to its nearest
#'     medoid for the chosen number of types, after BUILD and after SWAP,
#'     as `cluster::pam()` reports it.
#'   * `ch`: how the CH index was computed, `"medoids"` or `"centroids"`.
#'   * `method`: `"enterotype"`.
#'
#'   `plot()` draws the types' profiles or the selection (see
#'   [plot.ecodive_types()]); its title is "Enterotypes of" the distance's
#'   title.
#'
#' @section Same names in other packages:
#' None known.
#'
#' @references
#' Arumugam, M., Raes, J., Pelletier, E., et al. (2011). Enterotypes of the
#' human gut microbiome. *Nature*, 473, 174-180. \doi{10.1038/nature09944}
#'
#' Costea, P. I., Hildebrand, F., Arumugam, M., et al. (2018). Enterotypes
#' in the landscape of gut microbial community composition. *Nature
#' Microbiology*, 3, 8-16. \doi{10.1038/s41564-017-0072-8}
#'
#' Kaufman, L., & Rousseeuw, P. J. (1990). *Finding Groups in Data: An
#' Introduction to Cluster Analysis*. Wiley. \doi{10.1002/9780470316801}
#'
#' Calinski, T., & Harabasz, J. (1974). A dendrite method for cluster
#' analysis. *Communications in Statistics*, 3, 1-27.
#' \doi{10.1080/03610927408827101}
#'
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the
#' interpretation and validation of cluster analysis. *Journal of
#' Computational and Applied Mathematics*, 20, 53-65.
#' \doi{10.1016/0377-0427(87)90125-7}
#'
#' @family community types
#' @seealso [ecodive_types], [plot.ecodive_types()], [jensen()]
#' @export
#' @examples
#'     # Types of iris flowers, from the proportions of their four
#'     # measurements (a stand-in for a table of genus abundances).
#'     et <- enterotype(iris[, 1:4])
#'     et
#'     table(et$type, iris$Species)
#'
#'     # The types as groups in a test.
#'     permanova(jensen(iris[, 1:4]), et$type)
#'
#'     # Another distance, with the abundances for the profiles.
#'     enterotype(bray(iris[, 1:4]), n_types = 3, counts = iris[, 1:4])
enterotype <- function (
    x, n_types = 2:10, ch = c('medoids', 'centroids'), counts = NULL,
    margin = 1L, cpus = n_cpus() ) {

  ch <- match.arg(ch)
  validate_margin()
  validate_cpus()
  default_k <- missing(n_types)


  #--------------------------------------------------
  # Distances, and the samples' proportions.
  #--------------------------------------------------
  if (inherits(x, 'dist')) {

    dv     <- dist_vector(x)
    dm     <- dv$dm
    n      <- dv$n
    labels <- dv$labels
    title  <- dv$title

    p <- NULL
    if (!is.null(counts)) {
      p <- types_proportions(counts, margin)
      if (nrow(p) != n)
        stop('`counts` must have the same ', n, ' samples as the distances.', call. = FALSE)
      if (!is.null(labels) && !is.null(rownames(p))) {
        i <- match(labels, rownames(p))
        if (anyNA(i))
          stop('`counts` lacks samples in the distances: ',
               paste(utils::head(labels[is.na(i)], 5), collapse = ', '), call. = FALSE)
        p <- p[i, , drop = FALSE]
      }
    }
  }
  else {

    if (!is.null(counts))
      stop('Give `counts` only with distances as `x`.', call. = FALSE)

    p  <- types_proportions(x, margin)
    d  <- jensen(p, margin = 1L, cpus = cpus)
    dv <- dist_vector(d)
    dm     <- dv$dm
    n      <- dv$n
    labels <- rownames(p)
    title  <- dv$title
  }

  if (n < 3)
    stop('enterotype() needs at least 3 samples.', call. = FALSE)

  n_types <- validate_n_types(n_types, n, default_k)
  samples <- if (is.null(labels)) as.character(seq_len(n)) else as.character(labels)


  #--------------------------------------------------
  # PAM, CH and silhouettes for each number of types.
  #--------------------------------------------------
  fits <- lapply(n_types, function (k) {
    fit <- pam_fit(dm, n, k, cpus)
    sil <- .Call(C_silhouette, dm, n, fit$clustering, k, cpus)[[1]]
    fit$ch         <- .Call(C_calinski, dm, n, fit$clustering, k, as.integer(ch == 'medoids'), p)
    fit$silhouette <- sil
    fit
  })

  ch_v <- vapply(fits, `[[`, 0, 'ch')
  best <- which.max(replace(ch_v, is.nan(ch_v), -Inf))
  if (length(best) == 0) best <- 1L # # nocov
  fit  <- fits[[best]]
  k    <- n_types[[best]]

  selection <- data.frame(
    n_types    = n_types,
    ch         = ch_v,
    silhouette = vapply(fits, function (f) mean(f$silhouette), 0),
    objective  = vapply(fits, function (f) f$objective[[2]], 0),
    selected   = seq_along(n_types) == best )


  #--------------------------------------------------
  # The result.
  #--------------------------------------------------
  type <- factor(fit$clustering, levels = seq_len(k), labels = as.character(seq_len(k)))

  profiles <- if (!is.null(p)) types_profiles(p, type)

  types <- data.frame(
    stringsAsFactors = FALSE,
    type        = factor(levels(type), levels = levels(type)),
    n_samples   = tabulate(fit$clustering, k),
    top_feature = types_top_feature(profiles, k),
    medoid      = samples[fit$medoids],
    silhouette  = as.vector(tapply(fit$silhouette, type, mean)) )

  df <- data.frame(
    stringsAsFactors = FALSE,
    sample     = samples,
    type       = type,
    silhouette = fit$silhouette )

  res <- new_ecodive_types(
    df, 'enterotype', types, profiles, selection,
    objective = c(build = fit$objective[[1]], swap = fit$objective[[2]]),
    ch        = ch )

  set_display(res, title = title_of('Enterotypes', title))
}



# PAM (src/pam.c) with `k` medoids on the packed distances `dm` of `n`
# samples: list(medoids, clustering, objective, n_swaps). Medoids are
# 1-based and ascending; clusters are numbered in the medoids' order.
pam_fit <- function (dm, n, k, cpus = 1L) {
  res <- .Call(C_pam, dm, as.integer(n), as.integer(k), as.integer(cpus))
  names(res) <- c('medoids', 'clustering', 'objective', 'n_swaps')
  res
}


# The samples' proportions, samples x features, from a count table `x` of
# orientation `margin`. Every sample must have a positive total.
types_proportions <- function (x, margin) {
  counts <- x
  validate_counts()
  if (inherits(counts, 'simple_triplet_matrix'))
    counts <- as.matrix(counts)
  counts <- as.matrix(counts)
  if (margin == 2L) counts <- t(counts)
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)))
    stop('The count table must hold finite numbers, with no NA.', call. = FALSE)
  if (any(counts < 0))
    stop('The count table must not contain negative values.', call. = FALSE)
  tot <- rowSums(counts)
  if (any(tot <= 0))
    stop(
      'Every sample needs a positive total; these have none: ',
      paste(utils::head(if (is.null(rownames(counts))) which(tot <= 0) else rownames(counts)[tot <= 0], 5), collapse = ', '),
      call. = FALSE )
  counts / tot
}


# Checks `n_types` against `n` samples. The default is cut to 2 .. n - 1.
validate_n_types <- function (n_types, n, default = FALSE) {
  if (default) return (seq.int(2L, min(10L, n - 1L)))
  if (!is.numeric(n_types) || length(n_types) == 0 || anyNA(n_types) ||
      any(n_types %% 1 != 0))
    stop('`n_types` must be whole numbers.', call. = FALSE)
  n_types <- sort(unique(as.integer(n_types)))
  if (any(n_types < 2) || any(n_types > n - 1))
    stop('`n_types` must be between 2 and ', n - 1, ' (the number of samples minus 1).', call. = FALSE)
  n_types
}
