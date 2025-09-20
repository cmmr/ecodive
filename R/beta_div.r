# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Beta Diversity Wrapper Function
#' 
#' @inherit documentation
#' 
#' @param method   The name of a beta diversity metric. One of `c('aitchison',
#'   'bhattacharyya', 'bray', 'canberra', 'chebyshev', 'chord', 'clark',
#'   'divergence', 'euclidean', 'generalized_unifrac', 'gower', 'hamming',
#'   'hellinger', 'horn', 'jaccard', 'jensen', 'jsd', 'lorentzian', 'manhattan',
#'   'matusita', 'minkowski', 'morisita', 'motyka', 'normalized_unifrac',
#'   'ochiai', 'psym_chisq', 'soergel', 'sorensen', 'squared_chisq',
#'   'squared_chord', 'squared_euclidean', 'topsoe', 'unweighted_unifrac',
#'   'variance_adjusted_unifrac', 'wave_hedges', 'weighted_unifrac')`. Flexible
#'   matching is supported (see below). Programmatic access via
#'   `list_methods('beta')`.
#'   
#' @param ...  Additional options to pass through to the called function. I.e.
#'   `tree`, `pairs`, `alpha`, or `cpus`.
#' 
#' @return A numeric vector.
#' 
#' 
#' @details
#' 
#' **List of Beta Diversity Metrics**
#' 
#' | Option / Function Name      | Metric Name                                      |
#' | :-------------------------- | :----------------------------------------------- |
#' | `aitchison`                 | Aitchison distance                               |
#' | `bhattacharyya`             | Bhattacharyya distance                           |
#' | `bray`                      | Bray-Curtis dissimilarity                        |
#' | `canberra`                  | Canberra distance                                |
#' | `chebyshev`                 | Chebyshev distance                               |
#' | `chord`                     | Chord distance                                   |
#' | `clark`                     | Clark's divergence distance                      |
#' | `divergence`                | Divergence                                       |
#' | `euclidean`                 | Euclidean distance                               |
#' | `generalized_unifrac`       | Generalized UniFrac (GUniFrac)                   |
#' | `gower`                     | Gower distance                                   |
#' | `hamming`                   | Hamming distance                                 |
#' | `hellinger`                 | Hellinger distance                               |
#' | `horn`                      | Horn-Morisita dissimilarity                      |
#' | `jaccard`                   | Jaccard distance                                 |
#' | `jensen`                    | Jensen-Shannon distance                          |
#' | `jsd`                       | Jesen-Shannon divergence (JSD)                   |
#' | `lorentzian`                | Lorentzian distance                              |
#' | `manhattan`                 | Manhattan distance                               |
#' | `matusita`                  | Matusita distance                                |
#' | `minkowski`                 | Minkowski distance                               |
#' | `morisita`                  | Morisita dissimilarity                           |
#' | `motyka`                    | Motyka dissimilarity                             |
#' | `normalized_unifrac`        | Normalized Weighted UniFrac                      |
#' | `ochiai`                    | Otsuka-Ochiai dissimilarity                      |
#' | `psym_chisq`                | Probabilistic Symmetric Chi-Squared distance     |
#' | `soergel`                   | Soergel distance                                 |
#' | `sorensen`                  | Dice-Sorensen dissimilarity                      |
#' | `squared_chisq`             | Squared Chi-Squared distance                     |
#' | `squared_chord`             | Squared Chord distance                           |
#' | `squared_euclidean`         | Squared Euclidean distance                       |
#' | `topsoe`                    | Topsoe distance                                  |
#' | `unweighted_unifrac`        | Unweighted UniFrac                               |
#' | `variance_adjusted_unifrac` | Variance-Adjusted Weighted UniFrac (VAW-UniFrac) |
#' | `wave_hedges`               | Wave Hedges distance                             |
#' | `weighted_unifrac`          | Weighted UniFrac                                 |
#' 
#' 
#' 
#' 
#' 
#' 
#' **Flexible name matching**
#' 
#' Case insensitive and partial matching. Any runs of non-alpha characters are
#' converted to underscores. E.g. `metric = 'Weighted UniFrac` selects
#' `weighted_unifrac`.
#' 
#' UniFrac names can be shortened to the first letter plus "unifrac". E.g. 
#' `uunifrac`, `w_unifrac`, or `V UniFrac`. These also support partial matching.
#' 
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the metrics list.
#' 
#' 
#' @export
#' @examples
#'     # Example counts matrix
#'     ex_counts
#'     
#'     # Bray-Curtis distances
#'     beta_div(ex_counts, 'bray')
#'     
#'     # Generalized UniFrac distances
#'     beta_div(ex_counts, 'GUniFrac', tree = ex_tree)
#'     
beta_div <- function (x, method, ...) {
  match_method(method, div = 'beta')$func(x = x, ...)
}





