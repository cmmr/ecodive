# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' documentation
#' 
#' @name documentation
#' @keywords internal
#' 
#' @param alpha   How much weight to give to relative abundances; a value 
#'        between 0 and 1, inclusive. Setting `alpha=1` is equivalent to 
#'        `normalized_unifrac()`.
#'        
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features). 
#'        Typically contains absolute abundances (integer counts), though 
#'        proportions are also accepted.
#' 
#' @param cpus   How many parallel processing threads should be used. The
#'        default, `n_cpus()`, will use all logical CPU cores.
#' 
#' @param cutoff   The maximum number of observations to consider "rare".
#'        Default: `10`.
#' 
#' @param digits   Precision of the returned values, in number of decimal 
#'        places. E.g. the default `digits=3` could return `6.392`.
#' 
#' @param norm   Normalize the incoming counts. Options are:
#'   
#'   * `'none'`: No transformation.
#'   * `'percent'`: Relative abundance (sample abundances sum to 1).
#'   * `'binary'`: Unweighted presence/absence (each count is either 0 or 1).
#'   * `'clr'`: Centered log ratio.
#'   
#'   Default: `'none'`.
#' 
#' @param pairs   Which combinations of samples should distances be 
#'        calculated for? The default value (`NULL`) calculates all-vs-all. 
#'        Provide a numeric or logical vector specifying positions in the 
#'        distance matrix to calculate. See examples.
#' 
#' @param power   Scaling factor for the magnitude of differences between
#'        communities (\eqn{p}). Default: `1.5`
#' 
#' @param pseudocount   The value to add to all counts in `counts` to prevent 
#'        taking `log(0)` for unobserved features. The default, `NULL`, selects 
#'        the smallest non-zero value in `counts`.
#' 
#' @param margin  The margin containing samples. `1` if samples are rows, 
#'        `2` if samples are columns. Ignored when `counts` is a special object 
#'        class (e.g. `phyloseq`). Default: `1`
#' 
#' @param tree   A `phylo`-class object representing the phylogenetic tree for 
#'        the OTUs in `counts`. The OTU identifiers given by `colnames(counts)` 
#'        must be present in `tree`. Can be omitted if a tree is embedded with
#'        the `counts` object or as `attr(counts, 'tree')`.
#' 
#' @section Input Types:
#' 
#'   The `counts` parameter is designed to accept a simple numeric matrix, but 
#'   seamlessly supports objects from the following biological data packages:
#'   
#'   * `phyloseq`
#'   * `rbiom`
#'   * `SummarizedExperiment`
#'   * `TreeSummarizedExperiment`
#'   
#'   For large datasets, standard matrix operations may be slow. See 
#'   `vignette('performance')` for details on using optimized formats 
#'   (e.g. sparse matrices) and parallel processing.
#'   
NULL


#' documentation
#' 
#' @name adiv_assert_integer
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features) 
#'        or a numeric vector representing a single sample. Values must be 
#'        integers (non-integer counts will return an error).
#' 
NULL

#' documentation
#' 
#' @name adiv_percent_normalized
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features) 
#'        or a numeric vector representing a single sample. Raw counts are 
#'        automatically converted to relative abundances (proportions).
#' 
NULL

#' documentation
#' 
#' @name adiv_binary_normalized
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features) 
#'        or a numeric vector representing a single sample. Counts are 
#'        automatically converted to presence/absence (binary) values prior to 
#'        calculation.
#' 
NULL


#' documentation
#' 
#' @name bdiv_assert_integer
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features). 
#'        Values must be integers (non-integer counts will return an error).
#' 
NULL

#' documentation
#' 
#' @name bdiv_percent_normalized
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features). 
#'        Counts are automatically converted to relative abundances (proportions 
#'        summing to 1) prior to calculation.
#' 
NULL

#' documentation
#' 
#' @name bdiv_binary_normalized
#' @keywords internal
#' 
#' @param counts   A numeric matrix of count data (samples \eqn{\times} features). 
#'        Counts are automatically converted to presence/absence (binary) 
#'        values prior to calculation.
#' 
NULL
