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
#'        Also supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and 
#'        `TreeSummarizedExperiment` objects. See `vignette('performance')` for 
#'        optimizing large datasets.
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
#' @param norm   Normalize the incoming counts. Options are:
#'   
#'   * `'percent'`: Relative abundance (sample abundances sum to 1).
#'   * `'binary'`: Unweighted presence/absence (each count is either 0 or 1).
#'   * `'clr'`: Centered log ratio.
#'   * `'none'`: No transformation.
#'   
#'   Default: `'percent'`, which is the expected input for these formulas.
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
NULL
