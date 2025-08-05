# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' documentation
#' 
#' @name documentation
#' @keywords internal
#' 
#' @param counts   An OTU abundance matrix where each column is a sample, and 
#'        each row is an OTU. Any object coercible with `as.matrix()` can be 
#'        given here, as well as `phyloseq`, `rbiom`, `SummarizedExperiment`, 
#'        and `TreeSummarizedExperiment` objects. 
#' 
#' @param weighted   If `TRUE`, the algorithm takes relative abundances into 
#'        account. If `FALSE`, only presence/absence is considered.
#' 
#' @param normalized   For weighted UniFrac only, normalize distances by the 
#'        total branch length. Options: `TRUE` or `FALSE`.
#' 
#' @param alpha   How much weight to give to relative abundances; a value 
#'        between 0 and 1, inclusive. Setting `alpha=1` is equivalent to 
#'        `weighted_normalized_unifrac()`.
#' 
#' @param tree   A `phylo`-class object representing the phylogenetic tree for 
#'        the OTUs in `counts`. The OTU identifiers given by `colnames(counts)` 
#'        must be present in `tree`. Can be omitted if a tree is embedded with
#'        the `counts` object or as `attr(counts, 'tree')`.
#' 
#' @param pairs   Which combinations of samples should distances be 
#'        calculated for? The default value (`NULL`) calculates all-vs-all. 
#'        Provide a numeric or logical vector specifying positions in the 
#'        distance matrix to calculate. See examples.
#' 
#' @param cpus   How many parallel processing threads should be used. The
#'        default, `n_cpus()`, will use all logical CPU cores.
#' 
NULL
