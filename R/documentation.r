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
#'   * `'hellinger'`: Square root of relative abundance (the Hellinger
#'     transformation).
#'   * `'binary'`: Unweighted presence/absence (each count is either 0 or 1).
#'   * `'clr'`: Centered log ratio.
#'   * `'rclr'`: Robust centered log ratio.
#'   
#'   Default: `'none'`.
#' 
#' @param pairs   Which combinations of samples should distances be 
#'        calculated for? The default value (`NULL`) calculates all-vs-all. 
#'        Provide a numeric or logical vector specifying positions in the 
#'        distance matrix to calculate. See examples.
#' 
#' @param clust  Whether and how to cluster the samples (rows and columns)
#'        before returning the distance matrix: `NULL` (the default) or
#'        `FALSE` to keep them in their input order, or one of
#'        [stats::hclust()]'s linkage methods, `"average"` (UPGMA),
#'        `"complete"`, `"single"`, `"ward.D"`, `"ward.D2"`, `"mcquitty"`
#'        (WPGMA), `"median"` or `"centroid"`, or `"pcoa"` to order them
#'        along the first axis of [pcoa()]. Clustering **reorders the
#'        samples** of the result, with their labels (samples without names
#'        are labelled with their positions in the input): its `dist` record
#'        (see [ecodive_steps]) gives each one's position in the input as
#'        `order`, so that `meta[attr(d, "dist")$order, ]` lines metadata up
#'        with them. ecodive's ordinations and tests put the samples back in
#'        their input order before analysing them, so give them metadata in
#'        its usual order. Other packages see the reordered matrix: give
#'        them metadata reordered by `order`, or matched by name. The
#'        record's `clust` records `clust`
#'        (`FALSE` as `"none"`), which [plot.ecodive_dist()] follows: with `NULL` it clusters the samples
#'        for the heatmap (`"average"`), and with `"none"` it keeps their
#'        order. `clust` needs every distance, so it cannot be used with
#'        `pairs`. [rho_p()]'s similarities are clustered as \eqn{1 -
#'        \rho_p}.
#'
#' @param pairwise  Also test each pair of groups: run the same test on just
#'        the two groups' samples, for every pair, and adjust the pairs'
#'        p-values for multiple testing with `p_adjust`. See **Pairwise
#'        comparisons**. Default: `FALSE`.
#'
#' @param p_adjust  With `pairwise = TRUE`, how to adjust the pairs'
#'        p-values for multiple testing: one of `stats::p.adjust.methods`
#'        (`'holm'`, `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BH'`,
#'        `'BY'`, `'fdr'`, or `'none'`), matched as `stats::p.adjust()`
#'        matches them. Giving it without `pairwise = TRUE` is an error.
#'        Default: `'holm'`.
#'
#' @param power   Scaling factor for the magnitude of differences between
#'        communities (\eqn{p}). Default: `1.5`
#' 
#' @param pseudocount Value added to counts to handle zeros when 
#'        \code{norm = 'clr'}. Ignored for other normalization methods. See 
#'        **Pseudocount** section.
#' 
#' @param strata  Permute samples only within these groups (blocks), for
#'        designs such as repeated measures. With `data`, the name of one of
#'        its columns, or a vector with one value per row of `data`. Without
#'        `data`, a vector with one value per sample, matched to the samples
#'        by name when it has names, otherwise in the same order. Default:
#'        `NULL` (free permutation).
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
#' @name ordination_section
#' @keywords internal
#'
#' @section Ordination result:
#'   Every ordination function ([pcoa()], [pca()], [rpca()], [dbrda()], [capda()],
#'   [capscale()], [mds()], [nmds()], [isomap()], [tsne()], [umap()], and so
#'   [dist_ord()]) returns the same kind of
#'   object: a numeric matrix with one row per sample and one column per
#'   axis, so that it can go straight into a data frame. Rows are
#'   named after the samples, when the input has names. Columns are named
#'   after the method: `"PCoA1"`, `"PC1"`, `"RPCA1"`, `"dbRDA1"` or `"CAP1"`
#'   (then `"MDS1"` for unconstrained axes), `"MDS1"`, `"NMDS1"`,
#'   `"Isomap1"`, `"TSNE1"`, `"UMAP1"`, and so on.
#'
#'   Its attributes always begin with these three, in this order, each with
#'   one value per column:
#'
#'   * `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.
#'   * `eig`: Each axis's eigenvalue, on the method's own scale: [pcoa()]
#'     and [isomap()] give the eigenvalues of Gower's centered matrix (sums
#'     of squares); [pca()] and [rpca()] give variances, which are [pcoa()]'s
#'     divided by \eqn{n - 1} for \eqn{n} samples; and [dbrda()] and
#'     [capscale()] follow `vegan`, which also divides by \eqn{n - 1} when any
#'     distance is 4 or more. `NA` for [mds()], [nmds()], [tsne()], and
#'     [umap()], which have no eigenvalues.
#'   * `prop`: The proportion of the total variance (or inertia) that each
#'     axis explains. Unlike `eig`, it is comparable between methods, so use
#'     it to label axes. `NA` where it isn't defined.
#'
#'   Attributes particular to the method follow them, and last, the records
#'   of how it was made (see [ecodive_steps]): `dist`, the distance
#'   matrix's, for the methods that take distances, and `ord`, the
#'   ordination's own (its method, display name and settings). The matrix has class
#'   `"ecodive_ord"` (as well as `"matrix"`), so that `print()` lists the
#'   other attributes by name rather than in full, and `plot()` draws the
#'   samples on the first two axes, labelled from the attributes, with any
#'   loadings, constraints and centroids on top (see [plot.ecodive_ord()]).
#'   `plot(ord, "Species", iris)` colours the samples by a column of
#'   metadata. Subsetting with `[`, or `as.data.frame()`, gives a plain
#'   matrix or data frame without the attributes.
NULL


#' documentation
#'
#' @name pairwise_section
#' @keywords internal
#'
#' @section Pairwise comparisons:
#'   With `pairwise = TRUE`, each pair of groups is also tested on its own,
#'   as q2-diversity's `beta_group_significance(pairwise=True)` does: the
#'   distances are cut down to the two groups' samples, and the same test is
#'   run on them, with the same number of `permutations`, the same `strata`
#'   (cut down too), and the test's other arguments. The overall test is
#'   unchanged.
#'
#'   * **Permutations:** Each pair has random permutations of its own, drawn
#'     from `seed` and the pair's position in the table, never from the
#'     number of `cpus`. So a pair's p-value is reproducible, but differs
#'     (by sampling error) from that of the same test run by hand on the
#'     pair's samples with the same `seed`. A supplied `permutations` matrix
#'     is cut down to each pair's samples by following each permutation's
#'     cycles: each of the pair's samples goes to the next of the pair's
#'     samples on its cycle. That gives a permutation of the pair's samples,
#'     uniformly distributed when the full one is, which keeps any
#'     restriction the matrix has, such as permuting only within strata.
#'   * **Multiple testing:** The pairs' permutation p-values are adjusted
#'     together with `stats::p.adjust(method = p_adjust)`. The default,
#'     Holm's method, controls the family-wise error rate whatever the
#'     dependence between the pairs' tests. `'BH'` controls the false
#'     discovery rate instead, as q2-diversity's q-values do. The raw
#'     p-values are kept as well.
#'   * **Pairs that can't be tested:** A pair the test can't be run on, such
#'     as two groups of one sample each, gets `NA`s, with a warning, rather
#'     than an error. A test's other warnings for a pair, such as all of its
#'     distances being 0, are passed on, naming the pair.
#'
#'   The `pairwise` field of the result is a data frame with one row per pair
#'   of groups, named like `"a-b"`, in the order of the groups' factor levels
#'   (sorted values, for other vectors). Its columns are the test's
#'   statistics for the pair, `z`, `n` (the pair's number of samples),
#'   `n_perm`, `p_value` (unadjusted) and `p_adjusted`. The `p_adjust` field
#'   names the method, and `print()` shows both.
NULL


#' documentation
#'
#' @name pseudocount_section
#' @keywords internal
#' 
#' @section Pseudocount:
#'   The \code{pseudocount} parameter is only relevant when \code{norm = 'clr'}.
#'   
#'   Zeros are undefined in the centered log-ratio (CLR) transformation. If
#'   \code{norm = 'clr'}, \code{pseudocount} is \code{NULL} (the default), and 
#'   zeros are detected, the function uses half the minimum non-zero value 
#'   (\code{min(x[x>0]) / 2}) and issues a warning.
#'   
#'   To suppress the warning, provide an explicit value (e.g., \code{1}), or
#'   replace the zeros first with \code{\link{replace_zeros}()}, whose
#'   result needs no pseudocount.
#'   
#'   **Why this matters:** The choice of pseudocount is not neutral; it acts as
#'   a weighting factor that can significantly distort downstream results, 
#'   especially for sparse datasets. See Gloor et al. (2017) and Kaul et al. 
#'   (2017) for open-access discussions on the mathematical implications, or 
#'   Costea et al. (2014) for the impact on community clustering.
#'   
#'   See \code{\link{aitchison}} for references.
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
