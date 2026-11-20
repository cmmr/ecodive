# Phylogenetic Isometric Log-Ratio Transform (PhILR)

Silverman et al.'s PhILR coordinates: an isometric log-ratio transform
whose balances follow a phylogeny. Each internal node of the tree gives
one coordinate, the log ratio of its two child clades' (weighted)
geometric means, so every coordinate has an evolutionary meaning: "this
clade against its sister clade".

## Usage

``` r
philr(
  counts,
  tree = NULL,
  margin = 1L,
  pseudocount = NULL,
  part_weights = "uniform",
  ilr_weights = "uniform",
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of counts or proportions (samples \\\times\\
  features, see `margin`), or an object of one of the types in **Input
  Types**. Its features must have names, matched to the tips of `tree`.
  Counts, proportions and the output of
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  are all accepted, with one caveat: the `'gm.counts'` part weights
  assume counts.

- tree:

  A `phylo` object whose tips include every feature of `counts`, or
  `NULL` to use the tree inside `counts` (a phyloseq, rbiom or
  TreeSummarizedExperiment object, or `attr(counts, 'tree')`). See **The
  tree**.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Added to every value before taking logs. `NULL` (the default) adds
  nothing when `counts` has no zeros, and otherwise uses half the
  smallest non-zero value, with a warning. See **Pseudocount**.

- part_weights:

  The weights \\p_j\\ of the features (philr's `part.weights`), which
  set the reference measure of the simplex. One of

  - `'uniform'` (the default): all 1, the usual ilr;

  - `'gm.counts'`: each feature's geometric mean over the samples (of
    `counts + pseudocount`);

  - `'enorm'`: the Euclidean norm of each feature's values over the
    samples, after closing them to sum to 1 (so that a feature found in
    only a few samples gets a weight near 1, and one spread evenly over
    \\n\\ samples gets \\1/\sqrt{n}\\);

  - `'anorm'`: the Aitchison norm of each feature's values over the
    samples, \\\sqrt{\sum_k (\log x\_{kj} - \overline{\log x_j})^2}\\;

  - `'enorm.x.gm.counts'` and `'anorm.x.gm.counts'`: products of the
    two;

  - a vector of positive numbers, one per feature, named by feature (or
    in the order of the features of `counts`).

  Names may also be written with underscores (`'gm_counts'`). The
  weights computed from `counts` depend on every sample, so one sample's
  coordinates depend on the others.

- ilr_weights:

  The weights \\f_i\\ of the coordinates (philr's `ilr.weights`), which
  scale each balance. One of

  - `'uniform'` (the default): all 1;

  - `'blw'`: the sum of the node's two child branch lengths;

  - `'blw.sqrt'`: its square root;

  - `'mean.descendants'`: the sum of the two child branch lengths plus,
    for each child, its mean distance to its tips;

  - a vector of non-negative numbers, one per coordinate, named by
    coordinate (or in the order of the coordinates).

  For the branch-length weights, zero-length branches to tips first get
  the smallest positive branch length of the tree, as in philr. They
  need `tree` to have branch lengths.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A dense numeric `matrix` of the \\D - 1\\ coordinates, one per internal
node of the tree, in the orientation of `counts`: samples \\\times\\
coordinates for `margin = 1`, coordinates \\\times\\ samples for
`margin = 2`. Sample names are kept. Coordinates are named by the tree's
node labels when every node has a distinct one (nodes added to resolve a
polytomy at node `"X"` get `"X.1"`, `"X.2"`, ...), and otherwise `n1`,
`n2`, ..., in column order. Attributes:

- `method`: `'philr'`.

- `part_weights`: the part weights used, named by feature, in the order
  of the features of `counts`.

- `ilr_weights`: the ILR weights used, named by coordinate.

- `tree`: the tree used, as a `phylo` object (see **The tree**). Its
  internal nodes are numbered and labelled in column order (node \\D +
  k\\ is coordinate \\k\\), and for each node, its first child in
  `tree$edge` is the numerator ("+") clade.

- `pseudocount`: the pseudocount used (`0` for none).

These are what philr's `return.all = TRUE` gives as `p`, `ilr.weights`
and (through the tree) `sbp`. Its contrast matrix `V` is not kept: it
has \\D (D - 1)\\ entries.

## Details

For a sample \\x\\ (with the pseudocount added) and part weights \\p\\,
write \\y_j = x_j / p_j\\. At an internal node with clades \\A\\ (its
first child, the numerator) and \\B\\ (its second), let \\n_A = \sum\_{j
\in A} p_j\\ and \\\log g_A = \sum\_{j \in A} p_j \log y_j / n_A\\ (and
likewise for \\B\\). The node's coordinate is \$\$f \sqrt{\frac{n_A
n_B}{n_A + n_B}}\\\left(\log g_A - \log g_B\right),\$\$ with \\f\\ its
ILR weight (Silverman et al. 2017, Methods). It is computed in \\O(D)\\
per sample, by summing over each clade from the tips up, without forming
the \\D \times (D - 1)\\ basis.

With uniform weights the coordinates are an ordinary ilr: they equal
`ilr(counts, basis = sbp)` for the tree's sequential binary partition,
and `euclidean(philr(x, tree))` is `aitchison(x)` (with the same
pseudocount), as for any orthonormal basis. It is the weights that make
distances phylogenetic:
[`philr_distance()`](https://cmmr.github.io/ecodive/reference/philr_distance.md)
uses the paper's.

The coordinates feed straight into
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
[`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md),
which take negative values with `norm = 'none'`.

`philr(x, tree, pseudocount = c, part_weights = pw, ilr_weights = iw)`
equals philr 1.28.0's
`philr(x, tree, part.weights = pw, ilr.weights = iw, pseudocount = c)`
to rounding, for every weighting and for a binary tree whose tips are
the features of `x`.

## The tree

PhILR needs a rooted binary tree whose tips are the features. `tree` is
made one in three steps (the result is the `tree` attribute):

1.  **Tips that aren't features of `counts` are dropped**, and so are
    the internal nodes left without any. Every feature must be a tip.

2.  **Nodes left with one child are removed**, their branch length added
    to their child's (a root with one child is dropped), as
    [`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html) does.

3.  **Polytomies are resolved, with a warning:** a node with children
    \\c_1, \dots, c_k\\, \\k \> 2\\, keeps \\c_1\\ and gets a new child
    node, joined by a zero-length branch, that holds \\c_2, \dots, c_k\\
    (and so on), as `ape::multi2di(random = FALSE)` does. An unrooted
    tree (three children at the root) is resolved the same way. Results
    then depend on the order of the children: to choose the resolution,
    resolve the tree first.

Coordinates are in the order of the internal nodes' numbers in `tree`
(the root first), or, when a step changed the tree, in preorder of the
new tree (which is also the order of
[`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html)'s and
[`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html)'s
numbering). Branch lengths are needed only for the branch-length ILR
weights; they must be non-negative. A root edge is ignored.

**Same names in other packages**: `philr::philr()` computes the same
transform, with the same orientation and order. Whichever package was
attached last masks the other's `philr()`; write `ecodive::philr()` to
be sure.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## Pseudocount

Logs of zero are undefined, so zeros need replacing first. The
pseudocount is added to every value (not only to zeros), as in
`normalize(norm = 'clr')`. With `pseudocount = NULL` (the default), no
pseudocount is used when `counts` has no zeros; otherwise half the
smallest non-zero value is, with a warning. To avoid a pseudocount,
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md):
its result has none.

The choice of pseudocount is not neutral: it can distort the results,
especially for sparse tables. See
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
(2017). A phylogenetic transform enhances analysis of compositional
microbiota data. *eLife*, 6, e21887.
[doi:10.7554/eLife.21887](https://doi.org/10.7554/eLife.21887)

Egozcue, J. J., & Pawlowsky-Glahn, V. (2016). Changing the reference
measure in the simplex and its weighting effects. *Austrian Journal of
Statistics*, 45(4), 25–44.
[doi:10.17713/ajs.v45i4.126](https://doi.org/10.17713/ajs.v45i4.126)

## See also

[`philr_distance()`](https://cmmr.github.io/ecodive/reference/philr_distance.md),
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)

## Examples

``` r
    x <- philr(ex_counts, tree = ex_tree, pseudocount = 1)
    dim(x)
#> [1] 4 5
    x[1:3, 1:4]
#>               n1       n2        n3        n4
#> Saliva  1.965104 2.899056  2.200302  3.601825
#> Gums    1.186504 2.027911  3.338542  4.231282
#> Nose   -4.604698 0.000000 -1.569240 -1.649373

    # With uniform weights, distances are Aitchison distances
    all.equal(
      as.vector(euclidean(x)),
      as.vector(aitchison(ex_counts, pseudocount = 1)) )
#> [1] TRUE

    # The paper's weights
    x <- philr(
      ex_counts, tree = ex_tree, pseudocount = 1,
      part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt' )
    head(attr(x, 'ilr_weights'))
#>        n1        n2        n3        n4        n5 
#>  3.605551 10.583005  4.898979  5.385165  6.082763 
```
