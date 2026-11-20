# Non-metric Multidimensional Scaling (NMDS)

Ordinates samples from their distances, for plotting beta diversity, by
placing them so that their distances on the plot are in the same rank
order as the distances in `dm`, as nearly as possible. Uses the SMACOF
algorithm of scikit-learn's non-metric `MDS()`, ported to parallel C
that works from `dm` directly, with far less time and memory.

## Usage

``` r
nmds(
  dm,
  n_axes = 2L,
  n_starts = 4L,
  max_iters = 1000L,
  tol = 1e-07,
  ties = "primary",
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- dm:

  A `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  or [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), or a
  symmetric numeric matrix of distances. Must not contain `NA`s or
  negative values.

- n_axes:

  The number of dimensions. Default: `2`

- n_starts:

  The number of random starts to run after the start from the PCoA
  solution, from 0 to 1000. Default: `4`

- max_iters:

  The maximum number of SMACOF iterations for each start. Default:
  `1000`

- tol:

  Each start stops once an iteration lowers the raw stress by no more
  than `tol` times the sum of squared distances on the plot (i.e.,
  Kruskal's stress-1, squared, improves by less than about `tol`).
  Default: `1e-7`

- ties:

  How to treat tied distances in `dm`: `"primary"` (the default) lets
  them have different distances on the plot; `"secondary"` gives them
  equal disparities, as scikit-learn does. See "Ties and zeros" below.

- seed:

  An integer seed for the random starts. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A numeric matrix with one row per sample and one column per dimension,
named `"NMDS1"`, `"NMDS2"`, etc., as described in "Ordination result"
below. `method` is `"nmds"`; `eig` and `prop` are `NA`, as NMDS has no
eigenvalues. Then come:

- `stress`: Kruskal's stress-1 of the solution, from 0 (a perfect rank
  order) to 1.

- `converged`: `TRUE` if the best start stopped because a plain step met
  `tol`, `FALSE` if it stopped at `max_iters`.

- `n_iters`: The number of SMACOF iterations the best start ran.

- `n_starts`: The number of random starts run.

- `best_start`: Which start gave the solution: `0` for the PCoA start,
  else the random start's number.

- `n_best`: How many starts (of `n_starts + 1`) ended within 0.1% of the
  lowest stress. If it is 1, other starts might still find lower stress:
  try more.

- `ties`: The approach to ties used, `"primary"` or `"secondary"`.

With fewer than two samples there is nothing to fit: the result has no
columns, `stress` is 0, and no starts are run.

## Details

**Which ordination?**

[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md), and `nmds()`
all place the samples in `n_axes` dimensions from `dm`.
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) (classical
scaling) is exact and fast, and has eigenvalues and "variance
explained", but it fits the squared distances after double centering, so
large distances dominate, and non-Euclidean distances such as
Bray-Curtis lose part of their structure to negative eigenvalues.
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) (metric MDS)
fits the distances themselves, in their own units, so that distances on
the plot can be read as the data's, with `stress` saying how well.
`nmds()` fits only their rank order, which suits distances whose values
matter less than their order, such as the saturated Bray-Curtis values
of very different samples; it is the usual choice for community data.
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) and `nmds()`
iterate from several starts, the first of them the PCoA solution, and so
take longer than
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).

**Method**

NMDS minimizes Kruskal's stress-1, \\\sqrt{\sum(d\_{ij} -
\hat{d}\_{ij})^2 / \sum d\_{ij}^2}\\, where \\d\_{ij}\\ are the
distances between points on the plot, and the disparities
\\\hat{d}\_{ij}\\ are the closest values to them (in least squares) that
are in the same order as the distances in `dm`, found by isotonic
regression. It is at most 1, as constant disparities are allowed.

Each start minimizes stress with SMACOF (de Leeuw, 1977; Borg & Groenen,
2005), as scikit-learn's non-metric `smacof()` does: the first iteration
uses the distances in `dm` as disparities, and each later one fits them
by isotonic regression; disparities are scaled so that their squares sum
to the number of pairs; then the Guttman transform \\\Gamma\\ moves the
points. A start stops when its raw stress (the numerator above) falls by
no more than `tol` times the denominator in one iteration, or after
`max_iters` iterations. This rule is relative to the sum of squared
distances, not to the stress itself, so a nearly perfect fit (stress-1
below about 0.01) stops while stress can still fall noticeably in
relative terms; a smaller `tol` helps there. At the default `tol`, on
typical data, stress ends within about 0.01-0.1% of its converged value,
and points within about 1-5% of the configuration's spread of their
converged positions; `tol = 1e-8` cuts the stress error about tenfold
and the displacement about threefold, for about 1.5 times the
iterations.

Relaxed steps. To converge in fewer iterations, when `k > 1`, iterations
after the first move the points twice as far as the Guttman transform
would: \\X \leftarrow 2\Gamma(X) - X\\, the relaxed update of de Leeuw
and Heiser (1980). \\\Gamma\\ doesn't depend on the scale of \\X\\, so
relaxed steps don't converge in that one direction (they flip the scale
back and forth), nor at all when `n_axes = 1`, where \\\Gamma\\ only
depends on the order of the points. So a start takes relaxed steps only
until one makes no more than `tol` progress, then plain steps until a
plain step does (or until `max_iters`, whose last step is always plain),
and `n_axes = 1` uses plain steps throughout. The configuration returned
and measured thus always comes from a plain step: it is a fixed point of
\\\Gamma\\ to within `tol`, and its stress is what plain steps reach
from the same start. Relaxed steps typically save about half the
iterations, but can end in a different local minimum than plain steps
would.

In one dimension, SMACOF stops at the first order of the points that is
a fixed point, which is often not the best one: random starts
(`n_starts`) are the way to search for better ones.

**Ties and zeros**

Beta diversity distances are often tied, e.g. at 1 for every pair of
samples with no taxa in common. With `ties = "primary"` (Kruskal's
primary approach, as in `vegan::metaMDS()`), tied distances only have to
stay in order with the other distances: on the plot, they may differ
from each other. With `"secondary"` (as in scikit-learn), tied distances
get one shared disparity, the mean of their distances on the plot, so
the solution tries to make them equal too; this is a stronger
requirement, so its stress is higher. Heavily tied data (e.g. Jaccard or
Bray-Curtis on few taxa or rarefied counts) have many local minima of
similar stress, and benefit from more `n_starts`; see `n_best`.

Distances count as tied when they are within \\10^{-15}\\ times the
largest distance of the first of their group (sorted increasingly),
which absorbs rounding error in how they were computed.

Pairs with a distance of exactly 0, such as duplicate samples, have a
disparity of 0 with either approach, as in scikit-learn: they are pulled
towards the same point. (Under the primary approach, a zero would
otherwise just be the smallest tie, free to take any small distance on
the plot; but identical samples should coincide.) So when `dm` has no
zeros, `stress` with `ties = "primary"` equals `vegan::metaMDS()`'s
stress for the same configuration; when it has zeros, it is higher.

SMACOF only finds a local minimum of stress. So:

1.  The first start is the PCoA solution, from
    [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md). (An
    axis without a positive eigenvalue gets small random coordinates
    instead.)

2.  Then come `n_starts` random starts, with coordinates uniform on \[0,
    1).

3.  The start whose final configuration has the lowest stress-1 is kept;
    if several tie exactly, the earliest. Its configuration is centered,
    rotated to its principal axes, and scaled so that the distances
    between points have the same root mean square as the distances in
    `dm`. Each axis is flipped so that its coordinate with the largest
    absolute value is positive (the first such sample, if several are
    within a relative \\10^{-8}\\ of it).

**Differences from scikit-learn**

With plain steps only (an internal setting used by the tests) and
`ties = "secondary"`, a start follows scikit-learn's iterations from the
same starting configuration to within rounding error: the raw stress of
every iteration agrees to about \\10^{-13}\\, relatively. (Rarely more:
when points nearly coincide, scikit-learn's distance formula rounds
their tiny distances to 0, which can change an iteration count.)
Otherwise:

- scikit-learn takes plain Guttman steps; `nmds()` takes relaxed ones
  (see above).

- scikit-learn's defaults are `max_iter = 300` and `eps = 1e-6`; for a
  good fit (stress-1 of 0.05, squared 0.0025), `1e-6` stops while stress
  still falls by 0.04% per iteration. `nmds()` defaults to
  `max_iters = 1000` and `tol = 1e-7`.

- scikit-learn always uses the secondary approach to ties; `nmds()` uses
  the primary one unless `ties = "secondary"`.

- scikit-learn runs `n_init` random starts (`init = "random"`, its
  current default) or one start from classical MDS
  (`init = "classical_mds"`); `nmds()` runs both kinds.

- scikit-learn ranks starts by the stress of the last iteration,
  measured with the previous iteration's disparities; `nmds()` ranks
  them by stress-1 with disparities fitted to each start's final
  configuration, which is the value reported, and is never higher.

- scikit-learn returns the configuration as SMACOF leaves it; `nmds()`
  centers, rotates, scales, and flips it as described above. None of
  that changes the stress.

- scikit-learn stops when stress falls by *less* than `eps`; `nmds()`
  stops when it falls by no more than `tol`, so that `tol = 0` stops
  once stress no longer falls.

- scikit-learn's tie tolerance is an absolute \\10^{-15}\\; `nmds()`'s
  is relative to the largest distance, so that results don't depend on
  the units of `dm`.

- scikit-learn computes distances on the plot with the expansion
  \\\|x\|^2 - 2 x \cdot y + \|y\|^2\\, and replaces distances of 0 by
  \\10^{-5}\\ in the Guttman transform; `nmds()` computes them directly,
  and leaves out pairs at distance 0, whose contribution \\r\_{ij}(x_i -
  x_j)\\ is 0 either way.

- The random starts come from `seed`, not from R's or NumPy's random
  number generators, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect, and
  R's random number stream is left untouched.

**Speed and memory**

`dm`'s pairs are sorted by distance once, and the order is shared by all
starts: 4 bytes and 1 bit per pair. Each start that runs at the same
time needs 8 more bytes per pair, plus \\256 n k\\ bytes. With
`ties = "primary"`, it also needs 4 bytes per pair for its own order
within tie groups, and scratch space of 12 bytes per pair in the largest
tie group (more for several threads, but at most a quarter of all the
pairs' worth in all). Besides one copy of `dm` as a plain vector of
doubles, nothing of size \\n \times n\\ is formed, and distances on the
plot are recomputed from the coordinates as needed. Sorting the pairs at
the outset takes time proportional to \\n^2 \log n\\, and briefly needs
16 bytes per pair in its largest bucket (about \\N / 16384\\ of the
\\N\\ pairs, for large data) for each thread.

Each SMACOF iteration makes two passes over the pairs, and so takes time
proportional to \\n^2 k\\, plus the isotonic regression, a pass on one
thread. With `ties = "primary"`, each iteration also re-sorts the
distances within each group of tied ones, which is fast when they barely
change, but on heavily tied data can cost several times the rest of the
iteration (e.g. about 4 times, at 1 CPU, on Bray-Curtis of rarefied
counts with 78 distinct values); groups are sorted in parallel, in up to
64 parts. `ties = "secondary"` avoids that cost.

Several starts run at once, one per thread, while their memory totals at
most 1 GiB and one iteration takes under about \\2^{26}\\ pairs of work
(\\n\\ up to about 11,500). Otherwise, starts run one at a time, with
each iteration's passes over the pairs (and the re-sorting of tie
groups) split across `cpus` threads, for data with at least \\2^{18}\\
pairs (\\n\\ of about 725). Long passes run in parts between checks for
user interrupts.

`dm` may have at most 65,536 samples.

**Reproducibility**

The result, with all its attributes, is bit-identical from run to run,
in new R sessions, and for any value of `cpus`, given the same `dm`,
`n_axes`, `n_starts`, `max_iters`, `tol`, `ties`, and `seed`. Each start
splits its sums into a fixed number of parts that depends only on the
number of samples, and adds them in a fixed order; the sort order of the
distances is a total order (tied distances keep their order in `dm`),
and so is the order of each tie group's re-sort; and nothing depends on
the BLAS library. `dm` is read as `dm` times a power of 2 that brings
its largest value near 1, so `dm * 2^p` gives exactly the same result
times \\2^p\\, and neither tiny nor large scales overflow. (As for every
ordination, distances above \\10^{150}\\ are rejected.) Results can
still differ in the last bits between computers or compilers, as for any
floating point code.

`n_axes` larger than \\n - 1\\ is reduced to \\n - 1\\.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capda()`](https://cmmr.github.io/ecodive/reference/capda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md), `nmds()`,
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md))
returns the same kind of object: a numeric matrix with one row per
sample and one column per axis, so that it can go straight into a data
frame. Rows are named after the samples, when the input has names.
Columns are named after the method: `"PCoA1"`, `"PC1"`, `"RPCA1"`,
`"dbRDA1"` or `"CAP1"` (then `"MDS1"` for unconstrained axes), `"MDS1"`,
`"NMDS1"`, `"Isomap1"`, `"TSNE1"`, `"UMAP1"`, and so on.

Its attributes always begin with these three, in this order, each with
one value per column:

- `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.

- `eig`: Each axis's eigenvalue, on the method's own scale:
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) and
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) give
  the eigenvalues of Gower's centered matrix (sums of squares);
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) give
  variances, which are
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)'s divided
  by \\n - 1\\ for \\n\\ samples; and
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) and
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  follow `vegan`, which also divides by \\n - 1\\ when any distance is 4
  or more. `NA` for
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md), `nmds()`,
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), which
  have no eigenvalues.

- `prop`: The proportion of the total variance (or inertia) that each
  axis explains. Unlike `eig`, it is comparable between methods, so use
  it to label axes. `NA` where it isn't defined.

Attributes particular to the method follow them, and last `display`,
with the ordination's title for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`print()`](https://rdrr.io/r/base/print.html): the method's name, and
for the methods that take distances, the distance matrix's title, such
as "Principal Coordinates Analysis of Bray-Curtis Dissimilarity" (see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)).
The matrix has class `"ecodive_ord"` (as well as `"matrix"`), so that
[`print()`](https://rdrr.io/r/base/print.html) lists the other
attributes by name rather than in full, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
samples on the first two axes, labelled from the attributes, with any
loadings, constraints and centroids on top (see
[`plot.ecodive_ord()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md)).
`plot(ord, "Species", iris)` colours the samples by a column of
metadata. Subsetting with `[`, or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), gives a
plain matrix or data frame without the attributes.

## References

Kruskal, J. B. (1964). Nonmetric multidimensional scaling: A numerical
method. *Psychometrika*, 29(2), 115-129.
[doi:10.1007/BF02289694](https://doi.org/10.1007/BF02289694)

de Leeuw, J. (1977). Applications of convex analysis to multidimensional
scaling. In J. R. Barra, F. Brodeau, G. Romier, & B. van Cutsem (Eds.),
*Recent Developments in Statistics* (pp. 133-145). North Holland.

de Leeuw, J., & Heiser, W. J. (1980). Multidimensional scaling with
restrictions on the configurations. In P. R. Krishnaiah (Ed.),
*Multivariate Analysis V* (pp. 501-522). North Holland.

Borg, I., & Groenen, P. J. F. (2005). *Modern Multidimensional Scaling:
Theory and Applications* (2nd ed.). Springer.
[doi:10.1007/0-387-28981-X](https://doi.org/10.1007/0-387-28981-X)

Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
*Journal of Machine Learning Research*, 12, 2825-2830. The C code of
`nmds()` is ported from scikit-learn 1.9.1's `sklearn.manifold.smacof()`
and `sklearn.isotonic`, under scikit-learn's BSD-3-Clause license
(Copyright (c) 2007-2026 The scikit-learn developers).

## See also

[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)

## Examples

``` r
    dm  <- bray(iris[, 1:4])
    ord <- nmds(dm)
    plot(ord, 'Species', iris)

    ord
#> # A matrix: 150 samples × 2 axes
#> # Non-metric Multidimensional Scaling of Bray-Curtis Dissimilarity
#>    NMDS1    NMDS2
#>  1 0.17  -0.0231
#>  2 0.173  0.0162
#>  3 0.191  0.00692
#>  4 0.184  0.0171
#>  5 0.177 -0.025
#>  6 0.136 -0.0536
#>  7 0.19  -0.00604
#>  8 0.168 -0.0129
#>  9 0.199  0.0372
#> 10 0.173  0.011
#> # … 140 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method, eig, prop, stress, converged, n_iters, n_starts,
#> #   best_start, n_best, ties
    str(ord)
#>  'ecodive_ord' num [1:150, 1:2] 0.17 0.173 0.191 0.184 0.177 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "NMDS1" "NMDS2"
#>  - attr(*, "method")= chr "nmds"
#>  - attr(*, "eig")= num [1:2] NA NA
#>  - attr(*, "prop")= num [1:2] NA NA
#>  - attr(*, "stress")= num 0.0378
#>  - attr(*, "converged")= logi TRUE
#>  - attr(*, "n_iters")= int 18
#>  - attr(*, "n_starts")= int 4
#>  - attr(*, "best_start")= int 0
#>  - attr(*, "n_best")= int 1
#>  - attr(*, "ties")= chr "primary"
#>  - attr(*, "display")=List of 1
#>   ..$ title: chr "Non-metric Multidimensional Scaling of Bray-Curtis Dissimilarity"
```
