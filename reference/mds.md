# Metric Multidimensional Scaling (MDS)

Ordinates samples from their distances, for plotting beta diversity, by
placing them so that their distances on the plot match the distances in
`dm` as closely as possible, in least squares. Uses the SMACOF algorithm
of scikit-learn's metric `MDS()`, ported to parallel C that works from
`dm` directly, with far less time and memory.

## Usage

``` r
mds(
  dm,
  k = 2L,
  n_starts = 4L,
  max_iters = 1000L,
  tol = 1e-07,
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

- k:

  The number of dimensions. Default: `2`

- n_starts:

  The number of random starts to run after the start from the PCoA
  solution, from 0 to 1000. Default: `4`

- max_iters:

  The maximum number of SMACOF iterations for each start. Default:
  `1000`

- tol:

  Each start stops once an iteration lowers the raw stress by no more
  than `tol` times the sum of squared distances on the plot. Default:
  `1e-7`

- seed:

  An integer seed for the random starts. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A numeric matrix with one row per sample and one column per dimension,
named `"MDS1"`, `"MDS2"`, etc., as described in "Ordination result"
below. `method` is `"mds"`; `eig` and `prop` are `NA`, as metric MDS has
no eigenvalues. Then come:

- `stress`: The stress of the solution, from 0 (a perfect fit) to 1 (no
  better than putting every point in one place). See "Stress" below.

- `converged`: `TRUE` if the best start stopped because a plain step met
  `tol`, `FALSE` if it stopped at `max_iters`.

- `n_iters`: The number of SMACOF iterations the best start ran.

- `n_starts`: The number of random starts run.

- `best_start`: Which start gave the solution: `0` for the PCoA start,
  else the random start's number.

- `n_best`: How many starts (of `n_starts + 1`) ended within 0.1% of the
  lowest stress.

With fewer than two samples there is nothing to fit: the result has no
columns, `stress` is 0, and no starts are run.

## Details

**Which ordination?**
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md), `mds()`,
and [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) all fit
a configuration to `dm` in `k` dimensions.
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) (classical
scaling) is exact and fast, has eigenvalues and "variance explained",
but fits the squared distances after double centering, so large
distances dominate, and non-Euclidean distances such as Bray-Curtis lose
part of their structure to negative eigenvalues. `mds()` fits the
distances themselves, in their own units, so the plot's distances can be
read as the data's, with `stress` saying how well.
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) fits only
their rank order, which suits distances whose values matter less than
their order, such as the saturated Bray-Curtis values of very different
samples; it is the usual choice for community data. `mds()` and
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) iterate
from several starts and take longer than
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).

**Method**

Metric MDS minimizes the raw stress \\\sigma = \sum\_{i\<j}(d\_{ij} -
\delta\_{ij})^2\\, where \\d\_{ij}\\ are the distances between points on
the plot and \\\delta\_{ij}\\ the distances in `dm`. It does so by
SMACOF (de Leeuw, 1977; Borg & Groenen, 2005), as scikit-learn's metric
`smacof()` does: each iteration moves the points by the Guttman
transform, \\X \leftarrow \Gamma(X) = \frac{1}{n} B(X) X\\, which never
increases \\\sigma\\. A start stops when an iteration lowers \\\sigma\\
by no more than `tol` times \\\sum d\_{ij}^2\\, or after `max_iters`
iterations.

As in [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) (see
there), iterations take relaxed steps, \\X \leftarrow 2\Gamma(X) - X\\
(de Leeuw & Heiser, 1980), when `k > 1`, until one makes no more than
`tol` progress, and then plain steps until a plain step does; `k = 1`
uses plain steps only, and so does the last iteration allowed by
`max_iters`. So the returned configuration is a fixed point of
\\\Gamma\\ to within `tol`, at the right scale (\\\sum
d\_{ij}\delta\_{ij} = \sum d\_{ij}^2\\), typically in about half the
iterations that plain steps take. With plain steps only (an internal
setting used by the tests), a start follows scikit-learn's iterations
from the same starting configuration to within rounding error (about
\\10^{-13}\\, relatively, in every iteration's raw stress).

SMACOF only finds a local minimum of stress, so, as in
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md):

1.  The first start is the PCoA solution, from
    [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md). For
    distances that are Euclidean in `k` dimensions, PCoA recovers the
    configuration, and SMACOF stops at once. (An axis without a positive
    eigenvalue gets small random coordinates instead.)

2.  Then come `n_starts` random starts, with coordinates uniform on \[0,
    1). On beta diversity distances, a random start sometimes ends with
    a stress up to a few percent lower than the PCoA start's.

3.  The start with the lowest stress is kept; if several tie exactly,
    the earliest. Its configuration is centered, rotated to its
    principal axes, and each axis is flipped so that its coordinate with
    the largest absolute value is positive (the first such sample, if
    several are within a relative \\10^{-8}\\ of it). It is not
    rescaled: distances on the plot are in the units of `dm`.

**Stress**

The `stress` reported is \\\sqrt{\sigma / \sum \delta\_{ij}^2}\\: the
raw stress, normalized by the distances in `dm`. It only depends on
\\\sigma\\, so the lowest stress is the best fit, and it doesn't depend
on the units of `dm`. It is 1 when all points coincide, and at most 1
after any SMACOF iteration (`max_iters = 0` returns the starting
configuration, whose stress can be higher). Pairs with a distance of 0
in `dm`, such as duplicate samples, are pulled towards the same point.

scikit-learn reports \\\sigma\\ itself by default for metric MDS, or,
with `normalized_stress = True`, \\\sqrt{\sigma / \sum d\_{ij}^2}\\
(normalized by the distances on the plot instead). For a result `x`,
that is `sqrt(sum((dist(x) - dm)^2) / sum(dist(x)^2))`.

**Differences from scikit-learn**

- scikit-learn takes plain Guttman steps; `mds()` takes relaxed ones
  (see above).

- scikit-learn runs `n_init` random starts (`init = "random"`, its
  current default) or one start from classical MDS
  (`init = "classical_mds"`); `mds()` runs both kinds, and ranks them by
  the stress described above.

- scikit-learn's defaults are `max_iter = 300` and `eps = 1e-6`;
  `mds()`'s are `max_iters = 1000` and `tol = 1e-7`.

- scikit-learn returns the configuration as SMACOF leaves it; `mds()`
  centers, rotates, and flips it, which doesn't change the stress.

- scikit-learn stops when stress falls by *less* than `eps`; `mds()`
  stops when it falls by no more than `tol`, so that `tol = 0` stops
  once stress no longer falls.

- scikit-learn computes distances on the plot with the expansion
  \\\|x\|^2 - 2 x \cdot y + \|y\|^2\\, and replaces distances of 0 by
  \\10^{-5}\\ in the Guttman transform; `mds()` computes them directly,
  and leaves out pairs at distance 0, whose contribution \\r\_{ij}(x_i -
  x_j)\\ is 0 either way.

- The random starts come from `seed`, not from R's or NumPy's random
  number generators, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect, and
  R's random number stream is left untouched.

**Speed and memory**

`mds()` shares its C engine with
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), without
the sorting and the isotonic regression: `dm` is read in its own order,
once converted to a plain vector of doubles. Each start running at the
same time needs \\256 n k\\ bytes; nothing of size \\n \times n\\ is
formed, and there is no limit on the number of samples. Each SMACOF
iteration makes one pass over the pairs, which takes time proportional
to \\n^2 k\\. Several starts run at once, one per thread, while one
iteration takes under about \\2^{26}\\ pairs of work (\\n\\ up to about
11,500); otherwise starts run one at a time, with each iteration split
across `cpus` threads (for at least \\2^{18}\\ pairs), in parts between
checks for user interrupts.

**Reproducibility**

The result, with all its attributes, is bit-identical from run to run,
in new R sessions, and for any value of `cpus`, given the same `dm`,
`k`, `n_starts`, `max_iters`, `tol`, and `seed`, for the reasons given
in [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md).
`dm * 2^p` gives exactly the result for `dm`, times \\2^p\\. (As for
every ordination, distances above \\10^{150}\\ are rejected.) Results
can still differ in the last bits between computers or compilers, as for
any floating point code.

`k` larger than \\n - 1\\ is reduced to \\n - 1\\.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
`mds()`, [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
[`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md))
returns the same kind of object: a numeric matrix with one row per
sample and one column per axis, so that it can go straight into
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
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
  or more. `NA` for `mds()`,
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), which
  have no eigenvalues.

- `prop`: The proportion of the total variance (or inertia) that each
  axis explains. Unlike `eig`, it is comparable between methods, so use
  it to label axes. `NA` where it isn't defined.

Attributes particular to the method follow them. The matrix has class
`"ecodive_ord"` (as well as `"matrix"`), only so that
[`print()`](https://rdrr.io/r/base/print.html) lists the other
attributes by name rather than in full. Subsetting with `[`, or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), gives a
plain matrix or data frame without them.

## References

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
`mds()` is ported from scikit-learn 1.9.1's `sklearn.manifold.smacof()`,
under scikit-learn's BSD-3-Clause license (Copyright (c) 2007-2026 The
scikit-learn developers).

## See also

[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)

## Examples

``` r
    ord <- mds(bray(ex_counts))
    ord
#>              MDS1       MDS2
#> Saliva -0.4402325  0.2063831
#> Gums   -0.4036317 -0.2770339
#> Nose    0.3674288  0.5011040
#> Stool   0.4764354 -0.4304532
#> 
#> mds: 4 samples, 2 axes
#> Attributes: method, eig, prop, stress, converged, n_iters, n_starts, best_start, n_best

    attr(ord, 'stress')
#> [1] 0.1120873

    # Distances that are Euclidean in 2 dimensions are recovered exactly.
    pts <- cbind(x = c(0, 3, 0, 3, 1), y = c(0, 0, 4, 4, 2))
    ord <- mds(dist(pts))
    attr(ord, 'stress')
#> [1] 1.312143e-16
    all.equal(c(dist(ord)), c(dist(pts)))
#> [1] TRUE
```
