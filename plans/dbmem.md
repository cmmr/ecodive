# `dbmem()`: distance-based Moran's eigenvector maps

**Priority:** 78. **Tier:** clean room. Write it from Borcard & Legendre
(2002) and Dray, Legendre & Peres-Neto (2006). adespatial `dbmem()` (GPL)
is output-only. It uses `moran_i()` ([moran_i.md](moran_i.md)) to select
eigenvectors. Read [README.md](README.md) first.


## What it is

dbMEM (formerly PCNM) turns sample positions into a set of orthogonal
spatial variables, at scales from broad to fine. They are used as
explanatory variables in `rda()`, `varpart()` and `forward_select()`, to
model spatial structure. The steps are:

1. Truncate the geographic distance matrix at the minimum-spanning-tree
   threshold t (distances > t are set to 4t).
2. Do a PCoA.
3. Keep the eigenvectors with positive eigenvalues (or with positive
   Moran's I).


## Design

- **Signature:** `dbmem(space, threshold = NULL, moran_positive = TRUE, cpus = n_cpus())`.
  - `space` gives the samples' positions, named as in `moran_i()`: a dist
    of geographic distances, or a coordinate matrix (made Euclidean).
  - `threshold` defaults to the longest edge of the minimum spanning tree.
  - `moran_positive` keeps only the eigenvectors with positive Moran's I
    (adespatial's `MEM.autocor = "positive"`, its default).
- **Result:** a matrix of the spatial variables (samples × MEMs), columns
  `MEM1..`.
  - Attributes: `values` (eigenvalues), `moran_i` (each MEM's Moran's I),
    `threshold`.
  - It is **not** an `ecodive_ord` (it isn't an ordination for plotting;
    it is a set of explanatory variables), so it is a plain function.
  - Scale the vectors as adespatial does (by √n, or not — check this,
    since spacc's scaling is unverified).
- **Masking:** `adespatial::dbmem`.


## Reference implementations

- **No vetted permissive implementation.** spacc 0.11.0 (R, MIT) has
  `spatialEigenvectors(method = "dbmem")`, but it is a new single-author
  package, and its Moran's I uses a binary connectivity matrix that may
  select different vectors from adespatial's. Treat it as a hint at most,
  and don't port from it without checking its license file and results.
- **Write from the papers.** Reuse ecodive's code:
  - `pcoa()` (block Lanczos, no n × n matrix);
  - a minimum spanning tree over the distances of `space` (new and small: Prim's algorithm on
    the packed triangle);
  - `moran_i()`.
- **Output-only fixtures (GPL):** adespatial 0.3-30
  `dbmem(xyORdist, thresh, MEM.autocor = "positive")` and its attributes
  (`values`, `Moran's I`); vegan `pcnm()` as a cross-check of the
  truncation step.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Not an ordination, so not made by `ord_result()`. Its attributes are
  results: `values`, `moran_i` (one per MEM), and `threshold` (it
  qualifies the MEMs), plus `method = 'dbmem'`.
- A `plot()` method (MEMs mapped over coordinates, or their Moran's I)
  reads only those. Its title is the method's name, "Moran's Eigenvector
  Maps"; it names no input (positions have no title), so no `display` is
  needed.

## Tests

- Fixtures from adespatial on `mite.xy` and a regular grid. Compare the
  MEM vectors (sign-matched), eigenvalues, Moran's I and threshold at
  1e-10.
- The threshold equals the MST's longest edge, computed independently in
  the test.
- Duplicate positions; collinear points (a transect).
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page in the "Spatial" group with `moran_i()`.
- A vignette example: `forward_select()` over MEMs, then `varpart()` of
  environment against space.
- NEWS.


## References

- Borcard, D., & Legendre, P. (2002). All-scale spatial analysis of
  ecological data by means of principal coordinates of neighbour matrices.
  *Ecological Modelling*, 153, 51–68.
- Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a
  comprehensive framework for principal coordinate analysis of neighbour
  matrices (PCNM). *Ecological Modelling*, 196, 483–493.


## Decisions made while implementing

Asked of the user at the start (2026-11-26), each the recommended option:

- **`autocor`** replaces the plan's `moran_positive`: `'positive'` (the
  default), `'negative'`, `'non-null'` or `'all'`, adespatial's
  `MEM.autocor`. The signature is `dbmem(space, data = NULL, threshold =
  NULL, autocor = 'positive', cpus = n_cpus())`.
- **The result is a matrix of class `c('ecodive_dbmem', 'matrix',
  'array')`**, as `batch_correct()`'s `ecodive_batch`, not an
  `ecodive_ord` and not a data frame: a short `print()` preview, and
  `plot()` (a bar chart of each MEM's Moran's I, or, with `coords`, maps of
  `mems` over the samples' coordinates, given at plot time since the result
  doesn't record them). It works as is as `data` for `forward_select()`, in
  `cbind(meta, mem)` for `varpart()`, and with `as.data.frame()`; columns
  taken with `[` are a plain matrix.
- **A dense eigensolver in C, without the BLAS**: Householder reduction to
  tridiagonal form, accumulation of the reflections, and the QL algorithm
  with implicit shifts (JAMA's / EISPACK's tql2, as `lz_sym_eigen()`), with
  the O(n^3) steps in parallel. Block Lanczos was not used: dbmem keeps
  about a third of the n - 1 eigenvectors by default (all with `'all'`),
  where a dense solver is faster.
- **adespatial's scales:** MEMs centred with sum of squares n; `values` the
  eigenvalues of the doubly centred weights divided by n; `moran_i` each
  MEM's I with the weights as given, which is n^2 values / S0.
- **Neighbours are pairs within t (1 + 1e-6)**, as adespatial (found from
  its outputs: its neighbour list's `dnn` attribute is `t * (1 + 1e-6)`),
  with weights 1 - (d / 4t)^2 using t itself.
- **`data`**, so that `space` can name coordinate columns, as in
  `moran_i()`. A square weights matrix for `space` is an error.
- **Thresholds below the MST's longest edge** are allowed with a warning
  (some samples unconnected; how many have no neighbour at all), as
  adespatial allows them silently; no pair within the threshold is an
  error.

Made while implementing:

- **The ones vector is removed before the eigendecomposition**: the
  Householder reflection H that takes 1/sqrt(n) to -e1 gives H A H a zero
  first row and column, and the (n - 1) x (n - 1) trailing block is
  decomposed. So the n - 1 MEMs are exactly orthogonal to 1 (to rounding)
  however many eigenvalues are 0, with no constant vector to find and drop.
  B is built straight from the packed distances (row means of W, then the
  doubly centred entries less the reflection's rank-2 term), in parallel by
  column.
- **Parallel steps, bit-identical for any `cpus`:** the reduction's
  matrix-vector product and rank-2 update, and the accumulation of the
  reflections, are split by column; each QL sweep's plane rotations are
  recorded, then applied to the eigenvector matrix (stored by rows) in
  fixed blocks of 4 rows, interleaved, since each row's rotations are a
  chain of dependent steps. That cut the QL step from 18.8 s to 5.6 s at n
  = 2000 on one thread; blocks of 8 (as an array) were slower. Dot products
  use four partial sums, in a fixed order. Threads are woken only for
  steps of at least 5e5 flops per thread (`threads_for()`): of 2.5e5, 5e5
  and 1e6, the best balance between n = 500 (where 4 threads were slower
  than 1 at 2.5e5) and n = 1000 to 2000.
- **Null eigenvalues**: within sqrt(eps) (about 1.5e-8) of the largest
  absolute eigenvalue. Probing adespatial's outputs (a grid with exact
  zero eigenvalues, perturbed by varying amounts) put its cut between
  1.22e-8 and 1.54e-8 relative to the largest absolute eigenvalue,
  consistent with sqrt(eps). Null ones are left out of `'positive'`,
  `'negative'` and `'non-null'`, and have `values` of exactly 0 in
  `'all'`, as adespatial's.
- **Moran's I is computed directly** from the list of neighbour pairs,
  n / S0 * v'Wv / v'v, in parallel over MEMs (O(pairs) each), rather than
  with `moran_weights()` and `moran_stat()`: those need W as an n x n R
  matrix and run one variable per single-threaded pass. The formula is the
  same as `moran_i(style = 'none')`'s, which the tests check on three
  MEMs, and equals n^2 values / S0 to rounding.
- **`moran_space()`** (R/moran_i.r) gained `distances = FALSE`; with
  `TRUE` (dbmem), a square matrix is an error and zero distances are kept
  (they are errors for inverse-distance weights).
- **The threshold** defaults to the minimum spanning tree's longest edge,
  by Prim's algorithm on the packed triangle (O(n^2), in C). The warning
  for a threshold below it allows the same 1e-6.
- **Signs:** each MEM is made positive where its largest absolute value is
  (`lz_sign_index()`, as `pcoa()`).
- **n >= 3** is required. All samples at one place is an error (the
  default threshold would be 0). No MEM selected gives a 0-column matrix
  with a warning (e.g. `'positive'` for an equilateral triangle).
- **`values` and `moran_i` are named** by the MEMs.
- **vegan's `pcnm()` has no allowance for rounding at the threshold**:
  where computed distances fall just above a user threshold (mite.xy at 2
  or 0.8, where 2.0000000000000004 occurs), it truncates pairs that
  adespatial and ecodive keep, so the pcnm cross-check uses only cases
  without such distances (the generator checks this). pcnm's eigenvalues
  are 8 t^2 (n values + 1); with duplicate positions it returns fewer than
  n - 1 eigenvalues, so the tests compare its positive ones.
- **Fixtures** keep all n - 1 MEMs for mite, the 8 x 5 grid (default and
  threshold 1.5) and a transect, and only the positive ones for mite with
  threshold 2, mite's distances with threshold 0.8, and mite with
  duplicates, to keep the file at 128 KB. Where eigenvalues are equal (the
  grid), MEMs are compared as the projection onto their shared space.
- **The vignette example** is in `vignette("ordination")`, "Spatial
  structure", after `forward_select()`: simulated sites, `dbmem()`,
  `forward_select()` over the MEMs, then `varpart()` of pH against the
  selected MEMs (`attr(fs, 'selected')` is a formula, so `data` is
  `cbind(sites, mem)`).


## Outcome (from the implementing session's final report)

**Item 78, done 2026-11-26.** Commit `9325494` (function), plus the tick.

- Tests: 33,147 at home and 470 in CRAN mode (32,859 / 469 at 70dff27).
  `test-dbmem.r` has 287 tests.
- adespatial 0.3.30 `dbmem()` on mite.xy (default threshold, 2, 0.8 on
  its distances, and with 3 duplicate positions), an 8 x 5 grid (default
  and 1.5) and an irregular transect, each with `MEM.autocor` positive,
  negative, non-null and all: same number of MEMs and threshold; `values`
  within 4e-16; MEMs within 5e-11 up to sign where eigenvalues are distinct
  (projections onto shared spaces within 4e-13 on the grids, where they
  are equal); Moran's I (spdep 1.3.1 `moran()` with adespatial's weights)
  within 4e-15.
- vegan 2.6.4 `pcnm()`: eigenvalues 8 t^2 (n values + 1) within 1e-12 of
  the largest, its unit vectors within 1e-10, and the threshold equal to
  `spantree()`'s longest edge; but pcnm drops pairs that lie just above a
  user threshold by rounding (see the decisions).
- cpus 1, 2, 3, 4, 8 identical.
- Speed (uniform random points, 1 / 4 threads, reference BLAS): n = 500
  0.17 / 0.15 s (adespatial 0.56 s); n = 1000 1.5 / 1.0 s (2.6 s); n = 2000
  11.8 / 4.8 s (10.8 s).
- `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size).
