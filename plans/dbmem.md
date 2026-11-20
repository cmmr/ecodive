# `dbmem()`: distance-based Moran's eigenvector maps

**Priority:** 66. **Tier:** clean room. Write it from Borcard & Legendre
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


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- Not an ordination, so no `ord` record: a record of a new kind, attribute
  **`space`**: `ecodive_step("Moran's Eigenvector Maps", threshold =,
  moran_positive =, n_mem =)` (add the kind to `?ecodive_steps` when this
  lands). Functions that take `space` (`moran_i()`, constrained
  ordinations given MEMs as `data`) keep it when they can.
- Store what a plot needs: `values`, `moran_i` and `threshold` are
  already attributes; a `plot()` method (MEMs mapped over coordinates, or
  their Moran's I) titles from the record.

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
