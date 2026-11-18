# `moran_i()`: Moran's I test of spatial autocorrelation

**Priority:** 17. **Tier:** port from PySAL esda (BSD-3). Read
[README.md](README.md) first.

**Name:** `moran_i()`, after the statistic. A bare `moran()` would be easy
to confuse with `dbmem()`'s "Moran's eigenvector maps".


## What it is

Moran's I (Moran 1950) measures whether a numeric variable (e.g. a
taxon's abundance, or alpha diversity) is similar at nearby sites. It
needs spatial weights between samples. The test is analytic (under
normality or randomisation) or by permutation. `dbmem()` uses it to choose
which eigenvectors to keep ([dbmem.md](dbmem.md)).


## Design

- **Signature:** `moran_i(x, w, data = NULL, alternative = 'greater', permutations = 999L, seed = 0L, cpus = n_cpus())`.
  - `x` is a numeric vector, or a column name of `data`. `w` holds the
    spatial weights: a sample × sample matrix, a `dist` (converted to
    inverse-distance weights — decide and document), or coordinates (with
    a neighbour rule, e.g. k nearest).
  - Start with an explicit weight matrix, plus a documented helper for
    inverse distances. Row-standardise weights as esda does by default;
    `ape::Moran.I` doesn't. Document which, and offer an argument if both
    are needed for matching.
- **Result:** an `ecodive_test` of class `c('ecodive_moran_i', 'ecodive_test')`:
  - `statistic = c(I = ...)`, `p_value` (permutation), `z`, `n_perm`;
  - `expected` (−1/(n − 1)), `sd`, `p_value_norm` (analytic);
  - attribute `i_perm`.

  It is a plain function, not part of `dist_test()`, since it tests a
  variable against space, not distances against a variable.
- **Masking:** none (ape's is `Moran.I`; esda's is `Moran`).


## Reference implementations

- **Port from:** PySAL esda 2.10.0, `esda.Moran`, BSD-3,
  [repo](https://github.com/pysal/esda). Weights come from libpysal 4.15.1
  (BSD-3).
  - esda row-standardises weights by default (`transformation = 'r'`),
    and gives analytic (normality and randomisation) and permutation
    inference.
- **Output-only fixtures:**
  - ape `Moran.I()` (GPL-2|3), with an inverse-distance weight matrix:
    observed, expected, sd, p-value;
  - adespatial `moran.randtest()` for the permutation version.
- **Python fixture:** esda `Moran` (`data-raw/expected/moran_i.py`) with
  fixed weights, to check row standardisation and the analytic p-values.


## Tests

- I, E[I], sd and the analytic p-value match ape (raw weights) and esda
  (row-standardised weights) at 1e-10.
- With shared permutations, the permuted Is and p-value match exactly.
- Edge cases: a constant x (an error, or NA); all-zero weight rows.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page in a new "Spatial" group of the reference index, with
`dbmem()`; and NEWS.


## References

- Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
  *Biometrika*, 37, 17–23.
- Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models and
  Applications*. Pion.
