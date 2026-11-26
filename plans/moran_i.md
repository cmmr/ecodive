# `moran_i()`: Moran's I test of spatial autocorrelation

**Priority:** 59. **Tier:** port from PySAL esda (BSD-3). Read
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

- **Signature:** `moran_i(x, space, data = NULL, alternative = 'greater', permutations = 999L, seed = 0L, cpus = n_cpus())`.
  - `x` is a numeric vector, or a column name of `data`.
  - `space` gives the samples' positions, in any of three forms (the same
    argument name as `dbmem()`'s): spatial weights as a sample × sample
    matrix; geographic distances as a `dist` (converted to
    inverse-distance weights — decide and document); or a coordinate
    matrix (with a neighbour rule, e.g. k nearest). It isn't called
    `weights`, which means group weights in `mrpp()`.
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- An `ecodive_test`, so `new_ecodive_test()` sets the title. `x` is a
  numeric vector with no title, so it is the method's name alone:
  "Moran's I".
- Store whether `space` was coordinates, distances or weights as
  `space_type` in the result: it says how the weights behind I were made.
- `plot()` draws the permutation null, as for other tests, from `i_perm`.

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


## Decisions made while implementing

Asked of the user at the start (2026-11-25), each the recommended option:

- **`style`** chooses the transformation of the weights: `'row'` (the
  default, esda's and spdep's), `'binary'`, `'variance'` (variance-
  stabilising) and `'none'`, as libpysal's R, B, V and O. spdep calls it
  `style`; `standardize` and `transform` already mean other things
  (`linda()`, `maaslin3()`). esda's "D" (doubly standardised) is left out:
  a global rescaling gives the same I and moments as `'none'`. The
  signature is `moran_i(x, space, data = NULL, style = 'row', alternative =
  'greater', permutations = 999L, seed = 0L, cpus = n_cpus())`.
- **Distances and coordinates become weights `1 / d`** and nothing else: no
  `power`, no k-nearest rule, no exported helper. Other rules are a weight
  matrix (the docs' example builds k-nearest-neighbour weights). A zero
  distance between two samples is an error that names them.
- **Both analytic variances** are kept: `expected`, `sd_norm`,
  `p_value_norm`, `sd_rand`, `p_value_rand` (the plan had `sd` and
  `p_value_norm`). `p_value` and `z` stay the permutation ones.
- **The diagonal of a weight matrix is ignored** (no self-weights), so
  `1 / as.matrix(d)`, with `Inf` on its diagonal, works as is.

Made while implementing:

- **`space`'s form is read from its class:** a `dist` is distances; a
  square matrix (numeric or logical) is weights; a data frame or a
  non-square matrix is coordinates (Euclidean distances via
  `stats::dist()`); a character vector names coordinate columns of `data`.
  A square coordinate matrix (as many coordinates as samples) must be given
  as a data frame. Samples are matched to `x` and `data` by `space`'s names
  (`per_sample()` and `match_rows()`, with `space` as the reference).
- **Only the packed lower triangle of the symmetric part `(W + W') / 2`**
  is kept after the transformation, since `z'Wz = z'((W + W')/2)z`; the
  transformation and libpysal's `s0`, `s1`, `s2` are computed in one C pass
  over pairs (`C_moran_weights()`), from either a full matrix or packed
  distances, so a `dist` never becomes an n × n matrix. Each permutation
  is one fixed-order pass over the triangle (`C_moran_perm()`), blocks of
  8 permutations per task, `perm_source_get()` permutations (gather: sample
  i takes the value of sample p[i], as esda's `z[perm]`).
- **Weights must be finite and not negative** off the diagonal; weights all
  0 is an error.
- **Islands** (samples whose own weights are all 0) are kept, with a
  warning, and count in n and the mean, as in esda (and spdep with
  `adjust.n = FALSE`).
- **A constant `x`** gives NA statistics and p-values with a warning, as
  `mantel()` does for constant distances.
- **n ≥ 3 is required**; `sd_rand` is NA at n = 3 (its formula divides by
  n - 3).
- **A variance that is 0 to rounding** (below 1e-10 of its leading term,
  e.g. binary weights between every pair, where I is the same for every
  arrangement) gives NA sds and analytic p-values, where esda gives 0 and
  p = 0.
- **`alternative` sets the tail of the analytic p-values too.** esda's are
  two-tailed by default, or one-tailed in the direction of z; its `p_sim`
  is the smaller tail. The tests map each onto `alternative`.
- **ape's `Moran.I()` is not "raw weights, normality"** as the plan said:
  its outputs match `style = 'row'` with the randomisation sd (to about
  1e-15), so the docs and tests say so. spdep's `mat2listw(style = 'B')`
  keeps the matrix's own values, so the fixture gives it the binary matrix.
- **No `strata`** (the plan's signature has none).
- **adespatial's `moran.randtest()` was not run:** it is permutation-only,
  and with its own random permutations it can agree only to sampling error;
  the permuted statistics are checked exactly against esda's given the same
  permutations instead.
- **For `dbmem()`:** `moran_space()`, `moran_weights()`, `moran_stat()`
  (with `perm_matrix(0L, n)`, the observed I alone) and `moran_moments()`
  are internal helpers it can call; see the Outcome section.


## Outcome (from the implementing session's final report)

**Item 59, done 2026-11-26.** Commit `f548a26` (function), merged with
items 57 and 58 in `5bf8854`, plus the tick.

- Tests: 31,409 at home and 457 in CRAN mode before the merge (30,895 /
  450 at 15dc935); 32,351 and 467 on the merged tree. `test-moran_i.r`
  has 507 tests.
- esda 2.10.0 `Moran` (libpysal 4.15.1), 5 cases (inverse-distance,
  asymmetric k-nearest inverse-distance, an island, a negatively
  autocorrelated grid, n = 4) × 4 transformations (r, b, v, o): I within
  9e-16 (relative), sd under normality and randomisation within 2e-15,
  analytic p-values within 1.5e-14, s0, s1, s2 within 1.4e-15, and the 199
  permuted I per case within 4e-15 of esda's, given the same permutations;
  `p_sim` identical wherever no permuted I ties with the observed one.
- spdep 1.3.1 `moran.test()` (styles W, B, S, U; randomisation and
  normality; `adjust.n = FALSE`): I, variance and p-value within ~2e-14.
- ape 5.7.1 `Moran.I()`: matches `style = 'row'` with `sd_rand` (not raw
  weights with the normality sd, as this plan assumed), within ~1e-15.
- cpus 1, 2, 3, 4, 8 identical.
- Speed (n = 1000, 999 permutations, inverse-distance weights): ecodive
  0.38 s on 1 thread, 0.11 s on 4; esda 1.2 s; spdep `moran.mc()` 13 s.
- `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size).
- **For `dbmem()` (item 78):** internal helpers in `R/moran_i.r`:
  `moran_space(space, data)` (a `dist`, coordinates or weights matrix →
  raw weights, `n`, labels, type), `moran_weights(sp, style)` (the
  transformed weights' packed symmetric part `S` and `s0`, `s1`, `s2`, via
  `C_moran_weights()`), `moran_stat(x, wt, perm_matrix(0L, n))` (I alone;
  one O(n²) pass per variable) and `moran_moments(x, wt, n)`. dbmem can
  build its weights matrix (adespatial's choice of weights, to be checked
  against its outputs) and call `moran_weights()` once, then `moran_stat()`
  per MEM. `space`'s parsing (`moran_space()`) also gives dbmem its
  `dist`/coordinate handling.
