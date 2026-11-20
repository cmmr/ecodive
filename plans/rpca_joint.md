# `rpca_joint()`: Joint-RPCA across several tables

**Priority:** 30. **Tier:** port (gemelli 0.0.13,
`joint_rpca()`, BSD-3), on ecodive's OptSpace code from `rpca()`. Read
[README.md](README.md) first, and [compositional.md](compositional.md) for
the conventions the compositional plans share.


## What it is

Joint-RPCA (Martino et al., gemelli) ordinates the **same samples measured
by several 'omics** (16S, metagenomics, metabolomics, ...) at once: each
table is robust-clr transformed, and one low-rank model with **shared
sample loadings** and per-table feature loadings is fitted across all of
them by OptSpace-style gradient descent. Held-out test samples give a
cross-validated reconstruction error per iteration, which guides the
choice of rank. With one table it is gemelli's `rpca_with_cv()`.

**Name:** `rpca_joint()`, so it sorts with `rpca()` (rule 2); not
gemelli's `joint_rpca`. Docs aliases: "Joint-RPCA", "multi-omics RPCA".
**No wrapper:** it takes several tables (README keeps two-input methods as
plain functions).


## Design

### Signature

```r
rpca_joint(counts, margin = 1L, n_axes = 3L, max_iters = 5L,
           n_test_samples = 10L, test_samples = NULL, data = NULL,
           norm = 'rclr', min_sample_count = 0, min_feature_count = 0,
           min_feature_frequency = 0, distances = FALSE, cpus = n_cpus())
```

- `counts`: a **list** of count tables (any type `validate_counts()`
  accepts), named (the names label the tables in the output). All in the
  same orientation (`margin`). **Feature names must be unique across
  tables** (gemelli requires it); error otherwise, naming the clashes.
- `max_iters`: honoured here (gemelli's `joint_rpca` uses
  `max_iterations`, unlike `rpca`), with no early stopping, as gemelli.
- `n_test_samples`, `test_samples`, `data`: the held-out set.
  - `test_samples = NULL`: gemelli's automatic split (below).
  - Otherwise a character vector of sample names, a logical vector over
    the shared samples, or a column name of `data` holding `'train'` /
    `'test'` (gemelli's `train_test_column`).
- `norm`: `'rclr'` (default; gemelli's `rclr_transform_tables=True`) or
  `'none'` (tables already transformed: no filtering, NaN = missing).
  Reusing `norm` (rather than gemelli's boolean) keeps the glossary.
- Filters as `rpca()`, per table.

### Algorithm (gemelli 0.0.13, read from source)

1. Filter each table; keep the samples shared by all; filter again; rclr
   each.
2. **Split:** with `test_samples = NULL`, run plain RPCA on the **first**
   table, order samples by PC1, and take the samples at
   `round(linspace(0, n − 1, n_test_samples))` as the test set (spread
   along PC1). gemelli's **training-sample order** comes from Python set
   iteration (hash-randomised), so it isn't reproducible; ecodive uses the
   tables' sample order. Results should agree up to that ordering; fixtures
   use an explicit split.
3. **Joint solver:**
   - Test matrices are masked and centred (rows, then columns, over
     observed entries).
   - Training tables are stacked side by side and rescaled; the top k
     singular triplets (scipy `svds`, ascending, unsorted) initialise:
     \eqn{U \times \sqrt n}, \eqn{V \times \sqrt m}, S re-solved,
     \eqn{\rho = \epsilon n}. V is split into per-table blocks.
   - Each iteration, for each table: a gradient step and line search from
     the shared U on that table's **unrescaled** data; S re-solved with
     the **old** shared U. Then the test-set CV error: project the test
     data (test · V / diag(S)), reconstruct with S Vᵀ, RMS over observed
     entries.
   - Then \eqn{S_{shared}} = the top-k singular values of
     \eqn{\mathrm{mean}(U_i U_i^T)} (ascending, normalised to unit norm),
     \eqn{U_{shared} = \mathrm{mean}(U_i)}, each \eqn{V_i \leftarrow V_i
     S_{shared}}.
   - At the end, \eqn{S_{shared}} is divided by the rescaling factor and
     columns sorted descending.
   (These steps, including the "old U" and ascending orders, are gemelli's
   behaviour as read; port them as they are, and note them in a "Differences
   from the paper" comment, not in user docs.)
4. \eqn{X = U S [V_1; \dots; V_T]^T}, centred (columns, rows) and SVD'd:
   training samples get the scores; features get loadings stacked across
   tables; **test samples are projected** as
   (centred test rclr) · v / ‖s‖ (‖s‖ a single scalar, as gemelli).

### C

- Reuse `src/rpca.c`'s OptSpace kernels (sparse observed-entry layout,
  gradient, quartic line search, S solve).
- **Lanczos must return singular vectors** (ascending, as `svds`), not
  only the subspace it returns today for `rpca()`; the per-column scaling
  by \eqn{S_{shared}} depends on the vectors (not on their signs).
- \eqn{\mathrm{mean}(U_i U_i^T)} is n × n; never form it. Its non-zero
  eigenvalues are those of a small kT × kT matrix built from the
  \eqn{U_i}.
- Size estimate: 800–1,100 C lines, ~350 R lines.

### Result

`ecodive_ord` via `ord_result()` (`method = 'rpca_joint'`, prefix
`"RPCA"`), rows = **all** shared samples (training scores and projected
test samples), plus:

- `eig`, `prop` as `rpca()`.
- `loadings`: one row per feature across all tables (stacked), with
  attribute `table` (a factor naming each row's table), so
  `split(as.data.frame(loadings), attr(loadings, 'table'))` gives
  per-table loadings and `pcoa_biplot()`-style plotting works per table.
- `test`: logical, which rows were held out.
- `cv_error`: a data frame, one row per iteration: `iteration`, `mean`,
  `sd` (across tables), as gemelli's CV table.
- `n_iters`; `dm` when `distances = TRUE` (Euclidean over scores, which for
  gemelli includes the projected test samples and so differs from
  `rpca()`'s `dm`; say so).

Later, if asked: projecting **new** samples onto a fitted model (gemelli's
`transform()` / `rpca_transform()`), as a `predict()`-like function, and
gemelli's `feature_covariance_table()` (\eqn{V\,\mathrm{diag}(s)^2 V^T}).


## Reference implementations

- gemelli 0.0.13: `joint_rpca()`, `rpca_with_cv()` and their helpers in
  `gemelli/rpca.py`, `gemelli/optspace.py`, `gemelli/preprocessing.py`
  (BSD-3, already credited in `inst/COPYRIGHTS` for `rpca()`).
- **Fixtures** (`data-raw/expected/rpca_joint.py`): two and three random
  tables with partly shared samples, an **explicit** train/test column
  (gemelli's automatic split isn't reproducible), `n_axes` 2 and 3; and the
  one-table case (`rpca_with_cv`). Compare scores, loadings, eigenvalues
  and the CV table at about 1e-10 after `ord_flip()`.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The title comes from `ord_result()` automatically. Count tables have
  no title, so it is the method's name: "Joint Robust Aitchison PCA".
  The existing `plot.ecodive_ord()` uses it.
- A CV-curve plot type reads `cv_error`; loadings plots colour by the
  `table` attribute of `loadings` (both result attributes, above).
- Test: `display$title` for one call.

## Tests

- Fixtures as above.
- One table with no test samples behaves sensibly (CV error `NA`, or an
  error if gemelli requires test samples: match it).
- The automatic split picks the samples gemelli's rule picks given
  ecodive's ordering.
- Identical for `cpus` 1, 2 and 4.
- Errors: non-list `counts`; duplicate feature names across tables; no
  shared samples; `n_test_samples` ≥ number of shared samples; unknown
  `test_samples`; `norm` other than `'rclr'` / `'none'`.
- A smoke-test line.


## Docs

- Include the sequencing-depth caveat that [compositional_docs.md](compositional_docs.md)
  sets for the RPCA family (Schloss 2024, 2026), and cite the method's
  peer-reviewed paper: Cordazzo Vargas et al. (2026), *Molecular Systems
  Biology*, doi:10.1038/s44320-026-00236-3.
- Rd page under Ordination (a "Several tables or time series" subtitle,
  with `ctf()`), cross-linked with `rpca()`, `procrustes()` and the
  planned `coinertia()` (the classical alternatives for paired tables).
- NEWS.


## References

- Martino, C., et al. Joint-RPCA: gemelli documentation and tutorials
  (https://github.com/biocore/gemelli); cite the method's paper once
  published (check).
- Keshavan, R. H., Montanari, A., & Oh, S. (2010). Matrix completion from
  a few entries. *IEEE Transactions on Information Theory*, 56,
  2980–2998.
