# `alr()` and `ilr()`: log-ratio coordinates

**Priority:** 12. **Tier:** port (scikit-bio 0.7.4,
`skbio/stats/composition/_base.py`, BSD-3). Read [README.md](README.md)
first, and [compositional.md](compositional.md) for the conventions the
compositional plans share.


## What it is

Aitchison's log-ratio transforms map a composition of D parts to ordinary
real coordinates, where standard statistics apply. ecodive already has the
**clr** (and the robust **rclr**) as `norm` choices. These two complete
the set:

| Function | Transform | Output | Isometric? |
|:---------|:----------|:-------|:-----------|
| `alr()` | additive log ratio: \eqn{\log(x_j / x_{den})} for every part but the denominator | D − 1 coordinates | no |
| `ilr()` | isometric log ratio: \eqn{\mathrm{clr}(x)\,\Psi^T} for an orthonormal basis \eqn{\Psi} ((D − 1) × D, in clr coordinates) | D − 1 coordinates | yes |

- **ALR** is the easiest to interpret (everything relative to one chosen
  taxon) and is what multinomial regression models (e.g. Songbird) use.
  Distances on ALR coordinates aren't Aitchison distances.
- **ILR** coordinates are orthonormal, so Euclidean distance on them
  **equals the Aitchison distance** for any basis. `ilr()` is therefore
  not offered as a `norm` (it would repeat `aitchison()`); its value is
  feature-level: models, regression and `lda()` / `rda()` without the
  clr's singular covariance, and interpretable **balances** when the basis
  comes from a sequential binary partition (SBP).
- The phylogenetic ILR is `philr()` (see [philr.md](philr.md)); a single
  balance is `balance()` (see [log_ratio.md](log_ratio.md)).


## Names

- `alr()`, `ilr()`: the standard names in scikit-bio, compositions and
  robCompositions. Masking: `compositions::alr`, `compositions::ilr` and
  robCompositions' `addLR` / `pivotCoord` do the same thing with other
  defaults; say so in "Same names in other packages".
- No `clr()` (see compositional.md: it stays a `norm` choice).
- Inverses (`alr_inv`, `ilr_inv`) aren't planned. They matter for
  back-transforming model coefficients, which ecodive doesn't fit. Add them
  later as `alr_inv()` / `ilr_inv()` if asked.


## Design

### Signatures

```r
alr(counts, denominator, margin = 1L, pseudocount = NULL, cpus = n_cpus())
ilr(counts, basis = NULL, margin = 1L, pseudocount = NULL, cpus = n_cpus())
```

- **`denominator`** (alr): the reference part, as a feature name or a
  1-based index. **Required, no default.** scikit-bio defaults to the
  *first* part (`ref_idx=0`) and compositions to the *last*, so any default
  surprises someone, and the choice matters for interpretation. The error
  message says so and suggests a common, abundant, stable taxon. Same name
  as in `log_ratio()` and `rho_p()`.
- **`basis`** (ilr):
  - `NULL` (default): scikit-bio's `_gram_schmidt_basis(D)`, a Helmert
    basis. Row j (1-based) has \eqn{1/\sqrt{j(j+1)}} in columns 1..j,
    \eqn{-\sqrt{j/(j+1)}} in column j + 1, and 0 after: the balance of
    the first j parts against part j + 1. For D = 4:
    `[.7071068, -.7071068, 0, 0]`, `[.4082483, .4082483, -.8164966, 0]`,
    `[.2886751, .2886751, .2886751, -.8660254]`. It is the negative of the
    normalised transpose of `contr.helmert(D)`. Build it in skbio's
    arithmetic order (`(1/i) * sqrt(i/(i+1))`) to match to the last bit.
  - A **sequential binary partition**: a (D − 1) × D matrix of +1, −1
    and 0, one row per balance, +1 marking the numerator. Converted as
    scikit-bio's `sbp_basis()`: for a row with r entries of +1 and s of −1,
    +1 becomes \eqn{+\sqrt{s / (r(r+s))}} and −1 becomes
    \eqn{-\sqrt{r / (s(r+s))}}. Check that it is a valid SBP (each row has
    both signs; rows give an orthonormal basis), which scikit-bio doesn't.
  - An **orthonormal basis**: a numeric (D − 1) × D matrix in clr
    coordinates (scikit-bio's orientation: one row per coordinate), checked
    as scikit-bio does: \eqn{|\Psi\Psi^T - I| < 10^{-4} I + 10^{-6}}
    elementwise, else "Basis is not orthonormal.".
  - The two matrix forms can't be confused: an SBP has entries in
    {−1, 0, 1}, and no orthonormal clr basis does (a row of norm 1 with a
    single ±1 can't sum to 0). Still, check an SBP first, and document both.
  - Columns of `basis` are matched to the features **by name** when both
    have names (as `tree` tips are in `validate_tree()`), by position
    otherwise. Row names of `basis` name the output coordinates; otherwise
    they are `ilr1`, `ilr2`, ...
- **Zeros:** scikit-bio requires strictly positive input (`nozero=True`)
  and adds nothing. ecodive uses the shared `pseudocount` rule (see
  compositional.md), so `alr(x + pc)` here equals scikit-bio's
  `alr(x + pc)`. A `replace_zeros_*()` result needs none.
- Input needn't be closed: the logs of ratios are unchanged by scaling a
  sample.

### Computation

- One pass per sample: logs of (pseudocounted) values, then
  - alr: subtract the denominator's log; drop its column (the others keep
    their order, as scikit-bio);
  - ilr: clr, then multiply by \eqn{\Psi^T}. For the default Helmert basis
    use the O(D) recurrence (running sums of the logs), not a dense
    D × (D − 1) product. For a user basis, use BLAS (`dgemm`) on blocks of
    samples; SBP bases are sparse, so a sparse product is better for large
    D.
- Parallel over samples (blocks of samples per task; no shared
  accumulators), so results are identical for any `cpus`.
- Sparse input is read sparsely; the output is dense by nature.

### Result

- A dense base `matrix` in the **input's orientation** (`margin = 1`:
  samples × coordinates; `margin = 2`: coordinates × samples), with the
  sample names kept.
- Coordinate names: alr — the numerator features' names (scikit-bio keeps
  their order); ilr — `rownames(basis)`, or `ilr1`..`ilr<D-1>`.
- Attributes: `method` (`'alr'` / `'ilr'`), `denominator` (alr: the name
  or index), `basis` (ilr: the (D − 1) × D basis actually used, so the
  default is inspectable and `coords %*% basis` gives back the clr),
  `pseudocount`.
- These feed straight into `pca()`, `rda()`, `lda()`, `envfit()`, and
  `euclidean()` / `beta_div(norm = 'none')`, which all accept negative
  values when `norm = 'none'`.
- Warn above 1e8 values, as `normalize(norm = 'clr')` does.


## Reference implementations

- scikit-bio 0.7.4: `alr()`, `ilr()`, `_gram_schmidt_basis()`,
  `sbp_basis()`, `_check_basis()` (BSD-3; ecodive's `inst/COPYRIGHTS`
  already has scikit-bio stanzas for `pcoa_biplot` and `bioenv`; add the
  new files to the same form). The Python package is a directory since
  0.7: `skbio/stats/composition/_base.py`.
- Fixtures: scikit-bio's own test values, e.g. `alr([.1,.3,.4,.2])` =
  `[1.0986123, 1.3862944, 0.6931472]` (denominator 1),
  `ilr([.1,.3,.4,.2])` = `[-0.7768362, -0.68339802, 0.11704769]`,
  `ilr(closure([e,1,1]))` = `[0.70710678, 0.40824829]`; plus generated
  fixtures from `data-raw/expected/log_ratio_transforms.py` (random tables,
  default basis, an SBP, a user basis).
- Output-only cross-checks: `compositions::alr` (last part as
  denominator) and `compositions::ilr` (its default basis has the opposite
  sign, reportedly; check and document).


## Tests

- Fixtures at 1e-12.
- Identities: Euclidean distances on `ilr()` coordinates equal
  `aitchison()` (any basis, same pseudocount); `ilr()` %*% `basis` equals
  `normalize(x, 'clr')`; an SBP basis and its `sbp_basis()` matrix give the
  same result; `alr()` equals `log(x[, -k] / x[, k])`.
- Scale invariance per sample; `margin = 2` gives the transpose.
- Identical for `cpus` 1, 2 and 4; sparse and dense agree.
- Errors: missing `denominator`; unknown denominator name; a basis of the
  wrong shape, not orthonormal, or an invalid SBP; mismatched feature
  names; negative counts.
- Smoke-test lines.


## Docs

- One Rd page per function, cross-linked with `normalize()` (clr/rclr),
  `philr()`, `balance()`, `aitchison()`.
- "Compositional Data" group, subtitle "Log-ratio transforms".
- Say plainly: for distances use `aitchison()`; use `ilr()` for
  coordinate-wise analyses.
- NEWS.


## References

- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
  Barceló-Vidal, C. (2003). Isometric logratio transformations for
  compositional data analysis. *Mathematical Geology*, 35, 279–300.
- Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
  balances in compositional data analysis. *Mathematical Geology*, 37,
  795–828.
