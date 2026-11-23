# `alr()` and `ilr()`: log-ratio coordinates

**Priority:** 14. **Tier:** port (scikit-bio 0.7.4,
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


## Decisions made while implementing (2026-11-19)

Asked of the user at the start (all four answered with the recommended
option):

1. **Zeros with `pseudocount = NULL`:** the shared rule (half the smallest
   non-zero value, with a warning; none when there are no zeros). The rule
   now lives in `log_pseudocount(counts, pseudocount, what, constant_ok)` in
   `R/validate.r`; `validate_pseudocount()`'s clr branch calls it, with its
   messages unchanged. `alr()` / `ilr()` pass `constant_ok = TRUE`: a table
   of one repeated value gives coordinates of 0, where `norm = 'clr'` still
   refuses it. The pseudocount is added to every value, as `norm = 'clr'`
   does, so `alr(x, pseudocount = c)` is scikit-bio's `alr(x + c)`.
2. **Orientation:** the input's (`margin = 2` gives coordinates ×
   samples). A plain vector is one sample in a column (`validate_counts()`
   sets `margin = 2`), so it gives a (D − 1) × 1 matrix.
3. **Coordinate names:** alr, the numerator features' names, or `alr<k>`
   numbered by the feature's position when `counts` has no feature names
   (denominator 2 of 4 gives `alr1`, `alr3`, `alr4`); ilr,
   `rownames(basis)`, or `ilr1` .. `ilr<D-1>`.
4. **Trees:** `ilr()` takes no tree; that is `philr()`'s job.

Made while implementing:

5. **`denominator`** is a name, a 1-based index, or a logical vector with
   one `TRUE` (`feature_index()`). The `denominator` attribute is the name,
   or the integer index when the features have no names. No warning when
   the denominator has zeros (the pseudocount warning covers zeros).
6. **Basis checks.** A matrix whose entries are all in {−1, 0, 1} is an SBP
   (checked first); each SBP row needs both signs ("row k has no 1 / -1").
   Then, for any basis, scikit-bio's orthonormality test (|ΨΨᵀ − I| <
   1e-4 I + 1e-6), and **one check scikit-bio lacks**: every row must sum
   to 0 (|row sum| < 1e-6 √D), i.e. lie in clr space; without it the
   coordinates aren't isometric. An orthogonal but non-hierarchical
   partition (Hadamard-like) passes: it still gives a valid orthonormal
   basis. A failing SBP says "not a valid sequential binary partition: its
   balances are not orthogonal". A `data.frame` or Matrix-package basis is
   coerced with `as.matrix()`. A D × (D − 1) basis gets a "Transposed?"
   hint. Columns match features by name when both have names (same set, no
   duplicates, else an error naming up to 5 missing / extra), by position
   otherwise.
7. **The `basis` attribute** is the (D − 1) × D matrix used, columns in the
   features' order with their names, rows named as the coordinates. For the
   **default basis with D > 1000** (more than 1e6 values, 8 MB) it is the
   string `'helmert'` instead, so that a 5,000-taxon table doesn't carry a
   200 MB attribute. User bases and SBPs are always stored (the user's
   input was already that size).
8. **No BLAS.** `src/log_ratio_transforms.c` (`C_log_ratio_coords`) does one
   sample per task with a per-thread dense scratch vector: alr subtracts
   logs; ilr computes the clr (sequential sum, as `normalize()`), then
   either the Helmert recurrence (running sum of the clr values,
   `a_j * S_j + b_j * c_{j+1}`, with skbio's entries) or a loop over the
   basis's non-zeros by row (CSR, built in R by `ilr_basis()`), so an SBP
   costs its non-zeros and no BLAS call runs inside worker threads.
   Results are identical for any `cpus`. 2,000 × 5,000 sparse table: 0.2 s
   (vs 30 s for `normalize(x, 'clr') %*% t(basis)` on reference BLAS).
9. **Parity.** `helmert_basis(D)` and `sbp_basis()` reproduce
   `_gram_schmidt_basis(D)` and `sbp_basis()` **bit for bit** (same
   arithmetic: `(1 / i) * sqrt(i / (i + 1))`; `sbp * sqrt((s / r) ^ sbp /
   sum(abs(sbp)))`, where R's `^` and numpy's `**` both call `pow()`). `alr()`
   is **bit-identical** to scikit-bio's `alr()`. `ilr()` is not bit-identical
   (numpy's pairwise mean and BLAS dot products vs sequential sums): max
   absolute difference 1.3e-15 over the fixtures (default, SBP and random
   orthonormal bases, D up to 60). Fixtures use hexadecimal floats
   (`float.hex()`, read exactly by `as.numeric()` in R).
10. **compositions 2.0.9 (output only):** `compositions::alr()` = `alr(x,
    D)` (last part) and `compositions::ilr()` = `-ilr(x)` (the default
    basis is skbio's negated, **same coordinate order**, not reversed),
    both to 1.3e-15; stored in the fixture and tested.
11. **`sbp_basis()` stays internal** (README's function list doesn't have
    it); the `basis` attribute shows the converted basis.
12. **Large results:** warn above 1e8 output values, as `normalize(norm =
    'clr')` does.

For later compositional items, the reusable internals (all in
`R/log_ratio_transforms.r` unless noted):

- `log_pseudocount()` (`R/validate.r`): the shared zero rule; call it with
  `what` naming the method and `constant_ok = TRUE`.
- `log_ratio_input(what)`: validates a count table in the caller's frame
  (counts, margin, cpus, finiteness, ≥ 2 features, pseudocount, large-output
  warning) and returns `D`, `n`, `features`, `samples`.
- `feature_index(x, features, D, arg)`: one feature by name, index or
  logical. `log_ratio()` / `balance()` need a multi-feature version (sets
  of features); write it beside this one.
- `helmert_basis(D)`, `sbp_basis(sbp)` (any number of rows, so a one-row
  SBP gives one balance vector for `balance()`), and `ilr_basis(basis, D,
  features)` (checks, name matching, CSR) returning `basis`, `names`, `csr`,
  `sbp`.
- `C_log_ratio_coords(counts, margin, den, csr, pseudocount, cpus)`: alr
  when `den >= 0` (0-based), else ilr with a CSR basis (`list(p, i, x)`,
  rows of any count) or `NULL` for Helmert. `philr()` can feed it a
  tree-derived SBP's CSR (O(nnz) per sample, no dense basis) if it doesn't
  write its own O(D) recurrence; `balance()` can pass a one-row CSR.


## Outcome (from the implementing session's final report)

**Item 14, done 2026-11-19.** Commits `51fd7d0` (work), `2488e23` (tick),
`9396845` (merge).

- Tests: 9,367 at home, 217 in CRAN mode.
- scikit-bio 0.7.4: `alr()`, the default (Helmert) basis and `sbp_basis()`
  bit-identical; `ilr()` ≤ 1.3e-15. compositions 2.0.9: `alr(x, D)` and
  `-ilr(x)` equal its `alr()` / `ilr()` to 1.3e-15.
- robCompositions not installable (dependency left CRAN).
- Speed: 2,000 × 5,000 sparse in 0.2 s.
