# `log_ratio()` and `balance()`: one log ratio per sample

**Priority:** 16. **Tier:** port (Qurro, BSD-3, for
`log_ratio()`; scikit-bio / gneiss, BSD-3, for the balance formula).
Read [README.md](README.md) first, and [compositional.md](compositional.md)
for the conventions the compositional plans share.


## What it is

From relative abundances alone, the quantity that can be interpreted
without assumptions about total load is a **ratio between taxa**: if the
ratio of A to B goes up, it went up whatever happened to the total. This is
the "reference frames" argument of Morton et al. (2019), and the end point
of most current compositional advice: pick (or learn) a numerator and a
denominator, and analyse the log ratio as an ordinary per-sample variable.

There are two standard ways to combine several taxa into one ratio, so
there are two functions:

| Function | Ratio | Used by |
|:---------|:------|:--------|
| `log_ratio()` | \eqn{\log(\sum_{i \in N} x_i / \sum_{j \in D} x_j)}: log of the ratio of **summed** counts | Qurro (Fedarko et al. 2020), Morton et al. (2019); robust to zeros in single taxa |
| `balance()` | \eqn{\sqrt{rs/(r+s)}\,\log(g(x_N) / g(x_D))}: the scaled log ratio of **geometric means** (an ILR balance) | Egozcue & Pawlowsky-Glahn (2005), selbal, gneiss, coda4microbiome, PhILR |

where \eqn{N} (r taxa) and \eqn{D} (s taxa) are disjoint, and \eqn{g} is
the geometric mean. A balance is one coordinate of an `ilr()` with a
sequential binary partition (see [log_ratio_transforms.md](log_ratio_transforms.md));
`balance()` gives that one coordinate without building a full basis.

The README's rule (one function per named variant) applies: the two are
different quantities with different zero behaviour, and the literature
names them differently ("log ratio" vs "balance").


## Names

- `log_ratio()`: the generic name, as Qurro uses it.
- `balance()`: the CoDa term. Docs: "an isometric log-ratio balance, not a
  balanced design". The alias "ilr balance" is for the docs.
- Neither is a wrapper method.


## Design

### Signatures

```r
log_ratio(counts, numerator, denominator, margin = 1L, pseudocount = 0)
balance  (counts, numerator, denominator, margin = 1L, pseudocount = NULL)
```

- **`numerator`, `denominator`:** the features on each side: a character
  vector of feature names, a vector of feature indices (1-based), or a
  logical vector with one value per feature. Each must select at least one
  feature, and they must not overlap (an error names the shared features).
  The same names are used by `alr(denominator =)` (see
  [log_ratio_transforms.md](log_ratio_transforms.md)), so "denominator" means
  the same thing everywhere.
- No `cpus`: this is one pass over the selected features, cheap even for
  large tables. Keep sparse input sparse: sum the selected columns of a
  `dgCMatrix` directly rather than densifying.

### Zeros

- **`log_ratio()`:** zeros in single taxa are harmless as long as each
  sum is positive. A sample whose numerator or denominator sum is 0 gets
  `NA` (its ratio is undefined), with **one** warning naming how many (up
  to 5 names, via `name_list()` in `R/design.r`). Qurro drops such samples
  instead; ecodive keeps their place so the result lines up with `counts`.
  - `pseudocount = 0` by default. A positive `pseudocount` is added to
    every selected count before summing (so to each sum, `r * pseudocount`
    and `s * pseudocount`), as in `normalize(norm = 'clr')`.
  - The default differs from the other compositional functions' `NULL`
    because summing removes the need for one in almost every case; say so.
- **`balance()`:** a zero anywhere in N or D makes the geometric mean 0.
  `pseudocount = NULL` follows the shared rule of `validate_pseudocount()`
  (half the smallest non-zero value, with its warning), applied to the
  **whole** table, so a balance agrees with the matching `ilr()`
  coordinate. A zero-free table, or one passed through `replace_zeros()`,
  needs none.

### Result (both)

- A numeric vector, one value per sample, named after the samples (or
  unnamed if `counts` has none), in sample order. This drops straight into
  metadata (`data$ratio <- log_ratio(...)`) and into any test as `groups`
  (e.g. `permanova(dm, ratio)`, `envfit()`), with no conversion.
- Attributes: `method` (`'log_ratio'` or `'balance'`), `numerator` and
  `denominator` (the feature names, or indices if unnamed), and
  `pseudocount` (the value used). Natural log, as everywhere in ecodive.
- Orientation: `log(N / D)`, so a larger value means more of the
  numerator.


## Reference implementations

- **`log_ratio()`:** Qurro (biocore/qurro, BSD-3; verify the license file
  before reading). Its sample-plot log ratio is `log(sum(num) / sum(den))`
  and samples where either sum is 0 are excluded. No C needed.
- **`balance()`:** the scikit-bio / gneiss balance basis: for a partition
  with r taxa at +1 and s at −1, the basis row is
  \eqn{+\sqrt{s / (r(r+s))}} on the r taxa and
  \eqn{-\sqrt{r / (s(r+s))}} on the s taxa, and the coordinate is its
  inner product with the clr. That equals the closed form above.
  Fixtures: `skbio.stats.composition.ilr(clr_input, basis=sbp_basis(sbp))`
  with a one-row SBP padded to a full basis, or computed directly in
  numpy. scikit-bio 0.7.4's `sbp_basis()` uses exactly the row formula
  above, with +1 marking the numerator (see
  [log_ratio_transforms.md](log_ratio_transforms.md)).
- **Output-only cross-checks:** selbal and coda4microbiome balances (check
  their licenses; GPL means output-only).


## Tests

- `log_ratio()` equals `log(rowSums(x[, N]) / rowSums(x[, D]))` for
  dense, sparse and both margins; `NA` and one warning for a zero sum;
  names, indices and logical selectors agree.
- `balance()` equals the closed form; equals the matching coordinate of
  `ilr()` with an SBP containing that partition (at 1e-12); its sign flips
  when N and D are swapped; it is unchanged by multiplying a sample by a
  constant (scale invariance), as is `log_ratio()`.
- Errors: empty N or D, overlapping N and D, unknown names, out-of-range
  indices, a logical selector of the wrong length.
- Smoke-test lines.


## Docs

- One Rd page each, cross-linked, in the "Compositional Data" reference
  group, subtitle "Log ratios" (see [compositional.md](compositional.md)).
- Explain when to prefer which: sums are robust to zeros and easy to
  explain; balances are the coordinates the CoDa theory is built on, and
  have the right scale for comparing balances of different sizes.
- Point to Qurro for choosing ratios interactively, and to `vlr()` /
  `rho_p()` for finding taxa that move together.
- NEWS.


## References

- Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
  balances in compositional data analysis. *Mathematical Geology*, 37,
  795–828.
- Morton, J. T., Marotz, C., Washburne, A., et al. (2019). Establishing
  microbial composition measurement standards with reference frames.
  *Nature Communications*, 10, 2719.
- Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
  'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
  Bioinformatics*, 2, lqaa023.
- Rivera-Pinto, J., Egozcue, J. J., Pawlowsky-Glahn, V., et al. (2018).
  Balances: a new perspective for microbiome analysis. *mSystems*, 3,
  e00053-18.


## Decisions made while implementing (item 16, 2026-11-19)

Asked of the user at the start (all four answered with the recommended
option):

1. **`log_ratio(pseudocount = 0)`** stays the default, unlike the other
   compositional functions' `NULL`: summing removes the need in almost
   every case. A sample whose numerator or denominator sum is 0 gets `NA`
   with one warning naming up to 5 samples (indices when unnamed): "1
   sample ('s1') has a numerator or denominator that sums to 0, so its log
   ratio is NA. A positive `pseudocount` gives it a value." `pseudocount =
   NULL` is accepted too and applies the shared rule to the whole table.
2. **Mixed selector forms** are allowed: each side is resolved on its own
   (`numerator = c('a', 'b')`, `denominator = 5:7`).
3. **`margin = 2`** returns the same named vector (one value per sample,
   named after the columns): a vector has no orientation.
4. **`balance(pseudocount = 0)`** is accepted (no pseudocount): a sample
   with a zero in a selected feature gets `NA` with one warning ("... has a
   zero in a numerator or denominator feature, so its balance is NA. Use
   `pseudocount = NULL`, or replace_zeros() first, to give it a value.").
   Zeros only in other features don't matter: the clr is then undefined,
   but the balance is computed from the uncentred logs (the basis row sums
   to 0, so centring cancels). `log_pseudocount()` gained `zero_ok = TRUE`
   for this ("`pseudocount` must be a single non-negative number.");
   `alr()`, `ilr()` and `norm = 'clr'` still refuse 0.

Made while implementing:

5. **Selectors: `feature_set(x, features, D, arg)`** in `R/log_ratio.r`,
   beside `feature_index()` (unchanged, so `alr()` accepts exactly what it
   did). Names, 1-based whole-number indices or a logical vector of length
   D; factors become character. Errors: missing argument ("`numerator` is
   required"), empty or all-FALSE, NA, a non-atomic value, unknown names
   (up to 5 listed), indices outside 1..D (negative indices are not
   exclusions), duplicates, names when `counts` has none, a logical of the
   wrong length, and overlap ("must not share features, but they share 2
   features: 'f2', 'f3'."). The `numerator` / `denominator` attributes are
   the features **in the order of `counts`** (not the order given), names
   or integer indices when unnamed.
6. **C, not R, single-threaded, no `cpus`** (`src/log_ratio.c`:
   `C_log_ratio()`, `C_balance()`), so every table type goes through
   `new_ecomatrix()` and sparse tables stay sparse. `log_ratio_input()`
   validates `cpus` only when the caller has one, and takes `dense_out =
   FALSE` to skip the large-dense-result warning.
   - `log_ratio`: `log(sum_N(x + c) / sum_D(x + c))`, sums in feature
     order, Qurro's `log(Num_Sum / Denom_Sum)`.
   - `balance`: the SBP row from `sbp_basis()` (as `ilr_basis()` stores
     it), dotted with the clr of **all** D features computed exactly as
     `ilr()`'s worker does, so `balance()` is **bit-identical** to the
     matching `ilr()` coordinate (tested on every fixture case and after
     `replace_zeros()`), and its sign flip is exact.
   - Sparse 2,000 × 5,000 (10% dense), 1,500 selected features: 0.04 s
     (log_ratio), 0.07 s (balance); dense input ~0.3 s, mostly validation
     and conversion.
7. **References** (`data-raw/expected/log_ratio.py` / `.r`, 20 cases over
   5 tables: counts with zeros, scaled, sparse, wide (D = 60) and Dirichlet
   proportions; two features, small groups, every feature, random sizes;
   hex floats):
   - scikit-bio 0.7.4 `ilr(x + pc, sbp_basis(sbp))[:, 0]`, the SBP's first
     row being the balance and the rest completing it: max abs difference
     **5.6e-16**.
   - The closed form in numpy: 2.4e-15.
   - **Qurro 0.9.0 `qarcoal()`** (BSD-3, licence checked in the PyPI sdist;
     qiime2 stubbed, being only a type annotation): **bit-identical**,
     including `NA` exactly where Qurro drops samples.
   - `np.log(sum(x_N + 1) / sum(x_D + 1))`: 2.2e-16. Base R's
     `log(rowSums / rowSums)`: within 1e-14 (rowSums uses long doubles).
   - compositions 2.0.9 `balance(acomp(x), ~ (F1*F2)/(F3))` (GPL, output
     only): 1.3e-15, same sign; stored in the fixture and tested.
8. **Nothing is ported from Qurro** (its formula is one line), so Qurro
   gets no `Authors@R` entry, file header or `inst/COPYRIGHTS` stanza; the
   docs cite the paper. The balance's basis row is scikit-bio's
   `sbp_basis()`, already credited: `src/log_ratio.c` carries the BSD-3
   text, `R/log_ratio.r` points to it, the existing scikit-bio stanza now
   lists both files, and the `Authors@R` comment mentions balances.
9. **selbal and coda4microbiome were not run.** selbal's repository
   (UVic-omics/selbal) wasn't reachable from this session, so its licence
   is unchecked (treat it as output-only); selbal selects balances rather
   than computing a given one. coda4microbiome is MIT (CRAN DESCRIPTION)
   but needs Bioconductor's ComplexHeatmap, and its log-contrasts are
   weighted `1/k` sums, a different scale. compositions' `balance()` covers
   the cross-check instead.
10. **A vector** is one sample (as `validate_counts()` makes it): a
    length-1 result, unnamed.
11. **Docs:** `log_ratio.Rd` (Zeros section; when to prefer which),
    `balance.Rd` (Pseudocount section; "Same names": `compositions::
    balance()`); pkgdown "Compositional Data" → "Log ratios"; NEWS;
    `parameter_names.md` updated. `vlr()` / `rho_p()` are mentioned as
    plain code, not links, until they exist.

For later items:

- **`vlr()` / `rho_p()` / `phi_s()` (`denominator`, `NULL` = the clr's
  geometric mean):** use `feature_set()` for a multi-feature denominator,
  and `log_ratio_input(what, zero_ok = , dense_out = )` for validation
  (`cpus` is validated only when the function has one). Note that
  `log_pseudocount(zero_ok = TRUE)` returns 0 for an explicit 0 without
  checking for zeros: the caller must handle them.
- **`philr()`:** a balance is the dot product of one CSR row with the clr;
  `C_balance()` takes any (feature, coefficient) row, so a weighted PhILR
  row can be checked against it, and `C_log_ratio_coords()` remains the
  fast path for all D − 1 balances.


## Outcome (from the implementing session's final report)

**Item 16, done 2026-11-19.** Commits `c56ea2a` (work), `324a969` (tick).

- Tests: 9,942 at home, 228 in CRAN mode.
- scikit-bio 0.7.4 `ilr()` with a one-row SBP: 5.6e-16. Qurro 0.9.0
  `qarcoal()`: bit-identical, `NA` exactly where Qurro drops a sample.
  compositions 2.0.9 `balance()`: 1.3e-15. ecodive's `ilr()` with a
  matching SBP: bit-identical.
- Speed: 2,000 × 5,000 sparse: `log_ratio()` 0.04 s, `balance()` 0.07 s.
- Nothing ported from Qurro (one-line formula), so it isn't credited.
  selbal and coda4microbiome not run.
