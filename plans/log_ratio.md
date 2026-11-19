# `log_ratio()` and `balance()`: one log ratio per sample

**Priority:** 14. **Tier:** port (Qurro, BSD-3, for
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
