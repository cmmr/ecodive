# `mantel_correlogram()`: Mantel correlogram

**Priority:** 26. **Tier:** *derived*, built from `mantel()` against
distance-class indicator matrices; no permissive code exists. Read
[README.md](README.md) first. Part of the `mantel_*` family:
`mantel()`, `mantel_partial()`, `mantel_correlogram()` and
`mantel_regression()`.


## What it is

A Mantel correlogram (Oden & Sokal 1986; Legendre & Legendre 2012) splits
geographic (or other) distances into classes. For each class, it computes
the Mantel correlation between the community distances and a 0/1 matrix
marking the pairs in that class, with a permutation test. It shows how
community similarity decays with distance. vegan `mantel.correlog()`;
ecodist `mgram()`.


## Design

- **Signature:** `mantel_correlogram(dm, y, data = NULL, breaks = NULL, n_classes = NULL, cutoff = TRUE, method = 'pearson', permutations = 999L, p_adjust = 'holm', strata = NULL, alternative = 'greater', seed = 0L, cpus = n_cpus())`.
  - `y` holds the distances to classify (a dist, a matrix, coordinates, or
    a column of `data`), handled as `mantel()`'s `y` is.
  - `alternative` (`'greater'`, `'less'` or `'two.sided'`) applies to each
    class's test, as in every other Mantel function. Check how it maps
    onto vegan's `mantel.correlog()` p-values when writing the fixtures,
    and document it.
  - **Classes:** the default is Sturges' rule, as in vegan.
  - **`cutoff`:** as in vegan, test only classes up to half the largest
    distance.
  - **Correction:** `p_adjust` defaults to Holm, consistent with the
    pairwise tests ([pairwise_tests.md](pairwise_tests.md)).
    - vegan's own default is a *progressive* Holm correction, across
      classes in order of distance. Offer `p_adjust = 'progressive'` too,
      to match vegan, and decide which is the default with the user.
- **Result:** a plain function. A data frame with one row per distance
  class:
  - columns `class_from`, `class_to`, `class_mid`, `n_pairs`, `r`,
    `p_value`, `p_adjusted`;
  - class `c('ecodive_mantel_correlogram', 'data.frame')`.
- **Masking:** none (vegan's is `mantel.correlog`).


## Reference implementations

- **No permissive code.** jwcarr's `mantel` (Python, MIT) does plain
  Mantel tests only, and skbio's `mantel` / `pwmantel` (BSD-3) have no
  correlogram.
- **Build** on ecodive's `mantel()` C kernel: each class is a Mantel test
  against an indicator matrix. Share one permutation set across classes
  (as vegan does). Parallelise over permutations.
- **Output-only fixtures (GPL):**
  - vegan 2.7-x `mantel.correlog(D.eco, D.geo, ...)`: `$mantel.res`,
    holding the class index, n.dist, Mantel r, Pr(Mantel) and Pr(corrected);
  - ecodist `mgram()` as a second check.


## Tests

- Fixtures from vegan, written by `data-raw/expected/mantel_correlogram.r`
  on oribatid mites (`mite` with `mite.xy`), with default and custom
  breaks. Compare the classes, n, r and raw and corrected p-values (with
  shared permutations).
- A class with no pairs gives `NA`, as vegan does.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

An Rd page next to `mantel()` in the reference index, and NEWS.


## References

- Oden, N. L., & Sokal, R. R. (1986). Directional autocorrelation: an
  extension of spatial correlograms to two dimensions. *Systematic
  Zoology*, 35, 608–617.
- Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
  §13.1.
- Sturges, H. A. (1926). The choice of a class interval. *Journal of the
  American Statistical Association*, 21, 65–66.
