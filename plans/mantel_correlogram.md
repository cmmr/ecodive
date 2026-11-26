# `mantel_correlogram()`: Mantel correlogram

**Priority:** 58. **Tier:** *derived*, built from `mantel()` against
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The result is a data frame, not an `ecodive_test`, so no constructor
  titles it. Set `display$title` with `set_display()` and `title_of()`,
  naming both matrices as `mantel()` does when both have titles: "Mantel
  correlogram of Bray-Curtis Dissimilarity and Euclidean Distance"; with
  a title on `dm` only (`y` coordinates or a column of `data`), "Mantel
  correlogram of Bray-Curtis Dissimilarity".
- Store what the plot needs: the classes' breaks and midpoints, the
  per-class statistics and `p_adjusted` in plotting order, `alpha` for
  filled vs open points.

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


## Decisions made while implementing (2026-11-25)

vegan's `mantel.correlog()` was characterised from its documentation and
outputs only (vegan 2.7-8; its source was not read). What it does, as
`mantel_correlogram()` now does:

- **Classes:** Sturges' rule is `ceiling(log2(N) + 1)` for N pairs (R's
  `nclass.Sturges()`), with breaks `seq(min, max, length.out = K + 1)` over
  y's distances. A class holds the distances in `(from, to]`.
- **`cutoff`** is not "half the largest distance": in 300 random cases,
  vegan tests class k exactly when it has pairs and either
  `k <= floor(K / 2)` or its pairs include every sample. A class with all N
  pairs isn't tested either (its model matrix is constant).
- **Sign:** the model matrix is 0 for pairs in the class and 1 for the
  others, so a positive r means positive autocorrelation.
- **Permutations:** vegan runs `mantel()` with the model matrix first, so it
  permutes the classes (y's samples) and holds `D.eco` fixed;
  `mantel_correlogram()` does the same, which makes a supplied permutation
  matrix give vegan's p-values exactly. vegan has no argument for a
  permutation matrix (and draws new permutations for each class), so the
  fixture script makes vegan's internal `getPermuteMatrix()` return a fixed
  one via `assignInNamespace()`.
- **Direction:** each class is tested one-tailed in the direction of its r:
  `r > 0` upper tail, `r <= 0` lower tail, both with vegan's
  `sqrt(.Machine$double.eps)` tie tolerance (confirmed on tied data).
- **Progressive correction:** tested class j (in order of distance) is
  adjusted with the j - 1 before it, `p.adjust(p[1:j])[j]`; `'sidak'` is
  `1 - (1 - p)^j`; the first class is left exactly unadjusted.

The user's choices (asked at the start):

- **`p_adjust = 'holm'` plus `progressive = TRUE`**, two arguments as in
  vegan, whose defaults reproduce vegan's. `p_adjust` takes any
  `p.adjust()` method or vegan's `'sidak'`. (Instead of the plan's single
  `p_adjust` with a `'progressive'` value.)
- **`alternative = 'sign'`** (new value, the default): vegan's per-class
  direction. `'greater'`, `'less'` and `'two.sided'` work as in `mantel()`.
- **Gaps:** when an untested class lies between tested ones, vegan adjusts
  only the classes before the gap position (and gives tested classes after
  it an `NA` adjusted p-value; with `'sidak'`, non-progressive, it uses
  `n.tests` as the exponent). ecodive adjusts every tested class in order
  of distance, skipping untested ones. Identical to vegan whenever there is
  no gap; documented under "Differences".

My own decisions:

- **Not in `DIST_TESTS`.** It returns a table, not an `ecodive_test`, and
  needs class arguments, so it is a plain function, as the plan says and
  as the README's function table files it ("Other distance tests and
  models"). No aliases, so the registry check doesn't apply.
- **One C pass for all classes** (`C_mantel_classes()` in `src/mantel.c`),
  not one `C_mantel()` call per class: a class's correlation with dm
  depends on the permutation only through the sum of dm over the class's
  pairs (Pearson on standardized dm, Spearman on its standardized ranks,
  Kendall's tau-b on `2 * midrank - N - 1`, summed in int64, which is the
  Mann-Whitney form of tau's numerator). So each permutation is O(N) for
  every class and method, and Kendall is exact. It reuses `perm_source_*`,
  `perm_matrix()`, `mantel_std()`, `mantel_codes()`, `mantel_tied()` and
  `perm_p()`; `mantel_core()`'s input handling moved, unchanged, into
  `mantel_inputs()`, which the three functions share.
- **The first class includes its lower boundary.** vegan lowers the first
  break by `.Machine$double.eps`, which does nothing for distances of about
  2 or more, so their closest pairs fall in no class (e.g. mite × 1000
  loses 2 of 2,415 pairs). ecodive keeps them; documented.
- **`n_pairs` counts each pair once**; vegan's `n.dist` is twice that.
- **r within 1e-12 of 0 counts as 0** for the direction: with ties, a true
  0 can come out as ±1e-17, which flipped the tail against vegan's exact 0
  in 3 of 1,528 sweep cases.
- **`y` forms:** as `mantel()`'s, plus coordinates (a numeric matrix or data
  frame, one row per sample, that isn't a square symmetric zero-diagonal
  matrix; or several column names of `data`), whose Euclidean distances
  are used, as vegan's `XY`.
- **`breaks` and `n_classes` together are an error.** Distances outside
  `breaks` are in no class (they count as "other pairs" in every model
  matrix).
- **Result:** the plan's columns; attributes `breaks`, `method`,
  `alternative`, `p_adjust`, `progressive`, `n_perm`, `permutation`, and
  `r_perm` (classes × permutations, `NA` rows for untested classes). Title
  "Mantel correlogram of <dm> [and <y>]", as `mantel()` composes it.
- **Plot:** `plot(x, alpha = 0.05)`, vegan's look (line through filled
  squares for `p_adjusted <= alpha`, open otherwise, dashed zero line).
  `alpha` is a `plot()` argument, not stored on the result (nothing
  computed depends on it).
- **Docs:** `vignettes/stat_tests.Rmd` gains a short correlogram example;
  `mantel()` and `mantel_partial()` link to it.


## Outcome (from the implementing session's final report)

**Item 58, done 2026-11-26.** Commit `97ddfb1` (function), plus the tick.

- Tests: 31,134 at home and 454 in CRAN mode, before merging the other
  sessions' work (30,895 / 450 at 15dc935). `test-mantel_correlogram.r` has
  235 tests.
- vegan 2.7-8 `mantel.correlog()`, shared permutation matrix: the 5 stored
  mite cases (Sturges, custom breaks, n.class, cutoff = FALSE; Pearson,
  Spearman, Kendall; Holm, BH, Šidák, Bonferroni; progressive and not) and
  400 fresh random cases (n 5–60, tied and untied, random breaks, every
  `mult`): classes, n and tested classes all identical; r within 9.1e-16;
  1,528 raw p-values and 1,089 adjusted ones (cases without gaps) identical.
- ecodist 2.1.3 `mgram()`: Pearson and Spearman r within 2.2e-16 (its
  first class is half-open, so it drops pairs at the minimum distance).
- cpus 1, 2, 3, 4, 8 identical.
- Speed (n = 150, 999 permutations): vegan 2.5 s, ecodive 0.08 s on one
  thread.
- `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed size,
  the container's debug build flags).
