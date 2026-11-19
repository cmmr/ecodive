# `varpart()` and `forward_select()`: variation partitioning and forward selection

**Priority:** 12 (`varpart`) and 13 (`forward_select`). **Tier:**
*derived*. No permissive code exists, but both are thin layers over
`rda()` and `dbrda()` plus published formulas. **Requires
[rda.md](rda.md) first** (partial RDA). Read [README.md](README.md) first.


## What they are

- **Variation partitioning** (Borcard et al. 1992; Peres-Neto et al.
  2006) splits the variation in a community, explained by two to four
  sets of explanatory variables (e.g. environment, space, time), into
  fractions:
  - unique to each set;
  - shared between sets;
  - unexplained.

  It uses adjusted R² from a series of (partial) RDAs or db-RDAs. In high
  demand in ecology; vegan `varpart()`.
- **Forward selection** (Blanchet, Legendre & Borcard 2008) picks
  explanatory variables one at a time. At each step it adds the one that
  raises adjusted R² most, if its permutation test is significant, and
  stops when adjusted R² would exceed the global model's. It's usually run
  before varpart. vegan `ordiR2step()`; adespatial `forward.sel()`.


## Design

- **`varpart(x, ..., data = NULL, norm = 'none', pseudocount = NULL, margin = 1L, cpus = n_cpus())`.**
  - `x` is either a `dist` (db-RDA based, as vegan's `varpart(Y = dist)`)
    or a count table (RDA based, with `norm`, e.g. Hellinger).
  - `...` holds two to four explanatory sets, each a one-sided formula
    over `data` (`~ ph + temp`, `~ x + y`) or a data frame or matrix.
  - **Result:** a data frame of the fractions: rows [a], [b], [c], [d] for
    two sets (and [a]–[o] and residuals for four, as vegan names them).
    - Columns `df`, `r_squared`, `adj_r_squared`, `testable`.
    - Class `c('ecodive_varpart', 'data.frame')`, with a print method
      showing the Venn fractions as text.
  - Plotting is out of scope; the data frame is enough.
- **`forward_select(x, formula, data, norm = 'none', pseudocount = NULL, margin = 1L, alpha = 0.05, permutations = 999L, seed = 0L, cpus = n_cpus())`.**
  - `x` is a `dist` or a count table, as for `varpart()`. `formula` is a
    one-sided formula over `data` naming the candidate variables (vegan's
    `scope`), as for `bioenv()`.
  - **Result:** a data frame, one row per step: `variable`,
    `adj_r_squared`, `adj_r_squared_cum`, `F`, `p_value`.
    - Attribute `selected` holds the formula of the chosen variables,
      ready for `varpart()` or `rda()`.
    - Stop when a variable is non-significant (p > `alpha`), or when
      adjusted R² would exceed the global model's (Blanchet's double
      stopping rule, as in `ordiR2step()`).
  - Factors enter whole (all their dummy columns together), as vegan does.
- Both are plain functions, with no wrapper.
- **Masking:** `vegan::varpart`. `forward_select` clashes with nothing
  (vegan and adespatial use `ordiR2step` and `forward.sel`).


## Reference implementations

- **No permissive code.**
  - ecopy (MIT) reports adjusted R² from its `rda` using Ezekiel's
    formula, but its full dummy coding inflates the number of predictors;
    don't use it.
  - Build both on ecodive's `rda()` and `dbrda()` cores:
    `dbrda_r_squared()` in `R/dbrda.r` already computes R² and adjusted R²
    with `Condition()`, following Peres-Neto et al. (2006).
- **Output-only fixtures (GPL):**
  - vegan 2.7-x `varpart()` (for both `dist` and table Y, with
    `transfo = "hel"`), `ordiR2step()`, and `anova.cca` for the step
    tests;
  - adespatial `forward.sel()` as a second check of the selection order.
- **The step tests** are permutation F tests of the added variable,
  conditioned on those already selected (a partial RDA test). Reuse
  `perm_ss()` / `permanova()`'s machinery with a shared permutation
  matrix, so p-values match vegan's exactly.


## Tests

- Fixtures from vegan, written by `data-raw/expected/varpart.r`.
  - Data: mite (`mite`, `mite.env`, `mite.pcnm`) and varespec / varechem,
    with two, three and four sets, for both `dist` and Hellinger table Y.
  - Compare every fraction's adjusted R² at 1e-10.
- `forward_select()` selects the same variables in the same order as
  `ordiR2step()`, and its adjusted R² and F match at 1e-10. p-values match
  exactly with shared permutations.
- Edge cases: a set with one variable; collinear sets; a set that
  explains nothing (negative adjusted R², kept as vegan keeps it).
- Results are identical for any `cpus`; a smoke-test line for each
  function.


## Docs

- Rd pages, and a vignette section after constrained ordination: select
  variables, then partition.
- NEWS.


## References

- Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the
  spatial component of ecological variation. *Ecology*, 73, 1045–1055.
- Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006).
  Variation partitioning of species data matrices: estimation and
  comparison of fractions. *Ecology*, 87, 2614–2625.
- Blanchet, F. G., Legendre, P., & Borcard, D. (2008). Forward selection
  of explanatory variables. *Ecology*, 89, 2623–2632.


## Decisions made while implementing (varpart)

Item 12 (`varpart()` only; `forward_select()` is item 13 and is left open
here).

Asked of the user at the start (2026-11-19), and answered:

1. **CCA mode: yes,** as an argument: `varpart(counts, ..., chisquare =
   TRUE, permutations = 1000L, seed = 0L)`, vegan's `chisquare = TRUE`.
   Adjusted R² is estimated by permutation; the same permutations (a count,
   with PCG streams from `seed`, or a `shuffleSet()` matrix) serve every
   union. `permutations` and `seed` are ignored for RDA and db-RDA.
   `permutations` must be at least 1, and `norm` must be `'none'`.
2. **Sets may mix kinds:** each set in `...` is a one-sided formula over
   `data`, a data frame, a numeric (or logical) matrix, or a vector or
   factor (one variable).
3. **Data frames, matrices and named vectors are matched by row name** to
   `x`'s samples when both are named (extra rows ignored; `match_rows()`),
   else they must be in the samples' order. `data` is used only by formula
   sets. Data frames are dummy-coded as `model.matrix(~ .)` codes them
   (character and logical columns become factors); matrices are used as
   they are.
4. **Print:** a header (method, normalization, samples, sets), the total
   inertia (and, for RDA, the variance, as vegan prints it), one line per
   set (its label, df and adjusted R²) and for all sets together, a text
   Venn diagram of the adjusted fractions for two or three sets (none for
   four), then the fractions table. The union and conditional tables are
   kept as attributes but not printed.
5. **`testable` is kept,** as vegan's: `TRUE` for the unique fractions,
   `FALSE` for shared fractions and the residuals (and `TRUE` throughout
   the `unions` and `conditional` attributes).

Made while implementing:

6. **Result.** A data frame of class `c('ecodive_varpart', 'data.frame')`,
   rows named with vegan's letters (`[a]`–`[d]`, `[a]`–`[h]`, `[a]`–`[p]`,
   the last being the residuals; vegan's letter order for three and four
   sets is hard-coded in `varpart_fractions()`), columns `sets` (`"env"`,
   `"env & space"`, `"residuals"`), `df`, `r_squared`, `adj_r_squared`,
   `testable`. Attributes: `unions` (vegan's `fract`, rows named `X1`,
   `X1+X2`, ...; columns `fractions`, `df`, `r_squared`, `adj_r_squared`,
   `testable`), `conditional` (vegan's `contr1` and `contr2` together, rows
   named `"X1 | X3"`, `"X1 | X3+X4"`; `NULL` for two sets), `method`
   (`'rda'`, `'dbrda'`, `'cca'`), `norm` (RDA only), `total` (vegan's
   `SS.Y`), `n`, `set_labels` (the deparsed sets, named), and `n_perm`
   (CCA).
   - `r_squared` of the Venn fractions is the same inclusion–exclusion on
     unadjusted R² (vegan shows `NA`); documented as biased.
   - `df` is the rank a set adds to all the others for unique fractions,
     and `NA` (vegan: 0, or `NA` for residuals) for shared fractions and
     residuals.
   - Unnamed sets are `X1`..`X4`, as vegan's; named sets
     (`varpart(x, env = ~ ..., space = pcnm)`) use their names.
7. **One pass over the data (RDA and db-RDA).** New `C_cord_cross()` in
   `src/dbrda.c` (registered in `src/init.c`) centres the design of all
   the sets' columns, X, and returns `trace(B)` (unscaled) and
   `Xc' B Xc` (m × m), with B Gower's matrix of the distances (products by
   `gower_mult()`) or `A A'` of the centred, normalized table (via
   `pca_table_init()` / `pca_cross()`: `L = A' Xc`, cross = `L'L`, so sparse
   tables stay sparse). Every union's inertia then follows in R without
   another pass: `subset_inertia()` takes a pivoted `qr()` (tol 1e-7, as
   vegan) of the union's centred columns and computes
   `trace(R^-T P'Xc'BXcP R^-1)`. Deterministic for any `cpus` (products are,
   and cross products are summed on one thread in a fixed order).
8. **Shared pieces for item 13** (all in `R/varpart.r`):
   `varpart_response(x, norm, pseudocount, margin, chisquare)` (the
   response as `dbrda` / `rda` / `cca`, with labels), `response_cross(y,
   X, cpus)` (one C pass), `subset_inertia(model, cols)` (inertia, rank and
   kept columns of any subset of X's columns), and `adj_r_squared(r2, n,
   m)` (Ezekiel; `NA` without residual df). `forward_select()` can compute
   every candidate's adjusted R² from one `response_cross()` of all the
   candidate columns. Its permutation step tests are not covered by these
   (they need `perm_ss()` / `permanova()`-style machinery, as the plan
   says).
9. **CCA mode follows vegan's permutation scheme, found from outputs
   only.** vegan 2.7.8's `varpart(chisquare = TRUE)` does *not* give each
   union `RsquareAdj(cca(...))`'s value: with the same `shuffleSet()`
   matrix it permutes the rows of the chi-square matrix Q̄ against the
   fixed weighted design, even without conditions (the scheme `cca()`
   already uses for partial models), where `cca()`'s unconditioned
   adjusted R² permutes the design's rows. (mite, `~ SubsDens +
   WatrCont`, 99 permutations: 0.2593 against `RsquareAdj()`'s 0.2477;
   reproduced by hand in R before writing any code.) So `C_cca()` gained a
   10th argument, `sexp_fixed`, which sends unconditioned models to
   `cca_perm_inertia_fixed()` (with `rz = 0`). `cca()` and `ca()` pass
   `FALSE` and are unchanged. The unions' R² and the total inertia equal
   `cca()`'s exactly.
10. **vegan 2.7.8's four-set CCA is not matched:** with four sets and
    `chisquare = TRUE`, its `SS.Y` and union R² disagree with its own
    `cca()` (mite: `SS.Y` 13.8 where `cca()`'s total inertia, and vegan's
    own `SS.Y` for two or three sets, is 1.70). ecodive's four-set CCA is
    consistent with two and three sets and is tested against `cca()`'s R²
    for each union; the fixtures have CCA cases with two and three sets
    only.
11. **Collinearity:** columns aliased within a union are dropped by
    `qr()`'s rank decision, so shared columns move into the shared
    fractions (a set nested in another gets `[b] = 0`, df 0). When the
    sets' ranks sum to more than the rank of all sets, a warning says how
    many columns are shared (vegan warns "collinearity detected").
12. **Edge cases:** adjusted R² is `NA` when a union has no residual
    degrees of freedom (n − rank − 1 ≤ 0), and so are the fractions that
    use it; negative fractions are kept. Errors: fewer than 2 or more than
    4 sets; duplicate set names; fewer than 3 samples; `x` with no
    variation; a set with no variation (rank 0); `Condition()` in a set's
    formula; a two-sided set formula; a formula set without `data`; NA or
    infinite values; mismatched rows; `norm` or `chisquare` with
    distances; `chisquare` with `norm`; all-zero samples in CCA mode
    (`cca()`'s message).
13. **Not offered:** vegan's `scale` and `add` (no Lingoes / Cailliez for
    distances), and `transfo` beyond ecodive's `norm`. A plain matrix `x`
    is always a table; distances must be a `dist`.
14. **Docs:** `man/varpart.Rd` (method, speed, differences from vegan,
    masking), NEWS (2.4.0), the README's table (a "Model building" row with
    `varpart()` and `bioenv()`), the `ecodive` vignette's ordination
    section (one sentence), the `ordination` vignette (a "Partitioning
    variation" subsection with an example and the test of a unique
    fraction), and pkgdown's Ordination section (after `bioenv`).
15. **Fixtures use vegan 2.7.8 and permute 0.10.0** (scratch library);
    `data-raw/expected/varpart.r` writes
    `tests/tinytests/expected/varpart.rds` (58 kB, with the data, since
    tests can't load vegan's datasets).

Reference results:

- **vegan 2.7.8** `varpart()`, 23 cases: mite (`mite`, `mite.env`,
  `mite.pcnm`) and varespec / varechem; two, three and four sets; Hellinger
  table (`transfo = "hel"`), Bray-Curtis `dist` (`vegdist()` / `bray()`),
  and Euclidean distances of the Hellinger table; formula and data-frame
  sets (with factors); one variable per set; collinear and nested sets; a
  set of noise (negative fractions); CCA with two and three sets and stored
  99-permutation `shuffleSet()` matrices. Every Venn fraction's adjusted
  R², every union's df, R² and adjusted R², every conditional fraction,
  and `SS.Y` match at the 1e-10 tolerance. Largest absolute differences
  over all cases: Venn fractions 2.3e-15, unions' adjusted R² 2.3e-15,
  unions' R² 1.3e-15, conditional fractions 2.4e-15, `SS.Y` (relative)
  1.7e-15.
- **ecodive's own methods:** unique fractions equal partial `rda()` /
  `dbrda()` adjusted R² (`~ set + Condition(others)`), unions equal
  `rda()`'s, conditional fractions equal partial `rda()`'s, db-RDA of
  Euclidean distances of the Hellinger table equals RDA of the table, and
  CCA unions' R² equal `cca()`'s (all at 1e-12).
- **Speed:** 3,000 samples × 2,000 taxa (5% non-zero), sets `~ g + x`
  (4-level factor), `~ z` and 10 random columns: `varpart()` 0.13 s on 4
  threads for the Hellinger table and 0.12 s for Bray-Curtis distances;
  `vegan::varpart()` 3.5 s and 3.6 s. Fractions equal to rounding.
