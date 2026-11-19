# `varpart()` and `forward_select()`: variation partitioning and forward selection

**Priority:** 22 (`varpart`) and 23 (`forward_select`). **Tier:**
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
