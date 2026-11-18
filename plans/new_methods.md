# Plan: methods from the 2026-11-18 R package survey

This file plans items 1–18 of the survey's "methods ecodive lacks" list. For
each method it records:

- what the method mainly takes as input;
- what it returns;
- the name the ecodive function would have, and whether that name collides
  with anything;
- which wrapper, if any, the function belongs under.

Nothing here is implemented yet. Every method would still follow ecodive's
usual rules:

- match a reference implementation (named below) to rounding, checked against
  stored fixtures;
- multithreaded C, with results identical for any `cpus`;
- no new dependencies.


## 1. The methods

**Input** column:

- **counts** = a count table, as for `pca()`, `beta_div()` and the rest;
- **dm** = beta diversity distances;
- **meta** = sample metadata (`groups`, or `data` with a `formula`, as for
  `dbrda()` and `permanova()`);
- **ord** = an ordination result;
- **tree** = a phylogeny.

**Returns** column:

- **ord** = an ordination matrix (class `ecodive_ord`, with attributes
  `method`, `eig`, `prop`, ...);
- **test** = an `ecodive_test` object;
- **table** = a data frame.

| # | Method | Proposed function | Input | Returns | Reference to match |
|--:|:-------|:------------------|:------|:--------|:-------------------|
| 1 | Transformation-based RDA (constrained PCA) | `rda()` | counts (+ `norm`), meta | ord | `vegan::rda` |
| 2 | Canonical correspondence analysis | `cca()` | counts (raw abundances), meta | ord | `vegan::cca` |
| 3a | Correspondence analysis | `ca()` | counts | ord | `vegan::cca(X)`, `ade4::dudi.coa` |
| 3b | Detrended correspondence analysis | `dca()` | counts | ord (+ gradient lengths) | `vegan::decorana` |
| 4a | Variation partitioning | `varpart()` | dm or counts, + 2–4 sets of variables from meta | table (fractions [a]–[d], adjusted R²) | `vegan::varpart` |
| 4b | Forward selection of variables | `forward_select()` | dm or counts, meta, scope formula | table (steps) + final formula | `vegan::ordiR2step`, `adespatial::forward.sel` |
| 5a | Procrustes rotation | `procrustes()` | two ords (or two dms, via `pcoa()`) | ord (rotated) + m² | `vegan::procrustes` |
| 5b | PROTEST | `protest()` | two ords | test | `vegan::protest` |
| 6 | Pairwise PERMANOVA | a `pairwise` argument on `permanova()` (and `anosim()`, `mrpp()`) | dm, meta | test + pairwise table | `pairwiseAdonis`, `RVAideMemoire::pairwise.perm.manova` |
| 7 | SIMPER | `simper()` | counts, groups | table (contribution per taxon, per pair of groups) | `vegan::simper` |
| 8 | BIOENV | `bioenv()` | dm, data (numeric columns) | table (best subset per size) | `vegan::bioenv` |
| 9a | Multiple regression on distance matrices | `mrm()` | dm + explanatory distance matrices (or columns of data) | table (coefficients, R², p) | `ecodist::MRM` |
| 9b | Mantel correlogram | `mantel_correlogram()` | dm, + a distance matrix of distances (e.g. geographic) | table (r and p per distance class) | `vegan::mantel.correlog`, `ecodist::mgram` |
| 10 | Double principal coordinate analysis | `dpcoa()` | counts, tree | ord (+ taxon scores) | `ade4::dpcoa` |
| 11 | PLS-DA (sparse and orthogonal variants as options) | `plsda()` | counts (+ `norm`, usually CLR), groups | ord (+ loadings, VIP) | `mixOmics::plsda` / `splsda`, `ropls::opls` |
| 12 | MiRKAT kernel test | `mirkat()` | one or more dms, outcome, covariates (meta) | test | `MiRKAT::MiRKAT` |
| 13 | ~~Surface fitted over an ordination~~ (dropped; see §5) | ~~`ordisurf()`~~ | ord, one numeric variable | grid + deviance explained | `vegan::ordisurf` |
| 14 | Co-inertia analysis | `coinertia()` | two tables (counts, or counts and env) | ord for each table + RV coefficient and test | `ade4::coinertia`, `ade4::RV.rtest` |
| 15a | Discriminant canonical analysis of principal coordinates | `capdiscrim()` | dm, groups | ord + classification success + test | `BiodiversityR::CAPdiscrim` |
| 15b | Linear discriminant analysis | `lda()` | counts (or an ord, when features outnumber samples), groups | ord + classification | `MASS::lda` |
| 16 | Principal response curves | `prc()` | counts, time, treatment (meta) | table (coefficients per time × treatment) + taxon weights | `vegan::prc` |
| 17a | Distance-based Moran's eigenvector maps | `dbmem()` | sample coordinates, or a geographic dm | matrix of spatial variables (+ Moran's I of each) | `adespatial::dbmem` |
| 17b | Moran's I | `moran()` | a numeric variable + spatial weights or coordinates | test | `ape::Moran.I` |
| 18a | PHATE | `phate()` | dm | ord | `phateR` (Python `phate`) |
| 18b | Diffusion maps | `diffusion_map()` | dm | ord | `destiny::DiffusionMap` |
| 18c | PaCMAP | `pacmap()` | dm (via its kNN) | ord | Python `pacmap`; `manifoldsR` |
| 18d | LargeVis | `largevis()` | dm (via its kNN) | ord | `uwot::lvish` |

Items 1–18 hold 27 methods. Item 6 needs no new function: `permdisp()`
already has a `pairwise` argument, and the plan copies it to the other
group tests.


## 2. Name collisions

### Within ecodive

None of the proposed function names is an ecodive export, or the ID of a
metric, ordination or test. But two of them break the flexible name
matching of `dist_ord()` and `count_ord()`:

- **`ca` and `cca` against `capscale`.** `make_haystack()` (`R/metrics.r`)
  drops an alias that is a prefix of an alias already registered. The
  aliases `capscale` and `cap` are already there, so `ca` would be dropped.
  Checked on the current code: with a `ca` row added, the needle from
  `'CA'` prefix-matches only `capscale`, so `'CA'` would **silently run
  capscale**. Once `cca` is added too (alias
  `canonicalcorrespondenceanalysis`), it would stop as ambiguous instead.
  The lowercase ID `'ca'` works only because exact IDs are checked first.
  **Fix:** in `match_dist_ord()` (and the other `match_*()` functions),
  check whether the letters-only lowercase needle *equals* an alias before
  trying prefixes.
- **`cap` against discriminant CAP.** The `cap` alias already selects
  `capscale()`, and the documentation says so. But in Anderson's PRIMER
  usage, "CAP" is the discriminant method (15a). Keep `cap` for `capscale`
  (vegan's meaning, which R users know), and name 15a `capdiscrim()`, with
  aliases `capdiscriminant` and `discriminantcap`, not `cap`.

Some names are close to existing ones but don't collide:

- `rda` is not a prefix of `dbrda`, and `rpca` is a different ID.
- `dca` is distinct.
- `mrm` and `mirkat` are lowercase, while the references are `ecodist::MRM`
  and `MiRKAT::MiRKAT`.
- `hellinger` is already the name of a beta diversity function. For tb-RDA
  it would be an option *value* of `norm`, not a function name, so the two
  don't collide.

### With other packages (masking when both are attached)

ecodive already shares names with vegan, ape and uwot (`pcoa`, `pca`,
`dbrda`, `capscale`, `isomap`, `umap`, `anosim`, `mantel`, `mrpp`,
`envfit`), and documents this. The new names would add:

- vegan: `rda`, `cca`, `varpart`, `procrustes`, `protest`, `simper`,
  `bioenv`, `ordisurf`, `prc`
- ade4: `dpcoa`, `coinertia`
- `MASS::lda` (MASS is a recommended package and often attached)
- `ca::ca`
- `mixOmics::plsda`, `caret::plsda`
- `adespatial::dbmem`
- `phateR::phate`

That's the same trade-off as today: names users already know, with masking
that the documentation explains, and the wrappers or `ecodive::` to be sure.

The alternative is non-colliding names such as `dca()` instead of
`decorana()`. `dca()` already avoids vegan's name, since vegan calls it
`decorana`. Using the method's common abbreviation, as ecodive does now,
is the more consistent choice.


## 3. Is `dist_ord()` enough?

`dist_ord(dm, method, ...)` looks up the method, then passes `dm` and `...`
to its function. That's fine for any method whose first argument is
distances and whose result is an ordination matrix, constrained or not:
`groups`, `data` and `formula` go through `...`, as they do for `dbrda()`
now.

### Methods that fit `dist_ord()` as it is

Only registry rows are needed (plus the alias fix in §2):

- `phate`, `diffusion_map`, `pacmap`, `largevis` (18a–d)
- `capdiscrim` (15a), which takes `groups` like `dbrda`

### Methods that don't fit, and why

1. **Count-table ordinations:** `rda`, `cca`, `ca`, `dca`, `dpcoa`,
   `plsda`, `lda`. Their first argument is counts, not distances, and they
   take `norm`, `pseudocount` or `tree`. `pca()` and `rpca()` are already
   in this position: `match_dist_ord('pca')` stops with "call pca()
   directly". With 7 more such methods, 9 in all, a second wrapper earns
   its place.
2. **Methods that don't return a sample × axis matrix:**
   - `varpart`, `forward_select`, `bioenv`, `mrm`, `mantel_correlogram`,
     `simper`, `prc` and `dbmem` return tables or explanatory variables.
   - `ordisurf` returns a surface.
   - `procrustes` and `coinertia` take two inputs.

   A by-name wrapper is only useful when its methods take the same inputs
   and return the same shape, so that users can swap one name for another.
   These methods share neither, so wrapping them would add nothing.

### Recommendation

- **Keep `dist_ord()` as it is.** Add rows to its registry, and fix the
  alias matching.
- **Add one sibling wrapper, `count_ord(counts, method, ...)`,** with
  `list_count_ords()` and `match_count_ord()`, like the `dist_ord()` /
  `dist_test()` pair. It would run `pca`, `rpca`, `rda`, `cca`, `ca`,
  `dca`, `dpcoa`, `plsda` and `lda`.
  - The `ORDINATIONS` registry already has an `input` column (now always
    `"dist"`), and `list_dist_ords(input = )` already filters on it. So
    one registry can serve both wrappers: rows with `input = "counts"`
    go to `count_ord()`.
  - `match_dist_ord()`'s "call pca() directly" error becomes a pointer to
    `count_ord()`, and the reverse for `match_count_ord()`.
  - The name follows `beta_div(counts, ...)` and `alpha_div(counts, ...)`,
    whose first argument is also `counts`. `table_ord()` would also work.
- **Don't merge the two into one `ordinate()` that dispatches on its
  input.** `dist_ord()` accepts a symmetric matrix as distances, and a
  count table can be square, so the class of `x` can't reliably tell a
  distance matrix from a count table. Two explicit entry points avoid
  that.
- **`dist_test()` stays, and gains `mirkat`** (12), which takes `dm` and a
  numeric or binary outcome, plus covariates through `data` and `formula`.
  - Pairwise tests (6) go through `...`. `dist_test(dm, 'permanova', 'trt',
    meta, pairwise = TRUE)` would need no change to `dist_test()`.
  - `protest` (5b) and `moran` (17b) return `ecodive_test` objects, but
    they test two ordinations, or a variable against space, not `dm`
    against a variable, so they're called directly.
- **All other methods are plain functions,** grouped in the reference
  index rather than behind a wrapper:
  - Post-ordination, next to `envfit()` and `pcoa_biplot()`: `procrustes`,
    `protest`, `ordisurf`.
  - Model building: `varpart`, `forward_select`, `bioenv`.
  - Distance regression: `mrm`, `mantel_correlogram`.
  - Community composition: `simper`, `prc`.
  - Spatial: `dbmem`, `moran`.
  - Two-table: `coinertia`.


## 4. Shared groundwork

Done on 2026-11-18, before any new method:

- **`norm = 'hellinger'`**: the square root of each sample's proportions,
  available wherever `norm` is (beta metrics, `pca()`, `pcoa_biplot()`,
  `envfit()`, `normalize()`). It is `NORM_HELLINGER` (6) in
  `src/normalize.c` and `R/validate.r`. Negative values are an error.
- **Exact-alias matching**: `match_name()` in `R/utils.r`, used by every
  `match_*()`, tries an exact name, ID or alias before prefix matching, so
  `ca` and `CA` will select correspondence analysis once it is added, not
  `capscale`. It also fixed `match_metric('Jensen')`, which was ambiguous.
- **`count_ord()`**, with `list_count_ords()` and `match_count_ord()`
  (`R/count_ord.r`). It has its own registry, `COUNT_ORDS` (`pca`, `rpca`),
  separate from `ORDINATIONS`. Each wrapper's matcher points to the other
  wrapper, or to `dist_test()`, when given a method it doesn't run.

Still to do, but best done together with the first method that needs it:

- **A shared RDA core.** `rda`, `prc`, `varpart`, `forward_select` and
  `plsda` all need the constrained solver that `dbrda()` has (the QR of the
  design, `Condition()`, and the permutation F). Generalise `dbrda.c` to
  take a centred table as well as Gower's matrix, rather than writing a
  second solver. Do this with `rda()`.
- **A weighted (chi-square) variant of that core** for `cca`, `ca` and
  `dca`, with row and column weights. Do this with `cca()`.
- **A two-input result convention** for `procrustes` and `coinertia`.


## 5. Decisions (the user, 2026-11-18)

1. The wrapper for count-table ordinations is `count_ord()`.
2. Add `hellinger` to `norm` (done).
3. **Drop `ordisurf`** (item 13): it would need a GAM written in C.
4. Python reference implementations are fine for `phate`, `pacmap` and
   `diffusion_map`.
5. Pairwise tests default to the **Holm** p-value adjustment, and the user
   can choose another, or none (see §6).
6. Order of work:
   1. the alias fix, `count_ord()` and `hellinger` (done);
   2. `rda`, `cca`, `ca` and `dca` (they share the core);
   3. pairwise tests and `procrustes` / `protest`;
   4. `varpart` and `forward_select`;
   5. `simper` and `bioenv`;
   6. the rest.



## 6. Pairwise tests: design (the user, 2026-11-18)

Applies to `permanova()`, `anosim()`, `mrpp()` and `permdisp()`, each with
`pairwise = TRUE`.

- **Argument:** `p_adjust = 'holm'`. It accepts any method in
  `stats::p.adjust.methods` (`'holm'`, `'hochberg'`, `'hommel'`,
  `'bonferroni'`, `'BH'`, `'BY'`, `'fdr'`) or `'none'`, matched as
  `stats::p.adjust()` matches them. Giving `p_adjust` without
  `pairwise = TRUE` is an error. It is passed through `dist_test()`'s
  `...`, so `dist_test()` needs no change.
- **Output:** the pairwise table keeps the raw `p_value` column and adds
  `p_adjusted` beside it; the raw values are never overwritten. With
  `p_adjust = 'none'`, `p_adjusted` equals `p_value`. The method used is
  stored in the result (a `p_adjust` field), and `print()` names it.
- **Which p-values are adjusted:** each pairwise table's permutation
  p-values, across that table's pairs. The overall test is one p-value and
  isn't adjusted.
- **`permdisp()` changes too.** It already has `pairwise = TRUE`, which
  adjusts nothing (as vegan's `permutest()` doesn't). It gets the same
  `p_adjust` argument, defaulting to Holm like the others. Its table has two
  raw p-values: `p_value` (from the t distribution) and `p_value_perm` (from
  the permutations). Adjust each separately, giving `p_adjusted` and
  `p_adjusted_perm`. The vegan fixtures compare the raw columns, so they
  still apply. This changes `permdisp()`'s default output, so it needs a
  NEWS entry.
- **Tests:** `p_adjusted` equals `stats::p.adjust(p_value, method)` for
  every method and for `'none'`. The raw `p_value` columns still match
  vegan, `pairwiseAdonis` and `RVAideMemoire`. Invalid `p_adjust` values,
  and `p_adjust` without `pairwise = TRUE`, are errors.
