# `bioenv()`: the subset of environmental variables that best matches community distances

**Priority:** 9. **Tier:** port from scikit-bio (BSD-3). Read
[README.md](README.md) first.


## What it is

BIOENV (Clarke & Ainsworth 1993) searches subsets of numeric environmental
variables for the one whose Euclidean distances (on standardised variables)
have the highest rank correlation with the community distances. It's in
QIIME 2 (`diversity bioenv`) and vegan (`bioenv()`).


## Design

- **Signature:** `bioenv(dm, data, formula = NULL, method = 'spearman', metric = 'euclidean', max_vars = NULL, given = NULL, cpus = n_cpus())`.
  - `dm` is the community distances. `data` holds the sample metadata,
    matched to `dm` by row name. `formula` is a one-sided formula over
    `data` naming the candidate variables (`~ N + P + K + pH`), as
    variables come from `data` through `formula` everywhere else; the
    default `NULL` searches every numeric column of `data`.
  - `method`: the correlation (`'spearman'` by default, as in vegan; also
    `'pearson'` and `'kendall'`).
  - `metric`: the distances among the variables. vegan offers
    `'euclidean'`, `'mahalanobis'`, `'manhattan'` and `'gower'`. Start
    with `'euclidean'`, the default and the only one skbio has, and add the
    others if cheap (ecodive already has `manhattan` and `gower` code).
  - `max_vars`: the largest subset size (a cap, so `max_*`), as vegan's
    `upto`.
  - `given`: distances to partial out, as vegan's `partial`; the same name
    and concept as `mantel_partial()`'s `given`. This is a partial
    correlation; reuse `mantel_partial()`'s correlation code (see
    [mantel_partial.md](mantel_partial.md)).
- **Result:** a plain function, not part of `dist_test()` (it takes many
  variables and returns a table).
  - A data frame with one row per subset size: `size`, `variables` (a
    comma-joined string), `correlation`.
  - Class `c('ecodive_bioenv', 'data.frame')`, with attribute `best` (the
    best subset overall).
- **Cost:** 2^p − 1 subsets. Parallelise over subsets, with each subset's
  rank correlation computed against the precomputed ranks of `dm`. Warn
  when p > about 20.
- **Masking:** `vegan::bioenv`.


## Reference implementations

- **Port from:** scikit-bio 0.7.4 `skbio.stats.distance.bioenv()` (BSD-3).
  - It searches exhaustively, z-scores the variables, uses Euclidean
    distances and Spearman correlation, and breaks ties by taking the first
    subset, "to match vegan".
  - **It lacks** vegan's `upto`, `partial`, `metric` and `method` options
    (ecodive's `max_vars`, `given`, `metric` and `method`),
    and takes a distance matrix where vegan takes `comm` and `index`.
    ecodive takes `dm`, like skbio.
  - The degrees of freedom used for standardising can't change a Spearman
    result, but can change a Pearson one: use vegan's (sd with n − 1).
- **ecopy `bioenv` (MIT):** not inspected; secondary.
- **Output-only fixtures (GPL):** vegan 2.7-x `bioenv(comm, env)` and
  `bioenv(..., upto = , partial = , metric = , method = )`, using the same
  dissimilarity (`vegdist(comm, 'bray')` equals ecodive's `bray()`).
  vegan's `summary(bioenv())` gives the best subset per size.


## Tests

- Fixtures from vegan, written by `data-raw/expected/bioenv.r`: varespec
  against varechem, with every option. Compare the best subset per size
  (identical) and the correlations at 1e-10.
- skbio agreement on the default options. A Python fixture is optional,
  since vegan covers it.
- Ties between subsets follow the first-found rule, as skbio and vegan do;
  test a constructed tie.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page; add it to the reference index next to `mantel()` and
  `envfit()`.
- NEWS.


## References

- Clarke, K. R., & Ainsworth, M. (1993). A method of linking multivariate
  community structure to environmental variables. *Marine Ecology
  Progress Series*, 92, 205–219.


## Decisions made while implementing

- **`metric`: all four of vegan's** (the user's choice). `'euclidean'` on
  variables scaled by `scale()` (sd with n − 1); `'manhattan'` on variables
  scaled to unit range, as vegan's docs say; `'gower'` as `'manhattan'`
  divided by the number of variables (Gower's distance for numeric
  variables, so its correlations equal `'manhattan'`'s, as vegan's docs
  say); `'mahalanobis'` as \|L⁻¹(z_i − z_j)\| for L the Cholesky factor of the
  subset's covariance.
  - Mahalanobis is computed from each pair's difference vector, not by
    whitening each sample, so pairs with the same differences keep the same
    distance; for one variable, its ranks are Euclidean's exactly. vegan's
    whitening rounds some tied distances apart (254 distinct distances for
    varechem's Al, where `dist(scale(Al))` has 252), so for Spearman and
    Kendall the size-1 correlation differs from vegan's by about 5e-5;
    every best subset is the same, and the other sizes and Pearson match at
    1e-14. The tests check size 1 against Euclidean's on the same variable.
  - Linearly dependent variables (a Cholesky pivot ≤ 1e-10 of the variance)
    have no Mahalanobis distances: that subset's correlation is undefined
    and never best; a size with no defined subset has `NA`.
- **No warning for many variables** (the user's choice). The docs say the
  search doubles with each variable (about a million subsets for 20) and
  that `max_vars` limits it. The C code stops with an error only if one
  size has more than 2^52 subsets (beyond exact counting).
- **Non-numeric columns** (the user's choice): with `formula = NULL`, only
  numeric (integer or double) columns are searched and others are skipped
  silently (logical columns count as non-numeric); a formula that names a
  non-numeric variable is an error. Unlike vegan, factors aren't allowed
  even with `'gower'`.
- **`formula`**: one-sided (a left-hand side is ignored), each term one
  variable, evaluated with `model.frame()` so that `~ log(N) + P` and
  `~ . - K` work. Interactions and `Condition()` are errors (the latter
  points to `given`). A character formula is accepted, as in
  `model_design()`.
- **`given`**: a `dist` or symmetric matrix (matched by name, as in
  `mantel_partial()`), one value per row of `data`, or a column name of
  `data`; a column named by `given` is left out of the default variables.
  Partial correlation as `mantel_partial()`: the R code reuses
  `mantel_cor()` for r_xz and its standardizing (`mantel_std()`) and Kendall
  coding (`mantel_codes()`, `mantel_tied()`), now top-level helpers in
  `R/mantel.r`; the C code reuses mantel.c's Kendall merge sort, now
  exported as `kendall_sort_count()` and `kendall_tied_sorted()`. Kendall
  with `given` is allowed (vegan refuses it) and tested against a direct
  computation. `dm` and `given` perfectly correlated, or either constant, is
  an error.
- **`max_vars`** above the number of variables is capped silently; below 1,
  or not a whole number, is an error.
- **Result**: `size`, `variables` (joined by `", "`, as scikit-bio),
  `correlation`; attributes `best` (character vector; the smallest of
  equally good sizes, as vegan's `which.max`), `method`, `metric`,
  `partial` and `n_subsets`. No print method: it prints as a data frame.
- **Ties**: subsets of a size are visited in `combn()` order and the first
  strictly greatest wins. Threads take contiguous runs of subsets whose
  bests are combined in order, so this holds for any `cpus`.
- **Licensing**: the search follows scikit-bio's `bioenv()` (BSD-3), so
  `src/bioenv.c` carries its notice, `inst/COPYRIGHTS` has a stanza for
  `src/bioenv.c` and `R/bioenv.r`, and the existing scikit-bio `ctb`/`cph`
  entry in `DESCRIPTION` now also names BIOENV.
- **Speed**: 100 samples × 12 variables, Spearman: 1.2 s on one thread,
  0.45 s on four, where `vegan::bioenv()` takes 5.5 s; varechem's 14
  variables: about 0.1 s, where vegan takes 1.6–29 s depending on the
  options.


## Outcome (from the implementing session's final report)

**Item 9, done 2026-11-19.** Commits `b78f35f` (`counts` in sample-count
messages), `96381a6` (work), `5c0023d` (tick).

- Tests: 8,117 at home, 201 in CRAN mode.
- vegan 2.7.8, every method × metric × `partial` × `upto`: best subsets
  identical, correlations ≤ 1.1e-14; exception: Spearman / Kendall with
  Mahalanobis at size 1 differ ≤ 5.8e-5 (vegan's whitening breaks exact
  ties; ecodive keeps them).
- scikit-bio 0.7.4 defaults: identical subsets, ≤ 3.9e-16.
- Speed: 100 samples × 12 variables in 1.2 s (1 thread), 0.45 s (4); vegan
  5.5 s. Port of scikit-bio's search (credited).
