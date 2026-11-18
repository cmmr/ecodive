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
