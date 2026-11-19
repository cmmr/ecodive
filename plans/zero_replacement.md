# `replace_zeros()` and its methods: zero replacement for compositional data

**Priority:** 15 (the wrapper and `replace_zeros_mult()`) and 36 (the four `cmultRepl()` methods). **Tier:**
- `replace_zeros()` (wrapper, registry, `list_*` / `match_*`) and
  `replace_zeros_mult()`: port (scikit-bio 0.7.4 `multi_replace()`,
  BSD-3).
- `replace_zeros_czm()`, `replace_zeros_gbm()`, `replace_zeros_sq()`,
  `replace_zeros_bl()`: **clean room** (zCompositions `cmultRepl()`, GPL,
  is output-only; the defining paper is closed access, so some formulas
  must be confirmed from outputs; see "Unverified").

Read [README.md](README.md) first, and [compositional.md](compositional.md)
for the conventions the compositional plans share.


## What it is

Logs of zero are undefined, so every log-ratio method needs the zeros
dealt with. ecodive's `pseudocount` adds a constant to every cell, which
changes the ratios between non-zero parts and acts as an arbitrary
weighting (the docs already warn about it). **Multiplicative** replacement
(Martín-Fernández et al. 2003) instead gives each zero a small value and
shrinks the sample's non-zero parts by the same factor, so their ratios
are unchanged. **Bayesian-multiplicative** replacement
(Martín-Fernández et al. 2015) chooses the zero's value from a
Dirichlet-multinomial posterior, so it adapts to each sample's depth and
to the part's abundance elsewhere. zCompositions' `cmultRepl()` (default
`GBM`) is the de facto standard in microbiome CoDa workflows.

| Function | Method | Reference |
|:---------|:-------|:----------|
| `replace_zeros_mult()` | multiplicative replacement with a fixed value δ | scikit-bio `multi_replace()`; Martín-Fernández et al. (2003) |
| `replace_zeros_czm()` | count zero multiplicative: δ from one count's detection limit, per sample | zCompositions `cmultRepl(method = "CZM")` |
| `replace_zeros_gbm()` | geometric Bayesian-multiplicative | `cmultRepl(method = "GBM")` (its default) |
| `replace_zeros_sq()` | square-root Bayesian-multiplicative | `cmultRepl(method = "SQ")` |
| `replace_zeros_bl()` | Bayes-Laplace Bayesian-multiplicative | `cmultRepl(method = "BL")` |

They are separate named methods in the paper and in zCompositions'
documentation, so they are separate functions (README rule), joined by a
wrapper because they share one signature and one output, and because the
literature recommends checking that conclusions survive a change of zero
treatment, which is a loop over `list_zero_replacements()`.


## Names

- `replace_zeros_<code>`: the family sorts together (rule 2), and the code
  is zCompositions' method code, so `cmultRepl(method = "SQ")` users find
  `replace_zeros_sq()`. `mult` for the plain multiplicative method
  (zCompositions' `multRepl`, scikit-bio's `multi_replace`).
- Registry aliases (spelled out, lower case, no punctuation, as the other
  registries): `multiplicative`, `multreplacement`, `multirepl`,
  `countzeromultiplicative`, `geometricbayesianmultiplicative`, `cmultrepl`
  (→ gbm, the default of `cmultRepl`), `squareroot`, `sqrtbm`,
  `bayeslaplace`, `laplace`. Check each against the other registries
  (README snippet); `bl` / `sq` / `gbm` / `czm` / `mult` themselves must not
  prefix-match anything else.
- `replace_zeros()`'s selector is `method` (as `dist_ord()` /
  `count_ord()`), so no method may take an argument named `method` or a
  prefix of `method` / `counts` (README shadowing table, already updated).
  zCompositions' own `method` and `output` arguments therefore don't
  carry over: the method is the function, and the output is always
  proportions.


## Design

### Signatures

```r
replace_zeros(counts, method, ...)

replace_zeros_mult(counts, margin = 1L, delta = NULL)
replace_zeros_czm (counts, margin = 1L, frac = 0.65, threshold = 0.5)
replace_zeros_gbm (counts, margin = 1L, frac = 0.65, adjust = TRUE)
replace_zeros_sq  (counts, margin = 1L, frac = 0.65, adjust = TRUE)
replace_zeros_bl  (counts, margin = 1L, frac = 0.65, adjust = TRUE)

list_zero_replacements(...)       # as list_count_ords()
match_zero_replacement(method)    # as match_count_ord()
```

- No `cpus`: each is a few passes over the table. Add it only if
  benchmarks on large sparse tables call for it (then parallel over
  samples, as `normalize()`).
- `delta` (mult): the value every zero gets, as a proportion. `NULL`
  (default) is scikit-bio's \eqn{1/D^2}. Must be > 0 and small enough that
  no sample's non-zero parts go negative (scikit-bio's error:
  "Multiplicative replacement created negative proportions. Consider using
  a smaller `delta`.").
- `frac`, `threshold`, `adjust`: zCompositions' names and defaults.
- zCompositions' `z.warning = 0.8` / `z.delete = TRUE` (drop rows or
  columns with ≥ 80% zeros) are **not** carried over: ecodive never drops
  samples or features silently. Instead, warn once when any sample or
  feature is ≥ 80% zeros, naming how many (`name_list()`), and point to
  filtering. Fixtures must therefore use tables without such rows/columns,
  or pre-filter them.

### Input checks (all methods)

`validate_counts()`; finite, non-negative; no empty sample (all zeros).
The Bayesian methods need counts (they use \eqn{n_i}); warn if values
aren't whole numbers (proportions make \eqn{n_i = 1} and the prior
dominates). `mult` accepts proportions.

### Formulas

Notation: \eqn{c_{ij}} counts, \eqn{n_i = \sum_j c_{ij}},
\eqn{x_{ij} = c_{ij}/n_i}, \eqn{D} parts, \eqn{N} samples, \eqn{Z_i} the
zero parts of sample i.

**mult** (scikit-bio `multi_replace`, read from source):

\deqn{r_{ij} = \delta \text{ if } c_{ij} = 0, \qquad r_{ij} = x_{ij}\,(1 - |Z_i|\,\delta) \text{ otherwise.}}

Each sample sums to 1. scikit-bio fixture: `[[1,2,3,0,5],[1,0,0,4,5],[1,2,3,4,5]]`
→ `[[.087273,.174545,.261818,.04,.436364],[.092,.04,.04,.368,.46],[.066667,.133333,.2,.266667,.333333]]`.

**Bayesian-multiplicative** (Martín-Fernández et al. 2015; zCompositions
manual): the posterior estimate of part j in sample i is
\eqn{(c_{ij} + s_i t_{ij}) / (n_i + s_i)}, with prior "strength"
\eqn{s_i} and prior expectation \eqn{t_{ij}}. Only zeros are replaced:

\deqn{r_{ij} = \frac{s_i t_{ij}}{n_i + s_i} \text{ if } c_{ij} = 0, \qquad r_{ij} = x_{ij}\Bigl(1 - \sum_{k \in Z_i} \frac{s_i t_{ik}}{n_i + s_i}\Bigr) \text{ otherwise.}}

The methods differ in \eqn{s_i} and \eqn{t_{ij}}:

| Method | \eqn{t_{ij}} | \eqn{s_i} | Status |
|:-------|:-------------|:----------|:-------|
| BL | \eqn{1/D}? or the leave-one-out mean below? | \eqn{D} | **confirm** |
| SQ | as BL | \eqn{\sqrt{n_i}} | **confirm** |
| GBM | leave-one-out mean of part j's proportions over the other samples, \eqn{t_{ij} = \frac{1}{N-1}\sum_{k \ne i} x_{kj}} | \eqn{1 / g(t_{i\cdot})}, \eqn{g} the geometric mean over parts | **confirm** |

- **`adjust = TRUE`** (zCompositions manual): an imputed proportion above
  the smallest observed (non-zero) proportion of that part,
  \eqn{m_j = \min_{k: c_{kj} > 0} x_{kj}}, is set to \eqn{frac \cdot m_j};
  then the sample's non-zero parts are re-scaled multiplicatively so it
  sums to 1.
- **CZM:** every zero becomes \eqn{frac \cdot threshold / n_i}
  (0.325 / n_i by default: a fraction of the "half a count" detection
  limit), and the non-zero parts are scaled as in mult. Check whether
  `adjust` applies to CZM too.
- **Edge cases to establish from outputs:** a part that is non-zero in
  only one sample (its leave-one-out \eqn{t} is 0 there, making GBM's
  \eqn{g} zero and \eqn{s_i} infinite); a part that is zero everywhere
  (zCompositions likely errors or drops it; ecodive should error, naming
  it, since nothing can be estimated).

### Unverified: pin down before writing C

The 2015 paper is closed access and zCompositions is GPL, so the exact
\eqn{t} and \eqn{s} of each method, the edge cases and the adjust step
must be **established from documentation and outputs only**:

1. Get the paper (library access) and the zCompositions manual/vignette;
   record the formulas they state.
2. Black-box `cmultRepl()` (zCompositions 1.6.2, `output = "prop"`,
   `z.delete = FALSE`) on small probe tables designed to separate the
   candidates: 3 × 3 and 4 × 5 tables with one zero, two zeros in a sample,
   unequal depths, a part present in a single sample, and a zero whose
   imputed value would exceed \eqn{m_j} (to trigger `adjust`).
   `data-raw/expected/zero_replacement_probe.r` fits each candidate
   formula in R and reports which matches to 1e-12.
3. Write the confirmed formulas into this plan (replace the table above)
   before implementing. Don't read zCompositions' R code at any point.

### Result (all methods)

- A dense base `matrix` of **proportions** (each sample sums to 1, within
  rounding), in the input's orientation (`margin`), with its dimnames.
  Dense, because every zero gets a value. (zCompositions' `output =
  "p-counts"` isn't offered: every log-ratio method is scale invariant.)
- Attributes: `method` (`'mult'`, `'czm'`, `'gbm'`, `'sq'`, `'bl'`) and
  the parameters used (`delta` with its resolved value; `frac`,
  `threshold`, `adjust`), and `n_replaced` (the number of zeros
  replaced).
- `phyloseq`, `rbiom`, `SummarizedExperiment` input gives the matrix, as
  `normalize()` does.
- Warn above 1e8 values, as `normalize(norm = 'clr')`.
- The result has no zeros, so passing it to anything with
  `norm = 'clr'`, `aitchison()`, `alr()`, `ilr()`, `philr()`, `vlr()`,
  `ancom()` ... uses no pseudocount and gives no warning.

### The wrapper

- `ZERO_REPLACEMENTS` registry in `R/replace_zeros.r`, same columns as
  `COUNT_ORDS` (`name, id, alt_ids`):

  ```text
  name,                                    id,   alt_ids
  Multiplicative Replacement,              mult, multiplicative multreplacement multirepl
  Count Zero Multiplicative,               czm,  countzeromultiplicative
  Geometric Bayesian-Multiplicative,       gbm,  geometricbayesianmultiplicative cmultrepl
  Square-Root Bayesian-Multiplicative,     sq,   squareroot sqrtbm
  Bayes-Laplace Bayesian-Multiplicative,   bl,   bayeslaplace laplace
  ```

  The function for row `id` is `replace_zeros_<id>`. Copy
  `count_ord()` / `list_count_ords()` / `match_count_ord()` (and
  `count_ord_info()`), and add the registry to the alias-collision snippet
  in README.md (done) and to `test-wrapper_args.r`.
- `replace_zeros()` has **no default `method`**: the choice matters and
  should be visible. Its docs recommend `'gbm'` as the common default.
- Rows get the mult/czm functions in the first item; the Bayesian rows
  join in the second item (the registry is built so rows can be added).


## Reference implementations

- **mult:** scikit-bio 0.7.4 `multi_replace()` in
  `skbio/stats/composition/_base.py` (BSD-3; credit in `inst/COPYRIGHTS`
  as the existing scikit-bio stanzas). Fixtures from scikit-bio.
- **czm, gbm, sq, bl:** zCompositions 1.6.2 `cmultRepl()` (GPL-2 | GPL-3),
  **output-only**: `data-raw/expected/zero_replacement.r` with
  `output = "prop"`, `z.delete = FALSE`, on tables with no all-zero or
  ≥ 80%-zero rows/columns; compare at 1e-10.
- Papers: Martín-Fernández et al. (2003, 2015); Palarea-Albaladejo &
  Martín-Fernández (2015) for the zCompositions documentation.


## Tests

- Fixtures as above, both margins, dense and sparse input.
- Each sample sums to 1; non-zero parts keep their ratios exactly (mult,
  czm, and the Bayesian methods without adjustment); no zeros remain.
- `replace_zeros(x, 'gbm')` is `identical()` to `replace_zeros_gbm(x)`;
  aliases resolve; `list_zero_replacements()` lists all rows.
- Downstream: `aitchison(replace_zeros_gbm(x))` runs with no pseudocount
  warning.
- Errors: negative or non-finite values; an all-zero sample or part;
  `delta` too large (scikit-bio's message); `frac` outside (0, 1);
  unknown `method`.
- Smoke-test lines for every exported function.


## Docs

- `replace_zeros()` page explains the choice and links each method; one
  page per method (or one shared Rd with `@rdname` for the four
  `cmultRepl` methods, as the UniFrac family does).
- The pseudocount section (`R/documentation.r`, `pseudocount_section`)
  gains a sentence pointing to `replace_zeros()`.
- "Compositional Data" group, subtitle "Zeros".
- NEWS.


## References

- Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V.
  (2003). Dealing with zeros and missing values in compositional data
  sets using nonparametric imputation. *Mathematical Geology*, 35,
  253–278.
- Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
  Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of
  count zeros in compositional data sets. *Statistical Modelling*, 15,
  134–158.
- Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015).
  zCompositions — R package for multivariate imputation of left-censored
  data under a compositional approach. *Chemometrics and Intelligent
  Laboratory Systems*, 143, 85–96.
- Lubbe, S., Filzmoser, P., & Templ, M. (2021). Comparison of zero
  replacement strategies for compositional data with large numbers of
  zeros. *Chemometrics and Intelligent Laboratory Systems*, 210, 104248.
