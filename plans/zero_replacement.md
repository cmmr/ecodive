# `replace_zeros()` and its methods: zero replacement for compositional data

**Priority:** 15 (the wrapper and `replace_zeros_mult()`) and 55 (the four `cmultRepl()` methods and CZM). **Tier:**
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
replace_zeros_czm (counts, margin = 1L, frac = 0.65, threshold = 0.5, adjust = TRUE)  # adjust added in item 55
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

The methods differ in \eqn{s_i} and \eqn{t_{ij}}. **Confirmed on
2026-11-27 (item 55) by black-box probes of zCompositions 1.6.2
`cmultRepl()`** (see "Confirmed formulas" below; the table that stood here
had BL's and SQ's \eqn{t} as "1/D?" and GBM's \eqn{t} as the mean of the
other samples' proportions, both wrong for zCompositions):

| Method | \eqn{t_{ij}} | \eqn{s_i} |
|:-------|:-------------|:----------|
| BL | \eqn{T_{ij} = \sum_{k \ne i} c_{kj} \big/ \sum_{k \ne i} n_k} (pooled leave-one-out) | \eqn{D} |
| SQ | \eqn{T_{ij}} | \eqn{\sqrt{n_i}} |
| GBM | \eqn{T_{ij}} | \eqn{1 / g_i}, \eqn{g_i = \exp\bigl(\frac{1}{D}\sum_{j=1}^{D} \log T_{ij}\bigr)} over **all** D parts |

- **`adjust = TRUE`** (zCompositions manual): an imputed proportion above
  the smallest observed (non-zero) proportion of that part,
  \eqn{m_j = \min_{k: c_{kj} > 0} x_{kj}}, is set to \eqn{frac \cdot m_j};
  then the sample's non-zero parts are re-scaled multiplicatively so it
  sums to 1. **Confirmed**, for CZM too (see below).
- **CZM:** every zero becomes \eqn{frac \cdot threshold / n_i}
  (0.325 / n_i by default: a fraction of the "half a count" detection
  limit), and the non-zero parts are scaled as in mult. **Confirmed**, and
  `adjust` applies to CZM as well.

### Confirmed formulas (item 55, from documentation and outputs only)

Sources: the `cmultRepl()` help page and NEWS of zCompositions 1.6.2
(documentation), and its outputs. Its R source was never opened, and the
2015 paper was not available.

What the documentation states:

- BM methods: zeros "are imputed ... by posterior estimates of the
  multinomial probabilities generating the counts, assuming a Dirichlet
  prior distribution" with hyper-parameters \eqn{t} (prior estimates of
  the probabilities) and \eqn{s} (strength); `method = "user"` takes
  them as a matrix `t` and a vector `s`.
- `adjust`: imputations "falling above the lowest estimated probability of
  a multinomial part (c/n) ... [are] adjusted by using a fraction (frac)
  of the minimum c/n for that part. Lastly, the non-zero parts are
  multiplicatively adjusted."
- CZM "uses multiplicative simple replacement (multRepl) on the matrix of
  estimated probabilities", upper limit \eqn{threshold / n}, imputed
  value \eqn{frac} of it.
- NEWS 1.4.0-1: "error thrown when GBM method not applicable due to only
  one positive value in variable".

What the probes (`data-raw/expected/zero_replacement_probe.r`) established:

1. **The posterior form.** `method = "user"` with random `t` and `s`
   equals \eqn{r_{ij} = s_i t_{ij} / (n_i + s_i)} for zeros and
   \eqn{r_{ij} = x_{ij} (1 - \sum_{k \in Z_i} r_{ik})} otherwise
   (\eqn{2 \times 10^{-18}}).
2. **t and s.** With one or two zeros per row, each imputed value gives
   \eqn{s_i t_{ij}} directly. BL's values matched \eqn{s = D} only with
   the pooled leave-one-out \eqn{T_{ij}} above, not \eqn{1/D} and not the
   mean of the other samples' proportions (the probe tables have unequal
   depths, which separates the two means). SQ is the same \eqn{T} with
   \eqn{s_i = \sqrt{n_i}}; GBM the same \eqn{T} with \eqn{s_i = 1/g_i},
   the geometric mean taken over all D parts of \eqn{T_{i\cdot}}
   (including the sample's non-zero parts).
3. **The full candidate** (t, s, adjust, CZM) matches `cmultRepl()` on
   122 random tables (3–12 samples, 3–15 parts, 1/3 zeros, mixed depths),
   for every method, `adjust` TRUE and FALSE, `frac` 0.65 and 0.3 and
   `threshold` 0.5 and 0.7: largest relative difference
   \eqn{1.1 \times 10^{-15}}. The adjust step changed the result in 104
   (BL), 111 (SQ), 107 (GBM) and 83 (CZM) of the 122 tables, so it was
   exercised for every method.
4. **The adjust test is strict** (`>`): an imputed value exactly equal to
   \eqn{m_j} (1/16 against 1/16) is not adjusted.
5. **Edge cases.**
   - A part non-zero in only one sample: GBM errors ("not enough
     information to compute t hyper-parameter, probably there are columns
     with < 2 positive values"), because \eqn{T_{ij} = 0} for that sample
     makes \eqn{g_i = 0}. BL, SQ and CZM handle it (its \eqn{T} is 0 only
     in the sample where it is observed, which imputes nothing).
     zCompositions errors for GBM as soon as any part has fewer than two
     positive values, so with N = 2 and any zero, GBM always errors.
   - A part zero everywhere: GBM errors; BL and SQ impute 0 there (so
     zeros remain) with an R warning from `min()`; CZM imputes
     \eqn{frac \cdot threshold / n_i} (no adjust: no observed minimum).
   - One sample: every method errors ("at least two rows").
   - Non-integer values are accepted and used as they are (\eqn{n_i} is
     the row sum).
   - No zeros: an error ("Label 0 was not found").
   - `frac` and `threshold` are not checked (frac = 0 leaves zeros; frac =
     1.5 is accepted).

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
- Attributes: `method` (`'mult'`, `'czm'`, `'gbm'`, `'sq'`, `'bl'`), as
  `replace_zeros_mult()` already sets it, then results only: for mult
  `delta` (the value the zeros got, resolved from the data when not
  given), and `n_replaced` (the number of zeros replaced). The settings
  passed (`frac`, `threshold`, `adjust`) aren't echoed. Item 55 added
  `n_adjusted` (czm, gbm, sq, bl: the number of imputed values `adjust`
  lowered, which `cmultRepl()` prints), before `n_replaced`.
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- No plots. Every `replace_zeros_*()` function returns a plain table: no
  `display`, nothing saying how it was made, only the result attributes
  above. Later results computed from it don't mention the replacement.
- [da_sensitivity.md](da_sensitivity.md) labels its zero-treatment
  settings from the `ZERO_REPLACEMENTS` registry's `name` column, by the
  id it passed.

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


## Decisions made while implementing (item 15)

Asked of the user at the start:

1. **The 80%-zeros warning:** none for `replace_zeros_mult()` (scikit-bio
   gives none, and in microbiome tables most features are ≥ 80% zeros).
   Only the Bayesian methods (item 55) warn. The check is ready in
   `zero_replacement_input(mostly_zero = TRUE)`, covering samples and
   features: "1 sample ('a') and 1 feature ('f3') are at least 80% zeros,
   so most of their values will be imputed. Consider filtering them out
   first." (up to 5 names each, `name_list()`).
2. **`delta`:** one value, or one per sample (matched to the samples by
   name when both are named, else in order, through `per_sample()`). Each
   must be in (0, 1). The `delta` attribute is the scalar, or the
   per-sample vector in the samples' order, named after them.
3. **A feature that is zero in every sample:** allowed for mult (replaced
   like any zero, as scikit-bio does). The Bayesian methods should refuse
   it: `zero_replacement_input(empty_features_ok = FALSE)` errors, naming
   the features. A sample with no non-zero value is an error for every
   method, naming the samples.
4. **Registry and output, as planned:** `ZERO_REPLACEMENTS` has `name, id,
   alt_ids`, like `COUNT_ORDS`; the function of row `id` is
   `replace_zeros_<id>`. The result is a dense base `matrix` of proportions
   for every input type (sparse too), in the input's orientation with its
   dimnames; attributes `method`, the method's parameters, then
   `n_replaced` (a double: the number of zeros replaced). (Superseded on
   2026-11-20 by README.md § "Result attributes and plots": keep only
   `method`, `delta` and `n_replaced`; see "Result" above.)

Made while implementing:

5. **Bit-identical to scikit-bio 0.7.4 `multi_replace()`** on all 21
   fixture cases (D from 2 to 1037, counts and real values, a zero feature,
   a given delta, one delta per sample, one sample, no zeros), and on both
   margins. This needed:
   - each sample's total summed in **numpy's pairwise order**
     (`np_pairwise_sum()` in `src/replace_zeros.c`, following numpy's
     `pairwise_sum_DOUBLE()`: sequential below 8 values, 8 running sums up
     to 128, halves at a multiple of 8 above). NumPy (BSD-3) is credited in
     the C header, `inst/COPYRIGHTS` and `Authors@R` (`cph`);
   - scikit-bio's arithmetic step for step: closure first, zeros counted on
     the proportions, then `(1 - n_zeros * delta) * p`;
   - the **default delta from C's `pow(1.0 / D, 2.0)`**, as Python's `**`
     calls `pow()`, which differs by an ulp from `x * x` (R's `^2`) for 1,701
     of the D below 2e6 on glibc. The exponent is `volatile` so the
     compiler can't turn it into `x * x`. On another platform's libm the
     default may differ from that platform's Python by an ulp.
6. **zCompositions 1.6.2 `multRepl()`** (output only), given closed
   proportions, `frac = 1`, `dl = delta` for every part (a matrix for a
   per-sample delta), `z.warning = 1`, `z.delete = FALSE` (the argument is
   `z.delete`, not `delete`): within 3.3e-16 relative of scikit-bio (and so
   of ecodive) on the 16 cases it accepts (it refuses tables without zeros
   or with an all-zero part); stored and tested at 1e-14.
7. **A stricter `delta` check than scikit-bio's:** scikit-bio errors only
   when `1 - n_zeros * delta < 0`; ecodive also refuses `= 0`, which would
   leave the sample's non-zero parts at 0. The message keeps scikit-bio's
   words and names the samples: "Multiplicative replacement created
   negative proportions (or zeros) in 1 sample ('a'): `delta` times a
   sample's number of zeros must be below 1. Consider using a smaller
   `delta`."
8. **Aliases:** the plan's (`multiplicative`, `multreplacement`,
   `multirepl`) plus `multireplace`, scikit-bio's function name (`'multRepl'`
   and `'multi_replace'` both select mult). Checked with the README
   snippet: no alias or planned code (`mult`, `czm`, `gbm`, `sq`, `bl`) is
   an exact name in another registry, but **`sq` begins METRICS names**
   (`squared_euclidean`, `squared_chord`, `squared_chisq`, `squares`) and
   **`mult` begins DIST_TESTS names** (mrpp's "multi-response permutation
   procedure", permdisp's "multivariate dispersion"). That is harmless:
   each matcher searches only its own registry by prefix, and the other
   wrappers look at other registries' exact names only. No cross-wrapper
   redirects were added (`replace_zeros(x, 'pca')` says "does not match any
   option: mult"; `count_ord()` etc. are unchanged). `test-replace_zeros.r`
   builds the planned five-row registry and checks that every code and
   alias selects its own row, and that the codes are no other registry's
   exact names.
9. **No `cpus`:** one single-threaded C pass per sample (densified into a
   scratch vector). A sparse 2,000 × 5,000 table takes 0.2 s (0.45 s
   dense), vs 1.3 s for the vectorised R formula.
10. **`replace_zeros()` without `method`** errors, listing the ids, and
    passes `...` unchanged (no `seed` dropping, unlike `count_ord()`: no
    method has random parts). `replace_zeros(counts, method, ...)` is in
    `test-wrapper_args.r`.
11. **A vector** is one sample (as `validate_counts()` makes it): a D × 1
    result, without the vector's names (`validate_counts()` drops them).
    Integer, logical and data-frame tables are accepted; explicitly stored
    sparse zeros are zeros.
12. **Large results:** warn above 1e8 output values, as `alr()` / `ilr()`.
13. **Docs:** `replace_zeros.Rd` (choice, flexible matching, "Same names in
    other packages": none, zCompositions / scikit-bio named),
    `zero_replacements.Rd` (`list_*` / `match_*`), `replace_zeros_mult.Rd`.
    Both pseudocount sections (`pseudocount_section`,
    `log_ratio_pseudocount_section`) now point to `replace_zeros()`.
    pkgdown: "Compositional Data" → "Zeros".

For item 55 (the Bayesian methods and CZM):

- The plan's line "rows get the mult/czm functions in the first item" was
  superseded by TODO.md: **czm is item 55** too. The registry has only the
  `mult` row; add one row per method (the planned rows and aliases are in
  `test-replace_zeros.r`'s `planned` table).
- Write `replace_zeros_<id>(counts, margin = 1L, ...)`; start with
  `info <- zero_replacement_input(empty_features_ok = FALSE, whole = TRUE,
  mostly_zero = TRUE)` (sets `counts`, `margin`; returns `D`, `n`,
  `features`, `samples`, `sample_zeros`, `feature_nonzeros`) and end with
  `zero_replacement_result()`, passing no echoed settings (`frac`,
  `threshold`, `adjust` are not attributes; see "Result" above).
  `C_zero_counts()` and `dense_sample()` in `src/replace_zeros.c` are
  reusable. Whether `whole` / `mostly_zero` / `empty_features_ok`
  apply to CZM is open (CZM uses depths, so `whole = TRUE` likely).
- Update `replace_zeros()`'s docs (the `method` argument, the method table,
  flexible matching), the expectations of `df$id == 'mult'` in
  `test-replace_zeros.r`, `test-smoke.r`'s lines and the pkgdown "Zeros"
  list.
- The "confirm" formulas (t, s, the adjust step, CZM's adjust) remain
  open, as the plan says.


## Outcome (from the implementing session's final report)

### Item 15 (`replace_zeros()` wrapper, `replace_zeros_mult()`), done 2026-11-19

Commits `7d84f63` (work), `199eea8` (tick).

- Tests: 9,638 at home, 223 in CRAN mode.
- scikit-bio 0.7.4 `multi_replace()`: bit-identical on 21 cases (D 2–1,037;
  given and per-sample `delta`), after summing in NumPy's pairwise order
  and taking the default delta from C's `pow()`. zCompositions 1.6.2
  `multRepl()`: ≤ 3.3e-16 relative. NumPy credited.
- Open: the four Bayesian / CZM methods are a later item; the planned
  five-row registry is written into `test-replace_zeros.r`.


## Decisions made while implementing (item 55)

Asked of the user at the start (2026-11-27), after the probes:

1. **GBM and a feature non-zero in only one sample:** an **error** naming
   the features, as zCompositions errors ("not enough information to
   compute t hyper-parameter"). The message explains why (the feature's
   prior is zero in that sample, so the strength is infinite) and suggests
   filtering (features in at least 2 samples) or `sq` / `bl`. Not chosen:
   the \eqn{s \to \infty} limit, or a geometric mean over \eqn{t > 0} only
   (neither has a reference).
2. **No `cpus`** for the four methods (as the plan and
   `replace_zeros_mult()`): two O(N·D) passes. A sparse 2,000 × 5,000
   table (10% non-zero) takes 0.42 s for GBM, 0.26 s for CZM (0.21 s for
   mult); dense 0.9 s. The brief's "identical for cpus 1, 2 and 4" check is
   replaced by: identical for repeated calls, dense / `dgCMatrix` /
   `TsparseMatrix` / slam input, integer / double / data-frame input, and
   (transposed) margin 2.
3. **No 80%-zeros warning for any method** (overriding item 15's plan that
   the Bayesian methods would warn): the warning would fire on nearly every
   microbiome table. `zero_replacement_input(mostly_zero = TRUE)` stays,
   unused. The docs say that zCompositions warns about (or drops) such rows
   and columns and ecodive doesn't.

Made while implementing:

4. **The formulas** are zCompositions 1.6.2's, as confirmed (section
   "Confirmed formulas" above): BL's and SQ's prior is the pooled
   leave-one-out proportion, not the paper's 1/D (which I could not
   check; the paper is closed access). The docs state the formulas, so
   users see what is computed.
5. **`adjust` for CZM too** (`replace_zeros_czm(..., adjust = TRUE)`), as
   `cmultRepl()` applies it; the plan's signature lacked it. Strict `>`,
   as zCompositions.
6. **`n_adjusted`** attribute (a double) on all four results, between
   `method` and `n_replaced` (the wrapper's shape test keeps `n_replaced`
   last). It is a result (rule 1), and `cmultRepl()` prints it.
7. **Argument checks**, stricter than zCompositions (which checks none):
   `frac` a single number in (0, 1) exclusive, as the plan said;
   `threshold` in (0, 1] (a zero's upper limit can't be above one read);
   `adjust` a single `TRUE` / `FALSE`.
8. **Which inputs each method accepts:**
   - all four: no sample may be all zeros (shared check); a table without
     zeros gives its proportions (`cmultRepl()` errors); non-whole values
     give a warning (`whole = TRUE`; zCompositions accepts them silently)
     and are used as they are, as zCompositions does;
   - gbm, sq, bl: at least 2 samples (new `min_samples` option of
     `zero_replacement_input()`, checked before the empty-feature check so
     a single sample gets the right message); no feature zero everywhere
     (`empty_features_ok = FALSE`: BL / SQ in zCompositions would leave
     those zeros at 0);
   - czm: one sample is fine (a vector too), and a feature zero everywhere
     is replaced like any zero (`cmultRepl()` refuses one row; with an empty
     column it gives the same values, with a warning from `min()`);
   - czm refuses (naming the samples) a sample whose imputed values add up
     to 1 or more, where `cmultRepl()` returns negative proportions. That
     happened in 7 of 426 random CZM runs (213 tables, `adjust` on and
     off), all shallow samples with many zeros. The Bayesian methods can't reach it (each imputed value is
     below \eqn{t_{ij}\, s_i / (n_i + s_i)}).
9. **One C entry point**, `C_replace_zeros_cmult()` in
   `src/replace_zeros.c` (method code 0–3), two passes: depths, feature
   totals and each feature's smallest observed proportion; then each
   sample. Every sample is expanded densely before summing, so sparse and
   dense input give identical bits. \eqn{T_{ij}} is
   `(feature total - c_ij) / (grand total - n_i)`, exact for counts.
   The R wrappers share `cmult_replace()`.
10. **One Rd page** for the four, `replace_zeros_cmult.Rd` (`@rdname`;
    `?replace_zeros_gbm` etc. reach it), listed in pkgdown's "Zeros" as
    `replace_zeros_gbm`.
11. **Aliases as planned** (the `planned` table in `test-replace_zeros.r`
    is now asserted equal to the registry). The README snippet, run over
    all eight registries (adding `FEATURE_ASSOCS`; `BATCH_CORRECTIONS` has
    no aliases and its ids are `mmuphin`, `combat`): no new code or alias
    is an exact name elsewhere; only the known `sq` prefix of METRICS names
    (item 15). `'c'` alone is ambiguous (czm and `cmultrepl` → gbm), which
    the tests record. The collision test now also covers `DA_TESTS` and
    `FEATURE_ASSOCS`.
12. **`da_sensitivity(zeros = TRUE)`** tries every registry row, so it now
    runs five zero replacements instead of one; its lefse test was updated
    to expect them (NEWS says so).
13. **Fixtures:** `data-raw/expected/replace_zeros_cmult.r` →
    `tests/tinytests/expected/replace_zeros_cmult.rds` (86 kB): 8 tables
    (four probe tables, a 2-sample table, a non-integer table, `ex_counts`,
    and a sparse random 25 × 60 table with 72% zeros), every method ×
    `adjust` (GBM skipped where it refuses), plus non-default `frac` /
    `threshold`: 70 cases, with
    `z.warning = 1`, `z.delete = FALSE`. ecodive matches them to 1.9e-15
    (tested at 1e-13 relative), and `adjust` changes the result in at least
    one case of every method. The probes are kept in
    `data-raw/expected/zero_replacement_probe.r`.
14. **No new copyright stanza:** nothing was ported (clean room), so
    `inst/COPYRIGHTS` and `Authors@R` are unchanged; the C and R headers
    say how the formulas were established.
