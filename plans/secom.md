# `secom_pearson()`, `secom_distance()`: bias-corrected correlations (SECOM)

**Priority:** 57. **Tier:** port of R code from the Bioconductor package
ANCOMBC (**Artistic-2.0**: verified 2026-11-20 from the DESCRIPTION of the
`devel` branch, version 2.15.1, commit `fe21fa6`; author and copyright
holder Huang Lin). Artistic-2.0 code may be ported with its notices (README
"Licensing"), as `ancombc()` already does. SECOM's code calls three **GPL**
helpers, which are output-only and must be replaced by clean-room
equivalents written from their definitions: `Hmisc::rcorr()`,
`DescTools::Winsorize()` and `energy::dcor.test()` (all GPL ≥ 2, checked
on CRAN on 2026-11-20). Read
[README.md](README.md), then [sparcc.md](sparcc.md) (output structure),
[assoc_network.md](assoc_network.md), and [ancombc.md](ancombc.md) for the
notice pattern.


## What it is

SECOM (Lin, Eggesbø & Peddada 2022) corrects log abundances for
sample-specific biases (differing sampling fractions) and taxon-specific
ones (centring), then estimates sparse correlations. It has two
correlations, the paper's three named estimators:

- **SECOM (Pearson1):** Pearson, sparsified by a cross-validated threshold;
- **SECOM (Pearson2):** Pearson, sparsified by p-value filtering;
- **SECOM (Distance):** distance correlation (any dependence, not only
  linear), with permutation p-values.

It has 195 citations, rising from 19 in 2023 to 67 in 2026 (OpenAlex
2026-11-19). It is the most-cited permissively licensed newer network
estimator.


## Names

- `secom_pearson()` and `secom_distance()`: one function per correlation
  measure, the paper's "Pearson" and "Distance" names (rule 2: base `secom`,
  modifier). ANCOMBC calls them `secom_linear()` / `secom_dist()`. "Pearson1"
  and "Pearson2" are two sparsifications of the *same* estimate, so they are
  `assoc_network()` options, not functions:
  - Pearson1: `assoc_network(x)`, whose default threshold is SECOM's
    cross-validated one;
  - Pearson2: `assoc_network(x, alpha = 0.005, p_adjust = 'none')`.
- **Spearman** (ANCOMBC's `secom_linear(method = "spearman")`; the paper's
  Theorem 2) is not offered at first. **Ask the user**:
  - a separate `secom_spearman()` would follow the one-function-per-variant
    rule;
  - an argument would follow the `mantel(method)` precedent, but it can't be
    called `method` inside `feature_assoc()`.
- No registry alias starts with `sec` (checked 2026-11-20). Both join
  `feature_assoc()` if it is accepted.


## Design

### Signatures

```r
secom_pearson(counts, margin = 1L, pseudocount = 0, bias_counts = NULL,
              winsorize = c(0.05, 0.95), max_bias_corr = 0.5,
              n_thresholds = 100L, n_splits = 10L, clust = NULL,
              seed = 0L, cpus = n_cpus())

secom_distance(counts, margin = 1L, pseudocount = 0, bias_counts = NULL,
               winsorize = c(0.05, 0.95), max_bias_corr = 0.5,
               permutations = 1000L, clust = NULL, seed = 0L,
               cpus = n_cpus())
```

ANCOMBC → ecodive names, with ANCOMBC's defaults:

| ANCOMBC | ecodive |
|:--------|:--------|
| `pseudo = 0` | `pseudocount = 0`: zeros are **missing** (`log(0)` → `NA`, pairwise-complete), as `ancombc2(pseudocount = 0)` already treats them |
| `aggregate_data` vs `data` | `counts` vs `bias_counts`: the biases are estimated from `bias_counts` (default `counts`). ANCOMBC estimates them on the fine table and correlates the aggregated one. The paper wants many taxa (> 50, warned) for the bias estimate. **Name to confirm with the user.** |
| `wins_quant` | `winsorize` |
| `corr_cut` | `max_bias_corr` |
| `thresh_len` | `n_thresholds` |
| `n_cv` | `n_splits`: repeated random splits, not folds |
| `R` | `permutations` |
| `n_cl` | `cpus` |
| `thresh_hard`, `max_p` | `assoc_network(threshold, alpha)` |

Not offered (say so in the docs):

- `prv_cut = 0.5` / `lib_cut = 1000` filtering. Follow `ancombc()`'s
  precedent: the docs show the two lines. Keep **taxa with prevalence ≥
  0.5 first, then samples whose depth over the kept taxa is ≥ 1000**, the
  order of `.data_core()` in `R/ancombc_prep.R`. Fixtures apply the same
  filter.
- `soft = TRUE` (the paper's eq. 12 is soft thresholding; ANCOMBC's default
  is hard).
- `alpha_grid` (default 0).
- The `corr_reg` output (eigenvalue-regularised correlations).
- Several ecosystems in one call (`data = list(gut, tongue)`).
- `tax_level` aggregation.

### Algorithm (ANCOMBC `R/secom_prep.R`, `R/secom_linear.R`, `R/sparse_linear.R`, `R/secom_dist.R`, `R/sparse_dist.R`)

1. **Sample biases** (`.s_diff_est()`; paper eq. 7), on `bias_counts`:
   - \eqn{o = \log(\text{bias\_counts})}, with zeros → `NA` (no pseudocount
     here, whatever `pseudocount` is);
   - centre each taxon over samples (`na.rm`);
   - weight \eqn{w_j = 1/\mathrm{var}_j} (n − 1 divisor, `na.rm`), dropping
     taxa whose weight isn't finite;
   - \eqn{\hat s_i} is the weighted mean of sample i's observed centred
     values. R's `weighted.mean(na.rm = TRUE)` drops the weights of missing
     values too: Σ x w / Σ w over the observed taxa.
   - **Not `ancombc()`'s code.** ANCOM-BC estimates its bias by EM (`src/ancombc.c`);
     SECOM uses this closed-form weighted mean. Reuse `ancombc()`'s input
     handling (`log_ratio_input()`, the sparse per-feature lists of
     `src/proportionality.c`), not its estimator.
2. **Bias-corrected abundances:** \eqn{\hat y_{ij} = (o_{ij} -
   \bar o_{\cdot j}) - \hat s_i}, with \eqn{o = \log(\text{counts} +
   \text{pseudocount})}, zeros → `NA` when the pseudocount is 0.
3. **Winsorise** each taxon to its `winsorize` quantiles (R type 7, over
   non-`NA` values), clamping. Clean room: DescTools is GPL, and clamping is
   all it does here.
4. **Co-occurrence** \eqn{n_{lm}}: the number of samples where both are
   non-`NA`. ANCOMBC warns when any is < 10. Pairs with \eqn{n_{lm} < 2}
   get correlation 0 and p = 1.
5. **Pearson** (`secom_pearson`): pairwise-complete Pearson r on the
   winsorised values, and the two-sided t-test p-value with \eqn{n_{lm}-2}
   df. That is textbook, and it is what `Hmisc::rcorr()` reports.
   **Unverified:** whether rcorr computes in single precision. If its r
   differs from `cor()` beyond 1e-12, match `stats::cor(use =
   "pairwise.complete.obs")` (ANCOMBC uses that for thresholding) and
   compare rcorr's at its own precision.
6. **CV threshold** (`.sparse_linear()` and `.thresh_loss_grid()`; port
   them):
   - The grid is `seq(0, max|r|, length = n_thresholds)`, where the maximum
     is taken over entries ≠ 1.
   - Each of `n_splits` splits takes a random training set of
     \eqn{n_1 = n - \lfloor n/\log n\rfloor} samples, and computes the
     pairwise-complete correlations of both halves.
   - The loss is ½‖hard-threshold(corr₁) − corr₂‖_F over non-`NA` entries,
     including the diagonal, with suffix sums as ANCOMBC does.
   - The chosen threshold is the one with the smallest mean loss (`which.min`:
     the first on ties).
   - Splits come from PCG, one stream per split. ANCOMBC uses doRNG, which
     can't be reproduced.
7. **Distance correlation** (`secom_distance`). Clean room from Székely,
   Rizzo & Bakirov (2007), since energy is GPL.
   - Per pair, on the rows where both are non-`NA` (ANCOMBC first drops
     rows missing x, then those missing y), compute the double-centred
     distance matrices A, B.
   - \eqn{\mathrm{dCov}^2_n = \frac1{m^2}\sum A_{ij}B_{ij}} and
     \eqn{\mathrm{dCor} = \sqrt{\mathrm{dCov}^2/\sqrt{\mathrm{dVar}^2_x
     \mathrm{dVar}^2_y}}}, 0 when a dVar is 0. The paper's Definition 1.
   - p-value: permute y's rows `permutations` times,
     \eqn{\sum A_{ij}B_{\pi i\,\pi j}}, and `perm_p()`. **Unverified:**
     energy's p-value formula (check that the fixtures' p-values are
     multiples of 1/(R+1)) and its index = 1 default.
   - A pair with a single complete row gets 0 / 1, as ANCOMBC.
   - The diagonal is 1.
8. **False-positive filter:** for each taxon, \eqn{\mathrm{cor}(\hat s,
   \hat y_{\cdot j})}, pairwise-complete on the *un-winsorised*
   \eqn{\hat y}. If it is > `max_bias_corr` for **both** taxa of a pair, that
   pair gets r = 0 and p = 1.
9. Taxa stay in input order. ANCOMBC sorts them alphabetically; fixtures
   are reordered to compare.

### C and threads

Steps 1–3 and 8 are O(nD) in R or C. Steps 5 and 7 are pairwise, in C:
- threaded over rows of the packed triangle, with the
  `src/proportionality.c` pattern (per-feature lists, fixed-order sums);
- for distance correlation, an O(m²) workspace per thread;
- PCG permutation streams per pair (pair index), so results are identical
  for any `cpus`.

The CV runs one task per split. Cost of `secom_distance()`: D²/2 × R × m².
Document it (D = 100, n = 200, R = 1000: ~2e11 operations, minutes on 4
cores). Huo & Székely's O(m log m) algorithm is a later optimisation; it
must give bit-identical permutation decisions, or the tolerance in
`perm_p()` must absorb the difference.

### Result (feature-correlation structure from sparcc.md)

An `ecodive_dist` over features, a **similarity**: ANCOMBC's `corr` /
`dcorr` after the false-positive filter. Its record holds:

- `metric = "secom_pearson"` / `"secom_distance"`;
- `similarity = TRUE`, `features = TRUE`;
- `n_samples`;
- the settings;
- for Pearson, `cv_threshold` (ANCOMBC's `thresh_opt`) with
  `threshold_strict = TRUE`. ANCOMBC zeroes |r| **≤** the threshold, so
  `assoc_network(threshold = NULL)` must keep |r| **>** it here, unlike its
  usual ≥. Note it in both plans' docs.

Parallel attributes (listed in `parallel`): `p_value` (`corr_p` / `dcorr_p`)
and `n_pairs` (`mat_cooccur`, used by `assoc_network_diff()`'s Fisher z).
Other attributes:
- `bias`: \eqn{\hat s}, named by sample;
- for Pearson: `cv_error` (mean loss per threshold) and `thresholds`.

The registry marks `fisher_z = TRUE` for `secom_pearson` only.


## Licensing (as `ancombc()`)

- `R/secom.r` and `src/secom.c` headers carry ANCOMBC's notice ("ANCOMBC:
  Copyright (c) Huang Lin. Licensed under the Artistic License 2.0") and
  the section-4 statement: "C/R port of `.s_diff_est()`, `.abn_est()`,
  `.sparse_linear()`, `.thresh_loss_grid()`, `.sparse_dist()`,
  `secom_linear()` and `secom_dist()` from ANCOMBC <version>; restructured
  for ecodive (threaded C, PCG splits and permutations, clean-room
  replacements for Hmisc, DescTools and energy); see plans/secom.md".
- The user-visible differences go in the Rd's "Differences from ANCOMBC"
  section.
- Add a `Files:` stanza to `inst/COPYRIGHTS` (the Artistic-2.0 text is
  already there) and extend Huang Lin's `Authors@R` comment ("... and
  SECOM").
- **Version:** the devel code (2.15.1) was read for this plan. Port from the
  version used for the fixtures (2.8.1 installs on the R that the `ancombc()`
  fixtures used), after diffing its `secom_*` files against devel. The
  `.thresh_loss_grid()` suffix-sum refactor may be devel-only; the losses
  must agree either way.
- Note a docs/code mismatch in the reference: `thresh_hard` is documented
  as defaulting to 0.3 but defaults to 0 in the code.


## Reference implementations and fixtures

- `data-raw/expected/secom.r` runs ANCOMBC `secom_linear()` /
  `secom_dist()` on two tables (D = 60 with zeros; D = 120, n = 80) with
  `prv_cut`/`lib_cut` matching the pre-filter. It saves `s_diff_hat`,
  `y_hat`, `corr`, `corr_p`, `mat_cooccur`, `cv_error`, `thresh_grid`,
  `thresh_opt`, `corr_th`, `corr_fl`, `dcorr`, `dcorr_p`.
- **Deterministic parity at 1e-10:** `s_diff_hat`, `y_hat`, `corr`,
  `corr_p` (rcorr caveat above), `mat_cooccur`, `thresh_grid` and `dcorr`.
- **CV:** an internal `secom_details()` accepts a split-index matrix. The
  generator calls ANCOMBC's internal `.thresh_loss_grid()` (Artistic, so
  callable and readable) on the same splits. Match `cv_error` at 1e-10 and
  the chosen threshold exactly; `corr_th` / `corr_fl` then follow through
  `assoc_network()`.
- `dcorr_p` matches in distribution only (permutations differ). Check
  exactly with a supplied permutation matrix against a base-R
  implementation of the permutation test.


## Tests

- Fixtures as above. Identities: dCor = 1 for exact linear pairs and ~0.5
  for \eqn{y = (x-5)^2} (paper, Supplementary Fig. 3, approximate); Pearson
  p-values equal `cor.test()`'s.
- The false-positive filter on a constructed taxon that tracks \eqn{\hat s}.
- Edge cases:
  - a pair with no overlap (0 / 1);
  - all-zero taxa (error, or dropped with a warning; decide);
  - n < 3;
  - D < 50 (warned, as ANCOMBC);
  - `bias_counts` with different samples (error).
- `identical()` for `cpus` 1, 2 and 4; same `seed`, same result.
- Smoke-test lines.


## Docs

One Rd page for both functions. It covers:
- the two biases, and why correcting them matters for correlations
  (the paper's Fig. 1);
- Pearson1/Pearson2/Distance and how to get each with `assoc_network()`;
- the filters to apply first;
- the cost of `secom_distance()`;
- "Differences from ANCOMBC";
- "Same names": none.

NEWS. Add `bias_counts`, `winsorize`, `max_bias_corr`, `n_thresholds` and
`n_splits` to parameter_names.md.


## References

- Lin, H., Eggesbø, M., & Peddada, S. D. (2022). Linear and nonlinear
  correlation estimators unveil undescribed taxa interactions in microbiome
  data. *Nature Communications*, 13, 4946.
  doi:10.1038/s41467-022-32243-x
- Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and
  testing dependence by correlation of distances. *Annals of Statistics*,
  35, 2769–2794. doi:10.1214/009053607000000505
- Huo, X., & Székely, G. J. (2016). Fast computing for distance covariance.
  *Technometrics*, 58, 435–447. doi:10.1080/00401706.2015.1054435
- Lin, H., & Peddada, S. D. (2020), for ANCOM-BC; see ancombc.md.
