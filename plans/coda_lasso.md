# `coda_lasso()` and `coda_lasso_pairwise()`: penalized log-contrast regression

**Priority:** 74. **Tier:** port. `coda_lasso()` ports zeroSum 2.0.7's
solver (R/C++, **MIT**, github.com/rehbergT/zeroSum). `coda_lasso_pairwise()`
ports coda4microbiome 0.2.4's model setup (CRAN, **MIT**) on ecodive's own
coordinate descent. glmnet (GPL-2) is output-only. Read
[README.md](README.md) first, and [compositional.md](compositional.md) for
the conventions the compositional plans share.


## What it is

This method predicts an outcome (a continuous measurement, or case vs control)
from a microbiome table and picks out the taxa that matter. It does so
in a way that respects compositionality. It is the "machine learning"
that compositional-data guidance recommends for microbiome tables, in
place of a plain lasso or random forest on relative abundances.

The **log-contrast model** (Aitchison & Bacon-Shone 1984; Lin, Shi, Feng
& Li 2014):

\deqn{g(E[y_i]) = \beta_0 + \sum_j \beta_j \log x_{ij}, \quad
\text{subject to } \sum_j \beta_j = 0.}

Because the coefficients sum to zero, a sample's total cancels out. The
fit is the same on counts, proportions or any rescaling of a sample, and
\eqn{\sum_j \beta_j \log x_{ij} = \sum_j \beta_j \,\mathrm{clr}(x_i)_j}.
Adding an elastic-net penalty to \eqn{\beta} selects a sparse signature: the
positive coefficients sit over the negative ones in a single log ratio.

There are two named forms in current use, so per the README there are two functions:

| Function | Model | Reference | Use |
|:---------|:------|:----------|:----|
| `coda_lasso()` | the constrained elastic net on \eqn{\log x}, \eqn{\sum\beta = 0} | Lin et al. (2014); zeroSum (Altenbuchinger et al. 2017) | the original formulation; scales to many taxa |
| `coda_lasso_pairwise()` | an elastic net on **all pairwise log ratios** \eqn{\log(x_j/x_k)}, \eqn{j<k}; the taxon coefficients \eqn{\beta_j = \sum_k \theta_{jk} - \sum_k \theta_{kj}} sum to zero by construction | coda4microbiome's `coda_glmnet()` (Calle et al. 2023) | the form most cited in applied microbiome papers in 2023–2026 |

**They are equivalent in one case,** derived by the research agent (2026-11-19) and to be
verified by a test. With a pure lasso (`l1_ratio = 1`), no standardisation and the
same loss and intercept, the pairwise fit at \eqn{\lambda} gives exactly
the \eqn{\beta} of the constrained fit at \eqn{\lambda / 2}. The reason
is that the cheapest set of pairwise \eqn{\theta} producing a given
\eqn{\beta} has \eqn{\|\theta\|_1 = \|\beta\|_1 / 2} (a minimum-cost flow
on the complete graph). \eqn{\theta} itself is not unique, but \eqn{\beta}
is. Standardisation or `l1_ratio < 1` breaks the equivalence, and
coda4microbiome's defaults (`alpha = 0.9`, probably glmnet's
`standardize = TRUE`) use both. So the two functions differ in practice.

**Names.** "coda" (compositional data analysis) is the field's prefix
(coda4microbiome, CoDA-lasso in Susin et al. 2020). The variant is
`coda_lasso_<modifier>` (README rule 2). Another possible name is `log_contrast()`,
but "contrast" means something else in `dirmult_ttest(contrast =)`.
**Ask the user** to confirm the names before building.

Not planned, but related:
- **selbal** (GPL-3; a forward search for one balance): it could become
  a clean-room `balance_select()` later.
- **trac** (GPL-3; tree-aggregated log contrasts through c-lasso): it
  would need `tree` and c-lasso's path algorithm. Revisit after
  `philr()`.
- **FLORAL** and **Compack** (GPL-3): augmented-Lagrangian solvers for
  the same model as `coda_lasso()`, used only as output cross-checks.
- **Cox and multinomial families** (zeroSum and coda4microbiome both have
  them): a later commit, if the user wants them.


## Design

### Signatures

```r
coda_lasso(counts, y, data = NULL, given = NULL,
           family = c('gaussian', 'binomial'), margin = 1L,
           pseudocount = NULL, l1_ratio = 1, lambda = NULL,
           n_lambdas = 100L, lambda_min_ratio = NULL, standardize = FALSE,
           n_folds = 10L, folds = NULL, select = c('min', '1se'),
           tol = 1e-8, max_iters = 100L, seed = 0L, cpus = n_cpus())

coda_lasso_pairwise(counts, y, data = NULL, given = NULL,
                    family = c('gaussian', 'binomial'), margin = 1L,
                    pseudocount = NULL, l1_ratio = 0.9, lambda = NULL,
                    n_lambdas = 100L, lambda_min_ratio = NULL,
                    standardize = <coda4microbiome's choice>,
                    n_folds = 10L, folds = NULL, select = c('1se', 'min'),
                    tol = 1e-7, max_iters = 1e5L, seed = 0L, cpus = n_cpus())
```

- **`y`:** the outcome. It can be a vector with one value per sample, or a column name of `data`
  (resolved with `match_rows()` and `per_sample()`, as `groups` is). Two
  levels (a factor, logical or character) imply `family = 'binomial'`, a
  numeric vector `'gaussian'`. An explicit `family` is checked against `y`.
  Name it `y`, the glossary's name for a second input of a different kind.
  It is not `groups`, because a numeric outcome is not a grouping.
- **`given`:** covariates to adjust for (age, sex, BMI). This is the glossary's
  "partialled out" meaning. It can be a one-sided formula over `data` or a data
  frame. coda4microbiome fits the covariates first and passes their
  linear predictor as an **offset**. zeroSum puts them in the model
  unpenalized and outside the constraint (`zeroSum.weights = 0`,
  `penalty.factor = 0`; check this in zeroSum's code). Each function
  follows its reference.
- **`family`:** gaussian and binomial. The binomial fit uses IRLS on the same
  constrained solver, as zeroSum does: \eqn{w = p(1-p)}, \eqn{z = \eta +
  (y - p)/w}, with weights below \eqn{100\,\epsilon} set to 0.
- **`pseudocount`:** the shared rule in compositional.md (`NULL` means none if
  there are no zeros, otherwise half the smallest non-zero value with a warning). The
  output of `replace_zeros_*()` is the recommended input.
  coda4microbiome's `impute_zeros()` differs (it adds 1 if min > 1, else
  min/2; read its code to see whether it adds to every cell or only to
  zeros). Fixtures pass the same pre-imputed table to both packages.
- **`l1_ratio`:** the elastic-net mixing parameter, using scikit-learn's name.
  **Not `alpha`**, which means a significance level in `ancom()` and
  `forward_select()` and UniFrac's α in `generalized_unifrac()`. The
  default is 1 (zeroSum) for `coda_lasso()` and 0.9 (coda4microbiome) for
  `coda_lasso_pairwise()`. Penalty:
  \eqn{\lambda[(1-r)/2\,\|\beta\|_2^2 + r\|\beta\|_1]}.
- **`lambda`:** `NULL` builds a log-spaced path of `n_lambdas` values from
  \eqn{\lambda_{\max}} down to `lambda_min_ratio` × \eqn{\lambda_{\max}}.
  `lambda_min_ratio = NULL` means 0.01 if n < p, otherwise 1e-4, as zeroSum
  and glmnet both do. \eqn{\lambda_{\max}} under the zero-sum constraint
  is \eqn{(\max_j g_j - \min_j g_j) / 2}, where \eqn{g = X^\top(y -
  \bar y)/n} (the research agent's derivation). **Read zeroSum's
  `lambdaMax` code and use its formula.** A user-supplied
  vector is used as given (sorted decreasing).
- **`standardize`:** zeroSum's default is `FALSE`. With `TRUE` the
  constraint becomes \eqn{\sum_j s_j \beta_j = 0} on the standardised scale.
  zeroSum handles this through its `zeroSum.weights`; port that.
- **CV:**
  - `n_folds` sets the number of folds. `folds` takes a per-sample fold
    vector instead, which overrides `seed`. Passing the same `folds`
    (glmnet's `foldid`) is how fixtures match exactly, just as a
    `permutations` matrix makes p-values match vegan's.
  - Folds are drawn with PCG, stratified by class for binomial if zeroSum
    or coda4microbiome stratifies (check).
  - `select` picks `lambda.min` or `lambda.1se`. Each function defaults
    to its reference: zeroSum uses min, coda4microbiome uses 1se. **Ask
    the user** whether both should share one default instead.
  - The CV error is MSE for gaussian and deviance for binomial, unless
    the reference uses something else (coda4microbiome's
    `type.measure`; zeroSum's CV loss). Read the code.
- **`tol`, `max_iters`:** convergence tolerance and cap on outer iterations,
  as in each reference (zeroSum `precision = 1e-8`, 100 outer steps;
  glmnet `thresh = 1e-7`, `maxit = 1e5`).

### `coda_lasso()` algorithm (port zeroSum)

zeroSum's files are `src/zeroSumMoves.cpp` and `src/zeroSumCoordinateDescent.cpp`
(also the IRLS and CV drivers). Its key steps:

1. **Pair moves:** coordinate descent over *pairs* \eqn{(k, s)}.
   Increasing \eqn{\beta_k} by \eqn{\delta} decreases \eqn{\beta_s} by
   \eqn{(u_k/u_s)\delta}, so the constraint always holds. Each move is a
   one-dimensional quadratic plus two absolute values. Try the four sign
   patterns as candidate soft-threshold solutions and keep the lowest
   cost.
2. **Partner search:** inactive coordinates are paired with active ones.
3. **Rotated 3-coordinate moves** when the success rate drops to 0.90 or
   below.
4. **A final "polish" local search** (simulated annealing) on the active
   set. This step is random. Port it on PCG if fixtures need it. The problem is convex,
   so a converged solution should match zeroSum to the convergence tolerance
   whether or not polish runs. Measure this first, and drop polish if it
   isn't needed.
5. **Stopping:** the inner loop stops when |Δcost| < `tol`. There are at most
   `max_iters` outer (IRLS) steps.
6. **Warm starts** along the λ path. zeroSum's `cvStop = 0.1` (stop the path
   early once the CV error rises) is an option worth keeping.
7. **Intercept** update \eqn{\sum w(y - \eta) / \sum w}.

For the binomial case, zeroSum does not clamp p. Keep its behaviour, and note
that a fully separated class gives the same divergence as in zeroSum.

### `coda_lasso_pairwise()` algorithm

coda4microbiome builds the n × p(p−1)/2 matrix of pairwise log ratios
and calls `cv.glmnet`. ecodive should **never build that matrix**:

- Column \eqn{(j,k)} is \eqn{L_{\cdot j} - L_{\cdot k}} with \eqn{L = \log x}.
  So the per-pair gradient is \eqn{g_j - g_k}, computed from the
  per-taxon gradient \eqn{g = L^\top r} in O(np + p²). The fitted values are
  \eqn{X\theta = L\beta}, using the \eqn{\beta} from \eqn{\theta}.
- The column scale for `standardize = TRUE` is
  \eqn{\mathrm{sd}(L_j - L_k) = \sqrt{v_j + v_k - 2c_{jk}}}. This comes from the
  p × p covariance of \eqn{L}, computed once per fold.
- Coordinate descent on \eqn{\theta}: write it from scikit-learn's
  `_cd_fast.pyx` (BSD-3). Take glmnet's conventions from its documentation and
  Friedman, Hastie & Tibshirani (2010, open access), not its code:
  - the objective's \eqn{1/(2n)} scaling;
  - how λ and the penalty relate to `standardize`;
  - the λ sequence;
  - the convergence check;
  - the binomial quadratic approximation.
- Screening is essential, because p = 1,000 already gives 499,500 pairs.
  Use the sequential strong rules (Tibshirani et al. 2012) on \eqn{|g_j
  - g_k|}, with a KKT check over all pairs after convergence, threaded
  over j. Store \eqn{\theta} sparsely (active pairs only).
- Results should match glmnet to about 1e-6 on \eqn{\beta}, given the same λ
  sequence and `folds`. Elastic net with `l1_ratio < 1` is strictly
  convex, so \eqn{\theta} is unique and matchable. With `l1_ratio = 1`
  only \eqn{\beta} is compared.

**First step when implementing:** read coda4microbiome's
`coda4microbiome_functions.R` (MIT, so reading it is allowed). It settles four
things the research agent could not:
- whether `standardize` is left at TRUE;
- `type.measure`;
- the exact zero imputation;
- how `covar` enters, including the offset and how CV handles it.

Record the answers in this plan.

### Threads and determinism

- **The parallel unit is the CV fold:** `n_folds` path fits plus the
  full-data fit, each on its own thread with its own buffers. Results are
  stored per fold and reduced in fold order, so the output is the same for any `cpus`.
- Inside a fit, coordinate descent is sequential. In the pairwise
  function, thread the gradient, strong-rule and KKT sweeps over taxa
  in fixed blocks.
- Randomness comes only from the fold assignment and zeroSum's polish step,
  both on PCG with `seed`. Use one stream per fold for the polish.

### Result

The shared per-feature structure (compositional.md):

- A `data.frame` of class `c('ecodive_coda_lasso', 'data.frame')` (or
  `ecodive_coda_lasso_pairwise`), one row per taxon in table order.
  Columns: `feature`, `coefficient` (\eqn{\beta_j} at the chosen λ),
  `selected` (`coefficient != 0`). It has no `p_value` column: penalized fits
  give none, and the docs say so.
- Attributes:
  - `intercept`, `lambda` (the chosen value), `select`, `family`,
    `l1_ratio`, `pseudocount`, `folds`;
  - `path`: a data frame with `lambda`, `cv_error`, `cv_se`, `n_selected`;
  - `coef_path`: a taxa × λ matrix of \eqn{\beta}, dense (5,000 × 100
    doubles is 4 MB);
  - `fitted`: the per-sample linear predictor, plus probabilities for
    binomial;
  - `performance`: apparent and cross-validated AUC for binomial (mean
    and sd over folds, as coda4microbiome reports), or R² and CV MSE for
    gaussian. AUC is the Mann–Whitney statistic, with no new dependency;
  - pairwise only: `pairs`, a data frame of the non-zero
    \eqn{\theta_{jk}} (`numerator`, `denominator`, `coefficient`).
- **`predict()` S3 methods** for both classes:
  `predict(fit, counts, data = NULL, type = c('link', 'response'))`. They apply the
  same pseudocount and `margin`, and use `given` from `data`. A `print()` method
  shows the selected taxa split into numerator (β > 0) and denominator
  (β < 0) groups, which is how coda4microbiome presents the signature.
- **No plots.** ecodive doesn't plot. The docs show how to draw
  coda4microbiome's signature plot from the result.


## Plot method: signature plot

(Added 2026-11-20 from the compositional-visualization review; it
supersedes the "No plots" bullet above and the Docs bullet on drawing the
signature plot by hand: ecodive has `plot()` methods now.)

```r
plot.ecodive_coda_lasso(x, n_features = 30, ...)   # and, same Rd page,
plot.ecodive_coda_lasso_pairwise(x, n_features = 30, ...)
```

coda4microbiome's two standard figures, side by side
(`layout(matrix(1:2, 1), widths = c(1, 1.2))`, `par()` restored on exit):

- **Left, the signature:** horizontal bars of the non-zero `coefficient`s,
  largest positive at the top, most negative at the bottom; positive
  (numerator) bars `PLOT_HIGHLIGHT`, negative (denominator) bars
  `'royalblue3'` (the colours of the rank plot in
  [log_ratio_plots.md](log_ratio_plots.md)); each value written at its bar
  end to 2 decimals; a line at 0. More than `n_features` selected taxa:
  the largest |β|, with "k of m selected taxa shown" in the subtitle.
  Coefficients are drawn as stored. coda4microbiome rescales them to
  \eqn{\sum|\beta| = 2} before plotting (positive ones sum to 1, negative
  to −1); the docs say this changes neither signs nor order.
- **Right, the signature score by outcome:** the per-sample linear
  predictor (`attr(x, 'fitted')`). Binomial: a box plot per class with the
  samples as points, spread deterministically as in
  `plot.ecodive_log_ratio()` (no RNG). Gaussian: score (x) against the
  observed outcome (y). No fitted line or correlation (plot methods compute
  nothing new); the subtitle gives the stored `performance`: apparent and
  cross-validated AUC (mean ± sd) or apparent R² and CV MSE.
- This needs the outcome in the result: add attribute **`y`** (the outcome
  as fitted, one value per sample; a factor for binomial) to the Result
  list above when implementing.
- No taxa selected at the chosen λ: the left panel says so in text and the
  right panel is left empty (coda4microbiome prints a message and skips
  both plots); no error.
- Returns `x` invisibly. Rd in "Plotting"; tests as in `test-plot.r`
  (gaussian and binomial fits, `n_features = 2`, no selection,
  `.Random.seed` and `par('mfrow')` unchanged, `n_features = 0` error).
- **Reference:** coda4microbiome 0.2.4 (CRAN tarball; MIT, `LICENSE`
  "YEAR: 2022, COPYRIGHT HOLDER: coda4microbiome authors"; checked
  2026-11-20), `R/coda4microbiome_functions.R`: `plot_signature()` (ggpubr
  bar chart, sorted ascending, horizontal, filled by sign, labels rounded to
  2 digits), `plot_prediction()` (binary: box plot by outcome plus a
  density plot; continuous: scatter with a regression line and Pearson r),
  and `coda_glmnet()` (calls both; the ±1 rescaling). Nothing is copied
  (its code is ggplot2), so no notice beyond the credit this plan already
  gives coda4microbiome. Calle et al. (2023) describe the figures; 91
  citations, rising from a small base.


## Reference implementations

- **Port from:**
  - zeroSum 2.0.7 (MIT; holder Thorsten Rehberg 2019): the solver,
    IRLS, λ path and CV driver. Credit it in `inst/COPYRIGHTS` and
    `Authors@R`.
  - coda4microbiome 0.2.4 (MIT; holder "coda4microbiome authors"): the
    model setup, zero handling, covariate offset and performance summary.
    Credit it the same way.
  - scikit-learn (BSD-3): the coordinate-descent update and duality-gap
    check for the pairwise solver.
- **Output-only:** glmnet (GPL-2), FLORAL and Compack (GPL-3), selbal
  (GPL-3).
- **Cross-check (MIT, code readable):** c-lasso 1.0.x (Python, `classo`).
  Its R1 formulation \eqn{\|X\beta - y\|^2 + \lambda\|\beta\|_1} subject to
  \eqn{C\beta = 0}, with \eqn{C = 1^\top} and `intercept=False`, is
  `coda_lasso(l1_ratio = 1)` with λ rescaled by \eqn{2n} and centred data.
  Its exact homotopy path (`path_alg.py`) gives breakpoints to compare
  against.
- **Fixtures:**
  - `data-raw/expected/coda_lasso.r`:
    - zeroSum fits on a small table (about 40 samples × 30 taxa,
      pre-imputed, plus a sparse 200 × 500 case), gaussian and binomial,
      at a fixed λ sequence and fixed `foldid`;
    - coda4microbiome / `cv.glmnet` on the same data and `foldid` for
      the pairwise function;
    - FLORAL as an extra check at one λ.
  - `data-raw/expected/coda_lasso.py`: c-lasso R1 path breakpoints.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Attribute `method` (`'coda_lasso'` or `'coda_lasso_pairwise'`). The
  input is a count table, which has no title, so the plot title is the
  display name the registry gives that `method` ("Compositional Lasso",
  "Compositional Lasso (pairwise log ratios)"); no `display` attribute.
- The chosen `lambda` and `select` (the rule that chose it), `family`
  and `l1_ratio` qualify the coefficients, so they stay top-level
  attributes, as in Result.
- The signature plot reads only the result: the selected coefficients,
  `fitted`, `performance` and the `y` attribute (added for the signature
  plot).

## Tests

- Match zeroSum's coefficients and intercept at about 1e-6 (an iterative
  solver: use the reference's convergence tolerance), and c-lasso's R1
  solution at the same λ. Also check the CV curve and the chosen λ given the same
  `folds`.
- Match coda4microbiome's coefficients, chosen λ, and apparent and CV AUC given the same
  `folds`.
- **The equivalence:** with `l1_ratio = 1` and `standardize = FALSE`,
  `coda_lasso_pairwise(lambda = 2 * l)` and `coda_lasso(lambda = l)` give
  the same \eqn{\beta} to about 1e-6.
- Invariants: \eqn{\sum\beta = 0} to rounding at every λ. Rescaling each
  sample (counts vs proportions) gives the same fit.
- Identical results for `cpus` 1, 2 and 4 and the same `seed`.
- Edge cases:
  - a constant `y`, which gives an all-zero path and a warning;
  - a binomial class with fewer samples than `n_folds`;
  - p = 2;
  - `folds` with an empty fold;
  - NA in `y` (error, or drop with a warning, as `envfit(na.rm)` does);
  - `given` columns that are collinear;
  - `y` with three levels (error).
- `predict()` on the training table reproduces `fitted`.
- One smoke-test line per function (gaussian, small table, `n_folds = 3`).


## Docs

- An Rd page in "Compositional Data", subtitle "Regression and
  prediction".
- Explain the log-contrast model, why the zero-sum constraint makes it
  scale invariant, and the equivalence above. Say when each function
  matches which reference, and that coefficients come with no p-values.
- Show how to draw the signature plot (a bar chart of the non-zero
  coefficients, coloured by sign).
- "Same names in other packages": none known.
- NEWS.


## New argument names (add to [parameter_names.md](parameter_names.md) §1)

| Name | Meaning |
|:-----|:--------|
| `y` | (existing) here, the outcome to predict |
| `family` | the response distribution: `'gaussian'` or `'binomial'` |
| `l1_ratio` | elastic-net mixing (1 = lasso), scikit-learn's name |
| `lambda`, `n_lambdas`, `lambda_min_ratio` | the penalty path |
| `n_folds`, `folds` | cross-validation: a fold count, or a per-sample fold vector |
| `select` | which λ the result reports: `'min'` or `'1se'` |


## References

- Aitchison, J., & Bacon-Shone, J. (1984). Log contrast models for
  experiments with mixtures. *Biometrika*, 71, 323–330.
- Lin, W., Shi, P., Feng, R., & Li, H. (2014). Variable selection in
  regression with compositional covariates. *Biometrika*, 101, 785–797.
- Altenbuchinger, M., Rehberg, T., Zacharias, H. U., et al. (2017).
  Reference point insensitive molecular data analysis. *Bioinformatics*,
  33, 219–226.
- Calle, M. L., Pujolassos, M., & Susin, A. (2023). coda4microbiome:
  compositional data analysis for microbiome cross-sectional and
  longitudinal studies. *BMC Bioinformatics*, 24, 82.
- Susin, A., Wang, Y., Lê Cao, K.-A., & Calle, M. L. (2020). Variable
  selection in microbiome compositional data analysis. *NAR Genomics and
  Bioinformatics*, 2, lqaa029.
- Simpson, L., Combettes, P. L., & Müller, C. L. (2021). c-lasso: a
  Python package for constrained sparse and robust regression and
  classification. *Journal of Open Source Software*, 6, 2844.
- Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths
  for generalized linear models via coordinate descent. *Journal of
  Statistical Software*, 33(1), 1–22.
- Tibshirani, R., Bien, J., Friedman, J., et al. (2012). Strong rules for
  discarding predictors in lasso-type problems. *JRSS B*, 74, 245–266.


## Decisions made while implementing (2026-11-27)

The user's answers (one round, at the start):

1. **Names:** `coda_lasso()` and `coda_lasso_pairwise()`, as planned.
2. **`select`:** each function defaults to its reference: `'min'`
   (zeroSum) for `coda_lasso()`, `'1se'` (coda4microbiome) for
   `coda_lasso_pairwise()`.
3. **Folds:** stratified by class for binomial `y` (random_forest()'s
   `rf_folds()`, StratifiedKFold on PCG), plain KFold for gaussian. Neither
   reference stratifies (both shuffle `rep(1:K)`); fixtures pass `folds`.
4. **Zeros in `coda_lasso_pairwise()`:** ecodive's shared rule
   (`log_pseudocount()`), not coda4microbiome's `impute_zeros()`.

What the reference code settled (the plan's open questions):

- **coda4microbiome 0.2.4's `coda_glmnet()`:** `standardize` is glmnet's
  default (`TRUE`); `type.measure = 'auc'` (binomial) or `'deviance'` (= MSE,
  gaussian); `impute_zeros()` adds 1 to **every** cell if the smallest
  non-zero value is >= 1, else half of it, and only when the table has a
  zero; `covar` enters as an **offset**: `lm(y ~ ., covar)` or
  `glm(..., binomial)` fitted once on all samples, its linear predictor
  passed to `cv.glmnet(offset =)` (so the folds reuse the full-data
  offset) and to `predict(newoffset =)`. It rescales the reported
  coefficients to sum(|beta|) = 2 (we don't), reports glmnet's `cvsd` as
  "sd cv-AUC" (it is a standard error), sets the apparent AUC to 0.5 / R^2
  to 0 when nothing is selected, and reports the deviance as "mean cv-AUC"
  when glmnet falls back to it (we always report the AUC).
- **zeroSum 2.0.7:**
  - `standardize = TRUE` only multiplies the penalty factors `v` by each
    column's sd (x is not scaled; the constraint stays sum(beta) = 0 on the
    original scale), so it is **not** scale invariant: the sds depend on
    the samples' totals. ecodive therefore logs the values as given
    (`log(x + pc)`), not log proportions; without `standardize` the fit is
    invariant either way.
  - The gaussian fit multiplies lambda by sd(y) and divides `v` by it, which
    leaves the lasso term alone and divides the ridge term by sd(y)
    (glmnet's behaviour); ported.
  - `lambdaMax` loops over pairs `k = 2..P-1, s = 1..k-1` (0-based), so it
    **skips the first feature**. ecodive uses every pair; the fixture's data
    are chosen so that the paths agree (and they do, to 4e-15).
  - Its '1se' is not glmnet's: the lambda no smaller than lambda.min whose
    CV error is **closest to** min + SD. Ported as is for `coda_lasso()`.
  - CV error: -2 loglik per fold with the held-out samples weighted
    1/fold size (MSE or mean deviance), averaged over folds; SD
    `sqrt((mean(e^2) - mean(e)^2) / (K - 1))`. Training weights are
    1/(n - fold size).
  - `rotatedUpdates` is off by default, and the polish (`lsSaMove` with
    delta = -beta_k) is **deterministic**: no RNG is needed, so `seed` only
    draws the folds.
  - Covariates: `zeroSum.weights = 0` and `penalty.factor = 0` put them in
    the model unpenalized and outside the constraint (updated by `cdMove`).
    ecodive does that; its lambda max uses the residuals of `y` on the
    covariates (zeroSum's uses the intercept-only model).
  - `cvStop = 0.1` (stop once the CV error has risen for 10% of the path)
    is **not** ported: the whole path is always fitted (the fixtures use
    `cvStop = FALSE`).
- **glmnet (outputs only, checked numerically):** lambda max
  `max |<x_c, r>| / (n a_c max(alpha, 1e-3))`; lambda.min.ratio 0.01 when n
  < number of columns (here p(p-1)/2 pairs); the default path **stops
  early** after lambda m (m >= 5) when dev.ratio > 0.999, or (gaussian)
  dev.ratio_m - dev.ratio_{m-1} < 1e-5 dev.ratio_m, or (binomial) the same
  difference < 1e-5 absolute; with a user lambda it doesn't stop; in
  cv.glmnet the fold fits stop by the same rule and predict from their last
  fit beyond. `cvm` is the fold-size-weighted mean of the per-fold MSE /
  deviance / AUC, `cvsd = sqrt(weighted var / (K - 1))`; binomial deviance
  clamps p to [1e-5, 1 - 1e-5]; `type.measure = 'auc'` falls back to
  deviance (with a warning) when n / K < 10 (ported, with our own warning);
  lambda.1se is the largest lambda within one SE (for AUC, of the max). All
  of these reproduce cv.glmnet's `lambda`, `cvm`, `cvsd` (exactly, given its
  predictions) and `lambda.1se`.

Decided while implementing (my judgement):

- **zeroSum's algorithm, restructured:** each training set (fold or full
  data) runs its whole lambda path on its own thread (`run_parallel()`, a
  fold claimed by an atomic counter; each writes only its own slots; the CV
  error is reduced in fold order), so results are `identical()` for any
  `cpus`. Two speed-ups that leave the sequence of moves alone: (1) a pair
  with a zero coefficient is first tested, from cached gradients, for
  whether its 1-d minimum is "no change" (the four sign patterns then all
  fail), skipping the O(n) update; (2) gaussian fits (fixed weights) do the
  convergence sweeps and the polish with **covariance updates**: a cached
  weighted Gram matrix of the features ever active, gradients updated in
  O(|active|) per move, residuals recomputed by the cost function after each
  sweep. Binomial fits keep zeroSum's direct updates (the weights change
  after every move). Sums use four interleaved partial sums (zeroSum uses
  AVX). On the 200 x 500 fixture (12 lambdas, 5 folds) this is 5.4 s on one
  thread against zeroSum's 7.2 s (and 23 s before the changes), 2.4 s on
  4 threads (zeroSum 3.7 s). With the default 100-lambda path and 10 folds
  on 100 x 100: gaussian 10.7 s (zeroSum 12.2 s), 4.2 s on 4 threads;
  binomial about 10.5 s (zeroSum 10.1 s), 4.0 s on 4 threads: two thirds of
  a binomial fit is zeroSum's refresh of the quadratic approximation (an
  exp() per sample after every successful move), which the port keeps;
  samples outside a fold's training set skip it.
- **zeroSum can't be matched beyond its convergence:** its stop rule is
  |change in cost| < `precision` (1e-8), and the solution is flat in some
  directions, so two correct implementations (AVX sums vs ours) end up to
  ~1e-3 apart in the coefficients while their costs agree to ~1e-11. The
  tests compare coefficients at 2e-3 (gaussian) / 5e-3 (binomial) of the
  largest, CV errors at 5e-3, the chosen lambda exactly, and the costs at
  1e-6 relative (observed ~1e-11). With `tol` much below 1e-10 zeroSum's
  "the cost rose, restore and stop" rule trips on rounding and both
  implementations can stall early (seen at 1e-15); documented by the
  default, not changed.
- **The pairwise solver** is written from Friedman et al. (2010) and the
  strong rules, not ported (scikit-learn's `_cd_fast.pyx` wasn't needed,
  so its credit wasn't extended). It keeps theta sparse (only the working
  set: strong-rule pairs plus every pair ever in it) and never forms the
  pair matrix: column (j,k) is built on the fly, gradients are G_j - G_k
  from G = L'Wr, and with `standardize` the pair sds come from the p x p
  weighted covariance of L (one per thread). The strong-rule and KKT
  scans are O(p^2) per lambda and run inside each fold's thread (the
  parallel unit is the fold, as in `coda_lasso()`), not threaded over
  taxa as the plan suggested. Binomial: IRLS around the coordinate
  descent, p clamped to [1e-5, 1 - 1e-5] (glmnet.control's `pmin`).
  Convergence: the largest xv_c * change^2 below `tol` x the null
  deviance (glmnet's `thresh` semantics).
- **The pairwise path matches glmnet's**, including its early stop, to
  1e-15; CV curves to the solvers' tolerance (both 1e-7); the coefficients
  agree with glmnet at `thresh = 1e-14` to ~6e-3 of the largest at the far
  end of the path (near-collinear pairs: the objective is flat in theta),
  while the objectives agree to ~1e-9. coda_glmnet()'s selected taxa,
  chosen lambda and binomial AUCs (apparent and CV mean / SE) are matched
  exactly; its predictions to ~2e-3.
- **`y`:** numeric means gaussian unless `family = 'binomial'` and `y` is
  0/1; a factor, character (radix-sorted) or logical `y` with two values is
  binomial, the second level being the event; three or more classes error.
  `NA` in `y` (and in `given`) is left out with a warning, as
  `random_forest()`.
- **`given`:** a one-sided formula over `data`, column names of `data`, or a
  data frame (rows matched by name); `model.matrix()` columns; aliased
  columns dropped with a warning (`qr()` of `cbind(1, Z)`). The terms and
  factor levels are kept (attribute `given_model`) for `predict()`.
- **No CV:** `n_folds` 0 or 1 is allowed only with a single `lambda`
  (choosing lambda needs CV); `cv_error` etc. are then `NA`.
- **Constant gaussian `y`:** a warning and an all-zero path (zeroSum would
  divide by sd(y) = 0). A single class is an error.
- **Performance:** gaussian `r_squared` = cor(fitted, y)^2 (coda4microbiome's
  apparent R^2), `cv_mse` and `cv_mse_se` from the path at the chosen lambda;
  binomial `auc` (Mann–Whitney, ties 1/2), `cv_auc` and `cv_auc_se` averaged
  over folds as cv.glmnet does, for both functions.
- **Result:** as planned, plus `cv_measure`, `cv_link` (each sample's
  cross-validated linear predictor at the chosen lambda), `covariates` and
  the internal `given_model`. The display names live in a small registry
  (`CODA_LASSO_NAMES`, keyed by `method`); no `display` title (the title
  doesn't name an input).
- **print():** its own short method (header, lambda and selection,
  performance, then the numerator and denominator features by |beta|, the
  attributes in a footer), since most rows are zeros; a row subset falls
  back to `print_table_result()`.
- **plot():** `type = 'signature'` (the plan's two panels), plus `'path'`
  (coefficients against log lambda) and `'cv'` (the CV curve with SE bars
  and the number selected on top), both marking the chosen lambda.
- **pkgdown:** in "Prediction" next to `random_forest()` (where the README
  lists it), not under "Compositional Data"; the Rd inherits the "Machine
  learning" section.
- **Not done:** c-lasso and FLORAL cross-checks (the objective checks
  against zeroSum and glmnet and the pairwise/constrained equivalence test
  cover the same ground); Cox and multinomial families.


## Outcome (from the implementing session's final report)

Implemented in `1104dfd` on 2026-11-27/28.

- **Files:** `src/coda_lasso.c` (the zeroSum port: pair moves, polish,
  binomial approximation, CV weights and error; gradient shortcut and
  gaussian covariance updates), `src/coda_lasso_pairwise.c` (coordinate
  descent on pairwise log ratios without the pair matrix: strong rules, KKT
  checks, IRLS, glmnet's early path stop), `R/coda_lasso.r` (both
  functions, `predict()`, `print()`, `plot()` with `type = 'signature'`,
  `'path'`, `'cv'`), `test-coda_lasso.r` (224 tests) and two smoke lines,
  fixtures from `data-raw/expected/coda_lasso.r` (261 KB `.rds`).
- **zeroSum 2.0.7** (6 cases on 40 x 30 plus 200 x 500 with 48% zeros;
  gaussian and binomial, lasso and elastic net, standardize, covariates;
  fixed lambda path and folds): costs agree to ~1e-11; coefficients to
  zeroSum's convergence (largest difference 1.5e-3, binomial elastic net;
  2e-4 gaussian; the objective is flat there); CV errors to ~1e-3; the
  chosen lambda (min and zeroSum's 1se) identical; the default lambda path
  to 4e-15. Without the constraint the same solver matches glmnet's lasso
  to 1e-5.
- **coda4microbiome 0.2.4 / cv.glmnet** (gaussian and binomial, with and
  without covariates): lambda path (including early stops) to 1e-15, chosen
  lambda (1se) identical, selected taxa identical, binomial apparent and CV
  AUC (mean, SE) identical, CV curves to ~1e-2 (both at thresh 1e-7),
  coefficients to ~6e-3 of the largest against glmnet at thresh 1e-14 at
  the far end of the path (objectives to ~1e-9), predictions to ~2e-3.
- **Equivalence:** `coda_lasso_pairwise(lambda = 2 l, l1_ratio = 1,
  standardize = FALSE)` and `coda_lasso(lambda = l)` reach the same cost to
  1e-7 and select the same taxa.
- **Determinism:** `identical()` for `cpus` 1-4; sparse (slam, Matrix),
  dense and `margin = 2` input identical. Valgrind: 0 errors.
- **Speed** (4-core container): zeroSum's 200 x 500 case 5.4 s on one
  thread vs zeroSum's 7.2 s (2.4 s vs 3.7 s on 4); 100 x 100, default path,
  10 folds: gaussian 10.7 s vs 12.2 s, binomial ~10.5 s vs 10.1 s (4.0-4.2
  s on 4 threads). Pairwise vs coda_glmnet() on 100 samples: 30 taxa 0.04
  s vs 0.21 s (gaussian), 0.19 s vs 2.1 s (binomial); 100 taxa (4,950 pairs)
  0.31 s vs 0.97 s and 0.52 s vs 2.0 s on one thread, 3x less on 4.
- **Tests:** 35,709 at home and 528 in CRAN mode (35,483 and 526 at
  b90c956). `R CMD check --no-manual --ignore-vignettes`: the existing NOTE
  (installed size) only.
- **Open issues:** c-lasso and FLORAL cross-checks not done; Cox and
  multinomial families not done; binomial `coda_lasso()` is only on par
  with zeroSum on one thread (its per-move exp() refresh, kept).
- **For item 81 (ConQuR):** its "logistic lasso with CV" for presence is
  `coda_lasso()`'s binomial solver without the constraint:
  `ecodive:::cl_zerosum(L, y, Z, fold, K, binom = TRUE, alpha, lambda,
  n_lambdas, lmr, standardize, select, tol, max_iters, cpus, zero_sum =
  FALSE)` (or `C_coda_lasso` directly, with `u` = 1 and `zero_sum = 0`,
  which uses zeroSum's single-coordinate `cdMove`; `v` = 0 leaves
  covariates such as batch unpenalized). It returns the coefficient path,
  intercepts, CV error/SD (deviance) and held-out linear predictors per
  lambda, threaded over folds; folds from `rf_folds()`. A quantile lasso
  would need a different loss (not here).
