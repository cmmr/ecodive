# `coda_lasso()` and `coda_lasso_pairwise()`: penalized log-contrast regression

**Priority:** 36. **Tier:** port. `coda_lasso()` ports zeroSum 2.0.7's
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
