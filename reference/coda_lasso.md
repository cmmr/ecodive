# Compositional lasso: penalized log-contrast regression

Predicts an outcome (a measurement, or one of two classes) from a
feature table and picks out the features that matter, in a way that
respects compositionality: the **log-contrast model** (Aitchison &
Bacon-Shone 1984; Lin et al. 2014), fitted with an elastic-net penalty
and the penalty chosen by cross-validation. `coda_lasso()` constrains
the features' coefficients to sum to zero, as zeroSum does
(Altenbuchinger et al. 2017); `coda_lasso_pairwise()` puts the penalty
on all pairwise log ratios, as coda4microbiome's `coda_glmnet()` does
(Calle et al. 2023).

## Usage

``` r
coda_lasso(
  counts,
  y,
  data = NULL,
  given = NULL,
  family = NULL,
  margin = 1L,
  pseudocount = NULL,
  l1_ratio = 1,
  lambda = NULL,
  n_lambdas = 100L,
  lambda_min_ratio = NULL,
  standardize = FALSE,
  n_folds = 10L,
  folds = NULL,
  select = c("min", "1se"),
  tol = 1e-08,
  max_iters = 100L,
  seed = 0L,
  cpus = n_cpus()
)

coda_lasso_pairwise(
  counts,
  y,
  data = NULL,
  given = NULL,
  family = NULL,
  margin = 1L,
  pseudocount = NULL,
  l1_ratio = 0.9,
  lambda = NULL,
  n_lambdas = 100L,
  lambda_min_ratio = NULL,
  standardize = TRUE,
  n_folds = 10L,
  folds = NULL,
  select = c("1se", "min"),
  tol = 1e-07,
  max_iters = 100000L,
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of feature values (samples \\\times\\ features), or
  any of the table types listed under "Input Types". Typically counts;
  any positive scale (counts, proportions, the output of
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md))
  gives the same fit.

- y:

  What to predict: a vector with one value per sample, or the name of a
  column of `data`. Numeric `y` is a measurement
  (`family = 'gaussian'`); a factor, logical or character `y` with two
  values is a class (`'binomial'`; the second level, or `TRUE`, is the
  "event" whose log odds are modelled). Samples whose `y` is `NA` are
  left out, with a warning.

- data:

  A data frame of sample metadata, with one row per sample (matched to
  `counts` by row names when both have them), when `y`, `given` or
  `folds` name its columns. Default: `NULL`

- given:

  Covariates to adjust for (age, sex, BMI): a one-sided formula over the
  columns of `data` (such as `~ age + sex`), the names of columns of
  `data`, or a data frame with one row per sample. `coda_lasso()` fits
  them alongside the features, unpenalized and outside the zero-sum
  constraint (zeroSum's `penalty.factor = 0` and `zeroSum.weights = 0`);
  `coda_lasso_pairwise()` first regresses `y` on them alone and passes
  that fit as an offset (coda4microbiome's `covar`). Samples with `NA`
  in `given` are left out, with a warning; covariates that are linear
  combinations of others are dropped, with a warning. Default: `NULL`

- family:

  `'gaussian'` (least squares) or `'binomial'` (logistic regression), or
  `NULL` to choose from `y` (see `y`). Default: `NULL`

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- l1_ratio:

  The elastic-net mixing parameter, between 0 (ridge) and 1 (lasso): the
  penalty is \\\lambda\[(1 - r)/2 \\\beta\\\_2^2 + r \\\beta\\\_1\]\\.
  glmnet and zeroSum call it `alpha`. Default: `1` for `coda_lasso()`
  (zeroSum's), `0.9` for `coda_lasso_pairwise()` (coda4microbiome's).

- lambda:

  The penalties to fit, or `NULL` for a path of `n_lambdas` values from
  the smallest that makes every coefficient zero down to
  `lambda_min_ratio` times it, evenly spaced on the log scale. Default:
  `NULL`

- n_lambdas:

  The length of the default path. Default: `100`

- lambda_min_ratio:

  The smallest penalty of the default path, as a fraction of the
  largest; `NULL` for 0.01 when there are fewer samples than
  coefficients, else 1e-4 (zeroSum's and glmnet's rule). Default: `NULL`

- standardize:

  Penalize each coefficient in proportion to the standard deviation of
  its column (log abundance, or log ratio), as if the columns were
  scaled to unit variance; the constraint and the reported coefficients
  stay on the original scale. Default: `FALSE` for `coda_lasso()`
  (zeroSum's), `TRUE` for `coda_lasso_pairwise()` (coda4microbiome's,
  through glmnet's default).

- n_folds:

  The number of cross-validation folds, which choose \\\lambda\\; `0` or
  `1` for none (then `lambda` must be a single value). Default: `10`

- folds:

  A fold for each sample instead (any labels; each label is a fold, as
  glmnet's and zeroSum's `foldid`), matched to the samples like `y`.
  Overrides `n_folds` and `seed`. Default: `NULL`

- select:

  Which \\\lambda\\ the result reports: `'min'`, the one with the best
  cross-validated error, or `'1se'`, a larger (sparser) one within a
  standard error of it (see Details). Default: `'min'` for
  `coda_lasso()` (zeroSum's), `'1se'` for `coda_lasso_pairwise()`
  (coda4microbiome's).

- tol:

  The convergence tolerance: zeroSum's `precision` (the change in the
  cost function between sweeps) for `coda_lasso()`; glmnet's `thresh`
  (the largest change of any coefficient's contribution, relative to the
  null deviance) for `coda_lasso_pairwise()`. Default: `1e-8` and
  `1e-7`.

- max_iters:

  The most iterations: zeroSum's outer (active-set) steps for each
  \\\lambda\\ in `coda_lasso()`; glmnet's `maxit`, the most passes over
  the coefficients in all, for `coda_lasso_pairwise()`. Default: `100`
  and `1e5`.

- seed:

  An integer seed for the folds. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A data frame of class `ecodive_coda_lasso` (or
`ecodive_coda_lasso_pairwise`), one row per feature in table order:
`feature`, `coefficient` (\\\beta_j\\ at the chosen \\\lambda\\) and
`selected` (`coefficient != 0`). Its attributes:

- `method`:

  `'coda_lasso'` or `'coda_lasso_pairwise'`.

- `intercept`, `lambda`, `select`, `family`, `l1_ratio`, `pseudocount`:

  The fit's intercept, the chosen \\\lambda\\ and the rule that chose
  it, and the settings that qualify the coefficients.

- `covariates`:

  With `given`: the covariates' coefficients (`coda_lasso()`), or those
  of the offset model (`coda_lasso_pairwise()`, with its intercept).

- `folds`:

  Each sample's fold.

- `path`:

  A data frame: `lambda`, `cv_error`, `cv_se` and `n_selected`, one row
  per \\\lambda\\; `cv_measure` names the error (`'mse'`, `'deviance'`
  or `'auc'`).

- `coef_path`:

  The features' coefficients along the path (features \\\times\\
  \\\lambda\\).

- `fitted`:

  A data frame: `sample`, `link` (the linear predictor, the "signature
  score") and, for binomial `y`, `probability`; `cv_link` holds each
  sample's cross-validated linear predictor.

- `performance`:

  Gaussian: `r_squared` (the squared correlation of the fitted and
  observed values), `cv_mse` and `cv_mse_se`. Binomial: `auc`
  (apparent), `cv_auc` and `cv_auc_se` (the folds' AUCs averaged as
  glmnet does, at the chosen \\\lambda\\).

- `y`:

  The outcome as fitted (a factor for binomial `y`).

- `pairs`:

  `coda_lasso_pairwise()` only: the non-zero \\\theta\_{jk}\\, a data
  frame of `numerator`, `denominator` and `coefficient`.

## Details

### The log-contrast model

For samples \\i\\ and features \\j\\, with \\g\\ the identity (gaussian)
or the logit (binomial), \$\$g(E\[y_i\]) = \beta_0 + \sum_j \beta_j \log
x\_{ij}, \quad \sum_j \beta_j = 0.\$\$ Because the coefficients sum to
zero, multiplying a sample by any constant changes nothing: the fit is
the same on counts, proportions or rarefied counts, and \\\sum_j \beta_j
\log x\_{ij} = \sum_j \beta_j \mathrm{clr}(x_i)\_j\\. The features with
positive coefficients form the numerator, and those with negative ones
the denominator, of a single weighted log ratio (the "signature"; Calle
et al. 2023).

`coda_lasso()` minimizes \$\$-\frac{1}{n}\ell(\beta_0, \beta) + \lambda
\sum_j \left\[\frac{1 - r}{2} v_j^2 \beta_j^2 + r v_j
\|\beta_j\|\right\]\$\$ (\\\ell\\ the log-likelihood, \\r\\ =
`l1_ratio`, \\v_j\\ the column's standard deviation with `standardize`,
else 1) subject to the constraint, by zeroSum's coordinate descent over
pairs of coefficients (raising one while lowering another keeps the sum
at zero). For gaussian `y` the ridge term is divided by the standard
deviation of `y`, as zeroSum and glmnet do.

`coda_lasso_pairwise()` fits the same model through the log ratios of
every pair of features, \\\log(x_j / x_k)\\, \\j \< k\\, as
coda4microbiome does: \\\beta_j = \sum_k \theta\_{jk} - \sum_k
\theta\_{kj}\\ sums to zero by construction, and the penalty falls on
the \\\theta\\. It never forms the \\n \times p(p-1)/2\\ table of log
ratios. With a pure lasso (`l1_ratio = 1`), no standardization and no
covariates, its fit at \\2\lambda\\ is `coda_lasso()`'s at \\\lambda\\,
since the cheapest \\\theta\\ giving a \\\beta\\ has \\\\\theta\\\_1 =
\\\beta\\\_1 / 2\\. With the defaults (standardized, `l1_ratio = 0.9`)
they differ.

### Choosing lambda

Each fold's samples are predicted from a fit to the others along the
whole path. `coda_lasso()` follows zeroSum: the error is the mean
squared error (gaussian) or the binomial deviance per sample, averaged
over folds, with standard error \\\mathrm{sd}/\sqrt{K}\\ over the \\K\\
folds; `'1se'` is the \\\lambda\\ (no smaller than the best) whose error
is closest to the best error plus its standard error.
`coda_lasso_pairwise()` follows `cv.glmnet()` as coda4microbiome calls
it: the mean squared error (gaussian) or the AUC (binomial; the binomial
deviance, with a warning, when there are fewer than 10 samples per fold,
as glmnet does), averaged over folds weighted by their size; `'1se'` is
the largest \\\lambda\\ whose error is within a standard error of the
best. The default path stops early, as glmnet's does, once the fit
explains 99.9% of the deviance or stops improving.

Folds are drawn at random with `seed` (for binomial `y`, within each
class, so that each fold has both). Penalized fits come with no
p-values: a feature's selection is not a test of it.

### Zeros and covariates

Zeros have no logarithm: replace them first (see
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)),
or a `pseudocount` is added (see "Pseudocount" below).

## Differences from zeroSum and coda4microbiome

- **zeroSum:** the largest penalty of the default path includes every
  pair of features (zeroSum's `lambdaMax` skips the first feature), and
  with `given` it is computed from the residuals of `y` on the
  covariates; the path is always fitted to its end (zeroSum's
  `cvStop = 0.1` stops it once the CV error rises); folds are drawn with
  ecodive's random numbers (and stratified by class); `max_iters` caps
  zeroSum's fixed 100 outer steps. Coefficients agree with zeroSum's to
  its convergence tolerance; both reach the same cost to about 1e-10.

- **coda4microbiome:** zeros get ecodive's shared pseudocount rule
  (coda4microbiome adds 1, or half the smallest value below 1, to every
  value); `coefficient` is \\\beta\\ as fitted, not rescaled to
  \\\sum\|\beta\| = 2\\; the apparent \\R^2\\ and AUC are computed
  whether or not features are selected (coda4microbiome reports 0 and
  0.5 when none are); the reported CV AUC is the AUC even when the
  deviance chose \\\lambda\\. glmnet's solution is matched to its
  tolerance, not its iterations.

## Same names in other packages

None known.

## Machine learning

Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
and compositional geometry; don't expect them to improve prediction.
With random forests predicting 65 metadata variables in four datasets,
Yerke et al. (2024) found that alr, clr and ilr (PhILR) "performed
generally slightly worse or only as well as compositionally naïve
transformations", and relative abundances did best, by a small margin.
Compare with `norm = 'percent'` (see
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md))
before choosing a log-ratio predictor table.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## Pseudocount

The `pseudocount` parameter is only relevant when `norm = 'clr'` or
`norm = 'log'`.

Zeros are undefined in the centered log-ratio (CLR) transformation and
in the log of proportions. If `norm` is `'clr'` or `'log'`,
`pseudocount` is `NULL` (the default), and zeros are detected, the
function uses half the minimum non-zero value (`min(x[x>0]) / 2`) and
issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`), or
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
whose result needs no pseudocount.

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

See [`aitchison`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Aitchison, J., & Bacon-Shone, J. (1984). Log contrast models for
experiments with mixtures. *Biometrika*, 71, 323-330.
[doi:10.1093/biomet/71.2.323](https://doi.org/10.1093/biomet/71.2.323)

Lin, W., Shi, P., Feng, R., & Li, H. (2014). Variable selection in
regression with compositional covariates. *Biometrika*, 101, 785-797.
[doi:10.1093/biomet/asu031](https://doi.org/10.1093/biomet/asu031)

Altenbuchinger, M., Rehberg, T., Zacharias, H. U., et al. (2017).
Reference point insensitive molecular data analysis. *Bioinformatics*,
33, 219-226.
[doi:10.1093/bioinformatics/btw598](https://doi.org/10.1093/bioinformatics/btw598)

Calle, M. L., Pujolassos, M., & Susin, A. (2023). coda4microbiome:
compositional data analysis for microbiome cross-sectional and
longitudinal studies. *BMC Bioinformatics*, 24, 82.
[doi:10.1186/s12859-023-05205-3](https://doi.org/10.1186/s12859-023-05205-3)

Susin, A., Wang, Y., Lê Cao, K.-A., & Calle, M. L. (2020). Variable
selection in microbiome compositional data analysis. *NAR Genomics and
Bioinformatics*, 2, lqaa029.
[doi:10.1093/nargab/lqaa029](https://doi.org/10.1093/nargab/lqaa029)

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. *Journal of
Statistical Software*, 33(1), 1-22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01)

Tibshirani, R., Bien, J., Friedman, J., et al. (2012). Strong rules for
discarding predictors in lasso-type problems. *Journal of the Royal
Statistical Society B*, 74, 245-266.
[doi:10.1111/j.1467-9868.2011.01004.x](https://doi.org/10.1111/j.1467-9868.2011.01004.x)

## See also

[`predict.ecodive_coda_lasso()`](https://cmmr.github.io/ecodive/reference/predict.ecodive_coda_lasso.md),
[`plot.ecodive_coda_lasso()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_coda_lasso.md),
[`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md),
[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)

## Examples

``` r
    # Which features predict a measurement? 30 samples, 8 features.
    set.seed(1)
    counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
    y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
    fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
    fit
#> Compositional Lasso, gaussian y: a data.frame 8 × 3
#> lambda = 1.25e-06 (min), l1_ratio = 1: 8 of 8 features selected
#> Apparent R^2 0.828; cross-validated MSE 0.0109 +/- 0.00255
#> Numerator (coefficient > 0): A (1.29), G (0.137), E (0.131), D (0.0703)
#> Denominator (coefficient < 0): B (-1.2), F (-0.348), C (-0.0675), H
#>   (-0.0156)
#> # Attributes: method, intercept, lambda, select, family, l1_ratio,
#> #   pseudocount, folds, path, cv_measure, coef_path, fitted, cv_link,
#> #   performance, y

    # The same with all pairwise log ratios.
    coda_lasso_pairwise(counts, y, n_folds = 5, cpus = 1)
#> Compositional Lasso (pairwise log ratios), gaussian y: a data.frame 8 × 3
#> lambda = 0.0417 (1se), l1_ratio = 0.9: 3 of 8 features selected
#> Apparent R^2 0.787; cross-validated MSE 0.0132 +/- 0.00323
#> Numerator (coefficient > 0): A (0.874)
#> Denominator (coefficient < 0): B (-0.821), F (-0.0533)
#> # Attributes: method, intercept, lambda, select, family, l1_ratio,
#> #   pseudocount, folds, path, cv_measure, coef_path, fitted, cv_link,
#> #   performance, y, pairs
```
