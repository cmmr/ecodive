# Volatility of a Variable Over Time

The numbers behind q2-longitudinal's volatility chart (a control chart
of a variable over time), and the chart itself with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html): each group's
mean at each state with a bootstrap confidence interval, global control
limits (the mean \\\pm\\ 2 and 3 standard deviations over all the
samples), and each subject's trajectory.

## Usage

``` r
volatility(
  x,
  subjects,
  states,
  groups = NULL,
  data = NULL,
  conf_level = 0.95,
  n_resamples = 1000L,
  seed = 0L
)
```

## Arguments

- x:

  The variable: a numeric vector with one value per sample (such as an
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  result), or the name of a numeric column of `data`. `NA`s are left out
  of the means and limits. For the features of a count table, call
  `volatility()` for each.

- subjects:

  Optional: the subject of each sample, for the trajectories (q2's
  `individual_id_column`), given as `x` is. No `NA`s. `NULL` for none.

- states:

  The time point of each sample, as a number (days, weeks, or visit
  numbers): given as `subjects` is (q2's `state_column`). At least two
  distinct values; no `NA`s.

- groups:

  Optional: the group of each sample (q2's `default_group_column`),
  given as `x` is: a mean line per group. No `NA`s. Default: `NULL`, one
  group of every sample.

- data:

  Optional sample metadata: a data frame with one row per sample,
  holding the columns named in `x`, `subjects`, `states` and `groups`.
  When `x` names one of its columns, its rows are the samples; otherwise
  they are matched to the names of `x` by row name, or else taken in
  order. Default: `NULL`.

- conf_level:

  The confidence level of the intervals: a number strictly between 0
  and 1. Default: `0.95`, as q2's chart.

- n_resamples:

  The number of bootstrap resamples for each interval. Default: `1000`,
  as q2's chart.

- seed:

  An integer seed for the bootstrap resamples. Default: `0`.

## Value

A list of class `ecodive_volatility`:

- `means`:

  A data frame, one row per group and state (in the groups' order, then
  by state): `group`, `state`, `mean` (of the values that are not `NA`),
  `lower` and `upper` (the bootstrap interval of the mean) and `n` (the
  number of values).

- `global`:

  A named numeric vector over all the samples' values: `mean`, `sd` (the
  sample standard deviation, with \\n - 1\\; 0 for a single value, as
  Vega's), `lower_3sd`, `lower_2sd`, `upper_2sd` and `upper_3sd` (the
  mean \\\pm\\ 2 and 3 SDs).

- `trajectories`:

  With `subjects`, a data frame of each sample (`sample`, `subject`,
  `group`, `state`, `value`), ordered by subject and state; else `NULL`.

Attributes `conf_level` and `n_resamples` describe the intervals.

## Details

**The interval** is Vega's `ci0` and `ci1` aggregates, which q2's chart
draws (Vega's `bootstrapCI()`): `n_resamples` times, the mean of as many
values as the group has, drawn from them with replacement; then the
\\\alpha/2\\ and \\1 - \alpha/2\\ quantiles of those means (R's default
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) type 7),
with \\\alpha = 1 -\\ `conf_level`. Vega draws with JavaScript's
unseeded `Math.random()`, so its intervals change each time the chart is
drawn; ecodive's come from its own seeded generator (PCG, one stream per
group and state), so the same `seed` gives the same intervals, which
differ from Vega's as one draw does from another. A group of one value
has the interval of width 0.

**Differences from q2-longitudinal**: its chart lets the viewer pick the
variable and the grouping interactively; here, one of each per call. Its
tooltips count every sample of a group at a state (Vega's `count`),
those with a missing value too; `n` counts the values used.

**Same names in other packages**: `TTR::volatility()` (of financial
prices) is unrelated; whichever package was attached last masks the
other's, so write `ecodive::volatility()` to be sure.

## References

Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
longitudinal and paired-sample analyses of microbiome data. *mSystems*,
3, e00219-18.
[doi:10.1128/mSystems.00219-18](https://doi.org/10.1128/mSystems.00219-18)

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

Satyanarayan, A., Russell, R., Hoffswell, J., & Heer, J. (2016).
Reactive Vega: a streaming dataflow architecture for declarative
interactive visualization. *IEEE Transactions on Visualization and
Computer Graphics*, 22, 659–668.
[doi:10.1109/TVCG.2015.2467091](https://doi.org/10.1109/TVCG.2015.2467091)

## See also

Other longitudinal:
[`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md),
[`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md),
[`nmit()`](https://cmmr.github.io/ecodive/reference/nmit.md)

## Examples

``` r
    meta <- data.frame(
      subject = rep(c('a', 'b', 'c', 'd'), each = 3),
      day     = rep(c(0, 7, 14), 4),
      diet    = rep(c('x', 'y'), each = 6),
      shannon = c(2.1, 2.4, 2.2, 1.8, 2.0, 2.3, 2.5, 2.2, 2.6, 1.9, 2.1, 2.0) )

    v <- volatility('shannon', 'subject', 'day', 'diet', data = meta)
    v
#> Volatility
#> Global mean 2.18, SD 0.242; limits 1.69 to 2.66 (2 SD), 1.45 to 2.9 (3 SD)
#> Means with 95% bootstrap intervals (1000 resamples):
#> # A data.frame: 6 × 6
#>   group state  mean lower upper     n
#>   <fct> <dbl> <dbl> <dbl> <dbl> <int>
#> 1 x         0  1.95   1.8   2.1     2
#> 2 x         7  2.2    2     2.4     2
#> 3 x        14  2.25   2.2   2.3     2
#> 4 y         0  2.2    1.9   2.5     2
#> 5 y         7  2.15   2.1   2.2     2
#> 6 y        14  2.3    2     2.6     2
#> Trajectories of 4 subjects ($trajectories)
    plot(v)

```
