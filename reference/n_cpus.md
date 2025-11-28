# Number of CPU Cores

A thin wrapper around `parallely::availableCores()`. If the `parallely`
package is not installed, then it falls back to
`parallel::detectCores(all.tests = TRUE, logical = TRUE)`. Returns `1`
if `pthread` support is unavailable or when the number of cpus cannot be
determined.

## Usage

``` r
n_cpus()
```

## Value

A scalar integer, guaranteed to be at least `1`.

## Examples

``` r
    n_cpus()
#> [1] 4
```
