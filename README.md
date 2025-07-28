# fastbiom

<!-- badges: start -->

[![cran](https://www.r-pkg.org/badges/version/fastbiom)](https://CRAN.R-project.org/package=fastbiom)
[![conda](https://anaconda.org/conda-forge/r-fastbiom/badges/version.svg)](https://anaconda.org/conda-forge/r-fastbiom)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/fastbiom)](https://cranlogs.r-pkg.org/)
[![dev](https://github.com/cmmr/fastbiom/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cmmr/fastbiom/actions/workflows/R-CMD-check.yaml)
[![covr](https://codecov.io/gh/cmmr/fastbiom/graph/badge.svg)](https://app.codecov.io/gh/cmmr/fastbiom)
<!-- badges: end -->

`fastbiom` provides the fastest implementations of:

* 5 alpha diversity metrics: Shannon, Simpson, Inverse Simpson, Chao1, and Faith's Phylogenetic Diversity.
* 12 beta diversity metrics: Bray Curtis, Canberra, Euclidean, Gower, Jaccard, Kulczynski, Manhattan, Unweighted UniFrac, Weighted UniFrac, Normalized UniFrac, Generalized UniFrac, and Variance Adjusted UniFrac.
* Count matrix rarefaction
* Newick file parsing


## Installation

The latest stable version can be installed from CRAN.

``` r
install.packages('fastbiom')
```

The development version is available on GitHub.

``` r
install.packages('pak')
pak::pak('cmmr/fastbiom')
```

## Usage

#### Calculate alpha diverity

``` r
ex_counts
#>      A B C  D
#> OTU1 4 0 0  0
#> OTU2 0 8 9 10
#> OTU3 3 0 0  0
#> OTU4 2 0 0  0
#> OTU5 6 5 7  1

shannon(ex_counts)
#>         A         B         C         D 
#> 1.3095258 0.6662784 0.6853142 0.3046361 

faith(ex_counts, tree = ex_tree)
#>   A   B   C   D 
#> 3.4 2.2 2.2 2.2
```

#### Calculate beta diverity

``` r
bray_curtis(ex_counts)
#>           A         B         C
#> B 0.6428571                    
#> C 0.6129032 0.1034483          
#> D 0.9230769 0.2500000 0.2592593

generalized_unifrac(ex_counts, tree = ex_tree, alpha = 0.5)
#>            A          B          C
#> B 0.61036006                      
#> C 0.60260471 0.04873043           
#> D 0.75764452 0.25262174 0.29851111
```


## Documentation

The online manual for `fastbiom` is available at
<https://cmmr.github.io/fastbiom/>. It includes a getting started guide,
articles that explore specific use cases, and reference pages for each
function.


## Automated tests

The following commands will check if `fastbiom` passes the bundled testing
suite.

``` r
install.packages('testthat')
testthat::test_check('fastbiom')
```

## Community guidelines

### Support

Bug reports, feature requests, and general questions can be submitted at
<https://github.com/cmmr/fastbiom/issues>.

### Contributing

Pull requests are welcome. Please ensure contributed code is covered by
tests and documentation (add additional tests and documentation as
needed) and passes all automated tests.

