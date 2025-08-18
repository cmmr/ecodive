# ecodive <img src="man/figures/logo.png" align="right" width="174" height="200" alt="ecodive logo" />

<!-- badges: start -->
[![dev](https://github.com/cmmr/ecodive/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cmmr/ecodive/actions/workflows/R-CMD-check.yaml)
[![covr](https://codecov.io/gh/cmmr/ecodive/graph/badge.svg)](https://app.codecov.io/gh/cmmr/ecodive)
<!-- badges: end -->

`ecodive` is an R package for calculating ecological diversity metrics in a
parallelized and memory-efficient manner. It is designed to handle large
datasets, such as those common in microbiome research, with significant speed
and memory improvements over existing tools.


## Why ecodive?

Analyzing ecological diversity is often a computational bottleneck, especially
with large datasets. `ecodive` addresses this by providing:

* **High Performance:** `ecodive` is written in C and parallelized using pthreads, making it dramatically faster than other R packages. [Benchmarks](https://cmmr.github.io/ecodive/articles/benchmark.html) show it can be up to **43,000x faster** and use up to **33,000x less memory**.

* **Zero Dependencies:** The package has no external R dependencies, making it lightweight, stable, and easy to install. This also makes it an ideal and secure backend for other R packages.

* **Comprehensive Metrics:** It implements a wide range of common alpha and beta diversity metrics, including classic and phylogenetic-aware methods like **Faith's PD** and the **complete UniFrac family**.

* **Ease of Use:** The API is simple and integrates seamlessly with popular bioinformatics packages like [`phyloseq`](http://joey711.github.io/phyloseq/) and [`rbiom`](https://cmmr.github.io/rbiom/).



## Installation

The latest stable version can be installed from CRAN.

``` r
install.packages('ecodive')
```

The development version is available on GitHub. Please note that this method
requires a compiler - see http://www.rstudio.com/ide/docs/packages/prerequisites
if the installation does not succeed on the first try.

``` r
install.packages('pak')
pak::pak('cmmr/ecodive')
```


## Usage

`ecodive` functions are straightforward to use. Here are a few examples.


### With `phyloseq` objects

The easiest way to use `ecodive` is with a `phyloseq` or `rbiom` object, which
conveniently bundle the count data and phylogenetic tree.

``` r
library(ecodive)
data(esophagus, package = 'phyloseq')
data(hmp50,     package = 'rbiom')

# Calculate weighted UniFrac distance
w_unifrac <- weighted_unifrac(esophagus)
print(w_unifrac)
#>           B         C
#> C 0.1050480          
#> D 0.1401124 0.1422409

# Calculate Faith's Phylogenetic Diversity
faith_pd <- faith(hmp50)
print(faith_pd[1:4])
#>   HMP01   HMP02   HMP03   HMP04 
#> 6.22296 8.59432 8.93375 9.86597 
```


### With basic R objects

You can also provide the count data and phylogenetic tree as separate objects.

The `ex_counts` and `ex_tree` objects are included with `ecodive`.

``` r
## Example Data ----------------------

counts <- rarefy(ex_counts)
counts
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  309    6     1
#> Bacteroides            2    2    0   341
#> Corynebacterium        0    0  171     1
#> Haemophilus          180   34    0     1
#> Propionibacterium      1    0   82     0
#> Staphylococcus         0    0   86     1


## Alpha Diversity -------------------

shannon(counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.35692121 1.10615349 0.07927797 

faith(counts, tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    155    101    202 


## Beta Diversity --------------------

bray_curtis(counts)
#>          Saliva      Gums      Nose
#> Gums  0.4260870                    
#> Nose  0.9797101 0.9826087          
#> Stool 0.9884058 0.9884058 0.9913043

weighted_unifrac(counts, tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums   36.97681                    
#> Nose   67.23768  55.97101          
#> Stool 109.77971 109.44058 110.00870
```



## Documentation

The online manual for ecodive is available at
<https://cmmr.github.io/ecodive/>. It includes a getting started guide,
articles on alpha/beta diversity, and reference pages for each function.



## Automated tests

The following commands will check if ecodive passes the bundled testing
suite.

``` r
install.packages('testthat')
testthat::test_check('ecodive')
```



## Community guidelines


### Support

Bug reports, feature requests, and general questions can be submitted at
<https://github.com/cmmr/ecodive/issues>.


### Contributing

Pull requests are welcome. Please ensure contributed code is covered by
tests and documentation (add additional tests and documentation as
needed) and passes all automated tests.

New functions must leverage C and pthreads to minimize memory and CPU time.

Please note that the ecodive project is released with a [Contributor Code of
Conduct](https://cmmr.github.io/ecodive/CODE_OF_CONDUCT.html). By contributing
to this project, you agree to abide by its terms.
