# ecodive <img src="man/figures/logo.png" align="right" width="174" height="200" alt="ecodive logo" />

<!-- badges: start -->
[![cran](https://www.r-pkg.org/badges/version/ecodive)](https://CRAN.R-project.org/package=ecodive)
[![conda](https://anaconda.org/conda-forge/r-ecodive/badges/version.svg)](https://anaconda.org/conda-forge/r-ecodive)
[![covr](https://codecov.io/gh/cmmr/ecodive/graph/badge.svg)](https://app.codecov.io/gh/cmmr/ecodive)
[![joss](https://joss.theoj.org/papers/a93132f1a403729a2973a8fcc2be3685/status.svg)](https://joss.theoj.org/papers/a93132f1a403729a2973a8fcc2be3685)
<!-- badges: end -->

`ecodive` is an R package for calculating ecological diversity metrics in a
parallelized and memory-efficient manner. It is designed to handle large
datasets, such as those common in microbiome research.


## Why ecodive?

Analyzing ecological diversity is often a computational bottleneck, especially
with large datasets. `ecodive` addresses this by providing:

* **High Performance:** `ecodive` is written in C and parallelized using pthreads, making it dramatically faster than other R packages. [Benchmarks](https://cmmr.github.io/ecodive/articles/benchmark.html) show it can be up to **10,000x faster** and use up to **90,000x less memory**.

* **Zero Dependencies:** The package has no external R dependencies, making it lightweight, stable, and easy to install. This also makes it an ideal and secure backend for other R packages.

* **Comprehensive Metrics:** It implements a wide range of common alpha and beta diversity metrics, including classic and phylogenetic-aware methods like **Faith's PD** and the **complete UniFrac family**.

* **Ease of Use:** The API is simple and integrates seamlessly with popular bioinformatics packages like [`phyloseq`](http://joey711.github.io/phyloseq/) and [`rbiom`](https://cmmr.github.io/rbiom/).



## Installation

The latest stable version can be installed from CRAN.

``` r
install.packages('ecodive')
```

The development version is available on GitHub.

``` r
install.packages('pak')
pak::pak('cmmr/ecodive')
```


## Usage

`ecodive` functions are straightforward to use. Here are a few examples.


### With `phyloseq` or `rbiom` objects

The easiest way to use `ecodive` is with a `phyloseq` or `rbiom` object. These
objects conveniently bundle the count data and phylogenetic tree.

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
counts[,1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva           162           2               0         180
#> Gums             309           2               0          34
#> Nose               6           0             171           0
#> Stool              1         341               1           1


## Alpha Diversity -------------------
shannon(counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.35692121 1.10615349 0.07927797 

faith(counts, tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    155    101    202 


## Beta Diversity --------------------
bray(counts)
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


## Available Methods

Use `list_metrics()` to browse the metrics available for calculating diversity.

``` r
# Alpha Diversity
list_metrics('alpha', 'id')
#>  [1] "ace"         "berger"      "brillouin"   "chao1"       "faith"      
#>  [6] "fisher"      "simpson"     "inv_simpson" "margalef"    "mcintosh"   
#> [11] "menhinick"   "observed"    "shannon"     "squares"    

# Beta Diversity
list_metrics('beta', 'id')
#>  [1] "aitchison"                 "bhattacharyya"             "bray"                     
#>  [4] "canberra"                  "chebyshev"                 "chord"                    
#>  [7] "clark"                     "sorensen"                  "divergence"               
#> [10] "euclidean"                 "generalized_unifrac"       "gower"                    
#> [13] "hamming"                   "hellinger"                 "horn"                     
#> [16] "jaccard"                   "jensen"                    "jsd"                      
#> [19] "lorentzian"                "manhattan"                 "matusita"                 
#> [22] "minkowski"                 "morisita"                  "motyka"                   
#> [25] "normalized_unifrac"        "ochiai"                    "psym_chisq"               
#> [28] "soergel"                   "squared_chisq"             "squared_chord"            
#> [31] "squared_euclidean"         "topsoe"                    "unweighted_unifrac"       
#> [34] "variance_adjusted_unifrac" "wave_hedges"               "weighted_unifrac"   
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
