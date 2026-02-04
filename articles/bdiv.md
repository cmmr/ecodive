# Beta Diversity

## Introduction

Beta diversity is a measure of how different two samples are. The
different metrics described in this vignette quantify that difference,
referred to as the “distance” or “dissimilarity” between a pair of
samples. The distance is typically `0` for identical samples and `1` for
completely different samples.

### Input Data

We will use the `ex_counts` feature table included with ecodive. It
contains the number of observations of each bacterial genera in each
sample. In the text below, you can substitute the word ‘genera’ for the
feature of interest in your own data.

``` r
library(ecodive)

counts <- rarefy(ex_counts)
t(counts)
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  309    6     1
#> Bacteroides            2    2    0   341
#> Corynebacterium        0    0  171     1
#> Haemophilus          180   34    0     1
#> Propionibacterium      1    0   82     0
#> Staphylococcus         0    0   86     1
```

Looking at the matrix above, you can see that saliva and gums are
similar, while saliva and stool are quite different.

## Weighted vs Unweighted Metrics

Before selecting a formula, it is important to distinguish between
**weighted** and **unweighted** metrics.

- **Weighted metrics** take relative abundances into account.
- **Unweighted metrics** only consider presence/absence (binary data).

You can consult
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)
to see which category a specific metric falls into or to list all
available options programmatically.

``` r
list_metrics('beta', 'id', weighted = FALSE)
#> [1] "sorensen"  "hamming"  "jaccard"  "ochiai"  "unweighted_unifrac"

list_metrics('beta', 'id', weighted = TRUE)
#>  [1] "aitchison"            "bhattacharyya"              "bray"                     
#>  [4] "canberra"             "chebyshev"                  "chord"                    
#>  [7] "clark"                "divergence"                 "euclidean"                
#> [10] "generalized_unifrac"  "gower"                      "hellinger"                
#> [13] "horn"                 "jensen"                     "jsd"                      
#> [16] "lorentzian"           "manhattan"                  "matusita"                 
#> [19] "minkowski"            "morisita"                   "motyka"                   
#> [22] "normalized_unifrac"   "psym_chisq"                 "soergel"                  
#> [25] "squared_chisq"        "squared_chord"              "squared_euclidean"        
#> [28] "topsoe"               "variance_adjusted_unifrac"  "wave_hedges"
#> [31] "weighted_unifrac"
```

## Formulas

The following tables detail the mathematical definitions for the metrics
available in `ecodive`.

### Abundance-Based (Weighted)

Given:

- n : The number of features.
- X_i, Y_i : Absolute counts for the i-th feature in samples X and Y.
- X_T, Y_T : Total counts in each sample. X_T = \sum\_{i=1}^{n} X_i
- P_i, Q_i : Proportional abundances of X_i and Y_i. P_i = X_i / X_T
- X_L, Y_L : Mean log of abundances. X_L = \frac{1}{n}\sum\_{i=1}^{n}
  \ln{X_i}
- R_i : The range of the i-th feature across all samples (max - min).

[TABLE]

### Presence / Absence (Unweighted)

Given:

- A, B : Number of features in each sample.
- J : Number of features in common.

[TABLE]

> **Note:**
> [`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)
> is equivalent to `bray(norm = 'binary')`, and
> [`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md) is
> equivalent to `soergel(norm = 'binary')`.

### Phylogenetic

Given n branches with lengths L and a pair of samples’ binary (A and B)
or proportional abundances (P and Q) on each of those branches.

[TABLE]

See <https://cmmr.github.io/ecodive/articles/unifrac.html> for detailed
example UniFrac calculations.

## Partial Calculation

The default value of `pairs=NULL` in ecodive’s beta diversity functions
results in the returned all-vs-all distance matrix being completely
filled in.

``` r
bray(counts)
#>          Saliva      Gums      Nose
#> Gums  0.4260870                    
#> Nose  0.9797101 0.9826087          
#> Stool 0.9884058 0.9884058 0.9913043
```

If you are doing a reference-vs-all comparison, you can use the `pairs`
parameter to skip unwanted calculations and save some CPU time. The
larger the dataset, the more noticeable the improvement will be.

``` r
bray(counts, pairs = 1:3)
#>          Saliva      Gums      Nose
#> Gums  0.4260870                    
#> Nose  0.9797101        NA          
#> Stool 0.9884058        NA        NA
```

The `pairs` argument can be:

- A numeric vector, giving the positions in the result to calculate.
- A logical vector, indicating whether to calculate a position in the
  result.
- A `function(i,j)` that returns whether rows `i` and `j` should be
  compared.

Therefore, all of the following are equivalent:

``` r
bray(counts, pairs = 1:3)
bray(counts, pairs = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
bray(counts, pairs = function (i, j) i == 1)
```

The ordering of `pairs` follows the pairings produced by
[`combn()`](https://rdrr.io/r/utils/combn.html).

``` r
# Column index pairings
combn(nrow(counts), 2)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    1    2    2    3
#> [2,]    2    3    4    3    4    4

# Sample name pairings
combn(rownames(counts), 2)
#>      [,1]     [,2]     [,3]     [,4]   [,5]    [,6]   
#> [1,] "Saliva" "Saliva" "Saliva" "Gums" "Gums"  "Nose" 
#> [2,] "Gums"   "Nose"   "Stool"  "Nose" "Stool" "Stool"
```

So, for instance, to use gums as the reference sample:

``` r
my_combn <- combn(rownames(counts), 2)
my_pairs <- my_combn[1,] == 'Gums' | my_combn[2,] == 'Gums'

my_pairs
#> [1]  TRUE FALSE FALSE  TRUE  TRUE FALSE

bray(counts, pairs = my_pairs)
#>          Saliva      Gums      Nose
#> Gums  0.4260870                    
#> Nose         NA 0.9826087          
#> Stool        NA 0.9884058        NA
```

## References

Levy, A., Shalom, B. R., & Chalamish, M. (2024). A guide to similarity
measures. *arXiv*.

Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures
between probability density functions. *International Journal of
Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.
