# Beta Diversity

## Input Matrix

Here we’ll use the `ex_counts` feature table included with ecodive. It
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

## Beta Diversity

Beta diversity is a measure of how different two samples are.

Looking at the `counts` matrix above, you can easily see that saliva and
gums are similar, while saliva and stool are different. The different
metrics described below quantify that difference, referred to as the
“distance” or “dissimilarity” between a pair of samples. The distance is
`0` for identical samples and `1` for completely different samples.

### Weighted vs Unweighted

Weighted metrics take relative abundances into account, whereas
unweighted metrics only consider presence/absence. To determine which
metrics are weighted or unweighted, consult
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md).

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

Note that `sorsensen()` is equivalent to `bray(norm = 'binary')`, and
[`jaccard()`](https://cmmr.github.io/ecodive/reference/bdiv_functions.md)
is equivalent to `soergel(norm = 'binary')`.

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
