# VALENCIA's reference centroids

The 13 reference centroids of VALENCIA (France et al. 2020), against
which
[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)
classifies vaginal microbiome samples: the mean relative abundances of
199 taxa in the reference samples of each sub-CST.

## Usage

``` r
valencia_centroids
```

## Format

A numeric matrix of 13 sub-CSTs (rows: `"I-A"`, `"I-B"`, `"II"`,
`"III-A"`, `"III-B"`, `"IV-A"`, `"IV-B"`, `"IV-C0"` to `"IV-C4"`, `"V"`)
\\\times\\ 199 taxa (columns, named as VALENCIA names them; see
[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)).
Each row sums to 1, to rounding.

## Source

`CST_centroids_012920.csv` from VALENCIA, commit 8559d45 of
<https://github.com/ravel-lab/VALENCIA>, built by
`data-raw/valencia_centroids.r` in ecodive's source repository. Released
under the MIT licence:

Copyright (c) 2023 ravel-lab

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Details

The centroids were computed from 13,160 taxonomic profiles of 1,975
North American women, from 16S rRNA gene amplicon sequencing; only these
means are included, no individual data. A revised set of centroids (for
instance with all *Prevotella* merged into `"g_Prevotella"`, and the
samples' taxa merged the same way) can be given to
[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md) as
its `centroids`.

## References

France, M. T., Ma, B., Gajer, P., Brown, S., Humphrys, M. S., Holm, J.
B., Waetjen, L. E., Brotman, R. M., & Ravel, J. (2020). VALENCIA: a
nearest centroid classification method for vaginal microbial communities
based on composition. *Microbiome*, 8, 166.
[doi:10.1186/s40168-020-00934-6](https://doi.org/10.1186/s40168-020-00934-6)

## See also

[`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)

## Examples

``` r
    dim(valencia_centroids)
#> [1]  13 199
    round(valencia_centroids[, 1:5], 3)
#>       Lactobacillus_iners Lactobacillus_crispatus Gardnerella_vaginalis
#> I-A                 0.005                   0.958                 0.001
#> I-B                 0.062                   0.716                 0.013
#> II                  0.063                   0.006                 0.076
#> III-A               0.908                   0.003                 0.010
#> III-B               0.609                   0.123                 0.090
#> IV-A                0.056                   0.001                 0.187
#> IV-B                0.113                   0.011                 0.338
#> IV-C0               0.069                   0.040                 0.025
#> IV-C1               0.025                   0.009                 0.022
#> IV-C2               0.016                   0.037                 0.013
#> IV-C3               0.001                   0.004                 0.002
#> IV-C4               0.134                   0.008                 0.040
#> V                   0.252                   0.034                 0.067
#>       Lactobacillus_jensenii Lactobacillus_gasseri
#> I-A                    0.010                 0.001
#> I-B                    0.042                 0.021
#> II                     0.004                 0.658
#> III-A                  0.030                 0.006
#> III-B                  0.022                 0.005
#> IV-A                   0.000                 0.000
#> IV-B                   0.004                 0.003
#> IV-C0                  0.007                 0.016
#> IV-C1                  0.002                 0.015
#> IV-C2                  0.004                 0.019
#> IV-C3                  0.002                 0.007
#> IV-C4                  0.007                 0.007
#> V                      0.526                 0.031
```
