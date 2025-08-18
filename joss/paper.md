---
title: 'ecodive: Parallel and Memory-Efficient R Package for Ecological Diversity Analysis'
authors:
- name: Daniel P Smith
  orcid: "0000-0002-2479-2044"
  corresponding: true
  affiliation: "1, 2"
- name: Sara J Javornik Cregeen
  orcid: "0009-0000-2698-6478"
  affiliation: "1, 2"
- name: Joseph F Petrosino
  orcid: "0000-0002-4046-6898"
  affiliation: "1, 2"
date: "15 August 2025"
output:
  word_document: default
  pdf_document: default
tags:
- R
- microbiome
- unifrac
- bioinformatics
affiliations:
- name: The Alkek Center for Metagenomics and Microbiome Research, Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX 77030, USA
  index: 1
- name: Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX, USA
  index: 2
  ror: 02pttbw34
bibliography: paper.bib
---



# Summary

In ecology, diversity measures the composition of communities and is the first
step toward understanding the role communities play within their environment.
The most common measures of diversity in microbiome research are alpha-diversity
and beta-diversity. While alpha-diversity aims to describe the richness and
evenness of features within a single sample, beta-diversity assesses the
dissimilarities between two or more communities. Diversity calculations may
include the number of species or other features present, relative abundances,
evolutionary relationships, or a combination thereof.

Applying these metrics to large collections of communities, such as thousands of
gut microbiome samples, offers insights into predicting or diagnosing disease
states through ecological "fingerprints," a computationally intensive task.



# Statement of Need

Processing diversity metrics for thousands of communities is computationally
intensive. The speed and memory footprint of these calculations often become a
bottleneck for analysis, limiting the scope and scale of research studies. This
is especially true for Faith's PD [@Faith1992] and UniFrac [@Lozupone2005],
which require complex integration of species counts with evolutionary distances
by traversing phylogenetic trees. A faster and more memory-efficient
implementation enables researchers to analyze a greater number of samples,
leading to more robust and comprehensive insights. The `ecodive` R package
addresses these challenges by employing a compiled C library with pthreads
parallelization to efficiently compute these metrics, offering significant
performance gains.



## Related Works

There are currently ten other R packages that can calculate alpha and beta
diversity metrics: `abdiv` [@abdiv], `adiv` [@adiv], `ampvis2` [@ampvis2],
`ecodist` [@ecodist], `entropart` [@entropart], `GUniFrac` [@GUniFrac],
`phyloregion` [@phyloregion], `phyloseq` [@phyloseq], `picante` [@picante], and
`vegan` [@vegan]. While several R packages offer diversity metric calculations,
`ecodive` distinguishes itself by providing an implementation that is both
significantly faster and more memory efficient. This superior performance,
across various diversity metrics, is demonstrated in Figures 1-3 through
comprehensive benchmarking.

The `bench` R package [@bench] was used to compare `abdiv`, `adiv`, `ampvis2`,
`ecodist`, `ecodive`, `entropart`, `GUniFrac`, `phyloregion`, `phyloseq`,
`picante`, and `vegan`. The benchmarking runs are detailed in the benchmark
vignette, which is available from within R with `vignette('benchmark')` and
online at <https://cmmr.github.io/ecodive/articles/benchmark.html>. Note that
not all R packages offer all diversity metrics.


![Figure 1: UniFrac benchmarks. `ecodive` demonstrates substantial performance gains for UniFrac, being 2 to 3,900x faster and using 50 - 32,000x less memory, which helps overcome computational bottlenecks in large-scale analyses.](../man/figures/unifrac-benchmark.svg)


![Figure 2: Classic beta diversity benchmarks. `ecodive` is 6 to 2,300x faster and uses 1 to 1,800x less memory, enabling more efficient analysis of community dissimilarities.](../man/figures/bdiv-benchmark.svg)


![Figure 3: Alpha diversity benchmarks. `ecodive` is 2 to 43,000x faster and uses 1 to 33,000x less memory, significantly accelerating the analysis of diversity within single samples.](../man/figures/adiv-benchmark.svg)



# Algorithms

The full list of alpha and beta diversity metrics currently implemented by
`ecodive` is provided below. This set of metrics is subject to expansion as new
functionalities are developed. Refer to `ecodive`'s official documentation at
<https://cmmr.github.io/ecodive/reference/index.html> for the most up-to-date list
and detailed descriptions.


### Classic Alpha Diversity

* Shannon Index [@Shannon1948]
* Simpson Index [@Simpson1949; @Gini1912]
* Inverse Simpson Index [@Simpson1949]
* Chao1 [@Chao1984]


### Phylogenetic Alpha Diversity

* Faith's Phylogenetic Diversity [@Faith1992]


### Classic Beta Diversity

* Bray-Curtis Index [@Bray1957; @Sorenson1948]
* Canberra [@Godfrey1967]
* Euclidean [@Gower1986; @Legendre2013]
* Gower [@Gower1971; @Gower1986]
* Jaccard [@Jaccard1908]
* Kulczynski [@Kulczynski1927]
* Manhattan [@Kaufman1990]


### Phylogenetic Beta Diversity

* Unweighted UniFrac [@Lozupone2005]
* Weighted UniFrac [@Lozupone2007]
* Normalized Weighted UniFrac [@Lozupone2007]
* Generalized UniFrac [@Chen2012]
* Variance Adjusted Weighted UniFrac [@Chang2011]


# Usage

Users can easily compute alpha and beta diversity metrics using `ecodive`. For
example, to calculate weighted UniFrac distances with a `phyloseq` object:

``` r
library(phyloseq)
data(esophagus)

ecodive::weighted_unifrac(esophagus)
#>           B         C
#> C 0.1050480          
#> D 0.1401124 0.1422409
```



# Acknowledgements

This study was supported by NIH/NIAD (Grant number U19 AI44297), and Baylor
College of Medicine and Alkek Foundation Seed.

The authors acknowledge usage of Google AI's Gemini app for refinements of the
final draft of this manuscript.


# References
