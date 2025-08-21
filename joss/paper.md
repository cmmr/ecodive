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
date: "18 August 2025"
output:
  word_document: default
  pdf_document: default
tags:
- R
- ecology
- bioinformatics
- microbiome
- diversity
- high-performance computing
affiliations:
- name: The Alkek Center for Metagenomics and Microbiome Research, Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX 77030, USA
  index: 1
- name: Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX, USA
  index: 2
  ror: 02pttbw34
bibliography: paper.bib
---



# Summary

Characterizing the composition of biological communities is a fundamental task
in ecology, but the calculations involved can be computationally prohibitive.
`ecodive` is an R package that addresses this challenge by providing a highly
optimized implementation of common ecological diversity metrics, including
alpha-diversity (within-sample richness and evenness) and beta-diversity
(between-sample dissimilarity). These metrics can incorporate species counts,
relative abundances, and evolutionary relationships, providing a multi-faceted
view of ecological structure. By leveraging a compiled C library with pthreads
for parallelization, `ecodive` delivers substantial performance gains in both
speed and memory usage, enabling researchers to analyze larger datasets more
efficiently.


# Statement of Need

The analysis of ecological diversity in large-scale studies is often hampered by
the computational demands of calculating metrics across thousands of
communities, a common requirement in modern microbiome research. This is
particularly true for phylogenetic metrics like Faith's PD [@Faith1992] and the
UniFrac distance family [@Lozupone2005], which integrate species abundance with
evolutionary data from phylogenetic trees. The resulting high demand on
processing time and memory can limit the scope and scale of scientific inquiry.

`ecodive` overcomes these limitations by offering a significantly faster and
more memory-efficient solution. This allows researchers to analyze more samples,
explore more complex questions, and obtain more robust insights from their data.
By providing a high-performance, parallelized engine for these calculations,
`ecodive` empowers researchers to push the boundaries of large-scale ecological
analysis.



## Comparison to Existing Packages

While numerous R packages can calculate diversity metrics, our comparison
focuses on those that provide their own implementations: `abdiv` [@abdiv],
`adiv` [@adiv], `ampvis2` [@ampvis2], `ecodist` [@ecodist], `entropart`
[@entropart], `GUniFrac` [@GUniFrac], `phyloregion` [@phyloregion], `phyloseq`
[@phyloseq], `picante` [@picante], and `vegan` [@vegan]. `ecodive` sets itself
apart from these packages through its superior performance. Furthermore,
`ecodive` has zero external R dependencies. This makes it a lightweight, stable,
and secure computational backend, minimizing installation conflicts and
simplifying long-term maintenance for developers who build upon it.

Comprehensive benchmarks, conducted using the `bench` package [@bench],
demonstrate these advantages across a range of metrics (Figures 1-3). The
complete benchmark code and results are available in the package vignette
(`vignette('benchmark')`) and online.


![Figure 1: UniFrac benchmarks. `ecodive` demonstrates substantial performance gains for UniFrac, being 2 to 3,900x faster and using 50 - 32,000x less memory, which helps overcome computational bottlenecks in large-scale analyses.](../man/figures/unifrac-benchmark.svg)


![Figure 2: Classic beta diversity benchmarks. `ecodive` is 6 to 2,300x faster and uses 1 to 1,800x less memory, enabling more efficient analysis of community dissimilarities.](../man/figures/bdiv-benchmark.svg)


![Figure 3: Alpha diversity benchmarks. `ecodive` is 2 to 43,000x faster and uses 1 to 33,000x less memory, significantly accelerating the analysis of diversity within single samples.](../man/figures/adiv-benchmark.svg)



# Implemented Metrics

`ecodive` provides a comprehensive suite of alpha and beta diversity metrics.
The current implementation includes:


### Alpha Diversity

* Classic: Shannon Index [@Shannon1948], Simpson Index [@Simpson1949; @Gini1912], Inverse Simpson Index [@Simpson1949], and Chao1 [@Chao1984].

* Phylogenetic: Faith's Phylogenetic Diversity [@Faith1992].


### Beta Diversity

* Classic: Bray-Curtis [@Bray1957; @Sorenson1948], Canberra [@Godfrey1967], Euclidean [@Gower1986; @Legendre2013], Gower [@Gower1971; @Gower1986], Jaccard [@Jaccard1908], Kulczynski [@Kulczynski1927], and Manhattan [@Kaufman1990].

* Phylogenetic: Unweighted UniFrac [@Lozupone2005], Weighted UniFrac [@Lozupone2007], Normalized Weighted UniFrac [@Lozupone2007], Generalized UniFrac [@Chen2012], and Variance Adjusted Weighted UniFrac [@Chang2011].


For the most up-to-date list and detailed descriptions, please refer to the
official ecodive documentation at
<https://cmmr.github.io/ecodive/reference/index.html>.



# Example Usage

ecodive is designed for ease of use and integrates seamlessly with existing
bioinformatics workflows, such as those using phyloseq objects. For example,
calculating weighted UniFrac distances is straightforward:

``` r
library(phyloseq)
data(esophagus)

ecodive::weighted_unifrac(esophagus)
#>           B         C
#> C 0.1050480          
#> D 0.1401124 0.1422409
```



# Acknowledgements

This study was supported by NIH/NIAD (Grant number U19 AI144297), and Baylor
College of Medicine and Alkek Foundation Seed. The authors also acknowledge the
use of Google's Gemini for assistance in refining this manuscript.


# References
