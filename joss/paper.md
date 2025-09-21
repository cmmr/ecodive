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
date: "19 September 2025"
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
in ecology, but the calculations involved can be computationally prohibitive
when applied to large studies. `ecodive` is an R package that addresses this
challenge by providing highly optimized implementations of common ecological
diversity metrics, including alpha-diversity (within-sample richness and
evenness) and beta-diversity (between-sample dissimilarity). These metrics can
incorporate species counts, relative abundances, and evolutionary relationships,
providing a multi-faceted view of ecological structure. By leveraging a compiled
C library with pthreads for parallelization, `ecodive` delivers substantial
performance gains in both speed and memory usage, enabling researchers to
analyze large datasets quickly and efficiently.



# Statement of Need

A primary challenge in large-scale ecological analysis is the computational
complexity of beta diversity calculations. These algorithms exhibit $O(n^2)$
complexity, meaning their computational cost scales quadratically with the
number of samples ($n$). As microbiome and ecological studies grow to include
thousands of samples, this quadratic scaling creates a significant bottleneck,
demanding immense processing time and memory.

A second challenge is the fragmentation of diversity metrics across numerous R
packages. Researchers often need to install and manage a suite of dependencies
to access the full range of metrics required for a comprehensive analysis,
leading to potential version conflicts and a disjointed workflow.

`ecodive` addresses both of these critical needs. First, it provides a highly
optimized, parallelized C-based engine that dramatically reduces the time and
memory required by these algorithms, enabling the analysis of much larger
datasets. Second, it consolidates a vast collection of alpha and beta diversity
metrics into a single, dependency-free package. By solving the dual problems of
computational inefficiency and methodological fragmentation, `ecodive` empowers
researchers to push the boundaries of large-scale ecological analysis.



# Comparison to Existing Packages

To evaluate its performance, ecodive was benchmarked against numerous R packages
that provide their own implementations of diversity metrics, including `abdiv`
[@abdiv], `adiv` [@adiv], `ampvis2` [@ampvis2], `ecodist` [@ecodist],
`entropart` [@entropart], `GUniFrac` [@GUniFrac], `labdsv` [@labdsv],
`parallelDist` [@parallelDist], `philentropy` [@philentropy], `phyloregion`
[@phyloregion], `phyloseq` [@phyloseq], `picante` [@picante], `tabula`
[@tabula], and `vegan` [@vegan]. The results, conducted using the `bench` R
package, are summarized in Figure 1 and demonstrate `ecodive`'s superior speed
and memory efficiency for each of the metrics tested.

![Figure 1: `ecodive` performance benchmarks. Each point represents an R package, plotted by median calculation time (x-axis) and memory consumption (y-axis) from ten trials. (A) Benchmarks for Shannon Diversity Index, Bray-Curtis Dissimilarity, and Faith's Phylogenetic Diversity. (B) Benchmarks for the UniFrac family of metrics, with different variants distinguished by point shape. Not all packages implement every metric, but `ecodive` is consistently the fastest and most memory-efficient across all tested metrics, often by several orders of magnitude.](figures/fig1.svg){width="100%"}

Crucially, the benchmark suite confirms these performance gains do not come at
the cost of accuracy, as `ecodive` produces numerically identical output to the
other packages. Beyond its computational advantages, `ecodive` has zero external
R dependencies. This makes it a lightweight, stable, and secure backend,
minimizing installation conflicts and simplifying long-term maintenance for
developers who build upon it. The complete benchmark code and results are
available in the package vignette (`vignette('benchmark')`) and online.


# Implemented Metrics

`ecodive` stands out by offering an extensive and diverse collection of over 50
metrics for both alpha and beta diversity analysis, making it a uniquely
comprehensive tool. It provides researchers with a wide array of both
traditional and phylogeny-aware algorithms in a single, high-performance
package. The suite of alpha diversity metrics includes staples like the Shannon
Diversity Index [@Shannon1948] and Chao1 [@Chao1984], important estimators such
as the Abundance-based Coverage Estimator (ACE) [@Chao1992] and Fisher's Alpha
[@Fisher1943], and key phylogeny-aware metrics like Faith's Phylogenetic
Diversity [@Faith1992], offering robust ways to assess within-sample richness
and evenness. For assessing between-sample dissimilarity, `ecodive` implements
essential beta diversity metrics widely used in microbial ecology, including
Bray-Curtis Dissimilarity [@Bray1957; @Sorenson1948], the complete UniFrac
family [@Lozupone2005; @Lozupone2007; @Chen2012; @Chang2011], and the Aitchison
distance [@Aitchison1982] for compositional data analysis. This extensive
collection allows for a thorough and multi-faceted analysis of community
structure.

For the most up-to-date list and detailed descriptions, please refer to the
official `ecodive` documentation at <https://cmmr.github.io/ecodive>.



# Programmatic Use and API

Beyond interactive analysis, `ecodive` is engineered for programmatic use,
making it an ideal backend for applications like R Shiny web apps [@shiny]. The
package includes a `list_metrics()` function that allows developers to
dynamically filter and present available diversity metrics based on specific
criteria. For instance, metrics can be programmatically selected if they are
phylogeny-aware, abundance-weighted, capable of handling non-integer counts, or
are "true metrics" that satisfy the triangle inequality. This powerful API
simplifies the integration of `ecodive` into other software, enabling developers
to build sophisticated tools that offer users tailored diversity analysis
options based on their dataset and analytical needs.



# Example Usage

`ecodive` is designed for ease of use and integrates seamlessly with existing
bioinformatics workflows, such as those using `phyloseq` objects. For example,
calculating weighted UniFrac distances is straightforward:

``` r
data(esophagus, package = 'phyloseq')
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
