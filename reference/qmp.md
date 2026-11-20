# Quantitative Microbiome Profiling: Absolute Abundances from a Measured Load

Turn relative abundances into estimates of absolute abundance, using
each sample's measured total microbial **load** (flow cytometry cell
counts, total 16S qPCR, a spike-in, DNA yield, ...).

## Usage

``` r
qmp_scale(counts, load, data = NULL, margin = 1L, cpus = n_cpus())

qmp_rarefy(counts, load, data = NULL, margin = 1L, seed = 0L, cpus = n_cpus())
```

## Arguments

- counts:

  A count table (samples \\\times\\ features); see **Input Types**.
  `qmp_rarefy()` needs whole numbers (integer counts).

- load:

  The total microbial load of each sample, in any absolute unit (cells
  per gram, 16S copies per µl, ...): finite and greater than 0. A
  numeric vector with one value per sample, matched to the samples by
  name when both are named, else in the samples' order. Or, when `data`
  is given, the name of one of its columns (or a vector with one value
  per row of `data`).

- data:

  Optional sample metadata: a data frame with one row per sample,
  holding `load` as a column. Rows are matched to the samples in
  `counts` by row name, and rows for other samples are ignored; without
  row names, they must be in the same order as the samples. Default:
  `NULL`.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

- seed:

  An integer seed for the random rarefaction. The same seed gives the
  same result, for any `cpus`. Default: `0`.

## Value

A table of absolute abundances (doubles), in the form of `counts`, as
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
returns it: the same orientation and dimnames, a base `matrix` for a
base `matrix`, the same Matrix class (`dgCMatrix`, `dgTMatrix`,
`dgeMatrix`) or a slam `simple_triplet_matrix` for those, with zeros not
stored; and a `phyloseq`, `rbiom`, or `SummarizedExperiment` object
gives its count matrix. Each sample's values sum to its `load`.

`qmp_rarefy()`'s result also has attributes `sampling_depth` (each
sample's \\s_i\\) and `target` (its \\t_i\\, the number of reads it was
rarefied to), named by sample. They are attached to every class, sparse
ones too, but subsetting or converting a Matrix or slam object drops
them, so read them first.

## Details

- `qmp_scale()`: "relative abundance \\\times\\ load" (Jian et al. 2020;
  most qPCR and spike-in studies). Each sample's counts are divided by
  its total and multiplied by its load: \$\$q\_{ij} = x\_{ij} / n_i
  \times L_i\$\$

- `qmp_rarefy()`: quantitative microbiome profiling as published by
  Vandeputte et al. (2017). Each sample is first rarefied to the **same
  sampling depth relative to its load** (reads per cell), so that every
  sample has the same detection limit in absolute units, then scaled to
  its load. The steps follow the authors' script, `QMP.R`:

  1.  Each sample's sampling depth is \\s_i = n_i / L_i\\, its total
      number of reads over its load.

  2.  \\s\_{min} = \min_i s_i\\.

  3.  Sample \\i\\ is rarefied, without replacement, to \\t_i = \lfloor
      s\_{min} L_i \rfloor\\ reads (so \\t_i \le n_i\\; the sample with
      the smallest \\s_i\\ keeps all its reads).

  4.  The rarefied counts \\r\_{ij}\\ are scaled to the load: \\q\_{ij}
      = r\_{ij} / t_i \times L_i\\.

## Before and after

- **Neither method makes the data less noisy.** The error of the load
  measurement carries straight into every taxon (Vandeputte et al. 2017;
  Props et al. 2017): measure load in replicate.

- **16S copy-number correction** (dividing each taxon's counts by its
  16S gene copy number) is a separate, earlier step in Vandeputte's
  pipeline, and isn't part of either function. Do it with
  [`sweep()`](https://rdrr.io/r/base/sweep.html) (see the examples). The
  authors' script then rounds the corrected table **up** with
  [`ceiling()`](https://rdrr.io/r/base/Round.html); `qmp_rarefy()`
  requires whole numbers instead, so call `qmp_rarefy(ceiling(x), load)`
  to do as the script does.

- **Spike-ins:** the load is usually
  `spike_added * total_reads / spike_reads`. Drop the spike-in feature
  from `counts` before calling `qmp_*()`.

- **The results are absolute abundances, not compositions.** Don't
  follow them with a centred log ratio (`norm = 'clr'`), which removes
  the load again. Bray-Curtis on QMP data is the usual next step:
  [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) or
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  with `norm = 'none'`.

## Differences from QMP.R

`qmp_rarefy()` follows `rarefy_even_sampling_depth()` in the authors'
`QMP.R` (BSD-3, © 2017 Raes Lab), except:

- It doesn't apply [`ceiling()`](https://rdrr.io/r/base/Round.html) to
  `counts`; non-integer counts are an error.

- Targets \\s\_{min} L_i\\ are floored (as the script's call to
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  does), but a value within \\10^{-9}\\ (relative) of a whole number
  counts as that number, so that floating-point rounding never costs a
  read; the sample with the smallest sampling depth keeps all its reads.

- A sample whose target is 0 reads (a load so small that \\s\_{min} L_i
  \< 1\\) is an error, naming it, where the script returns `NaN`s. So is
  an empty sample.

- Rarefaction uses ecodive's own
  [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) code
  and `seed`, so the random draws differ from the script's
  (`rngseed = 711`), while the sampling depths, targets, and totals are
  the same.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## References

Vandeputte, D., Kathagen, G., D'hoe, K., et al. (2017). Quantitative
microbiome profiling links gut community variation to microbial load.
*Nature*, 551, 507–511.
[doi:10.1038/nature24460](https://doi.org/10.1038/nature24460)

Jian, C., Luukkonen, P., Yki-Järvinen, H., Salonen, A., & Korpela, K.
(2020). Quantitative PCR provides a simple and accessible method for
quantitative microbiota profiling. *PLoS ONE*, 15, e0227285.
[doi:10.1371/journal.pone.0227285](https://doi.org/10.1371/journal.pone.0227285)

Props, R., Kerckhof, F.-M., Rubbens, P., et al. (2017). Absolute
quantification of microbial taxon abundances. *The ISME Journal*, 11,
584–587.
[doi:10.1038/ismej.2016.117](https://doi.org/10.1038/ismej.2016.117)

## See also

[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md),
[`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md), and
the compositional methods, which recover only ratios when no load was
measured: [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md).

## Examples

``` r
    # Cell counts per gram for each sample of ex_counts.
    load <- c(Saliva = 2.0e10, Gums = 4.1e10, Nose = 1.0e10, Stool = 3.2e10)

    # Relative abundance x load
    q <- qmp_scale(ex_counts, load)
    q[1:3, 1:3]
#>        Streptococcus Bacteroides Corynebacterium
#> Saliva    9391304348   115942029               0
#> Gums     36696388262   185101580               0
#> Nose       217606330    19782394      4925816024
    all.equal(rowSums(q), load)
#> [1] TRUE

    # Vandeputte et al.'s QMP: rarefy to an even sampling depth first
    q <- qmp_rarefy(ex_counts, load)
    q[1:3, 1:3]
#>        Streptococcus Bacteroides Corynebacterium
#> Saliva    9391304348   115942029               0
#> Gums     36766619519   173974540               0
#> Nose       232558140           0      4941860465
    attr(q, 'sampling_depth')[1:3]  # reads per cell
#>       Saliva         Gums         Nose 
#> 1.725000e-08 2.160976e-08 1.011000e-07 
    attr(q, 'target')[1:3]          # reads kept
#> Saliva   Gums   Nose 
#>    345    707    172 

    # The load as a column of sample metadata
    meta <- data.frame(cells = load, row.names = names(load))
    q <- qmp_scale(ex_counts, 'cells', data = meta)

    # 16S copy-number correction comes first, as in the authors' script
    copies <- seq(1, 4, length.out = ncol(ex_counts))  # per taxon
    x <- sweep(ex_counts, 2, copies, '/')
    q <- qmp_rarefy(ceiling(x), load)

    # Usually followed by Bray-Curtis distances
    dm <- bray(qmp_rarefy(ex_counts, load))
```
