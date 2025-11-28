# Read a newick formatted phylogenetic tree.

A phylogenetic tree is required for computing UniFrac distance matrices.
You can load a tree from a file or by providing the tree string
directly. This tree must be in Newick format, also known as parenthetic
format and New Hampshire format.

## Usage

``` r
read_tree(newick, underscores = FALSE)
```

## Arguments

- newick:

  Input data as either a file path, URL, or Newick string. Compressed
  (gzip or bzip2) files are also supported.

- underscores:

  If `TRUE`, underscores in unquoted names will remain underscores. If
  `FALSE`, underscores in unquoted named will be converted to spaces.

## Value

A `phylo` class object representing the tree.

## Examples

``` r
    tree <- read_tree("
        (A:0.99,((B:0.87,C:0.89):0.51,(((D:0.16,(E:0.83,F:0.96)
        :0.94):0.69,(G:0.92,(H:0.62,I:0.85):0.54):0.23):0.74,J:0.1
        2):0.43):0.67);")
    class(tree)
#> [1] "phylo"
```
