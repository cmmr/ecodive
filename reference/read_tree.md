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

  If `TRUE`, underscores in unquoted labels will remain underscores. If
  `FALSE`, underscores in unquoted labels will be converted to spaces,
  as specified by the Newick standard. Quoted labels are never altered,
  regardless of this setting.

## Value

A `phylo` class object representing the tree.

## Details

Parsing follows the Newick standard described at
<https://evolution.genetics.washington.edu/phylip/newick_doc.html>.

**Quoted labels.** Any tip, internal node, or root label may be enclosed
in single quotes. Everything inside the quotes is taken literally,
including parentheses, square brackets, colons, semicolons, commas,
underscores, blanks, and newlines. The one exception is that two
adjacent single quotes stand for one literal apostrophe: `'it''s'` is
read as `it's`, and `''` is an empty label.

**Comments.** Outside of quotes, any text in square brackets (for
example `[&rate=0.5]`) is a comment and is ignored. Inside quotes,
square brackets are ordinary characters.

**Whitespace.** Outside of quotes, blanks around labels and branch
lengths are ignored. A line break, together with any blanks touching it,
is removed entirely, so that a label or branch length wrapped across two
lines is rejoined.

Only the first tree (up to the first `;`) is read. A branch length on
the root is accepted but discarded. Malformed input, such as an unclosed
quote, comment, or parenthesis, raises an error giving the position of
the problem.

## Examples

``` r
    tree <- read_tree("
        (A:0.99,((B:0.87,C:0.89):0.51,(((D:0.16,(E:0.83,F:0.96)
        :0.94):0.69,(G:0.92,(H:0.62,I:0.85):0.54):0.23):0.74,J:0.1
        2):0.43):0.67);")
    class(tree)
#> [1] "phylo"

    # Quoted labels can contain any character.
    tree <- read_tree("('it''s':1,'x[1]':2,b_c[a comment]:3);")
    tree$tip.label
#> [1] "it's" "x[1]" "b c" 
```
