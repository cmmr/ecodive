# documentation

documentation

## Pseudocount

Logs of zero are undefined, so zeros need replacing first. The
pseudocount is added to every value (not only to zeros), as in
`normalize(norm = 'clr')`. With `pseudocount = NULL` (the default), no
pseudocount is used when `counts` has no zeros; otherwise half the
smallest non-zero value is, with a warning.

The choice of pseudocount is not neutral: it can distort the results,
especially for sparse tables. See
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.
