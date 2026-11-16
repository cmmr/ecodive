# Threads of the BLAS Library

Get or set how many threads R's BLAS/LAPACK library uses, while R is
running.

## Usage

``` r
blas_threads(n = NULL)
```

## Arguments

- n:

  `NULL` to only read the count, or a whole number of at least 1 to set
  it to. Any other value sets nothing.

## Value

The thread count before the call: an integer of at least 1, or `NA` when
it can't be read. Invisibly when `n` is given.

## Details

R computes matrix operations such as `%*%`,
[`crossprod()`](https://rdrr.io/r/base/crossprod.html),
[`solve()`](https://rdrr.io/r/base/solve.html),
[`eigen()`](https://rdrr.io/r/base/eigen.html), and
[`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) with the
BLAS/LAPACK library it is linked to. Optimized libraries run each large
operation on several threads. Environment variables such as
`OPENBLAS_NUM_THREADS` set how many before the library loads;
`blas_threads()` changes it afterward.

The count can be read and set with OpenBLAS, Intel MKL, BLIS, and
FlexiBLAS (which passes it on to its backend). It can't with R's
reference BLAS (which is single-threaded), Apple's Accelerate, or any
other library: there `blas_threads()` returns `NA` and setting does
nothing. It never raises an error, warning, or message.

Setting is limited to what is known to be safe. The count is only read,
not set:

- With OpenBLAS's OpenMP builds before version 0.3.34, which take their
  count from OpenMP's, and OpenBLAS builds that don't report whether
  they use OpenMP, such as Ubuntu 20.04's.

- With FlexiBLAS, unless its backend is OpenBLAS (but not an OpenMP
  build), BLIS, or MKL, and only while that backend is in use.

- In a process forked from the one that loaded ecodive, e.g. by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).

Otherwise, the count can be lowered and raised again, but not above the
highest count read from the library, so that it never has to start
threads it didn't already have: a larger `n` is reduced to that. When no
count was read, only `n = 1` is set. OpenMP's thread count is left as it
is.

On Windows, the library is the one in `Rblas.dll` and `Rlapack.dll`,
which can be replaced with copies of an optimized library. If both are,
each copy has a count of its own: both are set, and the larger is
returned.

ecodive's own functions don't use the BLAS; their threads are set by
`cpus`.

## See also

[`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md)

## Examples

``` r
    blas_threads()
#> [1] 4

    # Use one thread for a matrix product, then restore the count.
    # Restoring NA (count unknown) sets nothing.
    old <- blas_threads(1)
    x   <- crossprod(matrix(runif(1e4), 100))
    blas_threads(old)
```
