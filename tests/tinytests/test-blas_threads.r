#test_that("blas_threads", {

  # blas_threads() reads and sets the thread count of R's BLAS/LAPACK
  # library. It is NA where the count can't be read, e.g. with R's reference
  # BLAS (as on CRAN), and never raises an error, warning, or message.

  # TRUE when `expr` runs without an error, warning, or message.
  is_silent <- function (expr) {
    n <- 0L
    withCallingHandlers(
      expr    = tryCatch(expr, error = function (e) n <<- n + 1L),
      warning = function (w) { n <<- n + 1L; invokeRestart('muffleWarning') },
      message = function (m) { n <<- n + 1L; invokeRestart('muffleMessage') })
    n == 0L
  }


  # -----------------------------------------------------------------------
  # Reading
  # -----------------------------------------------------------------------

  expect_true(is_silent(old <- blas_threads()))
  expect_true(is.integer(old) && length(old) == 1L)
  expect_true(is.na(old) || old >= 1L)
  expect_true(withVisible(blas_threads())$visible)
  expect_identical(blas_threads(NULL), old)


  # -----------------------------------------------------------------------
  # Values other than a whole number >= 1 set nothing.
  # -----------------------------------------------------------------------

  bad_n <- list(
    0, -1L, 1.5, NA, NA_integer_, NaN, Inf, 2^31, '2', TRUE, c(1L, 2L),
    integer(0), list(2L), factor(2), as.difftime(2, units = 'secs') )

  for (n in bad_n) {
    expect_true(is_silent(res <- withVisible(blas_threads(n))))
    expect_identical(res$value, old)
    expect_false(res$visible)
    expect_identical(blas_threads(), old)
  }


  # -----------------------------------------------------------------------
  # Setting
  # -----------------------------------------------------------------------

  # Whether the library found, if any, has a count that blas_threads() sets.
  # Some are only read, e.g. OpenBLAS's OpenMP builds before 0.3.34.
  settable <- .Call(C_blas_settable)
  expect_true(isTRUE(settable) || identical(settable, FALSE))

  if (is.na(old)) {

    # No count read (e.g. R's reference BLAS): only 1 could be set.
    expect_true(is_silent(res <- withVisible(blas_threads(2L))))
    expect_identical(res$value, NA_integer_)
    expect_false(res$visible)
    expect_identical(blas_threads(), NA_integer_)

  } else {

    # Counts are only lowered, or raised back to `old`, never above it.
    expect_true(is_silent(res <- withVisible(blas_threads(1L))))
    expect_identical(res$value, old)
    expect_false(res$visible)

    if (settable) {

      expect_identical(blas_threads(), 1L)

      # The count holds through BLAS and LAPACK calls.
      x <- matrix(seq_len(300 * 300) %% 7, 300)
      s <- crossprod(x %*% x)
      expect_true(isSymmetric(s))
      expect_identical(length(eigen(s, symmetric = TRUE, only.values = TRUE)$values), 300L)
      expect_identical(blas_threads(), 1L)

      # Above the starting count is reduced to it. Doubles are accepted.
      expect_identical(blas_threads(old + 1), 1L)
      expect_identical(blas_threads(), old)

      if (old >= 2L) {
        expect_identical(blas_threads(2L), old)
        expect_identical(blas_threads(), 2L)
      }

    } else {
      expect_identical(blas_threads(), old)
    }

    # Restore the count.
    blas_threads(old)
    expect_identical(blas_threads(), old)

    # In a forked child, the count is only read.
    if (settable && .Platform$OS.type == 'unix') {
      # Stop ecodive's thread pool first (it restarts when next needed): a
      # forked child inherits the memory of the pool's threads but not the
      # threads, and leaves by mcexit() without finalizers, so valgrind would
      # report that memory as possibly lost in the child.
      .Call(C_stop_thread_pool)
      job <- parallel::mcparallel({ blas_threads(1L); blas_threads() })
      expect_identical(parallel::mccollect(job)[[1]], old)
      expect_identical(blas_threads(), old)
    }
  }

#})
