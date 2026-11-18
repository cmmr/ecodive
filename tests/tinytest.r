if (requireNamespace("tinytest", quietly = TRUE)) {

  local({
    # Attach the namespace
    env_name <- "all_ecodive"
    attach(loadNamespace("ecodive"), name = env_name)

    # Guarantee detachment when this local block finishes or fails
    on.exit(detach(env_name, character.only = TRUE))

    # Two tiers.
    #  - At home (locally and on GitHub CI, which set NOT_CRAN=true, as
    #    devtools does): every test, including comparisons with stored
    #    reference results and with references computed by R.
    #  - On CRAN: only the smoke tests in test-smoke.r, which check that
    #    every function runs and returns results of the right shape. They
    #    compare nothing with stored results or with R's BLAS/LAPACK, whose
    #    builds differ between CRAN's machines.
    at_home <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

    # Off CRAN's machines without --as-cran, cap threads at 2 as CRAN does.
    if (!at_home) {
      old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
      on.exit(
        if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
        else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit),
        add = TRUE )
    }

    # Run the tests
    test_dir <- ifelse(dir.exists("tinytests"), "tinytests", "tests/tinytests")
    source(file = file.path(test_dir, 'helper.r'))

    results <- if (at_home) {
      tinytest::run_test_dir(test_dir, at_home = TRUE)
    } else {
      tinytest::run_test_file(file.path(test_dir, "test-smoke.r"), at_home = FALSE)
    }
    passed <- as.logical(results)

    print(results)

    if (!all(passed)) {
      stop("One or more tinytest tests failed.")
    }
  })

}
