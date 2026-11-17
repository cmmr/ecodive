if (requireNamespace("tinytest", quietly = TRUE)) {

  local({
    # Attach the namespace
    env_name <- "all_ecodive"
    attach(loadNamespace("ecodive"), name = env_name)

    # Guarantee detachment when this local block finishes or fails
    on.exit(detach(env_name, character.only = TRUE))

    # Tests that compare with stored reference results, or that depend on
    # the platform, run only "at home": locally and in CI, which set
    # NOT_CRAN=true (as devtools does). On CRAN, only the tests of the API,
    # input checks, edge cases, and identities within the package run.
    at_home <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

    # Off CRAN's machines without --as-cran, cap threads at 2 as CRAN does,
    # so that tests asking for more threads stay within CRAN's limit.
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

    results <- tinytest::run_test_dir(test_dir, at_home = at_home)
    passed  <- as.logical(results)

    print(results)

    if (!all(passed)) {
      stop("One or more tinytest tests failed.")
    }
  })

}
