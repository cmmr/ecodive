if (requireNamespace("tinytest", quietly = TRUE)) {
  
  local({
    # Attach the namespace
    env_name <- "all_ecodive"
    attach(loadNamespace("ecodive"), name = env_name)
    
    # Guarantee detachment when this local block finishes or fails
    on.exit(detach(env_name, character.only = TRUE))
    
    # Run the tests
    test_dir <- ifelse(dir.exists("tinytests"), "tinytests", "tests/tinytests")
    source(file = file.path(test_dir, 'helper.r'))
    
    results <- tinytest::run_test_dir(test_dir)
    passed  <- as.logical(results)
    
    print(results)
    
    if (!all(passed)) {
      stop("One or more tinytest tests failed.")
    }
  })
  
}
