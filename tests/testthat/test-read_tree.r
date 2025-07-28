test_that("read tree", {
  
  tree_file <- tempfile()
  on.exit({ if (file.exists(tree_file)) unlink(tree_file) })
  writeLines(text = tree_str, con = tree_file)
  
  x <- expect_silent(read_tree(tree_str))
  y <- expect_silent(read_tree(tree_file))
  
  expect_identical(x, y)
  
  expect_identical(
    object   = x, 
    expected = structure(
      list(
        edge = structure(
          .Data = as.integer(c(6,7,7,8,9,9,8,6,7,1,8,9,2,3,4,5)), 
          dim   = as.integer(c(8,2)) ), 
        Nnode       = 4L, 
        tip.label   = c("OTU4", "OTU2", "OTU3", "OTU1", "OTU5"),
        edge.length = c(0.056, 0.676, 0.276, 0.25, 0.548, 0.629, 0.751, 0.432) ), 
      class = "phylo", 
      order = "cladewise" ))
  
  newick <- "((OTU4,(('OTU2',OTU3),OTU1)),OTU5)ROOT;"
  expect_silent(read_tree(newick = newick, underscores = 0))
})
