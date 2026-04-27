#test_that("read tree", {
  
  tree_file <- tempfile()
  on.exit({ if (file.exists(tree_file)) unlink(tree_file) })
  writeLines(text = tree_str, con = tree_file)
  
  expect_silent(x <- read_tree(tree_str))
  expect_silent(y <- read_tree(tree_file))
  
  expect_identical(x, y)
  
  expect_identical(
    current = x, 
    target  = structure(
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
  
  
  tree <- "(t9,((t5,t2),(((t10,(t7,t4)),(t6,(t3,t1))),t8)));"
  expect_inherits(tree <- read_tree(tree), 'phylo')
#})
