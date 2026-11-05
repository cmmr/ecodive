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



#________________________________________________________
# Quoted labels follow the Newick standard.
# https://evolution.genetics.washington.edu/phylip/newick_doc.html
#________________________________________________________

# '' inside quotes is one literal apostrophe.
x <- read_tree("('it''s':1,b:2);")
expect_identical(x$tip.label,   c("it's", "b"))
expect_identical(x$edge.length, c(1, 2))

# Brackets inside quotes are literal; outside quotes they are comments.
x <- read_tree("('x[1]':1,b[a comment]:2);")
expect_identical(x$tip.label,   c("x[1]", "b"))
expect_identical(x$edge.length, c(1, 2))
expect_identical(read_tree("('x]':1,'[y':2);")$tip.label, c("x]", "[y"))

# Greengenes-style provisional names, e.g. 'o__[Saprospirales]'.
x <- read_tree("('o__[Saprospirales]':1,'c__[Cloacamonae]; o__[Cloacamonales]':2);")
expect_identical(x$tip.label, c("o__[Saprospirales]", "c__[Cloacamonae]; o__[Cloacamonales]"))

# Quoting applies to internal node and root labels, too.
x <- read_tree("((A,B)'node 1',C)'it''s';")
expect_identical(x$tip.label,  c("A", "B", "C"))
expect_identical(x$node.label, c("it's", "node 1"))
expect_identical(x$edge, matrix(c(4L, 5L, 5L, 4L, 5L, 1L, 2L, 3L), ncol = 2))

# Quotes inside a comment don't start a quoted label.
x <- read_tree("(A[it's],B);")
expect_identical(x$tip.label, c("A", "B"))
expect_identical(read_tree("[it's] (A[a [nested] 'comment],B)[x'y];")$tip.label, c("A", "B"))

# Branch lengths after a quoted label, with or without surrounding blanks.
x <- read_tree("('a b':0.5,c:1);")
expect_identical(x$tip.label,   c("a b", "c"))
expect_identical(x$edge.length, c(0.5, 1))
x <- read_tree("( 'a b' : 0.5 , c : 1 ) ;")
expect_identical(x$tip.label,   c("a b", "c"))
expect_identical(x$edge.length, c(0.5, 1))

x <- read_tree("(('a b','c d')'x y':0.5,e)'root':1.5;")
expect_identical(x$tip.label,   c("a b", "c d", "e"))
expect_identical(x$node.label,  c("root", "x y"))
expect_identical(x$edge.length, c(0.5, 0, 0, 0)) # root edge is discarded

# An empty quoted label is allowed.
x <- read_tree("('':1,b:2);")
expect_identical(x$tip.label,   c("", "b"))
expect_identical(x$edge.length, c(1, 2))

# Everything inside quotes is literal.
x <- read_tree("('(a,b):c;d':1,'  padded  ':2,'tab\there':3,'new\nline':4);")
expect_identical(x$tip.label, c("(a,b):c;d", "  padded  ", "tab\there", "new\nline"))
expect_identical(x$Nnode, 1L)

# Underscores: quoted labels are never altered.
x <- read_tree("('a_b':1,c_d:2)'e_f';", underscores = FALSE)
expect_identical(x$tip.label,  c("a_b", "c d"))
expect_identical(x$node.label, "e_f")
x <- read_tree("('a_b':1,c_d:2)'e_f';", underscores = TRUE)
expect_identical(x$tip.label,  c("a_b", "c_d"))
expect_identical(x$node.label, "e_f")

# Non-ASCII labels keep their encoding.
x <- read_tree("('naïve':1,café:2);")
expect_identical(x$tip.label, c("naïve", "café"))



#________________________________________________________
# Malformed input gives a clear R error, not a crash or hang.
#________________________________________________________
expect_error(read_tree("('abc:1,b:2);"),   "quoted label is never closed")
expect_error(read_tree("(it's,b);"),       "quoted label is never closed")
expect_error(read_tree("('abc:1,b:2);"),   "at character 2")
expect_error(read_tree("((A,B),C;"),       "unbalanced parentheses")
expect_error(read_tree("(A,B));"),         "unbalanced parentheses")
expect_error(read_tree("(A,B)C);"),        "unbalanced parentheses")
expect_error(read_tree("(A[abc,B);"),      "comment is never closed")
expect_error(read_tree("(A,B]);"),         "without a matching")
expect_error(read_tree("(A,B)C(D);"),      "unexpected text after the root")
expect_error(read_tree("(A,B),(C,D);"),    "outermost parentheses")
expect_error(read_tree("('a'b,c);"),       "unexpected text after a label")
expect_error(read_tree("(a'b',c);"),       "quotes must enclose an entire label")
expect_error(read_tree("(A:1:2,B);"),      "only one branch length")
expect_error(read_tree("(A:x,B);"),        "invalid branch length")
expect_error(read_tree("(A:,B);"),         "missing branch length")
expect_error(read_tree("[only a comment]"), "no tree found")



#________________________________________________________
# Deep trees parse without recursion, in linear time.
#________________________________________________________
n <- 20000
x <- read_tree(paste0(strrep("(", n - 1), "t1:1", paste0(",t", 2:n, ":1)", collapse = ""), ";"))
expect_identical(length(x$tip.label), as.integer(n))
expect_identical(x$Nnode,             as.integer(n - 1))
expect_identical(x$tip.label[c(1, n)], c("t1", paste0("t", n)))
rm(n)



#________________________________________________________
# A tree generated ahead of time with the 'ape' package
# (ape is not a dependency of ecodive):
#
#   set.seed(2026)
#   tr <- ape::rtree(20)
#   tr$node.label <- paste0("n", seq_len(tr$Nnode))
#   ape_newick <- ape::write.tree(tr)
#   ape_tree   <- ape::read.tree(text = ape_newick)
#   dput(ape_tree$edge); dput(ape_tree$edge.length); ...
#________________________________________________________

ape_newick <- paste0(
  "(t19:0.371201802,((((t15:0.03461930342,(t14:0.2945662972,t17:0.5192242875)",
  "n6:0.9286519818)n5:0.08195683337,(t12:0.9718098086,t10:0.6325115012)n7:",
  "0.2363964249)n4:0.4401367765,t6:0.9582691658)n3:0.4345973074,(((t2:",
  "0.7406427329,(t13:0.9577929138,(t5:0.2984141982,t8:0.8955726682)n12:",
  "0.2011366831)n11:0.04589632386)n10:0.7159156343,(((t11:0.08300503041,t20:",
  "0.6303507104)n15:0.3994683053,(t3:0.5723784592,t4:0.6622316933)n16:",
  "0.4650949899)n14:0.4555079846,((t7:0.7235743126,t9:0.8578579435)n18:",
  "0.3789369871,(t16:0.6326981222,t1:0.7181158245)n19:0.186255503)n17:",
  "0.8779082661)n13:0.566075772)n9:0.5231837637,t18:0.541806604)n8:",
  "0.8863539789)n2:0.8747784046)n1;" )

ape_tree <- structure(
  class = "phylo",
  order = "cladewise",
  .Data = list(
    edge = matrix(ncol = 2, c(
      21L, 21L, 22L, 23L, 24L, 25L, 25L, 26L, 26L, 24L, 27L, 27L, 23L, 22L,
      28L, 29L, 30L, 30L, 31L, 31L, 32L, 32L, 29L, 33L, 34L, 35L, 35L, 34L,
      36L, 36L, 33L, 37L, 38L, 38L, 37L, 39L, 39L, 28L,
      1L, 22L, 23L, 24L, 25L, 2L, 26L, 3L, 4L, 27L, 5L, 6L, 7L, 28L, 29L,
      30L, 8L, 31L, 9L, 32L, 10L, 11L, 33L, 34L, 35L, 12L, 13L, 36L, 14L,
      15L, 37L, 38L, 16L, 17L, 39L, 18L, 19L, 20L )),
    Nnode     = 19L,
    tip.label = c(
      "t19", "t15", "t14", "t17", "t12", "t10", "t6", "t2", "t13", "t5",
      "t8", "t11", "t20", "t3", "t4", "t7", "t9", "t16", "t1", "t18" ),
    edge.length = c(
      0.371201802, 0.8747784046, 0.4345973074, 0.4401367765, 0.08195683337,
      0.03461930342, 0.9286519818, 0.2945662972, 0.5192242875, 0.2363964249,
      0.9718098086, 0.6325115012, 0.9582691658, 0.8863539789, 0.5231837637,
      0.7159156343, 0.7406427329, 0.04589632386, 0.9577929138, 0.2011366831,
      0.2984141982, 0.8955726682, 0.566075772, 0.4555079846, 0.3994683053,
      0.08300503041, 0.6303507104, 0.4650949899, 0.5723784592, 0.6622316933,
      0.8779082661, 0.3789369871, 0.7235743126, 0.8578579435, 0.186255503,
      0.6326981222, 0.7181158245, 0.541806604 ),
    node.label = paste0("n", 1:19) ))

# ecodive reads the ape-written tree exactly as ape does.
x <- read_tree(ape_newick)
expect_inherits(x, "phylo")
expect_identical(x$edge,       ape_tree$edge)
expect_identical(x$Nnode,      ape_tree$Nnode)
expect_identical(x$tip.label,  ape_tree$tip.label)
expect_identical(x$node.label, ape_tree$node.label)
expect_equal(x$edge.length,    ape_tree$edge.length)



#________________________________________________________
# Round trip: hostile labels -> spec-quoted Newick -> read_tree().
#________________________________________________________

# Minimal Newick writer that quotes labels per the standard.
newick_quote <- function (x) {
  bare <- nzchar(x) & !grepl("[\\[\\]()':;,_\\s]", x, perl = TRUE)
  ifelse(bare, x, paste0("'", gsub("'", "''", x, fixed = TRUE), "'"))
}
write_newick <- function (tree) {
  n_tips <- length(tree$tip.label)
  labels <- newick_quote(c(tree$tip.label, tree$node.label))
  kids   <- split(
    x = seq_len(nrow(tree$edge)),
    f = factor(tree$edge[,1], levels = seq_len(n_tips + tree$Nnode)) )
  node_str <- function (node) {
    if (node <= n_tips) return (labels[[node]])
    children <- vapply(kids[[node]], FUN.VALUE = "", function (i) {
      paste0(node_str(tree$edge[i, 2]), ":", sprintf("%.17g", tree$edge.length[[i]]))
    })
    paste0("(", paste(children, collapse = ","), ")", labels[[node]])
  }
  paste0(node_str(n_tips + 1L), ";")
}

hostile <- ape_tree
hostile$tip.label <- c(
  "it's", "''", "'", "x[1]", "[not a comment]", "][", "f(x)", ")(",
  "a:b", "0.5:0.25", "semi;colon", ";", "comma,separated", ",",
  "two words", "  padded  ", "under_score", "_", "tab\there\nnewline",
  "all: 'of' [the] (above), _too_;" )
hostile$node.label <- c(
  "root's node", "", "(1)", "[2]", "3:4", "5;6", "7,8", "nine_ten", "  11  ",
  "naïve", "''''", "end'", "'start", "[[nested]]", "a b c", "_x_",
  "q", "r", "s" )

hostile_newick <- write_newick(hostile)
expect_true(grepl("'it''s':", hostile_newick, fixed = TRUE))
expect_true(grepl("'x[1]':",  hostile_newick, fixed = TRUE))

hostile_file <- tempfile(fileext = ".nwk")
writeLines(hostile_newick, hostile_file)

for (underscores in c(FALSE, TRUE)) {
  for (src in c(hostile_newick, hostile_file)) {
    x <- read_tree(src, underscores = underscores)
    expect_identical(x$tip.label,   hostile$tip.label)
    expect_identical(x$node.label,  hostile$node.label)
    expect_identical(x$edge,        hostile$edge)
    expect_identical(x$edge.length, hostile$edge.length)
    expect_identical(x$Nnode,       hostile$Nnode)
  }
}

unlink(hostile_file)
rm(src, underscores, hostile_file, hostile, hostile_newick, ape_tree, ape_newick, newick_quote, write_newick)
