#' Read a newick formatted phylogenetic tree.
#' 
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree from a file or by providing the tree string directly. 
#' This tree must be in Newick format, also known as parenthetic format and
#' New Hampshire format.
#'
#' @param newick   Input data as either a file path, URL, or Newick string. 
#'        Compressed (gzip or bzip2) files are also supported.
#'
#' @param underscores   If `TRUE`, underscores in unquoted names will remain
#'        underscores. If `FALSE`, underscores in unquoted named will be 
#'        converted to spaces.
#'        
#' @return A `phylo` class object representing the tree.
#' 
#' @export
#' @examples
#'     tree <- read_tree("
#'         (A:0.99,((B:0.87,C:0.89):0.51,(((D:0.16,(E:0.83,F:0.96)
#'         :0.94):0.69,(G:0.92,(H:0.62,I:0.85):0.54):0.23):0.74,J:0.1
#'         2):0.43):0.67);")
#'     class(tree)
#'
read_tree <- function (newick, underscores = FALSE) {
  
  validate_args()
  
  #________________________________________________________
  # Get the Newick data into a string
  #________________________________________________________
  if (!startsWith(newick, '(') && !startsWith(newick, '[')) {
    newick <- paste0(collapse = "", readLines(newick, warn = FALSE))
    stopifnot(is.character(newick), length(newick) == 1)
  }
  
  
  #________________________________________________________
  # Remove newlines, comments, and leading whitespace
  #________________________________________________________
  newick <- gsub("[\ \t]*[\r\n]+[\ \t]*", "", newick)
  newick <- gsub("\\[.*?\\]", "",             newick, perl=TRUE)
  newick <- sub("^[\ \t]+",  "",              newick)
  
  stopifnot(isTRUE(nchar(newick) >= 2))
  stopifnot(isTRUE(startsWith(newick, '(')))
  
  
  #________________________________________________________
  # Parse the Newick string into a phylo object
  #________________________________________________________
  tree        <- .Call(C_read_tree, newick, underscores)
  names(tree) <- c('edge', 'Nnode', 'tip.label', 'edge.length', 'node.label')
  attr(tree, 'class') <- 'phylo'
  attr(tree, 'order') <- 'cladewise'
  
  if (all(nchar(tree$node.label) == 0)) tree$node.label  <- NULL
  if (all(tree$edge.length == 0))       tree$edge.length <- NULL
  
  return (tree)
}
