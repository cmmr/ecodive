// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


static const char *tree_str;
static int         underscores;

static size_t n_leafs;
static size_t n_edges;

// Output arrays
static SEXP    sexp_node_label_vec;
static SEXP    sexp_leaf_label_vec;
static int    *edge_mtx;
static double *edge_length_vec;

// Track next open indices in output arrays
static size_t next_edge_index;
static size_t next_node_index;
static size_t next_leaf_index;




static SEXP extract_name(size_t x1, size_t x2) {
  
  char quoted = tree_str[x1] == '\'' && tree_str[x2] == '\'';
  
  // Quoted Name ==> Strip off quote marks
  if (quoted) {
    x1++;
    x2--;
  }
  
  char* node_name_str = calloc(x2 - x1 + 2, sizeof(char));
  strncpy(node_name_str, tree_str + x1, x2 - x1 + 1);
  node_name_str[x2 - x1 + 1] = '\0';
  
  // Unquoted Name ==> Replace underscores with spaces
  if (!quoted && !underscores) {
    for (size_t j = 0; j <= x2 - x1; j++) {
      if (node_name_str[j] == '_') node_name_str[j] = ' ';
    }
  }
  
  SEXP sexp_node_name = mkChar(node_name_str);
  free(node_name_str);
  
  return sexp_node_name;
}




static void recurse_tree(size_t x1, size_t x2, size_t parent) {
  
  
  size_t i, level;
  
  // Trim off whitespace from beginning and end of string section
  while ((tree_str[x1] == ' ' || tree_str[x1] == '\t') && x1 <= x2) x1++;
  while ((tree_str[x2] == ' ' || tree_str[x2] == '\t') && x1 <= x2) x2--;
  
  // Read backwards, extracting name and length if present
  for (i = x2; i >= x1; i--) {
    
    // Ignore special characters inside single quotes
    if (tree_str[i] == '\'') {
      do { i--; } while (tree_str[i] != '\'' && i >= x1);
      continue;
    }
    
    // Text after the colon is the branch length
    if (tree_str[i] == ':') {
      
      if (next_edge_index > 0) {
        char *junk_ptr;
        edge_length_vec[next_edge_index - 1] = strtod(tree_str + i + 1, &junk_ptr);
      }
      
      x2 = i - 1;
    }
    
    // Text after the end-parentheses is the node name
    if (tree_str[i] == ')') {
      
      if (i < x2)
        SET_STRING_ELT(sexp_node_label_vec, next_node_index, extract_name(i + 1, x2));
      
      if (next_edge_index > 0) {
        size_t mtx_row_i = next_edge_index - 1;
        edge_mtx[mtx_row_i + 0]       = parent + n_leafs;
        edge_mtx[mtx_row_i + n_edges] = next_node_index + n_leafs + 1;
      }
      next_node_index++;
      next_edge_index++;
      
      // Peel off parentheses
      x1++;
      x2 = i - 1;
      
      break;
    }
  }
  
  
  // No parens means we're at a leaf
  if (i <= x1) {
    
    if (x1 <= x2)
      SET_STRING_ELT(sexp_leaf_label_vec, next_leaf_index, extract_name(x1, x2));
    
    if (next_edge_index > 0) {
      size_t mtx_row_i = next_edge_index - 1;
      edge_mtx[mtx_row_i + 0]       = parent + n_leafs;
      edge_mtx[mtx_row_i + n_edges] = next_leaf_index + 1;
    }
    next_leaf_index++;
    next_edge_index++;
    
    return;
  }
  
  
  // Recurse into each of the subtrees
  level  = 0;
  parent = next_node_index;
  for (i = x1; i <= x2; i++) {
    
    // Ignore special characters inside single quotes
    if (tree_str[i] == '\'') {
      do { i++; } while (tree_str[i] != '\'' && i <= x2);
      continue;
    }
    
    // Find the other end of the current subclade
    if (tree_str[i] == '(') {
      level++;
    } else if (tree_str[i] == ')') {
      level--;
    } else if (tree_str[i] == ',' && level == 0) {
      recurse_tree(x1, i - 1, parent);
      x1 = i + 1;
    }
    
  }
  recurse_tree(x1, x2, parent);
  
  return;
}





SEXP C_read_tree(SEXP sexp_tree_str, SEXP sexp_underscores) {
  
  tree_str    = CHAR(asChar(sexp_tree_str));
  underscores = asLogical(sexp_underscores);
  
  // Start and End positions of the newick string
  size_t x1 = 0; 
  size_t x2 = strlen(tree_str) - 1;
  
  // Determine how many nodes are in this tree
  size_t n_nodes = 0;
         n_leafs = 1;
  
  for (size_t i = x1; i <= x2; i++) {
    
    // Ignore special characters inside single quotes
    if (tree_str[i] == '\'') {
      do { i++; } while (tree_str[i] != '\'' && i <= x2);
      continue;
    }
    
    // Only consider the first tree in the file
    if (tree_str[i] == ';') {
      x2 = i - 1;
      break;
    }
    
    if (tree_str[i] == '(') n_nodes++;
    if (tree_str[i] == ',') n_leafs++;
  }
  n_edges = n_nodes + n_leafs - 1;
  
  
  // Prepare a phylo-compatible data structure to return to R
  SEXP sexp_result          = PROTECT(allocVector(VECSXP, 5));
  SEXP sexp_edge_mtx        = PROTECT(allocMatrix(INTSXP,  n_edges, 2));
  SEXP sexp_edge_length_vec = PROTECT(allocVector(REALSXP, n_edges));
       sexp_node_label_vec  = PROTECT(allocVector(STRSXP,  n_nodes));
       sexp_leaf_label_vec  = PROTECT(allocVector(STRSXP,  n_leafs));
  
  SET_VECTOR_ELT(sexp_result, 0, sexp_edge_mtx);
  SET_VECTOR_ELT(sexp_result, 1, ScalarInteger(n_nodes));
  SET_VECTOR_ELT(sexp_result, 2, sexp_leaf_label_vec);
  SET_VECTOR_ELT(sexp_result, 3, sexp_edge_length_vec);
  SET_VECTOR_ELT(sexp_result, 4, sexp_node_label_vec);
  
  edge_mtx        = INTEGER(sexp_edge_mtx);
  edge_length_vec = REAL(sexp_edge_length_vec);
  memset(edge_length_vec, 0, n_edges * sizeof(double));
  
  next_edge_index = 0;
  next_node_index = 0;
  next_leaf_index = 0;
  
  // Start recursing at the highest level of parentheses; i.e. the whole tree
  recurse_tree(x1, x2, 0);
  
  UNPROTECT(5);
  
  return sexp_result;
}






