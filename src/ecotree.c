// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


/*
// Print ecotree_t for C-level debugging
static void print_ecotree(ecotree_t *et) {
  
  Rprintf("double edge_lengths_static[] = {\n");
  for (int i = 0; i < et->n_edges; i++) {
    Rprintf("  %g, // index %d\n", et->edge_lengths[i], i);
  }
  Rprintf("};\n\n");
  
  Rprintf("node_t node_vec_static[] = {\n");
  for (int i = 0; i < et->n_edges; i++) {
    node_t *node = et->node_vec + i;
    Rprintf("  { .edge = %d, .parent = %d, .length = %g }, // index %d\n", 
            node->edge, node->parent, node->length, i);
  }
  Rprintf("};\n\n");
  
  Rprintf("ecotree_t et_static = {\n");
  Rprintf("  .n_edges      = %d,\n", et->n_edges);
  Rprintf("  .edge_lengths = edge_lengths_static,\n");
  Rprintf("  .node_vec     = node_vec_static\n");
  Rprintf("};\n\n");
}
*/


// ecotree constructor
ecotree_t* new_ecotree(SEXP sexp_phylo_tree) {
  
  SEXP sexp_edge_mtx     = PROTECT(get(sexp_phylo_tree, "edge"));
  SEXP sexp_edge_lengths = PROTECT(get(sexp_phylo_tree, "edge.length"));
  SEXP sexp_nnode        = PROTECT(get(sexp_phylo_tree, "Nnode"));

  int    *edge_mtx      = INTEGER(sexp_edge_mtx);
  int     n_edges       = nrows(sexp_edge_mtx);
  double *edge_lengths  = REAL(sexp_edge_lengths);
  int     n_internal    = asInteger(sexp_nnode);
  
  ecotree_t *et       = (ecotree_t*) safe_malloc(sizeof(ecotree_t));
  node_t    *node_vec = (node_t*)    safe_malloc(sizeof(node_t) * n_edges);
  
  
  et->node_vec     = node_vec;
  et->n_edges      = n_edges;
  et->edge_lengths = edge_lengths;
  
  
  int n_otus = n_edges + 1 - n_internal;
  
  // sort edge data by child node
  for (int edge = 0; edge < n_edges; edge++) {
    
    int parent = edge_mtx[0 * n_edges + edge] - 2;
    int child  = edge_mtx[1 * n_edges + edge] - 1;
    
    if (child  > n_otus) child--;
    if (parent < n_otus) parent = -1;
    
    node_vec[child].edge   = edge;
    node_vec[child].parent = parent;
    node_vec[child].length = edge_lengths[edge];
  }
  
  UNPROTECT(3);
  return et;
}
