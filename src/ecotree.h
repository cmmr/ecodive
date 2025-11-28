// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#ifndef ECOTREE_H_INCLUDED
#define ECOTREE_H_INCLUDED


#include <R.h>
#include <Rinternals.h>


// ecotree data structure
typedef struct {
  int    edge;
  int    parent;
  double length;
} node_t;

typedef struct {
  int     n_edges;
  double *edge_lengths;
  node_t *node_vec;
} ecotree_t;


// ecotree constructor
ecotree_t* new_ecotree(SEXP sexp_phylo_tree);


#endif
