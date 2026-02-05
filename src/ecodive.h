// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#ifndef ECODIVE_H
#define ECODIVE_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h> // R_registerRoutines, R_useDynamicSymbols

#include <inttypes.h> // uint32_t, uint64_t
#include <math.h>     // exp, fabs, floor, log, lgamma, pow, round, sqrt
#include <string.h>   // memcpy, memset, strlen, strncpy
#include <stdlib.h>   // calloc, free, malloc, NULL, qsort, strtod
#include <stddef.h>   // size_t
#include "pcg_basic.h"


// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

// ecomatrix data structure
typedef struct {
  int     n_samples;
  int     n_otus;
  int     nnz;
  int    *sam_vec;
  int    *pos_vec;
  int    *otu_vec;
  double *val_vec;
  double *clr_vec;
  SEXP    sexp_sample_names;
} ecomatrix_t;

// ecotree data structures
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


/* --- ecomatrix.c --- */
ecomatrix_t* new_ecomatrix(SEXP sexp_matrix, SEXP sexp_margin);
double* rw_val_vec(ecomatrix_t *em);
double* rw_clr_vec(ecomatrix_t *em);


/* --- ecotree.c --- */
ecotree_t* new_ecotree(SEXP sexp_phylo_tree);


/* --- get.c --- */
SEXP get(SEXP, const char *);
void set(SEXP, const char *, SEXP);

/* --- memory.c --- */
void  init_n_ptrs(int n);
void* safe_malloc(size_t bytes);
int   is_safe_ptr(void *ptr);
void  free_all(void);
void* free_one(void *ptr);
void* maybe_free_one(void *ptr);

/* --- normalize.c --- */
void normalize(ecomatrix_t *em, SEXP norm, int n_threads_, SEXP sexp_pseudocount);


#endif
