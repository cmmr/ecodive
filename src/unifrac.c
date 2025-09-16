// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <math.h>   // fabs, pow
#include <stdlib.h> // calloc, free
#include <string.h> // memset
#include "get.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define UNWEIGHTED   1
#define WEIGHTED     2
#define NORMALIZED   3
#define GENERALIZED  4
#define VAR_ADJUSTED 5


typedef struct {
  int edge;
  int parent;
  double length;
} node_t;


//======================================================
// Variables shared between main and worker threads.
//======================================================
static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_edges;
static double *edge_lengths;
static int    *pairs_vec;
static int     n_pairs;
static char    all_pairs;
static double *weight_mtx;
static node_t *nodes;
static double *total_vec;
static int     n_threads;
static double *dist_vec;
static SEXP   *sexp_extra;



/*
 * START_SAMPLE_LOOP and END_SAMPLE_LOOP define a simple loop 
 * for iterating over all samples, ensuring that each is 
 * assigned to only a single thread. Provides the `sample` 
 * index and `weight_vec` for each sample.
 */

#define START_SAMPLE_LOOP                                                \
  int thread_i = *((int *) arg);                                         \
  for (int sample = thread_i; sample < n_samples; sample += n_threads) { \
    double *weight_vec = weight_mtx + (sample * n_edges);                \


#define END_SAMPLE_LOOP                                        \
  }                                                            \
  return NULL;                                                 \



/*
 * The START_PAIR_LOOP and END_PAIR_LOOP macros efficiently 
 * loop through all combinations of samples. Skips unwanted 
 * pairings and pairings not assigned to the current thread. 
 * Ensures all threads process the same number of pairs.
 * 
 * After calling START_PAIR_LOOP the code can expect 
 * `x_weight_vec` and `y_weight_vec` to point to the two 
 * samples' columns in `weight_mtx`.  The code should assign 
 * to `distance` before calling END_PAIR_LOOP.
 * 
 * START_PAIR_TOTAL_LOOP is identical to START_PAIR_LOOP, but
 * also assigns `x_total` and `y_total` from `total_vec`.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */


#define START_PAIR_LOOP                                        \
  int thread_i = *((int *) arg);                               \
  int pair_idx = 0;                                            \
  int dist_idx = 0;                                            \
  for (int i = 0; i < n_samples - 1; i++) {                    \
    double *x_weight_vec = weight_mtx + (i * n_edges);         \
    for (int j = i + 1; j < n_samples; j++) {                  \
      if (all_pairs || pairs_vec[pair_idx] == dist_idx) {      \
        if (pair_idx % n_threads == thread_i) {                \
          double *y_weight_vec = weight_mtx + (j * n_edges);   \
          double distance = 0;                                 \


#define START_PAIR_TOTAL_LOOP                                  \
  int thread_i = *((int *) arg);                               \
  int pair_idx = 0;                                            \
  int dist_idx = 0;                                            \
  for (int i = 0; i < n_samples - 1; i++) {                    \
    double  x_total      = total_vec[i];                       \
    double *x_weight_vec = weight_mtx + (i * n_edges);         \
    for (int j = i + 1; j < n_samples; j++) {                  \
      if (all_pairs || pairs_vec[pair_idx] == dist_idx) {      \
        if (pair_idx % n_threads == thread_i) {                \
          double  y_total      = total_vec[j];                 \
          double *y_weight_vec = weight_mtx + (j * n_edges);   \
          double  distance     = 0;                            \


#define END_PAIR_LOOP                                          \
          dist_vec[dist_idx] = distance;                       \
        }                                                      \
        pair_idx++;                                            \
        if (pair_idx == n_pairs) return NULL;                  \
      }                                                        \
      dist_idx++;                                              \
    }                                                          \
  }                                                            \
  return NULL;                                                 \





//======================================================
// Unweighted UniFrac.
//======================================================
static void *unweighted_mtx (void *arg) {
  START_SAMPLE_LOOP
    
  for (int otu = 0; otu < n_otus; otu++) {
    
    double abundance = otu_mtx[sample + otu * n_samples];
    if (abundance == 0) continue; // OTU not present in sample
    
    int node = otu;     // start at OTU tip/leaf in tree
    while (node > -1) { // traverse until we hit the tree's root
      int edge = nodes[node].edge;
      if (weight_vec[edge]) break; // already traversed
      weight_vec[edge] = 1;
      node = nodes[node].parent;   // proceed on up the tree
    }
    
  }
  
  END_SAMPLE_LOOP
}


static void *unweighted_dist (void *arg) {
  START_PAIR_LOOP
    
  double shared = 0;
  
  for (int edge = 0; edge < n_edges; edge++) {
    
    double x = x_weight_vec[edge];
    double y = y_weight_vec[edge];
    
    if (x || y) {
      if (x && y) { shared   += edge_lengths[edge]; }
      else        { distance += edge_lengths[edge]; }
    }
  }
  
  distance = distance / (distance + shared);
  
  END_PAIR_LOOP
}




//======================================================
// Weighted UniFrac.
//======================================================
static void *weighted_mtx (void *arg) {
  START_SAMPLE_LOOP
    
  double sample_depth = 0;
  for (int otu = 0; otu < n_otus; otu++)
    sample_depth += otu_mtx[sample + otu * n_samples];
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    double abundance = otu_mtx[sample + otu * n_samples];
    if (abundance == 0) continue; // OTU not present in sample
    
    int node = otu;     // start at OTU tip/leaf in tree
    while (node > -1) { // traverse until we hit the tree's root
      
      int    edge   = nodes[node].edge;
      double length = nodes[node].length;
      
      // relative abundance, weighted by branch length
      double relative_abundance = abundance / sample_depth;
      double weighted_abundance = length * relative_abundance;
      weight_vec[edge] += weighted_abundance;
      
      node = nodes[node].parent; // proceed on up the tree
    }
    
  }
  
  END_SAMPLE_LOOP
}

static void *weighted_dist (void *arg) {
  START_PAIR_LOOP
  
  for (int edge = 0; edge < n_edges; edge++) {
    
    double x = x_weight_vec[edge];
    double y = y_weight_vec[edge];
    
    if (x || y) {
      if (x > y) { distance += x - y; }
      else       { distance += y - x; }
    }
  }
  
  END_PAIR_LOOP
}




//======================================================
// Normalized Weighted UniFrac.
//======================================================
static void *normalized_mtx (void *arg) {
  START_SAMPLE_LOOP
    
  double sample_depth = 0;
  for (int otu = 0; otu < n_otus; otu++)
    sample_depth += otu_mtx[sample + otu * n_samples];
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    double abundance = otu_mtx[sample + otu * n_samples];
    if (abundance == 0) continue; // OTU not present in sample
    
    int node = otu;  // start at OTU tip/leaf in tree
    while (node > -1) { // traverse until we hit the tree's root
      
      int    edge   = nodes[node].edge;
      double length = nodes[node].length;
      
      // relative abundance, weighted by branch length
      double relative_abundance = abundance / sample_depth;
      double weighted_abundance = length * relative_abundance;
      weight_vec[edge]  += weighted_abundance;
      total_vec[sample] += weighted_abundance;
      
      node = nodes[node].parent; // proceed on up the tree
    }
    
  }
  
  END_SAMPLE_LOOP
}

static void *normalized_dist (void *arg) {
  START_PAIR_TOTAL_LOOP
    
  for (int edge = 0; edge < n_edges; edge++) {
    
    double x = x_weight_vec[edge];
    double y = y_weight_vec[edge];
    
    if (x || y) {
      if (x > y) { distance += x - y; }
      else       { distance += y - x; }
    }
  }
  
  distance /= x_total + y_total;
  
  END_PAIR_LOOP
}




//======================================================
// Generalized UniFrac.
//======================================================
static void *generalized_mtx (void *arg) {
  START_SAMPLE_LOOP
    
  double sample_depth = 0;
  for (int otu = 0; otu < n_otus; otu++)
    sample_depth += otu_mtx[sample + otu * n_samples];
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    double abundance = otu_mtx[sample + otu * n_samples];
    if (abundance == 0) continue; // OTU not present in sample
    
    int node = otu;     // start at OTU tip/leaf in tree
    while (node > -1) { // traverse until we hit the tree's root
      int edge = nodes[node].edge;
      weight_vec[edge] += abundance / sample_depth;
      node = nodes[node].parent; // proceed on up the tree
    }
    
  }
  
  END_SAMPLE_LOOP
}

static void *generalized_dist (void *arg) {
  
  double alpha = asReal(*sexp_extra);
  
  START_PAIR_LOOP
  
  double denominator = 0;
    
  for (int edge = 0; edge < n_edges; edge++) {
    
    double x = x_weight_vec[edge];
    double y = y_weight_vec[edge];
    
    if (x || y) {
      
      double sum  = x + y;
      double frac = fabs((x - y) / sum);
      double norm = edge_lengths[edge] * pow(sum, alpha);
      
      distance    += norm * frac;
      denominator += norm;
    }
  }
  
  distance = distance / denominator;

  END_PAIR_LOOP
}




//======================================================
// Variance Adjusted Weighted UniFrac.
//======================================================
static void *var_adjusted_mtx (void *arg) {
  START_SAMPLE_LOOP
  
  double sample_depth = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    double abundance = otu_mtx[sample + otu * n_samples];
    if (abundance == 0) continue; // OTU not present in sample
    sample_depth += abundance;
    
    int node = otu;     // start at OTU tip/leaf in tree
    while (node > -1) { // traverse until we hit the tree's root
      int edge = nodes[node].edge;
      weight_vec[edge] += abundance;
      node = nodes[node].parent; // proceed on up the tree
    }
  }
  
  total_vec[sample] = sample_depth;
  
  END_SAMPLE_LOOP
}

static void *var_adjusted_dist (void *arg) {
  START_PAIR_TOTAL_LOOP
  
  double denominator = 0;
  
  for (int edge = 0; edge < n_edges; edge++) {
    
    double x = x_weight_vec[edge];
    double y = y_weight_vec[edge];
    
    if (x || y) {
      
      double norm = (x + y) * (x_total + y_total - x - y);
             norm = edge_lengths[edge] / sqrt(norm);
      
      x /= x_total;
      y /= y_total;
      
      distance    += fabs(x - y) * norm;
      denominator +=     (x + y) * norm;
    }
  }
    
  distance = distance / denominator;
  
  END_PAIR_LOOP
}




//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_unifrac(
    SEXP sexp_algorithm, 
    SEXP sexp_otu_mtx,     SEXP sexp_phylo_tree, 
    SEXP sexp_pairs_vec,   SEXP sexp_n_threads, 
    SEXP sexp_result_dist, SEXP sexp_extra_args ) {
  
  algorithm     = asInteger(sexp_algorithm);
  otu_mtx       = REAL( sexp_otu_mtx);
  n_otus        = ncols(sexp_otu_mtx);
  n_samples     = nrows(sexp_otu_mtx);
  int *edge_mtx = INTEGER(get(sexp_phylo_tree, "edge"));
  n_edges       = nrows(  get(sexp_phylo_tree, "edge"));
  edge_lengths  = REAL(   get(sexp_phylo_tree, "edge.length"));
  n_threads     = asInteger(sexp_n_threads);
  dist_vec      = REAL(sexp_result_dist);
  sexp_extra    = &sexp_extra_args;
  
  
  // Avoid allocating pairs_vec for common all-vs-all case
  if (isNull(sexp_pairs_vec)) {
    all_pairs = 1;
    pairs_vec = NULL;
    n_pairs   = n_samples * (n_samples - 1) / 2;
  }
  else {
    all_pairs = 0;
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
  }
  
  
  // branch_weight/depth for each (sample,edge) combo.
  void * (*calc_weight_mtx)(void *) = NULL;
  
  // Calculate distance between pairs, using weight_mtx.
  void * (*calc_dist_vec)(void *) = NULL;
  
  switch (algorithm) {
    case UNWEIGHTED:
      calc_weight_mtx = unweighted_mtx;
      calc_dist_vec   = unweighted_dist;
      break;
    case WEIGHTED:
      calc_weight_mtx = weighted_mtx;
      calc_dist_vec   = weighted_dist;
      break;
    case NORMALIZED:
      calc_weight_mtx = normalized_mtx;
      calc_dist_vec   = normalized_dist;
      break;
    case GENERALIZED:
      calc_weight_mtx = generalized_mtx;
      calc_dist_vec   = generalized_dist;
      break;
    case VAR_ADJUSTED:
      calc_weight_mtx = var_adjusted_mtx;
      calc_dist_vec   = var_adjusted_dist;
      break;
  }
  
  if (calc_weight_mtx == NULL || calc_dist_vec == NULL) { // # nocov start
    error("Invalid UniFrac algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // intermediary values
  weight_mtx = calloc(n_samples * n_edges, sizeof(double));
  total_vec  = calloc(n_samples,           sizeof(double));
  nodes      = calloc(n_edges,             sizeof(node_t));
  
  if (weight_mtx == NULL || total_vec == NULL || nodes == NULL) { // # nocov start
    free(weight_mtx); free(total_vec); free(nodes);
    error("Insufficient memory for UniFrac calculation.");
    return R_NilValue;
  } // # nocov end
  
  memset(weight_mtx, 0, n_samples * n_edges * sizeof(double));
  memset(total_vec,  0, n_samples           * sizeof(double));
  
  
  // sort edge data by child node
  for (int edge = 0; edge < n_edges; edge++) {
    
    int parent = edge_mtx[0 * n_edges + edge] - 2;
    int child  = edge_mtx[1 * n_edges + edge] - 1;
    
    if (child  > n_otus) child--;
    if (parent < n_otus) parent = -1;
    
    nodes[child].edge   = edge;
    nodes[child].parent = parent;
    nodes[child].length = edge_lengths[edge];
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      if (tids == NULL || args == NULL) { // # nocov start
        free(tids); free(args);
        free(weight_mtx); free(total_vec); free(nodes);
        error("Insufficient memory for parallel UniFrac calculation.");
        return R_NilValue;
      } // # nocov end
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_weight_mtx, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free(tids); free(args);
      free(weight_mtx); free(total_vec); free(nodes);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_weight_mtx(&thread_i);
  calc_dist_vec(&thread_i);
  
  free(weight_mtx); free(total_vec); free(nodes);
  
  return sexp_result_dist;
}

