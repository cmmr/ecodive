#include <R.h>
#include <Rinternals.h>
#include <math.h>   // fabs, pow
#include <stdlib.h> // calloc, free
#include <string.h> // memset
#include <assert.h> // assert
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
  double *total_1;
  double *total_2;
  double *weight_vec_1;
  double *weight_vec_2;
  double *distance;
} pair_t;


typedef struct {
  int    edge;
  int    parent;
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
static pair_t *pair_vec;
static int     n_pairs;
static double  alpha;
static double *weight_mtx;
static node_t *nodes;
static double *total_vec;
static int     n_threads;




//======================================================
// Unweighted UniFrac.
//======================================================
static void *unweighted_mtx (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *otu_vec    = otu_mtx    + (sample * n_otus);
    double *weight_vec = weight_mtx + (sample * n_edges);
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      if (otu_vec[otu] == 0) continue; // OTU not present in sample
      
      int node = otu;     // start at OTU tip/leaf in tree
      while (node > -1) { // traverse until we hit the tree's root
        int edge = nodes[node].edge;
        if (weight_vec[edge]) break; // already traversed
        weight_vec[edge] = 1;
        node = nodes[node].parent;   // proceed on up the tree
      }
      
    }
  }
  
  return NULL;
}


static void *unweighted_dist (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in weight_mtx
    double *x_vec = pair->weight_vec_1;
    double *y_vec = pair->weight_vec_2;
    
    double distinct = 0, shared = 0;
    
    for (int edge = 0; edge < n_edges; edge++) {
      
      char x = x_vec[edge] > 0;
      char y = y_vec[edge] > 0;
      
      if      (x & y) { shared   += edge_lengths[edge]; }
      else if (x | y) { distinct += edge_lengths[edge]; }
    }
    
    // value to return
    *(pair->distance) = distinct / (distinct + shared);
  }
  
  return NULL;
}




//======================================================
// Weighted UniFrac.
//======================================================
static void *weighted_mtx (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *otu_vec    = otu_mtx    + (sample * n_otus);
    double *weight_vec = weight_mtx + (sample * n_edges);
    
    double sample_depth = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      sample_depth += otu_vec[otu];
    }
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double abundance = otu_vec[otu];
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
  }
  
  return NULL;
}

static void *weighted_dist (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in weight_mtx
    double *x_vec = pair->weight_vec_1;
    double *y_vec = pair->weight_vec_2;
    
    double distance = 0;
    
    for (int edge = 0; edge < n_edges; edge++) {
      double weight_1 = x_vec[edge];
      double weight_2 = y_vec[edge];
      distance += fabs(weight_1 - weight_2);
    }
    
    // value to return
    *(pair->distance) = distance;
  }
  
  return NULL;
}




//======================================================
// Normalized Weighted UniFrac.
//======================================================
static void *normalized_mtx (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *otu_vec    = otu_mtx    + (sample * n_otus);
    double *weight_vec = weight_mtx + (sample * n_edges);
    
    double sample_depth = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      sample_depth += otu_vec[otu];
    }
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double abundance = otu_vec[otu];
      if (abundance == 0) continue; // OTU not present in sample
      
      int node = otu;     // start at OTU tip/leaf in tree
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
  }
  
  return NULL;
}

static void *normalized_dist (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in weight_mtx
    double *x_vec = pair->weight_vec_1;
    double *y_vec = pair->weight_vec_2;
    
    double distance = 0;
    
    for (int edge = 0; edge < n_edges; edge++) {
      distance += fabs(x_vec[edge] - y_vec[edge]);
    }
    distance /= *(pair->total_1) + *(pair->total_2);
    
    // value to return
    *(pair->distance) = distance;
  }
  
  return NULL;
}




//======================================================
// Generalized UniFrac.
//======================================================
static void *generalized_mtx (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *otu_vec    = otu_mtx    + (sample * n_otus);
    double *weight_vec = weight_mtx + (sample * n_edges);
    
    double sample_depth = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      sample_depth += otu_vec[otu];
    }
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double abundance = otu_vec[otu];
      if (abundance == 0) continue; // OTU not present in sample
      
      int node = otu;     // start at OTU tip/leaf in tree
      while (node > -1) { // traverse until we hit the tree's root
        int edge = nodes[node].edge;
        weight_vec[edge] += abundance / sample_depth;
        node = nodes[node].parent; // proceed on up the tree
      }
      
    }
  }
  
  return NULL;
}

static void *generalized_dist (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in weight_mtx
    double *x_vec = pair->weight_vec_1;
    double *y_vec = pair->weight_vec_2;
    
    double distance    = 0;
    double denominator = 0;
    
    for (int edge = 0; edge < n_edges; edge++) {
      
      double x = x_vec[edge];
      double y = y_vec[edge];
      
      double sum = x + y;
      
      if (sum > 0) {
        double frac = fabs((x - y) / sum);
        double norm = edge_lengths[edge] * pow(sum, alpha);
        
        distance    += norm * frac;
        denominator += norm;
      }
    }
    
    // value to return
    *(pair->distance) = distance / denominator;
  }
  
  return NULL;
}




//======================================================
// Variance Adjusted Weighted UniFrac.
//======================================================
static void *var_adjusted_mtx (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double  sample_depth = 0;
    double *otu_vec      = otu_mtx    + (sample * n_otus);
    double *weight_vec   = weight_mtx + (sample * n_edges);
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double abundance = otu_vec[otu];
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
  }
  
  return NULL;
}

static void *var_adjusted_dist (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in weight_mtx
    double *x_vec = pair->weight_vec_1;
    double *y_vec = pair->weight_vec_2;
    
    double x_total = *(pair->total_1);
    double y_total = *(pair->total_2);
    
    double distance    = 0;
    double denominator = 0;
    
    for (int edge = 0; edge < n_edges; edge++) {
      
      double x = x_vec[edge];
      double y = y_vec[edge];
      
      double norm = (x + y) * (x_total + y_total - x - y);
      
      if (norm > 0) {
        norm = edge_lengths[edge] / sqrt(norm);
        x   /= x_total;
        y   /= y_total;
        
        distance    += fabs(x - y) * norm;
        denominator +=     (x + y) * norm;
      }
    }
    
    // value to return
    *(pair->distance) = distance / denominator;
  }
  
  return NULL;
}




//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_unifrac(
    SEXP sexp_algorithm, 
    SEXP sexp_otu_mtx,   SEXP sexp_phylo_tree, 
    SEXP sexp_alpha,     SEXP sexp_pair_idx_vec, 
    SEXP sexp_n_threads, SEXP sexp_result_dist ) {
  
  algorithm         = asInteger(sexp_algorithm);
  otu_mtx           = REAL( sexp_otu_mtx);
  n_otus            = nrows(sexp_otu_mtx);
  n_samples         = ncols(sexp_otu_mtx);
  int *edge_mtx     = INTEGER(get(sexp_phylo_tree, "edge"));
  n_edges           = nrows(  get(sexp_phylo_tree, "edge"));
  edge_lengths      = REAL(   get(sexp_phylo_tree, "edge.length"));
  alpha             = asReal(sexp_alpha);
  int *pair_idx_vec = INTEGER(sexp_pair_idx_vec);
  n_pairs           = LENGTH(sexp_pair_idx_vec);
  n_threads         = asInteger(sexp_n_threads);
  double *dist_vec  = REAL(sexp_result_dist);
  int n_dist        = LENGTH(sexp_result_dist);
  
  
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
    default:                         // # nocov start
      error("Invalid adiv metric.");
      return R_NilValue;             // # nocov end
  }
  
  
  // intermediary values
  weight_mtx = calloc(n_samples * n_edges, sizeof(double));
  nodes      = calloc(n_edges,             sizeof(node_t));
  pair_vec   = calloc(n_pairs,             sizeof(pair_t));
  total_vec  = calloc(n_samples,           sizeof(double));
  
  assert(weight_mtx != NULL);
  assert(nodes      != NULL);
  assert(pair_vec   != NULL);
  assert(total_vec  != NULL);
  
  if (weight_mtx == NULL || nodes == NULL || pair_vec == NULL || total_vec == NULL) { // # nocov start
    free(weight_mtx); free(nodes); free(pair_vec); free(total_vec);
    error("Unable to allocate memory for UniFrac calculation.");
    return R_NilValue;
  } // # nocov end
  
  memset(weight_mtx, 0, n_samples * n_edges * sizeof(double));
  
  
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
  
  // pointers to input and output for 
  // each pairwise comparison
  int pair_idx = 0;
  int dist_idx = 0;
  for (int i = 0; i < n_samples - 1; i++) {
    for (int j = i + 1; j < n_samples; j++) {
      if (pair_idx_vec[pair_idx] == dist_idx) {
        
        assert(pair_idx < n_pairs);
        assert(dist_idx < n_dist);
        
        pair_t *pair       = pair_vec   + pair_idx;
        pair->total_1      = total_vec  + i;
        pair->total_2      = total_vec  + j;
        pair->weight_vec_1 = weight_mtx + (i * n_edges);
        pair->weight_vec_2 = weight_mtx + (j * n_edges);
        pair->distance     = dist_vec   + dist_idx;
        pair_idx++;
        
        if (pair_idx == n_pairs) goto end_loops;
      }
      dist_idx++;
    }
  }
  end_loops:
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      assert(tids != NULL);
      assert(args != NULL);
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_weight_mtx, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      assert(tids       != NULL);
      assert(args       != NULL);
      assert(weight_mtx != NULL);
      assert(nodes      != NULL);
      assert(pair_vec   != NULL);
      assert(total_vec  != NULL);
      
      free(tids); free(args);
      free(weight_mtx); free(nodes); free(pair_vec); free(total_vec);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_weight_mtx(&thread_i);
  calc_dist_vec(&thread_i);
  
  assert(weight_mtx != NULL);
  assert(nodes      != NULL);
  assert(pair_vec   != NULL);
  assert(total_vec  != NULL);
  
  free(weight_mtx); free(nodes); free(pair_vec); free(total_vec);
  
  return sexp_result_dist;
}

