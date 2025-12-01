// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <math.h>     // fabs, pow
#include <string.h>   // memset
#include "ecomatrix.h"
#include "ecotree.h"
#include "memory.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define U_UNIFRAC 1
#define W_UNIFRAC 2
#define N_UNIFRAC 3
#define G_UNIFRAC 4
#define V_UNIFRAC 5


//======================================================
// Variables shared between main and worker threads.
//======================================================

static int     n_samples;
static int     n_otus;
static int     n_edges;
static int     n_pairs;
static int     n_dist;
static int     n_threads;
static int    *pos_vec;
static int    *otu_vec;
static double *val_vec;
static node_t *node_vec;
static double *edge_lengths;
static int    *pairs_vec;
static SEXP   *sexp_extra;
static double *weight_mtx;
static double *sample_norm_vec;
static double *dist_vec;


  
  
/*
 * Macro to loop over each sample based on current threading 
 * setup. Sets sam (sample index), weight_vec, offset and nnz.
 * `expression` may assign to *sample_norm.
 */
#define FOREACH_SAMPLE(expression)                             \
  do {                                                         \
    int sam = *((int *) arg);                                  \
    for (; sam < n_samples; sam += n_threads) {                \
      double *sample_norm = sample_norm_vec + sam;             \
      double *weight_vec  = weight_mtx + (sam * n_edges);      \
      int     offset      = pos_vec[sam];                      \
      int     nnz         = pos_vec[sam + 1] - offset;         \
                                                               \
      expression;                                              \
                                                               \
      (void)sample_norm;                                       \
    }                                                          \
  } while (0)



/*
 * Macro to loop over each value for the current sample.
 * Sets otu and val to point to each value in otu_vec and 
 * val_vec in turn, respectively.
 */
#define FOREACH_OTU_VAL(expression)                            \
  do {                                                         \
    int    *otu = otu_vec + offset;                            \
    double *val = val_vec + offset;                            \
    for (int i = 0; i < nnz; i++) {                            \
                                                               \
      expression;                                              \
                                                               \
      otu++;                                                   \
      val++;                                                   \
    }                                                          \
    (void)otu;                                                 \
    (void)val;                                                 \
  } while (0)



/*
 * Macro to traverse the phylogenetic tree from a leaf to the
 * root. It starts from the leaf node index given by *otu.
 * 
 * * In each iteration, it provides:
 * - `*node`:   A pointer to the current node_t struct.
 * - `*weight`: A pointer to the sample's weight for the 
 *              current node's edge.
 * 
 * * The expression is expected to read `*node` and assign a 
 *   value to `*weight`.
 */
#define FOREACH_NODE_WEIGHT(expression)                        \
  do {                                                         \
    int node_i = *otu;                                         \
    while (node_i > -1) {                                      \
      node_t *node   = node_vec + node_i;                      \
      double *weight = weight_vec + node->edge;                \
                                                               \
      expression;                                              \
                                                               \
      node_i = node->parent;                                   \
    }                                                          \
  } while (0)
  


/*
 * FOREACH_SAMPLE_PAIR runs `expression` on all unique sample 
 * pairs, dividing the workload evenly across multiple CPU 
 * threads. When `pairs_vec == NULL`, a simpler algorithm is used
 * to iterate all-vs-all pairs.
 * 
 * In all cases, FOREACH_SAMPLE_PAIR provides:
 *   - `*x_weight_vec` and `*y_weight_vec`   (from `weight_mtx`)
 *   - `*x_sample_norm` and `*y_sample_norm` (from `sample_norm_vec`)
 *   
 * And FOREACH_SAMPLE_PAIR expects `expression` to set:
 *   - `*distance`
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */

#define FOREACH_SAMPLE_PAIR(expression)                        \
  do {                                                         \
    int     thread_i = *((int *) arg);                         \
    int     dist_idx = 0;                                      \
    double *x_weight_vec, *x_sample_norm;                      \
    double *y_weight_vec, *y_sample_norm;                      \
                                                               \
    if (pairs_vec == NULL) { /* All vs All */                  \
                                                               \
      for (int i = 0; i < n_samples - 1; i++) {                \
        x_weight_vec  = weight_mtx + (i * n_edges);            \
        x_sample_norm = sample_norm_vec + i;                   \
                                                               \
        for (int j = i + 1; j < n_samples; j++) {              \
                                                               \
          if (dist_idx % n_threads == thread_i) {              \
            y_weight_vec  = weight_mtx + (j * n_edges);        \
            y_sample_norm = sample_norm_vec + j;               \
                                                               \
            double *distance = dist_vec + dist_idx;            \
            *distance = 0;                                     \
                                                               \
            expression;                                        \
          }                                                    \
          dist_idx++;                                          \
        }                                                      \
      }                                                        \
                                                               \
    } else { /* Specific Pairs of Samples */                   \
                                                               \
      int pair_idx = thread_i;                                 \
      for (; pair_idx < n_pairs; pair_idx += n_threads) {      \
                                                               \
        dist_idx = pairs_vec[pair_idx]; /* 1-based */          \
                                                               \
        int sam_i          = 0;                                \
        int sam_j          = dist_idx;                         \
        int pairs_in_block = n_samples - 1;                    \
                                                               \
        while (sam_j > pairs_in_block) {                       \
          sam_i++;                                             \
          sam_j -= pairs_in_block;                             \
          pairs_in_block--;                                    \
        }                                                      \
                                                               \
        sam_j += sam_i;                                        \
                                                               \
        x_sample_norm = sample_norm_vec + sam_i;               \
        y_sample_norm = sample_norm_vec + sam_j;               \
        x_weight_vec  = weight_mtx + (sam_i * n_edges);        \
        y_weight_vec  = weight_mtx + (sam_j * n_edges);        \
                                                               \
        double *distance = dist_vec + (dist_idx - 1);          \
        *distance = 0;                                         \
                                                               \
        expression;                                            \
      }                                                        \
    }                                                          \
                                                               \
    (void)x_sample_norm;                                       \
    (void)y_sample_norm;                                       \
                                                               \
  } while (0)



/*
 * FOREACH_WEIGHT_PAIR iterates through all edges, providing
 * the corresponding weights for the two samples being compared.
 * 
 * * In each iteration, it provides:
 * - `*x_weight` and `*y_weight` (from `weight_mtx`)
 * - `*edge_length` (from `edge_lengths`)
 * 
 * * `expression` is expected to accumulate its result into 
 *   `*distance` to satisfy the parent FOREACH_SAMPLE_PAIR macro.
 */

#define FOREACH_WEIGHT_PAIR(expression)                        \
  do {                                                         \
    double *x_weight    = x_weight_vec;                        \
    double *y_weight    = y_weight_vec;                        \
    double *edge_length = edge_lengths;                        \
                                                               \
    for (int edge = 0; edge < n_edges; edge++) {               \
                                                               \
      if (*x_weight || *y_weight) {                            \
        expression;                                            \
      }                                                        \
                                                               \
      x_weight++;                                              \
      y_weight++;                                              \
      edge_length++;                                           \
    }                                                          \
    (void)edge_length;                                         \
  } while (0)


  


//======================================================
// Unweighted UniFrac.
//======================================================
static void *unweighted_mtx (void *arg) {
  
  FOREACH_SAMPLE(
    FOREACH_OTU_VAL(
      FOREACH_NODE_WEIGHT(

        if (*weight) break; // already traversed
        *weight = 1;

      );
    );
  );
  return NULL;
}


static void *unweighted_dist (void *arg) {
  
  FOREACH_SAMPLE_PAIR(

    double shared = 0;

    FOREACH_WEIGHT_PAIR(

      if (*x_weight && *y_weight) {  shared   += *edge_length; }
      else                        { *distance += *edge_length; }

    );

    *distance /= *distance + shared;
  );
  return NULL;
}




//======================================================
// Weighted UniFrac.
//======================================================
static void *weighted_mtx (void *arg) {
  
  FOREACH_SAMPLE(
    
    double sample_depth = 0;
    
    FOREACH_OTU_VAL(sample_depth += *val);
    
    FOREACH_OTU_VAL(
      FOREACH_NODE_WEIGHT(
        
        // relative abundance, weighted by branch length
        *weight += node->length * (*val / sample_depth);
    
      );
    );
  );
  
  return NULL;
}

static void *weighted_dist (void *arg) {
  
  FOREACH_SAMPLE_PAIR(
    FOREACH_WEIGHT_PAIR(
      
      *distance += fabs(*x_weight - *y_weight);
      
    );
  );
  
  return NULL;
}




//======================================================
// Normalized Weighted UniFrac.
//======================================================
static void *normalized_mtx (void *arg) {
  
  FOREACH_SAMPLE(
    
    double sample_depth = 0;
    
    FOREACH_OTU_VAL(sample_depth += *val);
    
    FOREACH_OTU_VAL(
      FOREACH_NODE_WEIGHT(
        
        // relative abundance, weighted by branch length
        double abund = node->length * (*val / sample_depth);
        *weight      += abund;
        *sample_norm += abund;
        
      );
    );
  );
  
  return NULL;
}

static void *normalized_dist (void *arg) {
  
  
  FOREACH_SAMPLE_PAIR(
    FOREACH_WEIGHT_PAIR(
      
      *distance += fabs(*x_weight - *y_weight);
  
    );
    
    *distance /= *x_sample_norm + *y_sample_norm;
    
  );
  
  return NULL;
}




//======================================================
// Generalized UniFrac.
//======================================================
static void *generalized_mtx (void *arg) {
  
  FOREACH_SAMPLE(
    
    double sample_depth = 0;
  
    FOREACH_OTU_VAL(sample_depth += *val);
    
    FOREACH_OTU_VAL(
      FOREACH_NODE_WEIGHT(
        
        *weight += *val / sample_depth; // relative abundance
    
      );
    );
  );
  
  return NULL;
}

static void *generalized_dist (void *arg) {
  
  double alpha = asReal(*sexp_extra);
  
  FOREACH_SAMPLE_PAIR(
    
    double denominator = 0;
  
    FOREACH_WEIGHT_PAIR(
      
      double sum  = *x_weight + *y_weight;
      double frac = fabs((*x_weight - *y_weight) / sum);
      double norm = *edge_length * pow(sum, alpha);
      
      *distance   += norm * frac;
      denominator += norm;
  
    );
  
    *distance /= denominator;
  
  );
  
  return NULL;
}




//======================================================
// Variance Adjusted Weighted UniFrac.
//======================================================
static void *var_adjusted_mtx (void *arg) {
  
  FOREACH_SAMPLE(
    FOREACH_OTU_VAL(
      
      *sample_norm += *val;
      
      FOREACH_NODE_WEIGHT(
        *weight += *val; // absolute abundance
      );
    );
  );
  
  return NULL;
}

static void *var_adjusted_dist (void *arg) {
  
  FOREACH_SAMPLE_PAIR(
    
    double norm;
    double denominator = 0;
  
    FOREACH_WEIGHT_PAIR(
      
      norm  = *x_sample_norm + *y_sample_norm;
      norm -= *x_weight + *y_weight;
      norm *= *x_weight + *y_weight;
      
      // Check for div-by-zero
      // Contribution is 0 if variance is 0
      if (norm > 0) { norm = *edge_length / sqrt(norm); }
      else          { norm = 0;                         }
      
      double x = *x_weight / *x_sample_norm;
      double y = *y_weight / *y_sample_norm;
      
      *distance   += fabs(x - y) * norm;
      denominator +=     (x + y) * norm;
      
    );
    
    *distance /= denominator;
  
  );
  
  return NULL;
}




//======================================================
// R interface. Dispatches threads on unifrac variants.
//======================================================
SEXP C_unifrac(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx,   SEXP sexp_phylo_tree, 
    SEXP sexp_margin,    SEXP sexp_pairs_vec, SEXP sexp_n_threads,  
    SEXP sexp_extra_args ) {
  
  sexp_extra = &sexp_extra_args;
  n_threads  = asInteger(sexp_n_threads);
  init_n_ptrs(n_threads + 10);
  
  
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  ecotree_t   *et = new_ecotree(sexp_phylo_tree);
  
  n_samples    = em->n_samples;
  n_otus       = em->n_otus;
  pos_vec      = em->pos_vec;
  otu_vec      = em->otu_vec;
  val_vec      = em->val_vec;
  n_edges      = et->n_edges;
  edge_lengths = et->edge_lengths;
  node_vec     = et->node_vec;
  
  
  // branch_weight/depth for each (sample,edge) combo.
  void * (*calc_weight_mtx)(void *) = NULL;
  
  // Calculate distance between pairs, using weight_mtx.
  void * (*calc_dist_vec)(void *) = NULL;
  
  switch (asInteger(sexp_algorithm)) {
    case U_UNIFRAC:
      calc_weight_mtx = unweighted_mtx;
      calc_dist_vec   = unweighted_dist;
      break;
    case W_UNIFRAC:
      calc_weight_mtx = weighted_mtx;
      calc_dist_vec   = weighted_dist;
      break;
    case N_UNIFRAC:
      calc_weight_mtx = normalized_mtx;
      calc_dist_vec   = normalized_dist;
      break;
    case G_UNIFRAC:
      calc_weight_mtx = generalized_mtx;
      calc_dist_vec   = generalized_dist;
      break;
    case V_UNIFRAC:
      calc_weight_mtx = var_adjusted_mtx;
      calc_dist_vec   = var_adjusted_dist;
      break;
    default:                                // # nocov
      free_all();                           // # nocov
      error("Invalid UniFrac algorithm.");  // # nocov
  }
  
  
  // intermediary values
  weight_mtx      = (double *)safe_malloc(n_samples * n_edges * sizeof(double));
  sample_norm_vec = (double *)safe_malloc(n_samples           * sizeof(double));
  
  memset(weight_mtx,      0, n_samples * n_edges * sizeof(double));
  memset(sample_norm_vec, 0, n_samples           * sizeof(double));
  
  
  // Create the dist object to return
  n_dist                = n_samples * (n_samples - 1) / 2;
  SEXP sexp_result_dist = PROTECT(allocVector(REALSXP, n_dist));
  dist_vec              = REAL(sexp_result_dist);
  setAttrib(sexp_result_dist, R_ClassSymbol,      mkString("dist"));
  setAttrib(sexp_result_dist, mkString("Size"),   ScalarInteger(n_samples));
  setAttrib(sexp_result_dist, mkString("Diag"),   ScalarLogical(0));
  setAttrib(sexp_result_dist, mkString("Upper"),  ScalarLogical(0));
  setAttrib(sexp_result_dist, mkString("Labels"), em->sexp_sample_names);
  
  
  // Avoid allocating pairs_vec for common all-vs-all case
  if (isNull(sexp_pairs_vec)) {
    
    pairs_vec = NULL;
    n_pairs   = n_dist;
    
  } else {
    
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
    
    for (int i = 0; i < n_dist; i++)
      dist_vec[i] = R_NaReal;
    
    if (n_pairs == 0) {
      free_all();
      UNPROTECT(1);
      return sexp_result_dist;
    }
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = (pthread_t*) R_alloc(n_threads, sizeof(pthread_t));
      int       *args = (int*)       R_alloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_weight_mtx, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free_all();
      UNPROTECT(1);
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_weight_mtx(&thread_i);
  calc_dist_vec(&thread_i);
  
  free_all();
  UNPROTECT(1);
  return sexp_result_dist;
}

