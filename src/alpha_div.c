// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"

#define ADIV_ACE          1
#define ADIV_BERGER       2
#define ADIV_BRILLOUIN    3
#define ADIV_CHAO1        4
#define ADIV_FAITH        5
#define ADIV_FISHER       6
#define ADIV_INV_SIMPSON  7
#define ADIV_MARGALEF     8
#define ADIV_MCINTOSH     9
#define ADIV_MENHINICK   10
#define ADIV_OBSERVED    11
#define ADIV_SHANNON     12
#define ADIV_SIMPSON     13
#define ADIV_SQUARES     14

static int     n_samples;
static int    *pos_vec;
static int    *otu_vec;
static double *val_vec;
static SEXP   *sexp_extra;
static double *result_vec;

/*
 * FOREACH_SAMPLE iterates over all samples, ensuring that each 
 * is assigned to only a single thread. Provides `sample`, `nnz`, 
 * and `result` initialized to 0; expects `result` at the end.
 */
#define FOREACH_SAMPLE(expression)                             \
  do {                                                         \
    int sample    = ((worker_t *)arg)->i;                      \
    int n_threads = ((worker_t *)arg)->n;                      \
                                                               \
    for (; sample < n_samples; sample += n_threads) {          \
      double *val_begin = val_vec + pos_vec[sample];           \
      double *val_end   = val_vec + pos_vec[sample + 1];       \
      int     nnz       = val_end - val_begin;                 \
      double result     = 0;                                   \
      if (nnz) {                                               \
        expression;                                            \
      }                                                        \
      else {                                                   \
        result = NA_REAL;                                      \
      }                                                        \
      result_vec[sample] = result;                             \
    }                                                          \
  } while (0)

/*
 * FOREACH_VAL iterates over each value for the current sample.
 * Sets val to point to each value in val_vec in turn.
 */
#define FOREACH_VAL(expression)                                \
  do {                                                         \
    double *val = val_begin;                                   \
    for (; val != val_end; val++) {                            \
      expression;                                              \
    }                                                          \
  } while (0)



//======================================================
// Abundance-based Coverage Estimator (Chao & Lee 1992)
//======================================================

static int     ace_cutoff;
static double *ace_rare_nnz_k_mtx;

static void *ace(void *arg) {
  
  int     thread_i       = ((worker_t *)arg)->i;
  int     cutoff         = ace_cutoff;
  double *rare_nnz_k_vec = ace_rare_nnz_k_mtx + (thread_i * cutoff);
  
  FOREACH_SAMPLE(
    
    double abund_nnz = 0;
    double rare_sum  = 0;
    double rare_nnz  = 0;
    
    memset(rare_nnz_k_vec, 0, cutoff * sizeof(double));
    
    FOREACH_VAL(
      int x_int = (int)(ceil(*val));
      if (x_int < cutoff) {
        rare_sum += x_int;
        rare_nnz++;
        rare_nnz_k_vec[x_int]++;
      }
      else {
        abund_nnz++;
      }
    ); // FOREACH_VAL
    
    for (int k = 1; k < cutoff; k++)
      result += k * (k - 1) * rare_nnz_k_vec[k];
    
    double p = (1 - (rare_nnz_k_vec[1] / rare_sum));
    
    result = result * rare_nnz/(p * rare_sum * (rare_sum - 1)) - 1;
    if (result < 0) result = 0;
    result = abund_nnz + rare_nnz/p + result * rare_nnz_k_vec[1]/p;
    
  ); // FOREACH_SAMPLE
  
  return NULL;
}

static pthread_func_t ace_setup(int n_threads) {
  
  ace_cutoff         = asInteger(*sexp_extra) + 1;
  otu_vec            = maybe_free_one(otu_vec);
  ace_rare_nnz_k_mtx = safe_malloc(n_threads * ace_cutoff * sizeof(double));
  
  return ace;
}
  
  
//======================================================
// Berger-Parker (Berger & Parker 1970)
// max(x) / sum(x)
//======================================================
static void *berger(void *arg) {
  FOREACH_SAMPLE(
    FOREACH_VAL(
      if (*val > result) result = *val;
    );
  );
  
  return NULL;
}
  
  
//======================================================
// Brillouin (Brillouin 1956)
// (log(sum(x)!) - sum(log(x!))) / sum(x)
// note: lgamma(x + 1) == log(x!)
//======================================================
static void *brillouin(void *arg) {
  FOREACH_SAMPLE(
    
    double depth = 0;
  
    FOREACH_VAL(
      depth  += *val;
      result += lgamma(*val + 1);
    );
    
    result = (lgamma(depth + 1) - result) / depth;
  );
  
  return NULL;
}
  


//======================================================
// Chao1 (Chao 1984)
// sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))
//======================================================
static void *chao1(void *arg) {
  FOREACH_SAMPLE(
    
    double ones = 0;
    double twos = 0;
    
    FOREACH_VAL(
      if      (*val == 1) ones++;
      else if (*val == 2) twos++;
    );
    
    result = nnz + ((ones * ones) / (2 * twos));
  );
  
  return NULL;
}



//======================================================
// Faith's Phylogenetic Diversity
// Find the edges that an otu passes through.
//======================================================

static ecotree_t *faith_et;
static char      *faith_has_edge_mtx;

static void *faith(void *arg) {
  
  int     thread_i     = ((worker_t *)arg)->i;
  int     n_edges      = faith_et->n_edges;
  node_t *node_vec     = faith_et->node_vec;
  char   *has_edge_vec = faith_has_edge_mtx + (thread_i * n_edges);
  
  FOREACH_SAMPLE(
    
    memset(has_edge_vec, 0, n_edges * sizeof(char));
    
    int *otu_begin = otu_vec + pos_vec[sample];
    int *otu_end   = otu_vec + pos_vec[sample + 1];
    
    for (int *otu = otu_begin; otu != otu_end; otu++) {
      
      int node_i = *otu;    // start at OTU tip/leaf in tree
      while (node_i > -1) { // traverse until we hit the tree's root
        
        node_t *node   = node_vec + node_i;
        char *has_edge = has_edge_vec + node->edge;
        
        if (*has_edge) break; // already traversed
        *has_edge = 1;
        
        result += node->length;
        node_i  = node->parent;
      }
    }
    
  );
  
  return NULL;
}

static pthread_func_t faith_setup(int n_threads) {
  
  faith_et           = new_ecotree(*sexp_extra);
  faith_has_edge_mtx = safe_malloc(n_threads * faith_et->n_edges * sizeof(char));
  
  return faith;
}


//======================================================
// Fisher's diversity index (Fisher 1943)
// otus = fisher * log(1 + depth/fisher)
//======================================================
static void *fisher(void *arg) {
  
  int    digits = asInteger(*sexp_extra);
  double mult   = pow(10, digits);
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    FOREACH_VAL(depth += *val);
    
    if (depth == nnz) {
      result = R_PosInf;  // All singletons -> infinite loop
    }
    else {
      
      double alpha;
      double lo = 2;
      double hi = 16;
      
      // Sometimes the result will be less than 2 or greater than 16.
      while (lo * log(1 + depth/lo) > nnz) { hi = lo; lo /= 2; }
      while (hi * log(1 + depth/hi) < nnz) { lo = hi; hi *= 2; }
      
      // Check if range has converged to same value after rounding.
      while (round(hi * mult) != round(lo * mult)) {
        
        // This loop's guess for the alpha term.
        alpha = (lo + hi) / 2;
        
        // Update the range we need to examine.
        if (alpha * log(1 + depth/alpha) > nnz) { hi = alpha; }
        else                                    { lo = alpha; }
        
      }
      
      result = round(hi * mult) / mult;
    }
  );
  
  return NULL;
}


//======================================================
// Inverse Simpson
// p <- x / sum(x)
// 1 / sum(p ** 2)
//======================================================
static void *inv_simpson(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * *val);
  
    result = 1 / result;
  );
  
  return NULL;
}


//======================================================
// Margalef (Margalef 1958)
// (sum(x > 0) - 1) / log(sum(x))
//======================================================
static void *margalef(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    
    FOREACH_VAL(depth += *val);
    
    result = (nnz - 1) / log(depth);
  );
  
  return NULL;
}


//======================================================
// McIntosh (McIntosh 1967)
// (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
//======================================================
static void *mcintosh(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    
    FOREACH_VAL(
      depth  += *val;
      result += *val * *val;
    );
    
    result = (depth - sqrt(result)) / (depth - sqrt(depth));
  );
  
  return NULL;
}


//======================================================
// Menhinick (Menhinick 1964)
// sum(x > 0) / sqrt(sum(x))
//======================================================
static void *menhinick(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
  
    FOREACH_VAL(depth += *val);
    
    result = nnz / sqrt(depth);
  );
  
  return NULL;
}


//======================================================
// Observed Features
// sum(x > 0)
//======================================================
static void *observed(void *arg) {
  
  FOREACH_SAMPLE(result = nnz);
  
  return NULL;
}


//======================================================
// Shannon (Shannon 1948)
// p <- x / sum(x)
// -sum(p * log(p))
//======================================================
static void *shannon(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * log(*val));
  
    result *= -1;
  );
  
  return NULL;
}


//======================================================
// Simpson (Simpson 1949)
// p <- x / sum(x)
// 1 - sum(p ** 2)
//======================================================
static void *simpson(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * *val);
  
    result = 1 - result;
  );
  
  return NULL;
}


//======================================================
// Squares Estimator (Alroy 2018)
// N = sum(x)      # sampling depth
// S = sum(x > 0)  # number of non-zero OTUs
// F1 = sum(x == 1) # singletons
// ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S)) + S
//======================================================
static void *squares(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth      = 0;
    double singletons = 0;
    
    FOREACH_VAL(
      depth  += *val;
      result += *val * *val;
      if (*val == 1) singletons++;
    );
    
    double denominator = (depth * depth) - (singletons * nnz);
    
    if (denominator == 0) {
      result = R_PosInf; // All singletons
    } else {
      result *= (singletons * singletons);
      result /= denominator;
      result += nnz;
    }
  );
  
  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx, 
    SEXP sexp_margin,    SEXP sexp_norm, 
    SEXP sexp_n_threads, SEXP sexp_extra_args ) {
  
  init_n_ptrs(10);
  
  int norm       = asInteger(sexp_norm);
  int n_threads  = asInteger(sexp_n_threads);
  sexp_extra     = &sexp_extra_args;
  
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  if (norm) normalize(em, norm, n_threads, 0);
  
  n_samples = em->n_samples;
  pos_vec   = em->pos_vec;
  otu_vec   = em->otu_vec;
  val_vec   = em->val_vec;
  
  
  // function to run
  // void * (*adiv_func)(void *) = NULL;
  pthread_func_t adiv_func = NULL;
  switch (asInteger(sexp_algorithm)) {
    case ADIV_ACE:         adiv_func = ace_setup(n_threads);   break;
    case ADIV_BERGER:      adiv_func = berger;                 break;
    case ADIV_BRILLOUIN:   adiv_func = brillouin;              break;
    case ADIV_CHAO1:       adiv_func = chao1;                  break;
    case ADIV_FAITH:       adiv_func = faith_setup(n_threads); break;
    case ADIV_FISHER:      adiv_func = fisher;                 break;
    case ADIV_INV_SIMPSON: adiv_func = inv_simpson;            break;
    case ADIV_MARGALEF:    adiv_func = margalef;               break;
    case ADIV_MCINTOSH:    adiv_func = mcintosh;               break;
    case ADIV_MENHINICK:   adiv_func = menhinick;              break;
    case ADIV_OBSERVED:    adiv_func = observed;               break;
    case ADIV_SHANNON:     adiv_func = shannon;                break;
    case ADIV_SIMPSON:     adiv_func = simpson;                break;
    case ADIV_SQUARES:     adiv_func = squares;                break;
  }
  
  if (adiv_func == NULL) { // # nocov start
    error("Invalid alpha diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Create the diversity vector to return
  SEXP sexp_result_vec = PROTECT(allocVector(REALSXP, n_samples));
  result_vec           = REAL(sexp_result_vec);
  setAttrib(sexp_result_vec, R_NamesSymbol, em->sexp_sample_names);
  
  
  run_parallel(adiv_func, n_threads, n_samples);
  
  free_all();
  UNPROTECT(1);
  return sexp_result_vec;
}
