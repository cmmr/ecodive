// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <inttypes.h> // uint32_t, uint64_t
#include <math.h>     // fabs, floor
#include <string.h>   // memcpy, memset
#include <stdlib.h>   // qsort
#include "get.h"
#include "memory.h"
#include "pcg_basic.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

typedef void *(*pthread_func_t)(void *);

static uint32_t  target;
static uint64_t  seed;
static int       margin;
static int       n_threads;
static SEXP      sexp_val_mtx;
static SEXP      sexp_res_mtx;
static int       n_protected;
static int      *sam_vec;
static int      *pos_vec;
static double   *val_vec;
static double   *res_vec;
static int       n_sams;
static int       n_otus;
static int       n_vals;
static uint32_t *depth_vec;


typedef struct {
  uint32_t tried;
  uint32_t kept;
  pcg32_random_t rng;
} knuth_t;



/*
 * Functions for base matrix and dgeMatrix
 * 
 * Traversing a matrix by row (margin == 1) or by
 * column (margin == 2) is handled by the
 * rarefy_dense function by changing the otu_step
 * and sam_step variables.
 * 
 */

static void *rarefy_dense(void *arg) {
  
  int thread_i = *((int *) arg);
  
  int     otu_step    = (margin == 1) ? n_sams : 1;
  int     sam_step    = (margin == 1) ? 1 : n_otus;
  int     end_step    = otu_step * n_otus;
  int     thread_step = n_threads * sam_step;
  double *val_begin   = val_vec + thread_i * sam_step;
  double *res_begin   = res_vec + thread_i * sam_step;
  
  for (int sam = thread_i; sam < n_sams; sam += n_threads) {
    
    uint32_t  depth   = depth_vec[sam];
    
    // Sample can be be rarefied.
    if (depth > target) {
      
      // Seed the PRNG for this sample.
      pcg32_random_t rng;
      pcg32_srandom_r(&rng, seed, sam);
      
      double *val = val_begin; // Current # of observations
      double *res = res_begin; // Rarefied # of observations
      
      // Knuth algorithm for choosing target seqs from depth.
      uint32_t tried = 0, kept = 0;
      for (int otu = 0; otu < n_otus; otu++) {
        *res = 0;                
        
        for (int seq = 0; seq < *val && kept < target; seq++) {
          
          uint32_t not_tried  = depth - tried;
          uint32_t still_need = target - kept;
          uint32_t rand_int   = pcg32_random_r(&rng);
          
          if (rand_int % not_tried < still_need) {
            (*res)++; // retain this observation
            kept++;
          }
          
          tried++;
        }
        
        val += otu_step;
        res += otu_step;
      }
    }
    
    // Insufficient sequences - set all abundances to zero.
    else if (depth < target) {
      double *res_end = res_begin + end_step;
      for (double *res = res_begin; res < res_end; res += otu_step) {
        *res = 0;
      }
    }
    
    res_begin += thread_step;
    val_begin += thread_step;
  }
  
  return NULL;
}

static pthread_func_t setup_dense () {
  
  depth_vec = (uint32_t*) safe_malloc(n_sams * sizeof(uint32_t));
  memset(depth_vec, 0, n_sams * sizeof(uint32_t));
  
  if (margin == 1) {
    
    for (int sam = 0; sam < n_sams; sam++) {
      
      double  depth = 0;
      double *val   = val_vec + sam;
      
      for (int otu = 0; otu < n_otus; otu++) {
        depth += *val;
        val   += n_sams;
      }
      
      depth_vec[sam] = (uint32_t) depth;
    }
  }
  
  else {
    
    for (int sam = 0; sam < n_sams; sam++) {
      
      double  depth     = 0;
      double *val_begin = val_vec + (sam * n_otus);
      
      for (int otu = 0; otu < n_otus; otu++) {
        depth += val_begin[otu];
      }
      
      depth_vec[sam] = (uint32_t) depth;
    }
  }
  
  return rarefy_dense;
}




/*
 * Functions for slam, dgTMatrix, and margin 1 dgCMatrix
 * 
 */

static knuth_t *rarefy_triplet_knuth_vec;

static void *rarefy_triplet(void *arg) {
  
  int      thread_i = *((int *) arg);
  knuth_t *knuth_vec = rarefy_triplet_knuth_vec;
  
  // Initialize RNGs for this thread's samples.
  for (uint64_t sam = thread_i; sam < n_sams; sam += n_threads) {
    knuth_vec[sam].tried = 0;
    knuth_vec[sam].kept  = 0;
    pcg32_srandom_r(&(knuth_vec[sam].rng), seed, sam);
  }
  
  // Iterate over all tuples (sam,otu,val though otu is ignored).
  // Cannot assume any particular ordering.
  for (int i = 0; i < n_vals; i++) {
    
    int sam = sam_vec[i]; // Sample index
    
    // Only work on samples assigned to this thread.
    if (sam % n_threads == thread_i) {
      
      knuth_t   knuth = knuth_vec[sam]; // Tracks: tried, kept, and rng.
      uint32_t  depth = depth_vec[sam]; // Total observations in sample
      double    val   = val_vec[i];     // Current OTU # of observations
      double   *res   = res_vec + i;    // Rarefied OTU # of observations
      
      // Sample can be be rarefied.
      if (depth > target) {
        *res = 0;
        
        // Knuth algorithm for choosing target seqs from depth.
        for (int seq = 0; seq < val && knuth.kept < target; seq++) {
          
          uint32_t not_tried  = depth - knuth.tried;
          uint32_t still_need = target - knuth.kept;
          uint32_t rand_int   = pcg32_random_r(&(knuth.rng));
          
          if (rand_int % not_tried < still_need) {
            (*res)++; // retain this observation
            knuth.kept++;
          }
          
          knuth.tried++;
        }
      }
      
      // Insufficient sequences - set all abundances to zero.
      else if (depth < target) {
        *res = 0;
      }
      
    }
  }
  
  return NULL;
}

static pthread_func_t setup_triplet () {
  
  rarefy_triplet_knuth_vec = (knuth_t*)  safe_malloc(n_sams * sizeof(knuth_t));
  depth_vec                = (uint32_t*) safe_malloc(n_sams * sizeof(uint32_t));
  
  // Use a single pass to sum all samples' depths
  memset(depth_vec, 0, n_sams * sizeof(uint32_t));
  for (int i = 0; i < n_vals; i++) {
    depth_vec[sam_vec[i]] += (uint32_t) val_vec[i];
  }
  
  return rarefy_triplet;
}




/*
 * Base R dense `matrix`
 * 
 */


static pthread_func_t setup_matrix () {
  
  val_vec = REAL(sexp_val_mtx);
  res_vec = REAL(sexp_res_mtx);
  n_vals  = LENGTH(sexp_val_mtx);
  
  if (margin == 1) {
    n_sams  = nrows(sexp_val_mtx);
    n_otus  = ncols(sexp_val_mtx);
  }
  else {
    n_sams  = ncols(sexp_val_mtx);
    n_otus  = nrows(sexp_val_mtx);
  }
  
  return setup_dense();
}



/*
 * `Matrix` package's `dgeMatrix` dense matrix
 * 
 */

static pthread_func_t setup_dgeMatrix () {
  
  SEXP sexp_dim_vec = R_do_slot(sexp_res_mtx, install("Dim"));
  SEXP sexp_val_vec = R_do_slot(sexp_res_mtx, install("x"));
  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  n_protected++;
  R_do_slot_assign(sexp_res_mtx, install("x"), sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    n_sams  = INTEGER(sexp_dim_vec)[0];
    n_otus  = INTEGER(sexp_dim_vec)[1];
  }
  else {
    n_sams  = INTEGER(sexp_dim_vec)[1];
    n_otus  = INTEGER(sexp_dim_vec)[0];
  }
  
  return setup_dense();
}



/*
 * `slam` package's `simple_triplet_matrix`
 * 
 */

static pthread_func_t setup_slam () {
  
  SEXP sexp_val_vec = get(sexp_res_mtx, "v");
  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  n_protected++;
  set(sexp_res_mtx, "v", sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec = INTEGER(get(sexp_res_mtx, "i"));
    n_sams  = INTEGER(get(sexp_res_mtx, "nrow"))[0] + 1;
  }
  else {
    sam_vec = INTEGER(get(sexp_res_mtx, "j"));
    n_sams  = INTEGER(get(sexp_res_mtx, "ncol"))[0] + 1;
  }
  
  return setup_triplet();
}



/*
 * `Matrix` package's `dgTMatrix` triplet matrix
 * 
 */

static pthread_func_t setup_dgTMatrix () {
  
  SEXP sexp_val_vec = R_do_slot(sexp_res_mtx, install("x"));
  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  n_protected++;
  R_do_slot_assign(sexp_res_mtx, install("x"), sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec = INTEGER(R_do_slot(sexp_res_mtx, install("i")));
    n_sams  = INTEGER(R_do_slot(sexp_res_mtx, install("Dim")))[0];
  }
  else {
    sam_vec = INTEGER(R_do_slot(sexp_res_mtx, install("j")));
    n_sams  = INTEGER(R_do_slot(sexp_res_mtx, install("Dim")))[1];
  }
  
  return setup_triplet();
}



/*
 * `Matrix` package's `dgCMatrix` compressed sparse matrix
 * 
 */

static void *rarefy_compressed (void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sam = thread_i; sam < n_sams; sam += n_threads) {
    
    uint32_t depth     = depth_vec[sam];
    int      pos_begin = pos_vec[sam];
    int      pos_end   = pos_vec[sam + 1];
    int      sam_nnz   = pos_end - pos_begin;
    
    
    // Sample can be be rarefied.
    if (depth > target) {
      
      // Seed the PRNG for this sample.
      pcg32_random_t rng;
      pcg32_srandom_r(&rng, seed, sam);
      
      // Knuth algorithm for choosing target seqs from depth.
      uint32_t tried = 0, kept = 0;
      for (int pos = pos_begin; pos < pos_end; pos++) {
        
        double  val = val_vec[pos];  // Current # of observations
        double *res = res_vec + pos; // Rarefied # of observations
        
        *res = 0;
        for (int seq = 0; seq < val && kept < target; seq++) {
          
          uint32_t not_tried  = depth - tried;
          uint32_t still_need = target - kept;
          uint32_t rand_int   = pcg32_random_r(&rng);
          
          if (rand_int % not_tried < still_need) {
            (*res)++; // retain this observation
            kept++;
          }
          
          tried++;
        }
      }
    }
    
    // Insufficient sequences - set all abundances to zero.
    else if (depth < target) {
      memset(res_vec + pos_begin, 0, sam_nnz * sizeof(double));
    }
    
  }
  
  return NULL;
}

static pthread_func_t setup_dgCMatrix () {
  
  SEXP sexp_val_vec = R_do_slot(sexp_res_mtx, install("x"));
  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  n_protected++;
  R_do_slot_assign(sexp_res_mtx, install("x"), sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec = INTEGER(R_do_slot(sexp_res_mtx, install("i")));
    n_sams  = INTEGER(R_do_slot(sexp_res_mtx, install("Dim")))[0];
    return setup_triplet();
  }
  
  else {
    pos_vec = INTEGER(R_do_slot(sexp_res_mtx, install("p")));
    n_sams  = INTEGER(R_do_slot(sexp_res_mtx, install("Dim")))[1];
    
    depth_vec = (uint32_t*) safe_malloc(n_sams * sizeof(uint32_t));
    for (int sam = 0; sam < n_sams; sam++) {
      depth_vec[sam] = 0;
      int pos_begin = pos_vec[sam];
      int pos_end   = pos_vec[sam + 1];
      for (int i = pos_begin; i < pos_end; i++)
        depth_vec[sam] += (uint32_t) val_vec[i];
    }
    
    return rarefy_compressed;
  }
}


static void set_target (SEXP sexp_depth, SEXP sexp_n_samples) {
  
  double depth = asReal(sexp_depth);
  
  // Set target depth according to number/pct of samples to keep/drop.
  if (!isNull(sexp_n_samples)) {
    
    double n_samples = asReal(sexp_n_samples);
    
    if (n_samples == 0)      n_samples = n_sams;             // Keep all
    if (n_samples > n_sams)  n_samples = n_sams;             // Keep all
    if (fabs(n_samples) < 1) n_samples = n_sams * n_samples; // Keep/drop percentage
    if (n_samples <= -1)     n_samples = n_sams + n_samples; // Drop n_samples
    n_samples = (n_samples <= 1) ? 1.0 : floor(n_samples);   // Keep at least one
    
    int k = (int)(n_sams - n_samples);
    
    // Partially sort a copy of depth_vec in ascending order: only [k] is needed.
    int *quickselect_vec = (int*) safe_malloc(n_sams * sizeof(int));
    for (int sam = 0; sam < n_sams; sam++) {
      quickselect_vec[sam] = (int) (depth_vec[sam]);
    }
    iPsort(quickselect_vec, n_sams, k);
    
    target = (uint32_t) quickselect_vec[k];
    free_one(quickselect_vec);
  }
  
  // Depth is given as minimum percent of observations to keep.
  else if (depth < 1) {
    
    double depth_sum = 0;
    for (int sam = 0; sam < n_sams; sam++) {
      depth_sum += depth_vec[sam];
    }
    uint32_t min_depth_sum = (uint32_t) (depth_sum * depth);
    
    // Sort a copy of depth_vec in ascending order
    int *sorted_vec = (int*) safe_malloc(n_sams * sizeof(int));
    for (int sam = 0; sam < n_sams; sam++) {
      sorted_vec[sam] = (int) (depth_vec[sam]);
    }
    R_isort(sorted_vec, n_sams);
    
    for (int sam = 0; sam < n_sams; sam++) {
      target = (uint32_t) (sorted_vec[sam]);
      if (target * n_sams >= min_depth_sum) break;
    }
    free_one(sorted_vec);
  }
  
  // Depth is given as observations per sample to keep.
  else {
    target = (uint32_t) depth;
  }
  
}



//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_rarefy(
    SEXP sexp_otu_mtx,   SEXP sexp_depth,
    SEXP sexp_n_samples, SEXP sexp_seed,
    SEXP sexp_margin,    SEXP sexp_n_threads ) {
  
  init_n_ptrs(10);
  
  seed       = (uint64_t) asInteger(sexp_seed);
  margin     = asInteger(sexp_margin);
  n_threads  = asInteger(sexp_n_threads);
  
  
  /*
   * Shallow copy input matrix to result matrix.
   * For `VECSXP`s like slam and Matrix objects,
   * the underlying `REALSXP`s are still shared by
   * the input and result `SEXP`s at this point.
   */
  sexp_val_mtx = sexp_otu_mtx;
  sexp_res_mtx = PROTECT(duplicate(sexp_otu_mtx));
  n_protected  = 1;
  
  
  // function to run
  // void * (*rarefy_func)(void *) = NULL;
  pthread_func_t rarefy_func = NULL;
  
  
  // Select worker function and set *_vec and n_* variables.
  if (isMatrix(sexp_otu_mtx))                               { rarefy_func = setup_matrix();    }
  else if (inherits(sexp_otu_mtx, "simple_triplet_matrix")) { rarefy_func = setup_slam();      }
  else if (inherits(sexp_otu_mtx, "dgCMatrix"))             { rarefy_func = setup_dgCMatrix(); }
  else if (inherits(sexp_otu_mtx, "dgTMatrix"))             { rarefy_func = setup_dgTMatrix(); }
  else if (inherits(sexp_otu_mtx, "dgeMatrix"))             { rarefy_func = setup_dgeMatrix(); }
  else   { error("Unrecognized matrix format."); } // # nocov
  
  // Get the target rarefaction depth.
  // Note that `n_samples` isn't `n_sams`
  set_target(sexp_depth, sexp_n_samples);
  
  // free_all();
  // UNPROTECT(n_protected);
  // return ScalarInteger((int)target);
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_sams > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = (pthread_t*) R_alloc(n_threads, sizeof(pthread_t));
      int       *args = (int*)       R_alloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, rarefy_func, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free_all();
      UNPROTECT(n_protected);
      return sexp_res_mtx;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  rarefy_func(&thread_i);
  
  free_all();
  UNPROTECT(n_protected);
  return sexp_res_mtx;
}
