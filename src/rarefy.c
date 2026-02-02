// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


typedef void *(*pthread_func_t)(void *);

static uint32_t  target;
static uint64_t  seed;
static int       margin;
static int       n_threads;
static SEXP      sexp_val_mtx;
static SEXP      sexp_res_mtx;
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
  int     thread_step = n_threads * sam_step;
  double *val_begin   = val_vec + thread_i * sam_step;
  double *res_begin   = res_vec + thread_i * sam_step;
  
  for (int sam = thread_i; sam < n_sams; sam += n_threads) {
    
    uint32_t depth = depth_vec[sam];
    
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
    
    res_begin += thread_step;
    val_begin += thread_step;
  }
  
  return NULL;
}

static pthread_func_t setup_dense(void) {
  
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
      
      knuth_t  *knuth = &knuth_vec[sam]; // Tracks: tried, kept, and rng.
      uint32_t  depth = depth_vec[sam];  // Total observations in sample
      double    val   = val_vec[i];      // Current OTU # of observations
      double   *res   = res_vec + i;     // Rarefied OTU # of observations
      
      // Sample can be be rarefied.
      if (depth > target) {
        *res = 0;
        
        // Knuth algorithm for choosing target seqs from depth.
        for (int seq = 0; seq < val && knuth->kept < target; seq++) {
          
          uint32_t not_tried  = depth - knuth->tried;
          uint32_t still_need = target - knuth->kept;
          uint32_t rand_int   = pcg32_random_r(&(knuth->rng));
          
          if (rand_int % not_tried < still_need) {
            (*res)++; // retain this observation
            knuth->kept++;
          }
          
          knuth->tried++;
        }
      }
      
    }
  }
  
  return NULL;
}

static pthread_func_t setup_triplet(void) {
  
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


static pthread_func_t setup_matrix(void) {
  
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

static pthread_func_t setup_dgeMatrix(void) {
  
  SEXP sexp_dim_vec = PROTECT(R_do_slot(sexp_res_mtx, install("Dim")));
  SEXP sexp_val_vec = PROTECT(R_do_slot(sexp_res_mtx, install("x")));
  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
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
  
  pthread_func_t rarefy_func = setup_dense();

  UNPROTECT(3);
  return rarefy_func;
}



/*
 * `slam` package's `simple_triplet_matrix`
 * 
 */

static pthread_func_t setup_slam(void) {
  
  SEXP sexp_val_vec = PROTECT(get(sexp_res_mtx, "v"));
  SEXP sexp_i       = PROTECT(get(sexp_res_mtx, "i"));
  SEXP sexp_j       = PROTECT(get(sexp_res_mtx, "j"));
  SEXP sexp_nrow    = PROTECT(get(sexp_res_mtx, "nrow"));
  SEXP sexp_ncol    = PROTECT(get(sexp_res_mtx, "ncol"));

  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  set(sexp_res_mtx, "v", sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec = INTEGER(sexp_i);
    n_sams  = INTEGER(sexp_nrow)[0] + 1;
  }
  else {
    sam_vec = INTEGER(sexp_j);
    n_sams  = INTEGER(sexp_ncol)[0] + 1;
  }
  
  pthread_func_t rarefy_func = setup_triplet();
  
  UNPROTECT(6);
  return rarefy_func;
}



/*
 * `Matrix` package's `dgTMatrix` triplet matrix
 * 
 */

static pthread_func_t setup_dgTMatrix(void) {
  
  SEXP sexp_val_vec = PROTECT(R_do_slot(sexp_res_mtx, install("x")));
  SEXP sexp_i       = PROTECT(R_do_slot(sexp_res_mtx, install("i")));
  SEXP sexp_j       = PROTECT(R_do_slot(sexp_res_mtx, install("j")));
  SEXP sexp_dim     = PROTECT(R_do_slot(sexp_res_mtx, install("Dim")));

  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  R_do_slot_assign(sexp_res_mtx, install("x"), sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec = INTEGER(sexp_i);
    n_sams  = INTEGER(sexp_dim)[0];
  }
  else {
    sam_vec = INTEGER(sexp_j);
    n_sams  = INTEGER(sexp_dim)[1];
  }
  
  pthread_func_t rarefy_func = setup_triplet();
  
  UNPROTECT(5);
  return rarefy_func;
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
    
    
    // Sample can be be rarefied.
    if (depth > target) {
      
      // Seed the PRNG for this sample.
      pcg32_random_t rng;
      pcg32_srandom_r(&rng, seed, sam);
      
      // Knuth algorithm for choosing target seqs from depth.
      uint32_t tried = 0, kept = 0; // These are local to the sample
      for (int pos = pos_begin; pos < pos_end; pos++) {
        
        double  val = val_vec[pos];  // Current # of observations
        double *res = res_vec + pos; // Rarefied # of observations
        
        *res = 0;
        for (uint32_t seq = 0; seq < val && kept < target; seq++) {
          
          uint32_t not_tried  = depth - tried;
          uint32_t still_need = target - kept;
          uint32_t rand_int   = pcg32_random_r(&rng);
          
          if (rand_int % not_tried < still_need) {
            (*res)++;
            kept++;
          }
          
          tried++;
        }
      }
    }
    
  }
  
  return NULL;
}

static pthread_func_t setup_dgCMatrix(void) {
  
  pthread_func_t rarefy_func = NULL;

  SEXP sexp_val_vec = PROTECT(R_do_slot(sexp_res_mtx, install("x")));
  SEXP sexp_i       = PROTECT(R_do_slot(sexp_res_mtx, install("i")));
  SEXP sexp_p       = PROTECT(R_do_slot(sexp_res_mtx, install("p")));
  SEXP sexp_dim     = PROTECT(R_do_slot(sexp_res_mtx, install("Dim")));

  SEXP sexp_res_vec = PROTECT(duplicate(sexp_val_vec));
  R_do_slot_assign(sexp_res_mtx, install("x"), sexp_res_vec);
  
  val_vec = REAL(sexp_val_vec);
  res_vec = REAL(sexp_res_vec);
  n_vals  = LENGTH(sexp_val_vec);
  
  if (margin == 1) {
    sam_vec     = INTEGER(sexp_i);
    n_sams      = INTEGER(sexp_dim)[0];
    rarefy_func = setup_triplet();
  }
  
  else {
    pos_vec = INTEGER(sexp_p);
    n_sams  = INTEGER(sexp_dim)[1];
    
    depth_vec = (uint32_t*) safe_malloc(n_sams * sizeof(uint32_t));
    for (int sam = 0; sam < n_sams; sam++) {
      depth_vec[sam] = 0;
      int pos_begin = pos_vec[sam];
      int pos_end   = pos_vec[sam + 1];
      for (int i = pos_begin; i < pos_end; i++)
        depth_vec[sam] += (uint32_t) val_vec[i];
    }
    
    rarefy_func = rarefy_compressed;
  }

  UNPROTECT(5);
  return rarefy_func;
}


/*
 * Compacts a slam::simple_triplet_matrix (S3 object)
 * uses "v", "i", "j" components.
 */
static void compact_slam(SEXP sexp_mtx) {
    SEXP sexp_v = PROTECT(get(sexp_mtx, "v"));
    SEXP sexp_i = PROTECT(get(sexp_mtx, "i"));
    SEXP sexp_j = PROTECT(get(sexp_mtx, "j"));
    
    double *v     = REAL(sexp_v);
    int    *i     = INTEGER(sexp_i);
    int    *j     = INTEGER(sexp_j);
    int     n_all = LENGTH(sexp_v);
    
    // 1. Count new non-zeros
    int nnz_new = 0;
    for (int k = 0; k < n_all; k++) {
        if (v[k] != 0.0) nnz_new++;
    }
    
    if (nnz_new == n_all) { UNPROTECT(3); return; }
    
    // 2. Allocate new vectors
    SEXP sexp_new_v = PROTECT(allocVector(REALSXP, nnz_new));
    SEXP sexp_new_i = PROTECT(allocVector(INTSXP,  nnz_new));
    SEXP sexp_new_j = PROTECT(allocVector(INTSXP,  nnz_new));
    
    double *new_v = REAL(sexp_new_v);
    int    *new_i = INTEGER(sexp_new_i);
    int    *new_j = INTEGER(sexp_new_j);
    
    // 3. Filter data
    int idx = 0;
    for (int k = 0; k < n_all; k++) {
        if (v[k] != 0.0) {
            new_v[idx] = v[k];
            new_i[idx] = i[k];
            new_j[idx] = j[k];
            idx++;
        }
    }
    
    // 4. Update the list components
    set(sexp_mtx, "v", sexp_new_v);
    set(sexp_mtx, "i", sexp_new_i);
    set(sexp_mtx, "j", sexp_new_j);
    
    UNPROTECT(6);
}


/*
 * Compacts a dgCMatrix by removing explicit zeros (0.0) from the
 * x and i slots and updating the p slot.
 */
static void compact_dgCMatrix(SEXP sexp_mtx) {
    // 1. Access the slots
    SEXP sexp_p = PROTECT(R_do_slot(sexp_mtx, install("p")));
    SEXP sexp_i = PROTECT(R_do_slot(sexp_mtx, install("i")));
    SEXP sexp_x = PROTECT(R_do_slot(sexp_mtx, install("x")));

    int    *p    = INTEGER(sexp_p);
    int    *i    = INTEGER(sexp_i);
    double *x    = REAL(sexp_x);
    int     ncol = LENGTH(sexp_p) - 1;
    int     nnz_old = p[ncol];

    // 2. Pass 1: Count new non-zeros (nnz)
    int nnz_new = 0;
    for (int k = 0; k < nnz_old; k++) {
        if (x[k] != 0.0) nnz_new++;
    }

    // Optimization: If no zeros were created, do nothing.
    if (nnz_new == nnz_old) { UNPROTECT(3); return; }

    // 3. Allocate new x and i vectors of the correct size
    SEXP sexp_new_i = PROTECT(allocVector(INTSXP, nnz_new));
    SEXP sexp_new_x = PROTECT(allocVector(REALSXP, nnz_new));
    int    *new_i = INTEGER(sexp_new_i);
    double *new_x = REAL(sexp_new_x);

    // 4. Pass 2: Compact data and update 'p'
    // We must be careful updating 'p' in-place because it defines 
    // the boundaries of the data we are currently reading.
    
    int write_idx = 0;
    int col_start = p[0]; // Always 0
    
    for (int col = 0; col < ncol; col++) {
        int col_end = p[col + 1]; // Save old end before overwriting

        for (int k = col_start; k < col_end; k++) {
            if (x[k] != 0.0) {
                new_x[write_idx] = x[k];
                new_i[write_idx] = i[k];
                write_idx++;
            }
        }
        
        // Update p for the *next* column boundary
        p[col + 1] = write_idx; 
        
        // Move to next column
        col_start = col_end; 
    }

    // 5. Assign new slots back to the matrix
    R_do_slot_assign(sexp_mtx, install("i"), sexp_new_i);
    R_do_slot_assign(sexp_mtx, install("x"), sexp_new_x);

    UNPROTECT(5); // Unprotect new_i and new_x
}


/*
 * Compacts a Matrix::dgTMatrix (S4 object)
 * uses "x", "i", "j" slots.
 */
static void compact_dgTMatrix(SEXP sexp_mtx) {
    SEXP sexp_x = PROTECT(R_do_slot(sexp_mtx, install("x")));
    SEXP sexp_i = PROTECT(R_do_slot(sexp_mtx, install("i")));
    SEXP sexp_j = PROTECT(R_do_slot(sexp_mtx, install("j")));
    
    double *x     = REAL(sexp_x);
    int    *i     = INTEGER(sexp_i);
    int    *j     = INTEGER(sexp_j);
    int     n_all = LENGTH(sexp_x);
    
    // 1. Count new non-zeros
    int nnz_new = 0;
    for (int k = 0; k < n_all; k++) {
        if (x[k] != 0.0) nnz_new++;
    }
    
    if (nnz_new == n_all) { UNPROTECT(3); return; }
    
    // 2. Allocate new vectors
    SEXP sexp_new_x = PROTECT(allocVector(REALSXP, nnz_new));
    SEXP sexp_new_i = PROTECT(allocVector(INTSXP,  nnz_new));
    SEXP sexp_new_j = PROTECT(allocVector(INTSXP,  nnz_new));
    
    double *new_x = REAL(sexp_new_x);
    int    *new_i = INTEGER(sexp_new_i);
    int    *new_j = INTEGER(sexp_new_j);
    
    // 3. Filter data
    int idx = 0;
    for (int k = 0; k < n_all; k++) {
        if (x[k] != 0.0) {
            new_x[idx] = x[k];
            new_i[idx] = i[k];
            new_j[idx] = j[k];
            idx++;
        }
    }
    
    // 4. Update the S4 slots
    R_do_slot_assign(sexp_mtx, install("x"), sexp_new_x);
    R_do_slot_assign(sexp_mtx, install("i"), sexp_new_i);
    R_do_slot_assign(sexp_mtx, install("j"), sexp_new_j);
    
    UNPROTECT(6);
}



//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_rarefy(
    SEXP sexp_otu_mtx, SEXP sexp_depth,
    SEXP sexp_seed,    SEXP sexp_margin,
    SEXP sexp_n_threads ) {
  
  init_n_ptrs(10);
  
  target    = (uint32_t) asInteger(sexp_depth);
  seed      = (uint64_t) asInteger(sexp_seed);
  margin    = asInteger(sexp_margin);
  n_threads = asInteger(sexp_n_threads);
  
  
  /*
   * Shallow copy input matrix to result matrix.
   * For `VECSXP`s like slam and Matrix objects,
   * the underlying `REALSXP`s are still shared by
   * the input and result `SEXP`s at this point.
   */
  sexp_val_mtx = sexp_otu_mtx;
  sexp_res_mtx = PROTECT(duplicate(sexp_otu_mtx));
  
  
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
  
      // Post-process: Remove explicit zeros to restore sparsity
      if      (inherits(sexp_res_mtx, "simple_triplet_matrix")) { compact_slam(sexp_res_mtx);      }
      else if (inherits(sexp_res_mtx, "dgCMatrix"))             { compact_dgCMatrix(sexp_res_mtx); }
      else if (inherits(sexp_res_mtx, "dgTMatrix"))             { compact_dgTMatrix(sexp_res_mtx); }
      
      free_all();
      UNPROTECT(1);
      return sexp_res_mtx;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  rarefy_func(&thread_i);
  
  // Post-process: Remove explicit zeros to restore sparsity
  if      (inherits(sexp_res_mtx, "simple_triplet_matrix")) { compact_slam(sexp_res_mtx);      }
  else if (inherits(sexp_res_mtx, "dgCMatrix"))             { compact_dgCMatrix(sexp_res_mtx); }
  else if (inherits(sexp_res_mtx, "dgTMatrix"))             { compact_dgTMatrix(sexp_res_mtx); }
  
  free_all();
  UNPROTECT(1);
  return sexp_res_mtx;
}
