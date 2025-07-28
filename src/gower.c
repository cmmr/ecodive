#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // calloc, free
#include <string.h> // memset
#include <math.h>   // fabs

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  double *rescale_vec_1;
  double *rescale_vec_2;
  double *distance;
} pair_t;


//======================================================
// Variables shared between main and worker threads.
//======================================================
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static pair_t *pair_vec;
static int     n_pairs;
static int     weighted;
static int     n_threads;
static double *dist_vec;
static double *rescale_mtx;


static void *calc_rescale_mtx(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int otu = thread_i; otu < n_otus; otu += n_threads) {
    
    double min_count = otu_mtx[otu];
    double max_count = min_count;
    
    for (int sample = 1; sample < n_samples; sample++) {
      double count = otu_mtx[sample * n_otus + otu];
      if      (count < min_count) { min_count = count; }
      else if (count > max_count) { max_count = count; }
    }
    
    // Rescale to 0 - 1
    if (max_count > min_count) {
      double range = max_count - min_count;
      for (int sample = 0; sample < n_samples; sample++) {
        int idx = sample * n_otus + otu;
        rescale_mtx[idx] = (otu_mtx[idx] - min_count) / range;
      }
    }
    
  }
  
  return NULL;
}


static void *calc_dist_vec(void *arg) {
  
  int thread_i = *((int *) arg);
  
  //======================================================
  // Gower Weighted
  // sum(abs(x-y)) / n_otus
  //======================================================
  if (weighted) {
    
    for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
      
      pair_t *pair = pair_vec + pair_i;
      
      // pointers to each sample's column in rescale_mtx
      double *x_vec = pair->rescale_vec_1;
      double *y_vec = pair->rescale_vec_2;
      
      double distance = 0;
      
      for (int otu = 0; otu < n_otus; otu++) {
        
        // abundance of this OTU in the two samples
        double x = x_vec[otu];
        double y = y_vec[otu];
        
        // accumulate
        distance += fabs(x - y);
      }
      
      // value to return
      *(pair->distance) = distance / n_otus;
    }
  }
  
  
  //======================================================
  // Gower Unweighted
  // A = sum(x>0); B = sum(y>0); J = sum(x>0 & y>0)
  // (A + B - 2*J) / n_otus
  //======================================================
  else {
    
    for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
      
      pair_t *pair = pair_vec + pair_i;
      
      // pointers to each sample's column in rescale_mtx
      double *x_vec = pair->rescale_vec_1;
      double *y_vec = pair->rescale_vec_2;
      
      int A = 0, B = 0, J = 0;
      
      for (int otu = 0; otu < n_otus; otu++) {
        
        char a = x_vec[otu] > 0;
        char b = y_vec[otu] > 0;
        
        // accumulate if appropriate
        if (a && b) { A++; B++; J++; }
        else if (a) { A++; }
        else if (b) { B++; }
      }
      
      // value to return
      *(pair->distance) = (double)(A + B - 2*J) / n_otus;
    }
  }
  
  
  return NULL;
}



//======================================================
// R interface. Dispatches threads to bdiv algorithms.
//======================================================
SEXP C_gower(
    SEXP sexp_otu_mtx,      SEXP sexp_weighted, 
    SEXP sexp_pair_idx_vec, SEXP sexp_n_threads, 
    SEXP sexp_result_dist ) {
  
  otu_mtx   = REAL(sexp_otu_mtx);
  n_otus    = nrows(sexp_otu_mtx);
  n_samples = ncols(sexp_otu_mtx);
  
  weighted = asLogical(sexp_weighted);
  
  int *pair_idx_vec = INTEGER(sexp_pair_idx_vec);
  n_pairs           = LENGTH(sexp_pair_idx_vec);
  
  n_threads = asInteger(sexp_n_threads);
  dist_vec  = REAL(sexp_result_dist);
  
  
  pair_vec    = calloc(n_pairs,            sizeof(pair_t));
  rescale_mtx = calloc(n_otus * n_samples, sizeof(double));
  
  if (pair_vec == NULL || rescale_mtx == NULL) { // # nocov start
    free(pair_vec); free(rescale_mtx);
    error("Unable to allocate memory for Gower calculation.");
    return R_NilValue;
  } // # nocov end
  
  memset(rescale_mtx, 0, n_otus * n_samples * sizeof(double));
  
  
  // pointers to input and output for 
  // each pairwise comparison
  int pair_idx = 0;
  int dist_idx = 0;
  for (int i = 0; i < n_samples - 1; i++) {
    for (int j = i + 1; j < n_samples; j++) {
      if (pair_idx_vec[pair_idx] == dist_idx) {
        pair_t *pair        = pair_vec + pair_idx;
        pair->rescale_vec_1 = rescale_mtx + (i * n_otus);
        pair->rescale_vec_2 = rescale_mtx + (j * n_otus);
        pair->distance      = dist_vec + dist_idx;
        pair_idx++;
      }
      dist_idx++;
    }
  }
  
  
  // Run WITH multithreading
#ifdef HAVE_PTHREAD
  if (n_threads > 1) {
    
    // threads and their thread_i arguments
    pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
    int       *args = calloc(n_threads, sizeof(int));
    
    int i, n = n_threads;
    for (i = 0; i < n; i++) args[i] = i;
    for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_rescale_mtx, &args[i]);
    for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
    for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
    for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
    
    free(tids); free(args); free(pair_vec); free(rescale_mtx);
    
    return sexp_result_dist;
  }
#endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_rescale_mtx(&thread_i);
  calc_dist_vec(&thread_i);
  
  free(pair_vec); free(rescale_mtx);
  
  return sexp_result_dist;
}

