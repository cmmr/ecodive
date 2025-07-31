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
static int    *pairs_vec;
static int     n_pairs;
static int     weighted;
static int     n_threads;
static double *dist_vec;
static double *rescale_mtx;
static double *last_sample;


/*
 * The START_PAIR_LOOP and END_PAIR_LOOP macros efficiently 
 * loop through all combinations of samples. Skips unwanted 
 * pairings and pairings not assigned to the current thread. 
 * Ensures all threads process the same number of pairs.
 * 
 * After calling START_PAIR_LOOP the code can expect `x_vec` 
 * and `y_vec` to point to the two samples' columns in 
 * `rescale_mtx`. The code should assign to `distance` before 
 * calling END_PAIR_LOOP.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */

#define START_PAIR_LOOP                                        \
  int     thread_i = *((int *) arg);                           \
  double *x_vec    = rescale_mtx;                              \
  double *y_vec    = rescale_mtx + n_otus;                     \
  int     pair_idx = 0;                                        \
  int     dist_idx = 0;                                        \
  double distance;                                             \
  while (pair_idx < n_pairs) {                                 \
    if (pairs_vec[pair_idx] == dist_idx) {                     \
      if (pair_idx % n_threads == thread_i) {


#define END_PAIR_LOOP                                          \
        dist_vec[dist_idx] = distance;                         \
      }                                                        \
      pair_idx++;                                              \
    }                                                          \
    if (y_vec == last_sample) {                                \
      x_vec += n_otus;                                         \
      if (x_vec == last_sample) return NULL;                   \
      y_vec = x_vec + n_otus;                                  \
    } else {                                                   \
      y_vec += n_otus;                                         \
    }                                                          \
    dist_idx++;                                                \
  }                                                            \
  return NULL;                                                 \


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


//======================================================
// Gower
// sum(abs(x-y)) / n_otus
//======================================================
static void *gower_w(void *arg) {
  START_PAIR_LOOP
  
  distance = 0;
    
  for (int otu = 0; otu < n_otus; otu++)
    distance += fabs(x_vec[otu] - y_vec[otu]);
    
  distance = distance / n_otus;
  
  END_PAIR_LOOP
}


static void *gower_u(void *arg) {
  START_PAIR_LOOP
  
  int A = 0, B = 0, J = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    if (x_vec[otu] > 0) {
      A++;
      if (y_vec[otu] > 0) {
        B++; J++;
      }
    }
    else if (y_vec[otu] > 0) {
      B++;
    }
  }
  
  distance = (double)(A + B - 2*J) / n_otus;
  
  END_PAIR_LOOP
}



//======================================================
// R interface. Dispatches threads to bdiv algorithms.
//======================================================
SEXP C_gower(
    SEXP sexp_otu_mtx,   SEXP sexp_weighted, 
    SEXP sexp_pairs_vec, SEXP sexp_n_threads, 
    SEXP sexp_result_dist ) {
  
  otu_mtx   = REAL(sexp_otu_mtx);
  n_otus    = nrows(sexp_otu_mtx);
  n_samples = ncols(sexp_otu_mtx);
  weighted  = asLogical(sexp_weighted);
  pairs_vec = INTEGER(sexp_pairs_vec);
  n_pairs   = LENGTH(sexp_pairs_vec);
  n_threads = asInteger(sexp_n_threads);
  dist_vec  = REAL(sexp_result_dist);
  
  
  
  rescale_mtx = calloc(n_otus * n_samples, sizeof(double));
  
  if (rescale_mtx == NULL) { // # nocov start
    free(rescale_mtx);
    error("Insufficient memory for Gower calculation.");
    return R_NilValue;
  } // # nocov end
  
  memset(rescale_mtx, 0, n_otus * n_samples * sizeof(double));
  last_sample = rescale_mtx + ((n_samples - 1) * n_otus);
  
  
  // function to run
  void * (*calc_dist_vec)(void *) = NULL;
  if (weighted) { calc_dist_vec = gower_w; }
  else          { calc_dist_vec = gower_u; }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      if (tids == NULL || args == NULL) { // # nocov start
        free(tids); free(args);
        free(rescale_mtx);
        error("Insufficient memory for parallel Gower metric calculation.");
        return R_NilValue;
      } // # nocov end
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_rescale_mtx, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args);
      free(rescale_mtx);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_rescale_mtx(&thread_i);
  calc_dist_vec(&thread_i);
  
  free(rescale_mtx);
  
  return sexp_result_dist;
}

