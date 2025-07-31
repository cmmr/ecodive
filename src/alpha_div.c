#include <R.h>
#include <Rinternals.h>
#include <math.h>   // log
#include <stdlib.h> // calloc, free

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define CHAO1       1
#define INV_SIMPSON 2
#define SHANNON     3
#define SIMPSON     4

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_threads;
static double *result_vec;



/*
 * START_SAMPLE_LOOP and END_SAMPLE_LOOP define a simple loop 
 * for iterating over all samples, ensuring that each is 
 * assigned to only a single thread. Provides `otu_vec` and 
 * `result` initialized to 0; expects `result` at the end.
 */

#define START_SAMPLE_LOOP                                      \
  int thread_i = *((int *) arg);                               \
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {\
    double *otu_vec = otu_mtx    + (sample * n_otus);          \
    double  result  = 0;


#define END_SAMPLE_LOOP                                        \
    result_vec[sample] = result;                               \
  }                                                            \
  return NULL;


#define CALCULATE_DEPTH                                        \
  double depth = 0;                                            \
  for (int otu = 0; otu < n_otus; otu++) {                     \
    depth += otu_vec[otu];                                     \
  }
  


//======================================================
// Chao1
// sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))
//======================================================
static void *calc_chao1(void *arg) {
  START_SAMPLE_LOOP
    
  double nnz  = 0;
  double ones = 0;
  double twos = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      nnz++;
      if      (x <= 1) ones++;
      else if (x <= 2) twos++;
    }
  }
  
  result = nnz + ((ones * ones) / (2 * twos));
  
  END_SAMPLE_LOOP
}


//======================================================
// Inverse Simpson
// p <- x / sum(x)
// 1 / sum(p ** 2)
//======================================================
static void *calc_inv_simpson(void *arg) {
  START_SAMPLE_LOOP
  CALCULATE_DEPTH
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      double p = x / depth;
      result += p * p;
    }
  }
  
  if (result != 0) result = 1 / result;
    
  END_SAMPLE_LOOP
}


//======================================================
// Shannon
// p <- x / sum(x)
// -sum(p * log(p))
//======================================================
static void *calc_shannon(void *arg) {
  START_SAMPLE_LOOP
  CALCULATE_DEPTH
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      double p = x / depth;
      result += p * log(p);
    }
  }
  result *= -1;
  
  END_SAMPLE_LOOP
}


//======================================================
// Simpson
// p <- x / sum(x)
// 1 - sum(p ** 2)
//======================================================
static void *calc_simpson(void *arg) {
  START_SAMPLE_LOOP
  CALCULATE_DEPTH
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      double p = x / depth;
      result += p * p;
    }
  }
  result = 1 - result;
  
  END_SAMPLE_LOOP
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx, 
    SEXP sexp_n_threads, SEXP sexp_result_vec ) {
  
  algorithm  = asInteger(sexp_algorithm);
  otu_mtx    = REAL(sexp_otu_mtx);
  n_otus     = nrows(sexp_otu_mtx);
  n_samples  = ncols(sexp_otu_mtx);
  n_threads  = asInteger(sexp_n_threads);
  result_vec = REAL(sexp_result_vec);
  
  
  // function to run
  void * (*calc_adiv)(void *) = NULL;
  switch (algorithm) {
    case CHAO1:       calc_adiv = calc_chao1;       break;
    case INV_SIMPSON: calc_adiv = calc_inv_simpson; break;
    case SHANNON:     calc_adiv = calc_shannon;     break;
    case SIMPSON:     calc_adiv = calc_simpson;     break;
  }
  
  if (calc_adiv == NULL) { // # nocov start
    error("Invalid alpha diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_samples > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      if (tids == NULL || args == NULL) { // # nocov start
        free(tids); free(args);
        error("Insufficient memory for parallel alpha diversity calculation.");
        return R_NilValue;
      } // # nocov end
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_adiv, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args);
      
      return sexp_result_vec;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_adiv(&thread_i);
  
  return sexp_result_vec;
}
