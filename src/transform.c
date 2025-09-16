// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <math.h> // exp, log, sqrt

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <stdlib.h> // calloc, free
#    define HAVE_PTHREAD
#  endif
#endif

#define TRANSFORM_PCT    1
#define TRANSFORM_CLR    2
#define TRANSFORM_CHORD  3

static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     algorithm;
static int     n_threads;
static SEXP   *sexp_extra;
static double *result_mtx;
static int     mtx_len;


//======================================================
// Percent of Total
// x / sum(x)
//======================================================
static void *transform_pct(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    // Sum the counts for this sample
    double depth = 0;
    for (int i = sample; i < mtx_len; i += n_samples)
      depth += otu_mtx[i];
    
    for (int i = sample; i < mtx_len; i += n_samples)
      result_mtx[i] = otu_mtx[i] / depth;
  }
  
  return NULL;
}


//======================================================
// Centered Log Ratio
// log(x / exp(mean(log(x))))
//======================================================
static void *transform_clr(void *arg) {
  
  int    thread_i    = *((int *) arg);
  double pseudocount = asReal(*sexp_extra);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double norm = 0;
    for (int i = sample; i < mtx_len; i += n_samples)
      norm += log(otu_mtx[i] + pseudocount);
    
    norm = exp(norm / n_otus);
    
    for (int i = sample; i < mtx_len; i += n_samples)
      result_mtx[i] = log((otu_mtx[i] + pseudocount) / norm);
  }
  
  return NULL;
}


//======================================================
// Chord-Transformed Relative Abundance
// x / sqrt(sum(x ^ 2))
//======================================================
static void *transform_chord(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double sq_sum = 0;
    
    for (int i = sample; i < mtx_len; i += n_samples)
      sq_sum += otu_mtx[i] * otu_mtx[i];
    
    sq_sum = sqrt(sq_sum);
    
    for (int i = sample; i < mtx_len; i += n_samples)
      result_mtx[i] = otu_mtx[i] / sq_sum;
  }
  
  return NULL;
}




//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_transform(
    SEXP sexp_otu_mtx,   SEXP sexp_algorithm, 
    SEXP sexp_n_threads, SEXP sexp_extra_args ) {
  
  otu_mtx    = REAL(sexp_otu_mtx);
  n_otus     = ncols(sexp_otu_mtx);
  n_samples  = nrows(sexp_otu_mtx);
  algorithm  = asInteger(sexp_algorithm);
  n_threads  = asInteger(sexp_n_threads);
  sexp_extra = &sexp_extra_args;
  
  mtx_len = n_otus * n_samples;
  
  // Copy `otu_mtx` to a new object named `result_mtx`.
  SEXP sexp_result_mtx = duplicate(sexp_otu_mtx);
  PROTECT(sexp_result_mtx);
  result_mtx = REAL(sexp_result_mtx);
  
  
  // function to run
  void * (*transformer)(void *) = NULL;
  
  switch (algorithm) {
    case TRANSFORM_PCT:   transformer = transform_pct;   break;
    case TRANSFORM_CLR:   transformer = transform_clr;   break;
    case TRANSFORM_CHORD: transformer = transform_chord; break;
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_samples > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, transformer, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free(tids); free(args);
      
      UNPROTECT(1);
      return sexp_result_mtx;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  transformer(&thread_i);
  
  UNPROTECT(1);
  return sexp_result_mtx;
}

