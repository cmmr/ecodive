// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <math.h> // exp, log, sqrt
#include "ecomatrix.h"
#include "memory.h"


// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif


#define NORM_PERCENT 1
#define NORM_CLR     2
#define NORM_CHORD   3
#define NORM_BINARY  4

static int     n_threads;
static int     n_samples;
static int     n_otus;
static int    *pos_vec;
static double *val_vec;
static double *clr_vec;
static double  pseudocount;


// Macro to loop over each sample based on current threading setup.
// Sets sam (sample index), val_begin, and val_end (pointers to val_vec 
// for this sample's first value and the next sample's first value).
#define FOREACH_SAMPLE(block)                                       \
  do {                                                              \
    int thread_i = *((int *) arg);                                  \
    for (int sam = thread_i; sam < n_samples; sam += n_threads) {   \
      double *val_begin = val_vec + pos_vec[sam];                   \
      double *val_end   = val_vec + pos_vec[sam + 1];               \
      block                                                         \
    }                                                               \
  } while (0)

// Marco to loop over each value for the current sample.
// Sets val to point to each value in val_vec in turn.
#define FOREACH_VAL(block)                                          \
  do {                                                              \
    for (double *val = val_begin; val != val_end; val++) {          \
      block;                                                        \
    }                                                               \
  } while (0)



//======================================================
// Percent of Total
// x / sum(x)
//======================================================
static void *norm_percent(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    FOREACH_VAL(depth += *val);
    FOREACH_VAL(*val /= depth);
  );
  return NULL;
}


//======================================================
// Centered Log Ratio
// log(x / exp(mean(log(x))))
//======================================================
static void *norm_clr(void *arg) {
  FOREACH_SAMPLE(
    
    double logsum = 0;
    FOREACH_VAL(logsum += log(*val + pseudocount));
    
    // Add in all the zero abundances
    logsum += log(pseudocount) * (n_otus - (val_end - val_begin));
    
    double denominator = exp(logsum / n_otus);
    FOREACH_VAL(*val = log((*val + pseudocount) / denominator));
    
    // CLR value for zero abundance OTUs
    clr_vec[sam] = log(pseudocount / denominator);
    
  );
  return NULL;
}


//======================================================
// Chord-Transformed Relative Abundance
// x / sqrt(sum(x ^ 2))
//======================================================
static void *norm_chord(void *arg) {
  FOREACH_SAMPLE(
    
    double sq_sum = 0;
    FOREACH_VAL(sq_sum += *val * *val);
    
    sq_sum = sqrt(sq_sum);
    FOREACH_VAL(*val /= sq_sum);
    
  );
  return NULL;
}


//======================================================
// Set each abundance to 1 (present) or 0 (absent)
// x > 0
//======================================================
static void *norm_binary(void *arg) {
  FOREACH_SAMPLE(
    FOREACH_VAL(*val = 1);
  );
  return NULL;
}



void normalize(ecomatrix_t *em, SEXP sexp_norm, int n_threads_) {
  
  int norm  = asInteger(sexp_norm);
  n_threads = n_threads_;
  
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  val_vec   = rw_val_vec(em);
  
  
  // CLR's pseudocounts need special handling
  if (norm == NORM_CLR) {
    
    clr_vec = rw_clr_vec(em);
    
    SEXP sexp_pseudocount = getAttrib(sexp_norm, install("pseudocount"));
    
    // Default to the lowest nonzero value in val_vec
    if (isNull(sexp_pseudocount)) {
      pseudocount = val_vec[0];
      double *vp_end  = val_vec + em->nnz;
      for (double *vp = val_vec + 1; vp < vp_end; vp++)
        if (*vp < pseudocount)
          pseudocount = *vp;
    }
    else {
      pseudocount = asReal(sexp_pseudocount);
    }
  }
  
  
  // function to run
  void * (*norm_func)(void *) = NULL;
  switch (norm) {
    case NORM_PERCENT: norm_func = norm_percent; break;
    case NORM_CLR:     norm_func = norm_clr;     break;
    case NORM_CHORD:   norm_func = norm_chord;   break;
    case NORM_BINARY:  norm_func = norm_binary;  break;
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_samples > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = (pthread_t*) R_alloc(n_threads, sizeof(pthread_t));
      int       *args = (int*)       R_alloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, norm_func, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      return;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  norm_func(&thread_i);
}
