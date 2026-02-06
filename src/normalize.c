// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


#define NORM_NONE    0
#define NORM_PERCENT 1
#define NORM_CLR     2
#define NORM_CHORD   3
#define NORM_BINARY  4

static int     n_samples;
static int     n_otus;
static int    *pos_vec;
static double *val_vec;
static double *clr_vec;
static double  pseudocount;
static int     is_percent_normalized;


// Macro to loop over each sample based on current threading setup.
// Sets sam (sample index), val_begin, and val_end (pointers to val_vec 
// for this sample's first value and the next sample's first value).
#define FOREACH_SAMPLE(expression)                                  \
  do {                                                              \
    int thread_i  = ((worker_t *)arg)->i;                           \
    int n_threads = ((worker_t *)arg)->n;                           \
    for (int sam = thread_i; sam < n_samples; sam += n_threads) {   \
      double *val_begin = val_vec + pos_vec[sam];                   \
      double *val_end   = val_vec + pos_vec[sam + 1];               \
      expression                                                    \
    }                                                               \
  } while (0)

// Marco to loop over each value for the current sample.
// Sets val to point to each value in val_vec in turn.
#define FOREACH_VAL(expression)                                     \
  do {                                                              \
    for (double *val = val_begin; val != val_end; val++) {          \
      expression;                                                   \
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

static void *check_percent_normalized(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    FOREACH_VAL(
      if (*val < 0 || *val > 1) {
        is_percent_normalized = 0;
        return NULL;
      }
      depth += *val
    );
    if (fabs(depth - 1) > 1e-6) {
      is_percent_normalized = 0;
      return NULL;
    }
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



void normalize(ecomatrix_t *em, int norm, int n_threads, int pseudocount_) {
  
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  val_vec   = em->val_vec;
  
  if (norm == NORM_PERCENT) {
    // Check if it's already normalized to percent. If so, skip.
    if (*val_vec <= 1) {
      is_percent_normalized = 1;
      run_parallel(check_percent_normalized, n_threads, n_samples);
      if (is_percent_normalized) return;
    }
  }
  else if (norm == NORM_CLR) {
    clr_vec = rw_clr_vec(em);
  }
  
  // function to run
  void * (*norm_func)(void *) = NULL;
  switch (norm) {
    case NORM_PERCENT: norm_func = norm_percent; break;
    case NORM_CLR:     norm_func = norm_clr;     break;
    case NORM_CHORD:   norm_func = norm_chord;   break;
    case NORM_BINARY:  norm_func = norm_binary;  break;
  }
  
  pseudocount = pseudocount_;
  val_vec     = rw_val_vec(em);
  
  run_parallel(norm_func, n_threads, n_samples);
}
