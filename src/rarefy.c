#include <R.h>
#include <Rinternals.h>
#include <inttypes.h> // uint32_t, uint64_t
#include <string.h>   // memcpy

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <stdlib.h> // calloc, free
#    define HAVE_PTHREAD
#  endif
#endif

static int *otu_mtx;
static int  n_otus;
static int  n_samples;
static int  target;
static int  seed;
static int  n_threads;
static int *result_mtx;



// Begin vendored PCG code.
// https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.c

/*
 * PCG Random Number Generation for C.
 * Copyright 2014 Melissa O'Neill
 * Apache License 2.0
 * http://www.pcg-random.org
 */

typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_random_t;

static inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
  uint64_t oldstate   = rng->state;
  rng->state          = oldstate * 6364136223846793005ULL + rng->inc;
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot        = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static inline void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq) {
  rng->state = 0U;
  rng->inc   = (initseq << 1u) | 1u;
  pcg32_random_r(rng);
  rng->state += initstate;
  pcg32_random_r(rng);
}

// End vendored PCG code.





//======================================================
// Limits itself to samples based on modulo.
//======================================================
static void *rarefy_worker(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    int *otu_vec = otu_mtx + (sample * n_otus);
    
    // Sum the counts for this sample
    int  depth = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      depth += otu_vec[otu];
    }
    
    // Insufficient sequences - leave as all zeroes.
    if (depth < target) continue;
    
    // Already rarefied - copy input to output
    int *result_vec = result_mtx + (sample * n_otus);
    if (depth == target) {
      memcpy(result_vec, otu_vec, n_otus * sizeof(int));
      continue;
    }
    
    // Seed the PRNG for this thread.
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, seed, sample);
    
    // Knuth algorithm for choosing target seqs from depth.
    int tried = 0, kept = 0;
    for (int otu = 0; otu < n_otus && kept < target; otu++) {
      
      int n_seqs = otu_vec[otu];
      for (int seq = 0; seq < n_seqs && kept < target; seq++) {
        
        uint32_t not_tried  = depth - tried;
        uint32_t still_need = target - kept;
        uint32_t rand_int   = pcg32_random_r(&rng);
        
        if (rand_int % not_tried < still_need) {
          result_vec[otu]++; // retain this observation
          kept++;
        }
        
        tried++;
      }
    }
    
  }
  
  return NULL;
}



//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_rarefy(
    SEXP sexp_otu_mtx,   SEXP sexp_target, 
    SEXP sexp_seed,      SEXP sexp_n_threads, 
    SEXP sexp_result_mtx ) {
  
  otu_mtx    = INTEGER(sexp_otu_mtx);
  n_otus     = nrows(sexp_otu_mtx);
  n_samples  = ncols(sexp_otu_mtx);
  target     = asInteger(sexp_target);
  seed       = asInteger(sexp_seed);
  n_threads  = asInteger(sexp_n_threads);
  result_mtx = INTEGER(sexp_result_mtx);
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_samples > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, rarefy_worker, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free(tids); free(args);
      
      return sexp_result_mtx;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  rarefy_worker(&thread_i);
  
  return sexp_result_mtx;
}

