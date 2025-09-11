// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <math.h>   // log, lgamma, pow, round, sqrt
#include <stdlib.h> // calloc, free
#include <string.h> // memset

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define CHAO1        1
#define INV_SIMPSON  2
#define SHANNON      3
#define SIMPSON      4
#define BRILLOUIN    5
#define MCINTOSH     6
#define BERGER       7
#define MARGALEF     8
#define MENHINICK    9
#define SQUARES     10
#define ACE         11
#define FISHER      12

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_threads;
static SEXP   *sexp_extra;
static double *result_vec;

#define ace_rare_cutoff 10


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
// Chao1 (Chao 1984)
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
// Shannon (Shannon 1948)
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
// Simpson (Simpson 1949)
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
// Brillouin (Brillouin 1956)
// (log(sum(x)!) - sum(log(x!))) / sum(x)
// note: lgamma(x + 1) == log(x!)
//======================================================
static void *calc_brillouin(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth  += x;
    result += lgamma(x + 1);
  }
  result = (lgamma(depth + 1) - result) / depth;
  
  END_SAMPLE_LOOP
}


//======================================================
// McIntosh (McIntosh 1967)
// (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
//======================================================
static void *calc_mcintosh(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth  += x;
    result += x * x;
  }
  result = (depth - sqrt(result)) / (depth - sqrt(depth));
  
  END_SAMPLE_LOOP
}


//======================================================
// Berger-Parker (Berger & Parker 1970)
// max(x) / sum(x)
//======================================================
static void *calc_berger(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth += x;
    if (x > result) result = x;
  }
  result /= depth;
  
  END_SAMPLE_LOOP
}


//======================================================
// Margalef (Margalef 1958)
// (sum(x > 0) - 1) / log(sum(x))
//======================================================
static void *calc_margalef(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth += x;
    if (x > 0) result++;
  }
  result = (result - 1) / log(depth);
  
  END_SAMPLE_LOOP
}


//======================================================
// Menhinick (Menhinick 1964)
// sum(x > 0) / sqrt(sum(x))
//======================================================
static void *calc_menhinick(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth += x;
    if (x > 0) result++;
  }
  result /= sqrt(depth);
  
  END_SAMPLE_LOOP
}


//======================================================
// Squares Estimator (Alroy 2018)
// N = sum(x)      # sampling depth
// S = sum(x > 0)  # number of non-zero OTUs
// F1 = sum(x == 1) # singletons
// ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S)) + S
//======================================================
static void *calc_squares(void *arg) {
  START_SAMPLE_LOOP
  
  double depth      = 0;
  double singletons = 0;
  double nz_otus    = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    depth += x;
    if (x > 0) {
      nz_otus++;
      result += x * x;
      if (x == 1) singletons++;
    }
  }
  
  result *= (singletons * singletons);
  result /= (depth * depth) - singletons * nz_otus;
  result += nz_otus;
  
  END_SAMPLE_LOOP
}


//======================================================
// Abundance-based Coverage Estimator (Chao & Lee 1992)
//======================================================
static void *calc_ace(void *arg) {
  
  int cutoff         = asInteger(*sexp_extra) + 1;
  double *rare_nnz_k = calloc(cutoff, sizeof(double));
  
  START_SAMPLE_LOOP
  
  double abund_nnz = 0;
  double rare_sum  = 0;
  double rare_nnz  = 0;
  
  memset(rare_nnz_k, 0, cutoff * sizeof(double));
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      if (x < cutoff) {
        int x_int = (int)(ceil(x));
        rare_sum += x_int;
        rare_nnz++;
        rare_nnz_k[x_int]++;
      }
      else {
        abund_nnz++;
      }
    }
  }
  
  for (int k = 1; k < cutoff; k++)
    result += k * (k - 1) * rare_nnz_k[k];
  
  double p = (1 - (rare_nnz_k[1] / rare_sum));
  
  result = result * rare_nnz/(p * rare_sum * (rare_sum - 1)) - 1;
  if (result < 0) result = 0;
  result = abund_nnz + rare_nnz/p + result * rare_nnz_k[1]/p;
  
  END_SAMPLE_LOOP
  
  free(rare_nnz_k);
}


//======================================================
// Fisher's diversity index (Fisher 1943)
// otus = fisher * log(1 + depth/fisher)
//======================================================
static void *calc_fisher(void *arg) {
  
  int digits  = asInteger(*sexp_extra);
  double mult = pow(10, digits);
  
  START_SAMPLE_LOOP
  
  double nnz_otus = 0;
  double depth    = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    double x = otu_vec[otu];
    if (x > 0) {
      nnz_otus++;
      depth += x;
    }
  }
  
  if (nnz_otus > 0) {
    
    double alpha, lo = 2, hi = 16;
    
    // Sometimes the result will be less than 2 or greater than 16.
    while (lo * log(1 + depth/lo) > nnz_otus) { hi = lo; lo /= 2; }
    while (hi * log(1 + depth/hi) < nnz_otus) { lo = hi; hi *= 2; }
    
    // Check if range has converged to same value after rounding.
    while (round(hi * mult) != round(lo * mult)) {
      
      // This loop's guess for the alpha term.
      alpha = (lo + hi) / 2;
      
      // Update the range we need to examine.
      if (alpha * log(1 + depth/alpha) > nnz_otus) { hi = alpha; }
      else                                         { lo = alpha; }
      
    }
    
    result = round(hi * mult) / mult;
  }
  
  
  END_SAMPLE_LOOP
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx, 
    SEXP sexp_n_threads, SEXP sexp_result_vec,
    SEXP sexp_extra_args ) {
  
  algorithm  = asInteger(sexp_algorithm);
  otu_mtx    = REAL(sexp_otu_mtx);
  n_otus     = nrows(sexp_otu_mtx);
  n_samples  = ncols(sexp_otu_mtx);
  n_threads  = asInteger(sexp_n_threads);
  result_vec = REAL(sexp_result_vec);
  sexp_extra = &sexp_extra_args;
  
  
  // function to run
  void * (*calc_adiv)(void *) = NULL;
  switch (algorithm) {
    case CHAO1:       calc_adiv = calc_chao1;       break;
    case INV_SIMPSON: calc_adiv = calc_inv_simpson; break;
    case SHANNON:     calc_adiv = calc_shannon;     break;
    case SIMPSON:     calc_adiv = calc_simpson;     break;
    case BRILLOUIN:   calc_adiv = calc_brillouin;   break;
    case MCINTOSH:    calc_adiv = calc_mcintosh;    break;
    case BERGER:      calc_adiv = calc_berger;      break;
    case MARGALEF:    calc_adiv = calc_margalef;    break;
    case MENHINICK:   calc_adiv = calc_menhinick;   break;
    case SQUARES:     calc_adiv = calc_squares;     break;
    case ACE:         calc_adiv = calc_ace;         break;
    case FISHER:      calc_adiv = calc_fisher;      break;
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
