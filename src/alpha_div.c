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

#define ADIV_ACE          1
#define ADIV_BERGER       2
#define ADIV_BRILLOUIN    3
#define ADIV_CHAO1        4
#define ADIV_FISHER       5
#define ADIV_INV_SIMPSON  6
#define ADIV_MARGALEF     7
#define ADIV_MCINTOSH     8
#define ADIV_MENHINICK    9
#define ADIV_OBSERVED    10
#define ADIV_SHANNON     11
#define ADIV_SIMPSON     12
#define ADIV_SQUARES     13

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_threads;
static SEXP   *sexp_extra;
static double *result_vec;
static int     mtx_len;


/*
 * START_SAMPLE_LOOP and END_SAMPLE_LOOP define a simple loop 
 * for iterating over all samples, ensuring that each is 
 * assigned to only a single thread. Provides `sample`, and 
 * `result` initialized to 0; expects `result` at the end.
 */

#define START_SAMPLE_LOOP                                                 \
  int thread_i = *((int *) arg);                                          \
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {  \
    double result = 0;

#define END_SAMPLE_LOOP                                        \
    result_vec[sample] = result;                               \
  }


//======================================================
// Abundance-based Coverage Estimator (Chao & Lee 1992)
//======================================================
static void *ace(void *arg) {
  
  int cutoff         = asInteger(*sexp_extra) + 1;
  double *rare_nnz_k = calloc(cutoff, sizeof(double));
  
  START_SAMPLE_LOOP
  
  double abund_nnz = 0;
  double rare_sum  = 0;
  double rare_nnz  = 0;
  
  memset(rare_nnz_k, 0, cutoff * sizeof(double));
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
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
  return NULL;
}
  
  
//======================================================
// Berger-Parker (Berger & Parker 1970)
// max(x) / sum(x)
//======================================================
static void *berger(void *arg) {
  START_SAMPLE_LOOP
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x > result) result = x;
  }
  
  END_SAMPLE_LOOP
  return NULL;
}
  
  
//======================================================
// Brillouin (Brillouin 1956)
// (log(sum(x)!) - sum(log(x!))) / sum(x)
// note: lgamma(x + 1) == log(x!)
//======================================================
static void *brillouin(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      depth  += x;
      result += lgamma(x + 1);
    }
  }
  
  result = (lgamma(depth + 1) - result) / depth;
  
  END_SAMPLE_LOOP
  return NULL;
}
  


//======================================================
// Chao1 (Chao 1984)
// sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))
//======================================================
static void *chao1(void *arg) {
  START_SAMPLE_LOOP
    
  double nnz = 0, ones = 0, twos = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      nnz++;
      if      (x <= 1) ones++;
      else if (x <= 2) twos++;
    }
  }
  
  result = nnz + ((ones * ones) / (2 * twos));
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Fisher's diversity index (Fisher 1943)
// otus = fisher * log(1 + depth/fisher)
//======================================================
static void *fisher(void *arg) {
  
  int    digits = asInteger(*sexp_extra);
  double mult   = pow(10, digits);
  
  START_SAMPLE_LOOP
  
  double nnz = 0, depth = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      nnz++;
      depth += x;
    }
  }
  
  if (nnz) {
    
    double alpha, lo = 2, hi = 16;
    
    // Sometimes the result will be less than 2 or greater than 16.
    while (lo * log(1 + depth/lo) > nnz) { hi = lo; lo /= 2; }
    while (hi * log(1 + depth/hi) < nnz) { lo = hi; hi *= 2; }
    
    // Check if range has converged to same value after rounding.
    while (round(hi * mult) != round(lo * mult)) {
      
      // This loop's guess for the alpha term.
      alpha = (lo + hi) / 2;
      
      // Update the range we need to examine.
      if (alpha * log(1 + depth/alpha) > nnz) { hi = alpha; }
      else                                    { lo = alpha; }
      
    }
    
    result = round(hi * mult) / mult;
  }
  
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Inverse Simpson
// p <- x / sum(x)
// 1 / sum(p ** 2)
//======================================================
static void *inv_simpson(void *arg) {
  START_SAMPLE_LOOP
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) result += x * x;
  }
  
  if (result) result = 1 / result;
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Margalef (Margalef 1958)
// (sum(x > 0) - 1) / log(sum(x))
//======================================================
static void *margalef(void *arg) {
  START_SAMPLE_LOOP
  
  double nnz = 0, depth = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      nnz++;
      depth += x;
    }
  }
  
  result = (nnz - 1) / log(depth);
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// McIntosh (McIntosh 1967)
// (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
//======================================================
static void *mcintosh(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      depth  += x;
      result += x * x;
    }
  }
  
  result = (depth - sqrt(result)) / (depth - sqrt(depth));
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Menhinick (Menhinick 1964)
// sum(x > 0) / sqrt(sum(x))
//======================================================
static void *menhinick(void *arg) {
  START_SAMPLE_LOOP
  
  double depth = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      depth += x;
      result++;
    }
  }
  
  result /= sqrt(depth);
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Observed Features
// sum(x > 0)
//======================================================
static void *observed(void *arg) {
  START_SAMPLE_LOOP
  
  for (int i = sample; i < mtx_len; i += n_samples)
    if (otu_mtx[i])
      result++;
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Shannon (Shannon 1948)
// p <- x / sum(x)
// -sum(p * log(p))
//======================================================
static void *shannon(void *arg) {
  START_SAMPLE_LOOP
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) result += x * log(x);
  }
  
  result *= -1;
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Simpson (Simpson 1949)
// p <- x / sum(x)
// 1 - sum(p ** 2)
//======================================================
static void *simpson(void *arg) {
  START_SAMPLE_LOOP
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) result += x * x;
  }
  
  result = 1 - result;
  
  END_SAMPLE_LOOP
  return NULL;
}


//======================================================
// Squares Estimator (Alroy 2018)
// N = sum(x)      # sampling depth
// S = sum(x > 0)  # number of non-zero OTUs
// F1 = sum(x == 1) # singletons
// ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S)) + S
//======================================================
static void *squares(void *arg) {
  START_SAMPLE_LOOP
  
  double nnz = 0, depth = 0, singletons = 0;
  
  for (int i = sample; i < mtx_len; i += n_samples) {
    double x = otu_mtx[i];
    if (x) {
      nnz++;
      depth  += x;
      result += x * x;
      if (x == 1) singletons++;
    }
  }
  
  result *= (singletons * singletons);
  result /= (depth * depth) - singletons * nnz;
  result += nnz;
  
  END_SAMPLE_LOOP
  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx, 
    SEXP sexp_n_threads, SEXP sexp_extra_args ) {
  
  algorithm  = asInteger(sexp_algorithm);
  otu_mtx    = REAL(sexp_otu_mtx);
  n_otus     = ncols(sexp_otu_mtx);
  n_samples  = nrows(sexp_otu_mtx);
  n_threads  = asInteger(sexp_n_threads);
  sexp_extra = &sexp_extra_args;
  
  mtx_len = n_otus * n_samples;
  
  
  // function to run
  void * (*adiv_func)(void *) = NULL;
  switch (algorithm) {
    case ADIV_ACE:         adiv_func = ace;         break;
    case ADIV_BERGER:      adiv_func = berger;      break;
    case ADIV_BRILLOUIN:   adiv_func = brillouin;   break;
    case ADIV_CHAO1:       adiv_func = chao1;       break;
    case ADIV_FISHER:      adiv_func = fisher;      break;
    case ADIV_INV_SIMPSON: adiv_func = inv_simpson; break;
    case ADIV_MARGALEF:    adiv_func = margalef;    break;
    case ADIV_MCINTOSH:    adiv_func = mcintosh;    break;
    case ADIV_MENHINICK:   adiv_func = menhinick;   break;
    case ADIV_OBSERVED:    adiv_func = observed;    break;
    case ADIV_SHANNON:     adiv_func = shannon;     break;
    case ADIV_SIMPSON:     adiv_func = simpson;     break;
    case ADIV_SQUARES:     adiv_func = squares;     break;
  }
  
  if (adiv_func == NULL) { // # nocov start
    error("Invalid alpha diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Create the diversity vector to return
  SEXP sexp_result_vec = PROTECT(allocVector(REALSXP, n_samples));
  result_vec           = REAL(sexp_result_vec);
  SEXP sexp_mtx_dimnames = getAttrib(sexp_otu_mtx, R_DimNamesSymbol);
  if (sexp_mtx_dimnames != R_NilValue) {
    SEXP sexp_mtx_rownames = VECTOR_ELT(sexp_mtx_dimnames, 0);
    if (sexp_mtx_rownames != R_NilValue) {
      setAttrib(sexp_result_vec, R_NamesSymbol, sexp_mtx_rownames);
    }
  }
  
  
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
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, adiv_func, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args);
      
      UNPROTECT(1);
      return sexp_result_vec;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  adiv_func(&thread_i);
  
  UNPROTECT(1);
  return sexp_result_vec;
}
