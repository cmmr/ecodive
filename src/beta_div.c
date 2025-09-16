// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

// https://pdodds.w3.uvm.edu/research/papers/others/everything/cha2007a.pdf

#include <R.h>
#include <Rinternals.h>
#include <math.h>   // fabs, log, pow, sqrt
#include <stdlib.h> // calloc, free

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define BDIV_BHATTACHARYYA  1
#define BDIV_BRAY           2
#define BDIV_CANBERRA       3
#define BDIV_CHEBYSHEV      4
#define BDIV_CLARK          6
#define BDIV_DIVERGENCE     7
#define BDIV_EUCLIDEAN      8
#define BDIV_GOWER          9
#define BDIV_HAMMING       10
#define BDIV_HORN          11
#define BDIV_JACCARD       12
#define BDIV_JSD           13
#define BDIV_LORENTZIAN    14
#define BDIV_MANHATTAN     15
#define BDIV_MINKOWSKI     16
#define BDIV_MORISITA      17
#define BDIV_MOTYKA        18
#define BDIV_OCHIAI        19
#define BDIV_SOERGEL       20
#define BDIV_SORENSEN      21
#define BDIV_SQUARED_CHISQ 22
#define BDIV_SQUARED_CHORD 23
#define BDIV_WAVE_HEDGES   24

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_pairs;
static int     n_threads;
static int    *pairs_vec;
static char    all_pairs;
static double *dist_vec;
static SEXP   *sexp_extra;


/*
 * The START_PAIR_LOOP and END_PAIR_LOOP macros efficiently 
 * loop through all combinations of samples. Skips unwanted 
 * pairings and pairings not assigned to the current thread. 
 * Ensures all threads process the same number of pairs.
 * 
 * After calling START_PAIR_LOOP the code can expect `x1` 
 * and `y1` to be the two samples' first indices in 
 * `otu_mtx`. The code should assign to `distance` before 
 * calling END_PAIR_LOOP.
 * 
 * After calling START_OTU_LOOP the code can expect `x` 
 * and `y` to be the current OTU's values from `otu_mtx`.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */

#define START_PAIR_LOOP                                        \
int thread_i = *((int *) arg);                                 \
int pair_idx = 0; /* The pairs that we're asked to compute */  \
int dist_idx = 0; /* All combinations; length of return vec */ \
for (int x1 = 0; x1 < n_samples - 1; x1++) {                   \
  for (int y1 = x1 + 1; y1 < n_samples; y1++) {                \
    if (all_pairs || pairs_vec[pair_idx] == dist_idx) {        \
      if (pair_idx % n_threads == thread_i) {                  \
        double distance = 0;                                   \

        
#define START_OTU_LOOP                                         \
        for (int otu = 0; otu < n_otus; otu++) {               \
          double x = otu_mtx[x1 + otu * n_samples];            \
          double y = otu_mtx[y1 + otu * n_samples];            \


#define END_OTU_LOOP                                           \
        }                                                      \


#define END_PAIR_LOOP                                          \
        dist_vec[dist_idx] = distance;                         \
      }                                                        \
      pair_idx++;                                              \
      if (pair_idx == n_pairs) goto end_loops;                 \
    }                                                          \
    dist_idx++;                                                \
  }                                                            \
}                                                              \
end_loops:                                                     \





//======================================================
// Bhattacharyya
// -log(sum(sqrt(x * y)))
//======================================================
static void *bhattacharyya(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  distance += sqrt(x * y);
  END_OTU_LOOP
  
  distance = -1 * log(distance);
  
  END_PAIR_LOOP
  return NULL;
}



//======================================================
// Dice-Sorensen; Bray-Curtis
// sum(abs(x-y)) / sum(x+y)
//======================================================
static void *bray(void *arg) {
  START_PAIR_LOOP
  
  double diffs = 0;
  double sums  = 0;
  
  START_OTU_LOOP
  sums  += x + y;
  diffs += (x > y) ? x - y : y - x;
  END_OTU_LOOP
  
  distance = diffs / sums;
  
  END_PAIR_LOOP
  return NULL;
}



//======================================================
// Canberra
// nz = (x+y) > 0; x = x[nz]; y = y[nz]
// sum(abs(x-y) / (x + y)) / sum(nz)
//======================================================
static void *canberra(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) {
    if (x > y) { distance += (x - y) / (x + y); }
    else       { distance += (y - x) / (x + y); }
  }
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Chebyshev
// max(abs(x - y))
//======================================================
static void *chebyshev(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  double d = (x > y) ? x - y : y - x;
  if (d > distance) distance = d;
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Clark
// sqrt(sum((abs(x - y) / (x + y)) ^ 2))
//======================================================
static void *clark(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y)
    distance += ((x - y) / (x + y)) * ((x - y) / (x + y));
  END_OTU_LOOP
  
  distance = sqrt(distance);
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Divergence
// 2 * sum((x-y)^2 / (x+y)^2)
//======================================================
static void *divergence(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y)
    distance += ((x - y) * (x - y)) / ((x + y) * (x + y));
  END_OTU_LOOP
  
  distance = 2 * distance;
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Euclidean
// sqrt(sum((x-y)^2))
//======================================================
static void *euclidean(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y)
    distance += (x - y) * (x - y);
  END_OTU_LOOP
  
  distance = sqrt(distance);
  
  END_PAIR_LOOP
  return NULL;
}



//======================================================
// Gower
// sum(abs(x-y) / r) / n
//======================================================
static void *gower(void *arg) {
  
  double *range_vec = REAL(*sexp_extra);
  
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if      (x > y) { distance += (x - y) / range_vec[otu]; }
  else if (y > x) { distance += (y - x) / range_vec[otu]; }
  END_OTU_LOOP
  
  distance /= n_otus;
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Hamming
// sum(xor(x, y))
//======================================================
static void *hamming(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (!x ^ !y) distance++;
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Horn
// 
// z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
// 1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))
//======================================================
static void *horn(void *arg) {
  START_PAIR_LOOP
  
  double sum_x = 0, sum_x2 = 0;
  double sum_y = 0, sum_y2 = 0;
  
  START_OTU_LOOP
  if (x || y) {
    distance += x * y;
    sum_x += x; sum_x2 += x * x;
    sum_y += y; sum_y2 += y * y;
  }
  END_OTU_LOOP
  
  sum_x2 /= sum_x * sum_x;
  sum_y2 /= sum_y * sum_y;
  
  distance = 1 - (2 * distance) / ((sum_x2 + sum_y2) * sum_x * sum_y);
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Jaccard
// sum(xor(x, y)) / sum(x | y)
//======================================================
static void *jaccard(void *arg) {
  START_PAIR_LOOP
  
  double D = 0, U = 0;
  
  START_OTU_LOOP
  if      (x) { U++; if (!y) D++; }
  else if (y) { U++; D++; }
  END_OTU_LOOP
  
  distance = D / U;
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Jensen-Shannon Divergence (JSD)
// sum(x * log(2*x / (x+y)), y * log(2*y / (x+y))) / 2
//======================================================
static void *jsd(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x) distance += x * log(2 * x / (x + y));
  if (y) distance += y * log(2 * y / (x + y));
  END_OTU_LOOP
  
  distance /= 2;
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Lorentzian
// sum(log(1 + abs(x - y)))
//======================================================
static void *lorentzian(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) distance += log(1 + fabs(x - y));
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Manhattan
// sum(abs(x-y))
//======================================================
static void *manhattan(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) {
    if (x > y) { distance += x - y; }
    else       { distance += y - x; }
  }
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Minkowski
// sum(abs(x - y)^p) ^ (1/p)
//======================================================
static void *minkowski(void *arg) {
  
  double power     = asReal(*sexp_extra);
  double inv_power = 1 / power;
  
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) {
    if (x > y) { distance += pow(x - y, power); }
    else       { distance += pow(y - x, power); }
  }
  END_OTU_LOOP
  
  distance = pow(distance, inv_power);
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Morisita
// 
// simpson_x <- sum(x * (x - 1)) / (sum(x) * (sum(x) - 1))
// simpson_y <- sum(y * (y - 1)) / (sum(y) * (sum(y) - 1))
// 1 - ((2 * sum(x * y)) / ((simpson_x + simpson_y) * sum(x) * sum(y)))
//======================================================
static void *morisita(void *arg) {
  START_PAIR_LOOP
  
  double simpson_x = 0, sum_x = 0;
  double simpson_y = 0, sum_y = 0;
  
  START_OTU_LOOP
  if (x || y) {
    distance += x * y;
    sum_x += x; simpson_x += x * (x - 1);
    sum_y += y; simpson_y += y * (y - 1);
  }
  END_OTU_LOOP
  
  simpson_x /= sum_x * (sum_x - 1);
  simpson_y /= sum_y * (sum_y - 1);
  
  distance = 1 - (2 * distance) / ((simpson_x + simpson_y) * sum_x * sum_y);
  
  END_PAIR_LOOP
  return NULL;
}




//======================================================
// Motyka
// sum(pmax(x, y)) / sum(x, y)
//======================================================
static void *motyka(void *arg) {
  START_PAIR_LOOP
  
  double sums = 0;
  
  START_OTU_LOOP
  distance += (x > y) ? x : y;
  sums     += x + y;
  END_OTU_LOOP
  
  distance /= sums;
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Dice-Sorensen
// 2 * sum(x & y) / sum(x>0, y>0)
//======================================================
static void *sorensen(void *arg) {
  START_PAIR_LOOP
  
  double top = 0, bot = 0;
  
  START_OTU_LOOP
  if (x || y) {
    bot++;
    if (x && y) { top++; bot++; }
  }
  END_OTU_LOOP
  
  distance = 1 - 2 * top / bot;
    
  END_PAIR_LOOP
  return NULL;
}



//======================================================
// Ochiai
// sum((x & y)) / sqrt(sum(x > 0) * sum(y > 0))
//======================================================
static void *ochiai(void *arg) {
  START_PAIR_LOOP
  
  double A = 0, B = 0, J = 0;
  
  START_OTU_LOOP
  if      (x) { A++; if (y) { B++; J++; } } 
  else if (y) { B++; }
  END_OTU_LOOP
  
  distance = 1 - J / sqrt(A * B);
  
  END_PAIR_LOOP
  return NULL;
}
  
  
//======================================================
// Soergel
// 1 - sum(pmin(x, y)) / sum(pmax(x, y))
//======================================================
static void *soergel(void *arg) {
  START_PAIR_LOOP
  
  double min_sum = 0;
  double max_sum = 0;
  
  START_OTU_LOOP
  if (x < y) { min_sum += x; max_sum += y; } 
  else       { min_sum += y; max_sum += x; }
  END_OTU_LOOP
  
  distance = 1 - (min_sum / max_sum);
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Squared Ch-Squared
// sum((x - y) ^ 2 / (x + y))
//======================================================
static void *squared_chisq(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) distance += (x - y) * (x - y) / (x + y);
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Squared Chord
// sum((sqrt(x) - sqrt(y)) ^ 2)
//======================================================
static void *squared_chord(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) {
    double d = sqrt(x) - sqrt(y);
    distance += d * d;
  }
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}


//======================================================
// Wave Hedges
// sum(abs(x - y) / pmax(x, y))
//======================================================
static void *wave_hedges(void *arg) {
  START_PAIR_LOOP
  
  START_OTU_LOOP
  if (x || y) {
    if (x > y) { distance += (x - y) / x; }
    else       { distance += (y - x) / y; }
  }
  END_OTU_LOOP
  
  END_PAIR_LOOP
  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_beta_div(
    SEXP sexp_algorithm,   SEXP sexp_otu_mtx,   
    SEXP sexp_pairs_vec,   SEXP sexp_n_threads, 
    SEXP sexp_result_dist, SEXP sexp_extra_args ) {
  
  algorithm  = asInteger(sexp_algorithm);
  otu_mtx    = REAL(sexp_otu_mtx);
  n_otus     = ncols(sexp_otu_mtx);
  n_samples  = nrows(sexp_otu_mtx);
  n_threads  = asInteger(sexp_n_threads);
  dist_vec   = REAL(sexp_result_dist);
  sexp_extra = &sexp_extra_args;
  
  // Avoid allocating pairs_vec for common all-vs-all case
  if (isNull(sexp_pairs_vec)) {
    all_pairs = 1;
    pairs_vec = NULL;
    n_pairs   = n_samples * (n_samples - 1) / 2;
  }
  else {
    all_pairs = 0;
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
  }
  
  
  // function to run
  void * (*bdiv_func)(void *) = NULL;
  
  switch (algorithm) {
    case BDIV_BHATTACHARYYA: bdiv_func = bhattacharyya; break;
    case BDIV_BRAY:          bdiv_func = bray;          break;
    case BDIV_CANBERRA:      bdiv_func = canberra;      break;
    case BDIV_CHEBYSHEV:     bdiv_func = chebyshev;     break;
    case BDIV_CLARK:         bdiv_func = clark;         break;
    case BDIV_DIVERGENCE:    bdiv_func = divergence;    break;
    case BDIV_EUCLIDEAN:     bdiv_func = euclidean;     break;
    case BDIV_GOWER:         bdiv_func = gower;         break;
    case BDIV_HAMMING:       bdiv_func = hamming;       break;
    case BDIV_HORN:          bdiv_func = horn;          break;
    case BDIV_JACCARD:       bdiv_func = jaccard;       break;
    case BDIV_JSD:           bdiv_func = jsd;           break;
    case BDIV_LORENTZIAN:    bdiv_func = lorentzian;    break;
    case BDIV_MANHATTAN:     bdiv_func = manhattan;     break;
    case BDIV_MINKOWSKI:     bdiv_func = minkowski;     break;
    case BDIV_MORISITA:      bdiv_func = morisita;      break;
    case BDIV_MOTYKA:        bdiv_func = motyka;        break;
    case BDIV_OCHIAI:        bdiv_func = ochiai;        break;
    case BDIV_SOERGEL:       bdiv_func = soergel;       break;
    case BDIV_SORENSEN:      bdiv_func = sorensen;      break;
    case BDIV_SQUARED_CHISQ: bdiv_func = squared_chisq; break;
    case BDIV_SQUARED_CHORD: bdiv_func = squared_chord; break;
    case BDIV_WAVE_HEDGES:   bdiv_func = wave_hedges;   break;
  }
  
  if (bdiv_func == NULL) { // # nocov start
    error("Invalid beta diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      if (tids == NULL || args == NULL) { // # nocov start
        free(tids); free(args);
        error("Insufficient memory for parallel beta diversity calculation.");
        return R_NilValue;
      } // # nocov end
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, bdiv_func, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
      n_threads = 1;
  int thread_i  = 0;
  bdiv_func(&thread_i);
  
  return sexp_result_dist;
}



