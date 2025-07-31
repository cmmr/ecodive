#include <R.h>
#include <Rinternals.h>
#include <math.h>   // fabs, sqrt
#include <stdlib.h> // calloc, free

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif

#define BRAY_CURTIS 1
#define CANBERRA    2
#define EUCLIDEAN   3
#define JACCARD     4
#define KULCZYNSKI  5
#define MANHATTAN   6

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     weighted;
static int     n_pairs;
static int     n_threads;
static int    *pairs_vec;
static double *dist_vec;
static double *last_sample;


/*
 * The START_PAIR_LOOP and END_PAIR_LOOP macros efficiently 
 * loop through all combinations of samples. Skips unwanted 
 * pairings and pairings not assigned to the current thread. 
 * Ensures all threads process the same number of pairs.
 * 
 * After calling START_PAIR_LOOP the code can expect `x_vec` 
 * and `y_vec` to point to the two samples' columns in 
 * `otu_mtx`. The code should assign to `distance` before 
 * calling END_PAIR_LOOP.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */

#define START_PAIR_LOOP                                        \
  int     thread_i = *((int *) arg);                           \
  double *x_vec    = otu_mtx;                                  \
  double *y_vec    = otu_mtx + n_otus;                         \
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

  
/*
 * Calculations for unweighted beta diversity are derived 
 * from the same three values:
 *   A = sum(x_vec > 0)
 *   B = sum(y_vec > 0)
 *   J = sum(x_vec > 0 & y_vec > 0)
 * CALCULATE_ABJ will make these values available.
 */

#define CALCULATE_ABJ                                          \
  double A = 0, B = 0, J = 0;                                  \
  for (int otu = 0; otu < n_otus; otu++) {                     \
    if (x_vec[otu]) {                                          \
      A++;                                                     \
      if (y_vec[otu]) {                                        \
        B++;                                                   \
        J++;                                                   \
      }                                                        \
    }                                                          \
    else if (y_vec[otu]) {                                     \
      B++;                                                     \
    }                                                          \
  }
  
  
  


//======================================================
// Bray-Curtis
// sum(abs(x-y)) / sum(x+y)
//======================================================
static void *bray_curtis_w(void *arg) {
  START_PAIR_LOOP
  
  double diffs = 0;
  double sums  = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    // abundance of this OTU in the two samples
    double x = x_vec[otu];
    double y = y_vec[otu];
    
    // accumulate
    sums  += x + y;
    diffs += fabs(x - y);
  }
          
  distance = diffs / sums;
  
  END_PAIR_LOOP
}

static void *bray_curtis_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = (A + B - 2*J) / (A + B);
  
  END_PAIR_LOOP
}


//======================================================
// Canberra
// nz = (x+y) > 0; x = x[nz]; y = y[nz]
// sum(abs(x-y) / (x + y)) / sum(nz)
//======================================================
static void *canberra_w(void *arg) {
  START_PAIR_LOOP
  
  int nnz  = 0;
  distance = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    // abundance of this OTU in the two samples
    double x = x_vec[otu];
    double y = y_vec[otu];
    
    // accumulate if appropriate
    if (x > 0 || y > 0) {
      nnz++;
      distance += fabs(x - y) / (x + y);
    }
  }
  
  distance = distance / nnz;
  
  END_PAIR_LOOP
}

static void *canberra_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = (A + B - 2*J) / (A + B - J);
  
  END_PAIR_LOOP
}


//======================================================
// Euclidean
// sqrt(sum((x-y)^2))
//======================================================
static void *euclidean_w(void *arg) {
  START_PAIR_LOOP
  
  distance = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    // abundance of this OTU in the two samples
    double x = x_vec[otu];
    double y = y_vec[otu];
    
    // accumulate
    distance += (x - y) * (x - y);
  }
  
  distance = sqrt(distance);
  
  END_PAIR_LOOP
}

static void *euclidean_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = sqrt(A + B - 2*J);
  
  END_PAIR_LOOP
}


//======================================================
// Jaccard
// 2 * BrayCurtis_W / (1 + BrayCurtis_W)
//======================================================
static void *jaccard_w(void *arg) {
  START_PAIR_LOOP
  
  double diffs = 0;
  double sums  = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    // abundance of this OTU in the two samples
    double x = x_vec[otu];
    double y = y_vec[otu];
    
    // accumulate
    sums  += x + y;
    diffs += fabs(x - y);
  }
  
  distance = diffs / sums;
  distance = 2 * distance / (1 + distance);
  
  END_PAIR_LOOP
}

static void *jaccard_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = (A + B - 2*J) / (A + B);
  distance = 2 * distance / (1 + distance);
  
  END_PAIR_LOOP
}


//======================================================
// Kulczynski
// 1 - (sum(pmin(x,y)) / sum(x) + sum(pmin(x,y)) / sum(y)) / 2
//======================================================
static void *kulczynski_w(void *arg) {
  START_PAIR_LOOP
  
  double x_sum = 0;
  double y_sum = 0;
  double min_sum = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    
    // abundance of this OTU in the two samples
    double x = x_vec[otu];
    double y = y_vec[otu];
    
    x_sum += x;
    y_sum += y;
    min_sum += (x < y) ? x : y;
  }
  
  distance = 1 - (min_sum/x_sum + min_sum/y_sum) / 2;
  
  END_PAIR_LOOP
}

static void *kulczynski_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = 1 - (J/A + J/B) / 2;
  
  END_PAIR_LOOP
}


//======================================================
// Manhattan
// sum(abs(x-y))
//======================================================
static void *manhattan_w(void *arg) {
  START_PAIR_LOOP
  
  distance = 0;
  
  for (int otu = 0; otu < n_otus; otu++) {
    distance += fabs(x_vec[otu] - y_vec[otu]);
  }
  
  END_PAIR_LOOP
}

static void *manhattan_u(void *arg) {
  START_PAIR_LOOP
  CALCULATE_ABJ
  
  distance = A + B - 2 * J;
  
  END_PAIR_LOOP
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_beta_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx,   
    SEXP sexp_weighted,  SEXP sexp_pairs_vec, 
    SEXP sexp_n_threads, SEXP sexp_result_dist ) {
  
  algorithm = asInteger(sexp_algorithm);
  otu_mtx   = REAL(sexp_otu_mtx);
  n_otus    = nrows(sexp_otu_mtx);
  n_samples = ncols(sexp_otu_mtx);
  weighted  = asLogical(sexp_weighted);
  pairs_vec = INTEGER(sexp_pairs_vec);
  n_pairs   = LENGTH(sexp_pairs_vec);
  n_threads = asInteger(sexp_n_threads);
  dist_vec  = REAL(sexp_result_dist);
  
  last_sample = otu_mtx + ((n_samples - 1) * n_otus);
  
  
  // function to run
  void * (*calc_dist_vec)(void *) = NULL;
  
  if (weighted) {
    switch (algorithm) {
    case BRAY_CURTIS: calc_dist_vec = bray_curtis_w; break;
    case CANBERRA:    calc_dist_vec = canberra_w;    break;
    case EUCLIDEAN:   calc_dist_vec = euclidean_w;   break;
    case JACCARD:     calc_dist_vec = jaccard_w;     break;
    case KULCZYNSKI:  calc_dist_vec = kulczynski_w;  break;
    case MANHATTAN:   calc_dist_vec = manhattan_w;   break;
    }
  }
  else {
    switch (algorithm) {
    case BRAY_CURTIS: calc_dist_vec = bray_curtis_u; break;
    case CANBERRA:    calc_dist_vec = canberra_u;    break;
    case EUCLIDEAN:   calc_dist_vec = euclidean_u;   break;
    case JACCARD:     calc_dist_vec = jaccard_u;     break;
    case KULCZYNSKI:  calc_dist_vec = kulczynski_u;  break;
    case MANHATTAN:   calc_dist_vec = manhattan_u;   break;
    }
  }
  
  if (calc_dist_vec == NULL) { // # nocov start
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
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_dist_vec, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
      n_threads = 1;
  int thread_i  = 0;
  calc_dist_vec(&thread_i);
  
  return sexp_result_dist;
}



