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

typedef struct {
  double *otu_vec_1;
  double *otu_vec_2;
  double *distance;
} pair_t;

static int     algorithm;
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     weighted;
static pair_t *pair_vec;
static int     n_pairs;
static int     n_threads;


//======================================================
// Bray-Curtis Weighted
// sum(abs(x-y)) / sum(x+y)
//======================================================
static void *calc_bray_curtis(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
  
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
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
    
    // value to return
    *(pair->distance) = diffs / sums;
  }
  
  return NULL;
}


//======================================================
// Canberra Weighted
// nz = (x+y) > 0; x = x[nz]; y = y[nz]
// sum(abs(x-y) / (x + y)) / sum(nz)
//======================================================
static void *calc_canberra(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
    int    nnz      = 0;
    double distance = 0;
    
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
    
    // value to return
    *(pair->distance) = distance / nnz;
  }
  
  return NULL;
}


//======================================================
// Euclidean Weighted
// sqrt(sum((x-y)^2))
//======================================================
static void *calc_euclidean(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
    double distance = 0;
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = x_vec[otu];
      double y = y_vec[otu];
      
      // accumulate
      distance += (x - y) * (x - y);
    }
    
    // value to return
    *(pair->distance) = sqrt(distance);
  }
  
  return NULL;
}


//======================================================
// Jaccard Weighted
// 2 * BrayCurtis_W / (1 + BrayCurtis_W)
//======================================================
static void *calc_jaccard(void *arg) {
  
  int thread_i = *((int *) arg);
  
  calc_bray_curtis(&thread_i);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    double  bray = *(pair->distance);
    
    // value to return
    *(pair->distance) = 2 * bray / (1 + bray);
  }
  
  return NULL;
}


//======================================================
// Kulczynski Weighted
// 1 - (sum(pmin(x,y)) / sum(x) + sum(pmin(x,y)) / sum(y)) / 2
//======================================================
static void *calc_kulczynski(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
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
    
    // value to return
    *(pair->distance) = 1 - (min_sum/x_sum + min_sum/y_sum) / 2;
  }
  
  return NULL;
}


//======================================================
// Manhattan Weighted
// sum(abs(x-y))
//======================================================
static void *calc_manhattan(void *arg) {
  
  int thread_i = *((int *) arg);
  
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
    double distance = 0;
    
    for (int otu = 0; otu < n_otus; otu++) {
      distance += fabs(x_vec[otu] - y_vec[otu]);
    }
    
    // value to return
    *(pair->distance) = distance;
  }
  
  return NULL;
}



static void *calc_unweighted(void *arg) {
  
  int thread_i = *((int *) arg);
    
  for (int pair_i = thread_i; pair_i < n_pairs; pair_i += n_threads) {
    
    pair_t *pair = pair_vec + pair_i;
    
    // pointers to each sample's column in otu_mtx
    double *x_vec = pair->otu_vec_1;
    double *y_vec = pair->otu_vec_2;
    
    // accumulate if appropriate
    double A = 0; // sum(x>0)
    double B = 0; // sum(y>0)
    double J = 0; // sum(x>0 & y>0)
    
    for (int otu = 0; otu < n_otus; otu++) {
      if (x_vec[otu]) {
        A++;
        if (y_vec[otu]) {
          B++;
          J++;
        }
      }
      else if (y_vec[otu]) {
        B++;
      }
    }
    
    // where to save the result
    double *d = pair->distance;
    
    // shortcut equations from vegan::vegdist()
    switch (algorithm) {
      case BRAY_CURTIS: *d = (A + B - 2*J) / (A + B);     break;
      case CANBERRA:    *d = (A + B - 2*J) / (A + B - J); break;
      case EUCLIDEAN:   *d = sqrt(A + B - 2*J);           break;
      case KULCZYNSKI:  *d = 1 - (J/A + J/B) / 2;         break;
      case MANHATTAN:   *d = A + B - 2 * J;               break;
      case JACCARD:     *d = (A + B - 2*J) / (A + B); *d = 2 * *d / (1 + *d);
    }
  }
  
  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_beta_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx,   
    SEXP sexp_weighted,  SEXP sexp_pair_idx_vec, 
    SEXP sexp_n_threads, SEXP sexp_result_dist ) {
  
  algorithm         = asInteger(sexp_algorithm);
  otu_mtx           = REAL(sexp_otu_mtx);
  n_otus            = nrows(sexp_otu_mtx);
  n_samples         = ncols(sexp_otu_mtx);
  weighted          = asLogical(sexp_weighted);
  int *pair_idx_vec = INTEGER(sexp_pair_idx_vec);
  n_pairs           = LENGTH(sexp_pair_idx_vec);
  n_threads         = asInteger(sexp_n_threads);
  double *dist_vec  = REAL(sexp_result_dist);
  
  
  // function to run
  void * (*calc_bdiv)(void *) = NULL;
  
  if (weighted) {
    switch (algorithm) {
      case BRAY_CURTIS: calc_bdiv = calc_bray_curtis; break;
      case CANBERRA:    calc_bdiv = calc_canberra;    break;
      case EUCLIDEAN:   calc_bdiv = calc_euclidean;   break;
      case JACCARD:     calc_bdiv = calc_jaccard;     break;
      case KULCZYNSKI:  calc_bdiv = calc_kulczynski;  break;
      case MANHATTAN:   calc_bdiv = calc_manhattan;   break;
      default: error("Invalid bdiv metric."); return R_NilValue; // # nocov
    }
  }
  else {
    calc_bdiv = calc_unweighted;
  }
  
  
  // pointers to input and output for 
  // each pairwise comparison
  pair_vec = calloc(n_pairs, sizeof(pair_t));
  int pair_idx = 0;
  int dist_idx = 0;
  for (int i = 0; i < n_samples - 1; i++) {
    for (int j = i + 1; j < n_samples; j++) {
      if (pair_idx_vec[pair_idx] == dist_idx) {
        pair_t *pair    = pair_vec + pair_idx;
        pair->otu_vec_1 = otu_mtx  + (i * n_otus);
        pair->otu_vec_2 = otu_mtx  + (j * n_otus);
        pair->distance  = dist_vec + dist_idx;
        pair_idx++;
        
        if (pair_idx == n_pairs) goto end_loops;
      }
      dist_idx++;
    }
  }
  end_loops:
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_bdiv, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free(tids); free(args); free(pair_vec);
      
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  calc_bdiv(&thread_i);
  
  free(pair_vec);
  
  return sexp_result_dist;
}



