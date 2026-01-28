// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

// https://pdodds.w3.uvm.edu/research/papers/others/everything/cha2007a.pdf

#include "ecodive.h"


typedef void *(*pthread_func_t)(void *);

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

static int     n_samples;
static int     n_otus;
static int     n_dist;
static int     n_pairs;
static int     n_threads;
static int    *pos_vec;
static int    *otu_vec;
static double *val_vec;
static double *clr_vec;
static int    *pairs_vec;
static double *dist_vec;
static SEXP   *sexp_extra;


/*
 * The FOREACH_PAIR macro iterates through all combinations of 
 * samples. Skips unwanted pairings and pairings not assigned to 
 * the current thread. Ensures all threads process the same 
 * number of pairs. The code should assign to `distance`.
 * 
 * The FOREACH_OTU macro iterates through all OTU abundances for 
 * a given pair of samples, assigning the values to `x` and `y`.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */
#define FOREACH_PAIR(expression)                               \
  do {                                                         \
    int thread_i = *((int *) arg);                             \
    int dist_idx = 0;                                          \
                                                               \
    if (pairs_vec == NULL) { /* All vs All */                  \
                                                               \
      for (int sam_i = 0; sam_i < n_samples - 1; sam_i++) {    \
        for (int sam_j = sam_i + 1; sam_j < n_samples; sam_j++) {\
          if (dist_idx % n_threads == thread_i) {              \
                                                               \
            double distance = 0;                               \
                                                               \
            expression;                                        \
                                                               \
            dist_vec[dist_idx] = distance;                     \
          }                                                    \
          dist_idx++;                                          \
        }                                                      \
      }                                                        \
                                                               \
    } else { /* Specific Pairs of Samples */                   \
                                                               \
      int pair_idx = thread_i;                                 \
      for (; pair_idx < n_pairs; pair_idx += n_threads) {      \
                                                               \
        dist_idx = pairs_vec[pair_idx]; /* 1-based */          \
                                                               \
        int sam_i          = 0;                                \
        int sam_j          = dist_idx;                         \
        int pairs_in_block = n_samples - 1;                    \
                                                               \
        while (sam_j > pairs_in_block) {                       \
          sam_i++;                                             \
          sam_j -= pairs_in_block;                             \
          pairs_in_block--;                                    \
        }                                                      \
                                                               \
        sam_j += sam_i;                                        \
                                                               \
        double distance = 0;                                   \
                                                               \
        expression;                                            \
                                                               \
        dist_vec[dist_idx - 1] = distance;                     \
      }                                                        \
    }                                                          \
  } while (0)


#define FOREACH_OTU(expression)                                \
  do {                                                         \
    int    *i      = otu_vec + pos_vec[sam_i];                 \
    int    *j      = otu_vec + pos_vec[sam_j];                 \
    int    *i_end  = otu_vec + pos_vec[sam_i + 1];             \
    int    *j_end  = otu_vec + pos_vec[sam_j + 1];             \
    double *val_i  = val_vec + pos_vec[sam_i];                 \
    double *val_j  = val_vec + pos_vec[sam_j];                 \
    double  x_zero = clr_vec ? clr_vec[sam_i] : 0;             \
    double  y_zero = clr_vec ? clr_vec[sam_j] : 0;             \
    int     otu = 0, n_ops = 0;                                \
    double  x, y;                                              \
    while (1) {                                                \
      if (i != i_end && j != j_end) {                          \
        if (*i == *j) {                                        \
          otu = *i;                                            \
          x   = *val_i; i++; val_i++;                          \
          y   = *val_j; j++; val_j++;                          \
        }                                                      \
        else if (*i < *j) {                                    \
          otu = *i;                                            \
          x   = *val_i; i++; val_i++;                          \
          y   = y_zero;                                        \
        }                                                      \
        else {                                                 \
          otu = *j;                                            \
          x   = x_zero;                                        \
          y   = *val_j; j++; val_j++;                          \
        }                                                      \
      }                                                        \
      else if (i != i_end) {                                   \
        otu = *i;                                              \
        x   = *val_i; i++; val_i++;                            \
        y   = y_zero;                                          \
      }                                                        \
      else if (j != j_end) {                                   \
        otu = *j;                                              \
        x   = x_zero;                                          \
        y   = *val_j; j++; val_j++;                            \
      }                                                        \
      else {                                                   \
        break;                                                 \
      }                                                        \
      n_ops++;                                                 \
      expression;                                              \
    }                                                          \
    if (clr_vec) { /* Double zeros */                          \
      x = x_zero, y = y_zero;                                  \
      while (n_ops < n_otus) {                                 \
        n_ops++;                                               \
        expression;                                            \
      }                                                        \
    }                                                          \
    (void)otu;                                                 \
  } while (0)


#define WITH_ABJ(expression)                                   \
  do {                                                         \
    int *i     = otu_vec + pos_vec[sam_i];                     \
    int *j     = otu_vec + pos_vec[sam_j];                     \
    int *i_end = otu_vec + pos_vec[sam_i + 1];                 \
    int *j_end = otu_vec + pos_vec[sam_j + 1];                 \
    double A = 0, B = 0, J = 0;                                \
    while (i != i_end && j != j_end) {                         \
      if      (*i == *j) { A++; B++; J++; i++; j++; }          \
      else if (*i < *j)  { A++; i++; }                         \
      else               { B++; j++; }                         \
    }                                                          \
    A += i_end - i;                                            \
    B += j_end - j;                                            \
                                                               \
    expression;                                                \
                                                               \
  } while (0)








//======================================================
// Bhattacharyya
// -log(sum(sqrt(x * y)))
//======================================================
static void *bhattacharyya(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(distance += sqrt(x * y));
  
    distance = -1 * log(distance);
  );
  
  return NULL;
}



//======================================================
// Dice-Sorensen; Bray-Curtis
// sum(abs(x-y)) / sum(x+y)
//======================================================
static void *bray(void *arg) {
  FOREACH_PAIR(
  
    double diffs = 0;
    double sums  = 0;
    
    FOREACH_OTU(
      sums  += x + y;
      diffs += fabs(x - y);
    );
  
    distance = diffs / sums;
  );
  
  return NULL;
}



//======================================================
// Canberra
// nz = (x+y) > 0; x = x[nz]; y = y[nz]
// sum(abs(x-y) / (x + y)) / sum(nz)
//======================================================
static void *canberra(void *arg) {
  FOREACH_PAIR(
    FOREACH_OTU(distance += fabs(x - y) / (x + y));
  );
  
  return NULL;
}


//======================================================
// Chebyshev
// max(abs(x - y))
//======================================================
static void *chebyshev(void *arg) {
  FOREACH_PAIR(
    FOREACH_OTU(
      double d = fabs(x - y);
      if (d > distance) distance = d;
    );
  );
  
  return NULL;
}


//======================================================
// Clark
// sqrt(sum((abs(x - y) / (x + y)) ^ 2))
//======================================================
static void *clark(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double d = (x - y) / (x + y);
      distance += d * d;
    );
  
    distance = sqrt(distance);
  );
  
  return NULL;
}


//======================================================
// Divergence
// 2 * sum((x-y)^2 / (x+y)^2)
//======================================================
static void *divergence(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double diff = x - y;
      double sum  = x + y;
      distance += (diff * diff) / (sum * sum);
    );
  
    distance = 2 * distance;
  );
  
  return NULL;
}


//======================================================
// Euclidean
// sqrt(sum((x-y)^2))
//======================================================
static void *euclidean(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double d = x - y;
      distance += d * d;
    );
  
    distance = sqrt(distance);
  );
  
  return NULL;
}



//======================================================
// Gower
// sum(abs(x-y) / r) / n
//======================================================

static double *gower_range_vec;

static void *gower(void *arg) {
  
  FOREACH_PAIR(
    
    FOREACH_OTU(
      
      double range = gower_range_vec[otu];
      
      if (range) {
        distance += fabs(x - y) / range;
      }
    );
  
    distance /= n_otus;
  );
  
  return NULL;
}

static pthread_func_t gower_setup(void) {
  
  gower_range_vec = (double*) safe_malloc(n_otus * sizeof(double));
  double *min_vec = (double*) safe_malloc(n_otus * sizeof(double));
  double *max_vec = (double*) safe_malloc(n_otus * sizeof(double));
  int    *obs_vec = (int*)    safe_malloc(n_otus * sizeof(int));
  
  memset(min_vec, 0, n_otus * sizeof(double));
  memset(max_vec, 0, n_otus * sizeof(double));
  memset(obs_vec, 0, n_otus * sizeof(int));
  
  
  // Initialize min to first sample's values.
  for (int i = 0; i < pos_vec[1]; i++) {
    min_vec[otu_vec[i]] = val_vec[i];
  }
  
  // Search for min/max values across all samples.
  for (int i = 0; i < pos_vec[n_samples]; i++) {
    int    otu = otu_vec[i];
    double val = val_vec[i];
    obs_vec[otu]++;
    if (val < min_vec[otu]) min_vec[otu] = val;
    if (val > max_vec[otu]) max_vec[otu] = val;
  }
  
  // Assign final range values
  for (int otu = 0; otu < n_otus; otu++) {
    if (obs_vec[otu] < n_samples) {
      gower_range_vec[otu] = max_vec[otu];
    } else {
      gower_range_vec[otu] = max_vec[otu] - min_vec[otu];
    }
  }
  
  free_one(min_vec);
  free_one(max_vec);
  free_one(obs_vec);
  
  return gower;
}


//======================================================
// Hamming
// sum(xor(x, y))
//======================================================
static void *hamming(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = A + B - 2 * J);
  );
  
  return NULL;
}


//======================================================
// Horn
// z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
// 1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))
//======================================================
static void *horn(void *arg) {
  
  FOREACH_PAIR(
    
    double sum_x  = 0;
    double sum_y  = 0;
    double sum_x2 = 0;
    double sum_y2 = 0;
    
    FOREACH_OTU(
      distance += x * y;
      sum_x    += x;
      sum_y    += y;
      sum_x2   += x * x;
      sum_y2   += y * y;
    );
    
    sum_x2 /= sum_x * sum_x;
    sum_y2 /= sum_y * sum_y;
    
    distance = 1 - (2 * distance) / ((sum_x2 + sum_y2) * sum_x * sum_y);
  );
  
  return NULL;
}


//======================================================
// Jaccard
// sum(xor(x, y)) / sum(x | y)
//======================================================
static void *jaccard(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = (A + B - 2 * J) / (A + B - J));
  );
  
  return NULL;
}


//======================================================
// Jensen-Shannon Divergence (JSD)
// sum(x * log(2*x / (x+y)), y * log(2*y / (x+y))) / 2
//======================================================
static void *jsd(void *arg) {
  
  FOREACH_PAIR(
    
    FOREACH_OTU(
      if (x) distance += x * log(2 * x / (x + y));
      if (y) distance += y * log(2 * y / (x + y));
    );
    
    distance /= 2;
  );
  
  return NULL;
}


//======================================================
// Lorentzian
// sum(log(1 + abs(x - y)))
//======================================================
static void *lorentzian(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(distance += log(1 + fabs(x - y)));
  );
  
  return NULL;
}


//======================================================
// Manhattan
// sum(abs(x-y))
//======================================================
static void *manhattan(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(distance += fabs(x - y));
  );
  
  return NULL;
}


//======================================================
// Minkowski
// sum(abs(x - y)^p) ^ (1/p)
//======================================================
static void *minkowski(void *arg) {
  
  double power     = asReal(*sexp_extra);
  double inv_power = 1 / power;
  
  FOREACH_PAIR(
    
    FOREACH_OTU(distance += pow(fabs(x - y), power));
  
    distance = pow(distance, inv_power);
  );
  
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
  
  FOREACH_PAIR(
    
    double sum_x = 0;
    double sum_y = 0;
    double simpson_x = 0;
    double simpson_y = 0;
    
    FOREACH_OTU(
      distance += x * y;
      sum_x    += x; simpson_x += x * (x - 1);
      sum_y    += y; simpson_y += y * (y - 1);
    );
    
    simpson_x /= sum_x * (sum_x - 1);
    simpson_y /= sum_y * (sum_y - 1);
    
    distance = 1 - (2 * distance) / ((simpson_x + simpson_y) * sum_x * sum_y);
  );
  
  return NULL;
}




//======================================================
// Motyka
// sum(pmax(x, y)) / sum(x, y)
//======================================================
static void *motyka(void *arg) {
  
  FOREACH_PAIR(
    
    double sums = 0;
  
    FOREACH_OTU(
      distance += (x > y) ? x : y;
      sums     += x + y;
    );
    
    distance /= sums;
  );
  
  return NULL;
}


//======================================================
// Dice-Sorensen
// 2 * sum(x & y) / sum(x>0, y>0)
//======================================================
static void *sorensen(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = 1 - (2 * J) / (A + B));
  );
  
  return NULL;
}



//======================================================
// Ochiai
// sum((x & y)) / sqrt(sum(x > 0) * sum(y > 0))
//======================================================
static void *ochiai(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = 1 - J / sqrt(A * B));
  );
  
  return NULL;
}
  
  
//======================================================
// Soergel
// 1 - sum(pmin(x, y)) / sum(pmax(x, y))
//======================================================
static void *soergel(void *arg) {
  
  FOREACH_PAIR(
    
    double min_sum = 0;
    double max_sum = 0;
  
    FOREACH_OTU(
      if (x < y) { min_sum += x; max_sum += y; } 
      else       { min_sum += y; max_sum += x; }
    );
    
    distance = 1 - (min_sum / max_sum);
  );
  
  return NULL;
}


//======================================================
// Squared Ch-Squared
// sum((x - y) ^ 2 / (x + y))
//======================================================
static void *squared_chisq(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      double d  = (x - y);
      distance += (d * d) / (x + y);
    );
  );
  
  return NULL;
}


//======================================================
// Squared Chord
// sum((sqrt(x) - sqrt(y)) ^ 2)
//======================================================
static void *squared_chord(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      double d = sqrt(x) - sqrt(y);
      distance += d * d;
    );
  );
  
  return NULL;
}


//======================================================
// Wave Hedges
// sum(abs(x - y) / pmax(x, y))
//======================================================
static void *wave_hedges(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      if (x > y) { distance += (x - y) / x; }
      else       { distance += (y - x) / y; }
    );
  );
  
  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_beta_div(
    SEXP sexp_algorithm,   SEXP sexp_otu_mtx,   
    SEXP sexp_margin,      SEXP sexp_norm, 
    SEXP sexp_pairs_vec,   SEXP sexp_n_threads, 
    SEXP sexp_extra_args ) {
  
  sexp_extra = &sexp_extra_args;
  n_threads  = asInteger(sexp_n_threads);
  init_n_ptrs(10);
  
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  if (!isNull(sexp_norm)) normalize(em, sexp_norm, n_threads);
  
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  otu_vec   = em->otu_vec;
  val_vec   = em->val_vec;
  clr_vec   = em->clr_vec;
  
  
  // function to run
  // void * (*bdiv_func)(void *) = NULL;
  pthread_func_t bdiv_func = NULL;
  
  switch (asInteger(sexp_algorithm)) {
    case BDIV_BHATTACHARYYA: bdiv_func = bhattacharyya; break;
    case BDIV_BRAY:          bdiv_func = bray;          break;
    case BDIV_CANBERRA:      bdiv_func = canberra;      break;
    case BDIV_CHEBYSHEV:     bdiv_func = chebyshev;     break;
    case BDIV_CLARK:         bdiv_func = clark;         break;
    case BDIV_DIVERGENCE:    bdiv_func = divergence;    break;
    case BDIV_EUCLIDEAN:     bdiv_func = euclidean;     break;
    case BDIV_GOWER:         bdiv_func = gower_setup(); break;
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
  
  
  // Create the dist object to return
  n_dist                = n_samples * (n_samples - 1) / 2;
  SEXP sexp_result_dist = PROTECT(allocVector(REALSXP, n_dist));
  dist_vec              = REAL(sexp_result_dist);
  
  SEXP sexp_dist_class = PROTECT(mkString("dist"));
  SEXP sexp_size_val   = PROTECT(ScalarInteger(n_samples));
  SEXP sexp_diag_val   = PROTECT(ScalarLogical(0));
  SEXP sexp_upper_val  = PROTECT(ScalarLogical(0));
  
  setAttrib(sexp_result_dist, R_ClassSymbol,     sexp_dist_class);
  setAttrib(sexp_result_dist, install("Size"),   sexp_size_val);
  setAttrib(sexp_result_dist, install("Diag"),   sexp_diag_val);
  setAttrib(sexp_result_dist, install("Upper"),  sexp_upper_val);
  setAttrib(sexp_result_dist, install("Labels"), em->sexp_sample_names);
  
  
  // Avoid allocating pairs_vec for common all-vs-all case
  if (isNull(sexp_pairs_vec)) {
    
    pairs_vec = NULL;
    n_pairs   = n_dist;
    
  } else {
    
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
    
    for (int i = 0; i < n_dist; i++)
      dist_vec[i] = NA_REAL;
    
    if (n_pairs == 0) {
      free_all();
      UNPROTECT(5);
      return sexp_result_dist;
    }
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_pairs > 100) {
      
      // threads and their thread_i arguments
      pthread_t *tids = (pthread_t*) R_alloc(n_threads, sizeof(pthread_t));
      int       *args = (int*)       R_alloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, bdiv_func, &args[i]);
      for (i = 0; i < n; i++) pthread_join(tids[i], NULL);
      
      free_all();
      UNPROTECT(5);
      return sexp_result_dist;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads     = 1;
  int thread_i  = 0;
  bdiv_func(&thread_i);
  
  free_all();
  UNPROTECT(5);
  return sexp_result_dist;
}
