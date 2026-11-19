// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


#define NORM_NONE    0
#define NORM_PERCENT 1
#define NORM_CLR     2
#define NORM_CHORD   3
#define NORM_BINARY  4
#define NORM_RCLR    5
#define NORM_HELLINGER 6

static int     n_samples;
static int     n_otus;
static int    *pos_vec;
static double *val_vec;
static double *clr_vec;
static double  pseudocount;
static int     is_percent_normalized;

// Worker threads may clear the flag at the same time. They all
// store the same value, and it is read only after they join,
// but a plain store would still be a data race (TSan reports it).
#if defined(__GNUC__) || defined(__clang__)
  #define CLEAR_FLAG(f) __atomic_store_n(&(f), 0, __ATOMIC_RELAXED)
#else
  #define CLEAR_FLAG(f) ((f) = 0)
#endif


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
        CLEAR_FLAG(is_percent_normalized);
        return NULL;
      }
      depth += *val
    );
    if (fabs(depth - 1) > 1e-6) {
      CLEAR_FLAG(is_percent_normalized);
      return NULL;
    }
  );
  return NULL;
}


//======================================================
// Centered Log Ratio
// log(x) - mean(log(x)), where x = counts + pseudocount
//======================================================
static void *norm_clr(void *arg) {
  FOREACH_SAMPLE(

    int n_zeros = n_otus - (int)(val_end - val_begin);

    double logsum = 0;
    FOREACH_VAL(logsum += log(*val + pseudocount));

    // Add in all the zero abundances. Skip when there are none:
    // pseudocount is 0 for zero-free data, and log(0) * 0 is NaN.
    if (n_zeros > 0)
      logsum += log(pseudocount) * n_zeros;

    double log_mean = logsum / n_otus;
    FOREACH_VAL(*val = log(*val + pseudocount) - log_mean);

    // CLR value for zero abundance OTUs (unused when n_zeros is 0)
    clr_vec[sam] = (n_zeros > 0) ? log(pseudocount) - log_mean : 0;

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

    // Very small or large values: the sum of squares under- or
    // overflows. Scale by the largest value first (only then, so
    // that other results are unchanged).
    if (!(sq_sum >= DBL_MIN && sq_sum <= DBL_MAX)) {
      double mx = 0;
      FOREACH_VAL(if (fabs(*val) > mx) mx = fabs(*val));
      if (mx > 0 && mx <= DBL_MAX) {
        FOREACH_VAL(*val /= mx);
        sq_sum = 0;
        FOREACH_VAL(sq_sum += *val * *val);
      }
    }

    sq_sum = sqrt(sq_sum);
    FOREACH_VAL(*val /= sq_sum);
    
  );
  return NULL;
}


//======================================================
// Hellinger-Transformed Relative Abundance
// sqrt(x / sum(x))
//======================================================
static void *norm_hellinger(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    FOREACH_VAL(depth += *val);
    FOREACH_VAL(*val = sqrt(*val / depth));
  );
  return NULL;
}


//======================================================
// Set each abundance to 1 (present) or 0 (absent)
// x > 0
//======================================================
static void *norm_binary(void *arg) {
  FOREACH_SAMPLE(
    FOREACH_VAL(if (!ISNAN(*val)) *val = 1); // NA stays NA, not a presence
  );
  return NULL;
}


//======================================================
// Robust Centered Log Ratio (rclr)
// log(x) - mean(log(x[x>0]))
//======================================================
static void *norm_rclr(void *arg) {
  FOREACH_SAMPLE(
    
    int count = val_end - val_begin;
  
    if (count > 0) {
      double logsum = 0;
      
      // 1. Sum of logs for strictly positive elements
      FOREACH_VAL(logsum += log(*val));
      
      // 2. Log of the robust geometric mean
      double log_mean = logsum / count;
      
      // 3. Apply rclr transformation
      FOREACH_VAL(*val = log(*val) - log_mean);
    }
  
  );
  return NULL;
}


void normalize(ecomatrix_t *em, int norm, int n_threads, double pseudocount_) {
  
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  val_vec   = em->val_vec;
  
  if (norm == NORM_PERCENT) {
    // Check if it's already normalized to percent. If so, skip.
    if (em->nnz > 0 && *val_vec <= 1) {
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
    case NORM_RCLR:    norm_func = norm_rclr;    break;
    case NORM_HELLINGER: norm_func = norm_hellinger; break;
  }
  
  pseudocount = pseudocount_;
  val_vec     = rw_val_vec(em);
  
  run_parallel(norm_func, n_threads, n_samples);
}


//======================================================
// An ecomatrix as a table in the input's orientation
// (samples in rows for margin 1, in columns for margin 2),
// for normalize() and qmp_scale().
//
// dense = 1: an n_rows x n_cols double matrix, with zeros
//            (or, with z0, each sample's value for zero).
// dense = 0: list(i, j, x), the non-zero entries' 0-based
//            row and column, and value, in column-major
//            order.
//======================================================
static SEXP ecomatrix_table(
    ecomatrix_t *em, int margin, int dense, const double *z0 ) {

  int        ns  = em->n_samples;
  int        no  = em->n_otus;
  int        nnz = em->nnz;
  const int *pos = em->pos_vec;
  const int *otu = em->otu_vec;
  const double *val = em->val_vec;

  int n_rows = (margin == 1) ? ns : no;
  SEXP result;

  if (dense) {

    result = PROTECT(allocMatrix(REALSXP, n_rows, (margin == 1) ? no : ns));
    double *out = REAL(result);

    for (int s = 0; s < ns; s++) {

      double zero = z0 ? z0[s] : 0;

      if (margin == 1) {
        for (int o = 0; o < no; o++) out[s + (size_t) n_rows * o] = zero;
        for (int e = pos[s]; e < pos[s + 1]; e++)
          out[s + (size_t) n_rows * otu[e]] = val[e];
      }
      else {
        double *col = out + (size_t) n_rows * s;
        for (int o = 0; o < no; o++) col[o] = zero;
        for (int e = pos[s]; e < pos[s + 1]; e++) col[otu[e]] = val[e];
      }
    }
  }
  else {

    const char *names[] = { "i", "j", "x", "" };
    result = PROTECT(mkNamed(VECSXP, names));
    SEXP sexp_i = allocVector(INTSXP,  nnz); SET_VECTOR_ELT(result, 0, sexp_i);
    SEXP sexp_j = allocVector(INTSXP,  nnz); SET_VECTOR_ELT(result, 1, sexp_j);
    SEXP sexp_x = allocVector(REALSXP, nnz); SET_VECTOR_ELT(result, 2, sexp_x);
    int    *ri = INTEGER(sexp_i);
    int    *rj = INTEGER(sexp_j);
    double *rx = REAL(sexp_x);

    if (margin == 2) {
      // Samples are the columns: already column-major.
      for (int s = 0; s < ns; s++)
        for (int e = pos[s]; e < pos[s + 1]; e++) {
          ri[e] = otu[e];
          rj[e] = s;
          rx[e] = val[e];
        }
    }
    else {
      // Features are the columns: a stable counting sort
      // by feature keeps each column's samples in order.
      int *fill = (int *) R_alloc((size_t) no + 1, sizeof(int));
      memset(fill, 0, ((size_t) no + 1) * sizeof(int));
      for (int e = 0; e < nnz; e++) fill[otu[e] + 1]++;
      for (int o = 0; o < no; o++)  fill[o + 1] += fill[o];
      for (int s = 0; s < ns; s++)
        for (int e = pos[s]; e < pos[s + 1]; e++) {
          int d = fill[otu[e]]++;
          ri[d] = s;
          rj[d] = otu[e];
          rx[d] = val[e];
        }
    }
  }

  UNPROTECT(1);
  return result;
}


//======================================================
// R interface, for normalize(): the normalized table
// (see ecomatrix_table()).
//======================================================
SEXP C_normalize(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_dense, SEXP sexp_n_threads ) {

  int    margin = asInteger(sexp_margin);
  int    norm   = asInteger(sexp_norm);
  double pc     = asReal(sexp_pseudocount);
  int    dense  = asLogical(sexp_dense);
  int    n_thr  = asInteger(sexp_n_threads);

  init_n_ptrs(10);

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  if (norm) normalize(em, norm, n_thr, pc);

  const double *z0 = (norm == NORM_CLR) ? em->clr_vec : NULL;
  SEXP result = PROTECT(ecomatrix_table(em, margin, dense, z0));

  free_all();
  UNPROTECT(1);
  return result;
}


//======================================================
// Scaling to a measured load (qmp_scale(), qmp_rarefy()):
// x / sum(x) * load, per sample. Each sample's total is
// written to scale_total_vec; a sample with total 0 is
// left as it is (R reports it).
//======================================================
static const double *scale_load_vec;
static double       *scale_total_vec;

static void *norm_scale(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    FOREACH_VAL(depth += *val);
    scale_total_vec[sam] = depth;
    if (depth > 0) {
      double load = scale_load_vec[sam];
      FOREACH_VAL(*val = *val / depth * load);
    }
  );
  return NULL;
}


//======================================================
// R interface, for qmp_scale() and qmp_rarefy(): the
// scaled table (see ecomatrix_table()), with attribute
// "totals", each sample's total before scaling.
//======================================================
SEXP C_qmp_scale(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_load,
    SEXP sexp_dense, SEXP sexp_n_threads ) {

  int margin = asInteger(sexp_margin);
  int dense  = asLogical(sexp_dense);
  int n_thr  = asInteger(sexp_n_threads);

  init_n_ptrs(10);

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  if (LENGTH(sexp_load) != em->n_samples) {
    free_all();                                                    // # nocov
    error("`load` must have one value per sample.");               // # nocov
  }

  SEXP sexp_totals = PROTECT(allocVector(REALSXP, em->n_samples));

  n_samples       = em->n_samples;
  n_otus          = em->n_otus;
  pos_vec         = em->pos_vec;
  val_vec         = rw_val_vec(em);
  scale_load_vec  = REAL(sexp_load);
  scale_total_vec = REAL(sexp_totals);

  run_parallel(norm_scale, n_thr, n_samples);

  SEXP result = PROTECT(ecomatrix_table(em, margin, dense, NULL));
  setAttrib(result, install("totals"), sexp_totals);

  free_all();
  UNPROTECT(2);
  return result;
}
