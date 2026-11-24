// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The sums behind composition_barplot(): each feature's mean proportion,
// and the drawn segments' proportions (and the features' centred log
// ratios) of each sample. Read sparsely, through new_ecomatrix(). Samples
// are split over threads, each writing its own values only, and every sum
// runs in a fixed order, so the results are the same for any number of
// threads. See plans/composition_barplot.md.

#include "ecodive.h"


static int     cb_n_samples;
static int     cb_n_otus;
static int    *cb_pos_vec;
static int    *cb_otu_vec;
static double *cb_val_vec;
static double *cb_totals;  // each sample's total


// Each sample's total.
static void *cb_total_worker (void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  for (int s = thread_i; s < cb_n_samples; s += n_threads) {
    double total = 0;
    for (int i = cb_pos_vec[s]; i < cb_pos_vec[s + 1]; i++) total += cb_val_vec[i];
    cb_totals[s] = total;
  }
  return NULL;
}


// Reads the table and sums each sample. Returns the number of samples.
static int cb_read (SEXP sexp_counts, SEXP sexp_margin, int n_threads) {

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  cb_n_samples = em->n_samples;
  cb_n_otus    = em->n_otus;
  cb_pos_vec   = em->pos_vec;
  cb_otu_vec   = em->otu_vec;
  cb_val_vec   = em->val_vec;
  cb_totals    = (double *) safe_malloc(((size_t) cb_n_samples + 1) * sizeof(double));

  run_parallel(cb_total_worker, n_threads, cb_n_samples);
  return cb_n_samples;
}



//======================================================
// R interface: each sample's total, and each feature's
// mean proportion over the samples whose total is above
// 0 (NaN when there are none). The proportions are summed
// in the samples' order.
//
// Returns list(totals, means).
//======================================================
SEXP C_comp_bar_means (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_n_threads) {

  init_n_ptrs(20);

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  int n = cb_read(sexp_counts, sexp_margin, n_threads);
  int p = cb_n_otus;

  SEXP sexp_totals = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_means  = PROTECT(allocVector(REALSXP, p));
  double *means = REAL(sexp_means);
  for (int j = 0; j < p; j++) means[j] = 0;

  int n_kept = 0;
  for (int s = 0; s < n; s++) {
    double total = cb_totals[s];
    REAL(sexp_totals)[s] = total;
    if (!(total > 0)) continue;
    n_kept++;
    for (int i = cb_pos_vec[s]; i < cb_pos_vec[s + 1]; i++)
      means[cb_otu_vec[i]] += cb_val_vec[i] / total;
  }
  for (int j = 0; j < p; j++) means[j] = n_kept ? means[j] / n_kept : R_NaN;

  free_all();

  SEXP sexp_result = PROTECT(allocVector(VECSXP, 2));
  SEXP sexp_names  = PROTECT(allocVector(STRSXP, 2));
  SET_VECTOR_ELT(sexp_result, 0, sexp_totals); SET_STRING_ELT(sexp_names, 0, mkChar("totals"));
  SET_VECTOR_ELT(sexp_result, 1, sexp_means);  SET_STRING_ELT(sexp_names, 1, mkChar("means"));
  setAttrib(sexp_result, R_NamesSymbol, sexp_names);

  UNPROTECT(4);
  return sexp_result;
}



static int    *cb_seg;       // each feature's segment, or -1
static int     cb_n_seg;
static int     cb_other;     // the segment that holds the rest, or -1
static int    *cb_samples;   // the drawn samples, in the drawn order
static int     cb_m;         // their number
static int    *cb_row;       // each feature's heatmap row, or -1
static int     cb_n_rows;
static double  cb_pc;        // the pseudocount of the log ratios
static double *cb_prop;      // n_seg x m
static double *cb_clr;       // n_rows x m, or NULL


static void *cb_segment_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int K = cb_n_seg, R = cb_n_rows, p = cb_n_otus;
  double log_pc = cb_clr ? log(cb_pc) : 0;

  for (int c = thread_i; c < cb_m; c += n_threads) {

    int     s     = cb_samples[c];
    int     beg   = cb_pos_vec[s], end = cb_pos_vec[s + 1];
    double  total = cb_totals[s];
    double *prop  = cb_prop + (size_t) c * K;

    for (int k = 0; k < K; k++) prop[k] = 0;
    for (int i = beg; i < end; i++) {
      int k = cb_seg[cb_otu_vec[i]];
      if (k >= 0) prop[k] += cb_val_vec[i] / total;
    }

    // The rest: what the other segments leave of 1.
    if (cb_other >= 0) {
      double drawn = 0;
      for (int k = 0; k < K; k++) if (k != cb_other) drawn += prop[k];
      prop[cb_other] = drawn < 1 ? 1 - drawn : 0;
    }

    if (!cb_clr) continue;

    // clr: log(x + pc) minus the sample's mean of those logs, over every
    // feature (zeros give log(pc)).
    double *clr = cb_clr + (size_t) c * R;
    int     n_0 = p - (end - beg);   // zeros: none when the pseudocount is 0
    double  sum = n_0 ? (double) n_0 * log_pc : 0;
    for (int i = beg; i < end; i++) sum += log(cb_val_vec[i] + cb_pc);
    double mean = sum / p;
    for (int r = 0; r < R; r++) clr[r] = log_pc - mean;
    for (int i = beg; i < end; i++) {
      int r = cb_row[cb_otu_vec[i]];
      if (r >= 0) clr[r] = log(cb_val_vec[i] + cb_pc) - mean;
    }
  }

  return NULL;
}



//======================================================
// R interface: the segments' proportions in each drawn
// sample, and the clr of the heatmap's features.
//
// sexp_seg:     for each feature, its 0-based segment, or
//               -1 (in `other`, or not drawn).
// sexp_other:   the 0-based segment that is 1 minus the
//               others, or -1.
// sexp_samples: the 0-based samples to draw, in order;
//               each must have a total above 0.
// sexp_row:     for each feature, its 0-based heatmap row,
//               or -1; NULL for no heatmap.
//
// Returns list(proportions (n_seg x m), clr (n_rows x m,
// or NULL)).
//======================================================
SEXP C_comp_bar_segments (
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_seg, SEXP sexp_n_seg,
    SEXP sexp_other, SEXP sexp_samples, SEXP sexp_row, SEXP sexp_n_rows,
    SEXP sexp_pseudocount, SEXP sexp_n_threads ) {

  init_n_ptrs(20);

  int n_threads = asInteger(sexp_n_threads);
  if (n_threads < 1) n_threads = 1;

  cb_read(sexp_counts, sexp_margin, n_threads);

  cb_seg     = INTEGER(sexp_seg);
  cb_n_seg   = asInteger(sexp_n_seg);
  cb_other   = asInteger(sexp_other);
  cb_samples = INTEGER(sexp_samples);
  cb_m       = LENGTH(sexp_samples);
  cb_row     = isNull(sexp_row) ? NULL : INTEGER(sexp_row);
  cb_n_rows  = cb_row ? asInteger(sexp_n_rows) : 0;
  cb_pc      = asReal(sexp_pseudocount);

  SEXP sexp_prop = PROTECT(allocMatrix(REALSXP, cb_n_seg, cb_m));
  SEXP sexp_clr  = PROTECT(cb_row ? allocMatrix(REALSXP, cb_n_rows, cb_m) : R_NilValue);
  cb_prop = REAL(sexp_prop);
  cb_clr  = cb_row ? REAL(sexp_clr) : NULL;

  run_parallel(cb_segment_worker, n_threads, cb_m);

  free_all();

  SEXP sexp_result = PROTECT(allocVector(VECSXP, 2));
  SEXP sexp_names  = PROTECT(allocVector(STRSXP, 2));
  SET_VECTOR_ELT(sexp_result, 0, sexp_prop); SET_STRING_ELT(sexp_names, 0, mkChar("proportions"));
  SET_VECTOR_ELT(sexp_result, 1, sexp_clr);  SET_STRING_ELT(sexp_names, 1, mkChar("clr"));
  setAttrib(sexp_result, R_NamesSymbol, sexp_names);

  UNPROTECT(4);
  return sexp_result;
}
