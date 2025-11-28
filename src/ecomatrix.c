// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit



#include <R.h>
#include <Rinternals.h>
#include <string.h> // memcpy
#include "get.h"
#include "ecomatrix.h"
#include "normalize.h"
#include "memory.h"


// For C-level debugging
void debug_ecomatrix(ecomatrix_t *em, const char *msg) {
  
  Rprintf("--- DEBUG: %s ---\n", msg);
  Rprintf("n_samples: %d\n", em->n_samples);
  Rprintf("n_otus: %d\n", em->n_otus);
  Rprintf("nnz: %d\n", em->nnz);
  
  if (em->sam_vec) {
    Rprintf("sam_vec (safe=%d):", is_safe_ptr(em->sam_vec));
    for (int i = 0; i < em->nnz; i++) {
      Rprintf("  %d", em->sam_vec[i]);
    }
    Rprintf("\n");
  } else {
    Rprintf("sam_vec: NULL\n");
  }
  
  if (em->pos_vec) {
    Rprintf("pos_vec (safe=%d):", is_safe_ptr(em->pos_vec));
    for (int i = 0; i < em->n_samples + 1; i++) {
      Rprintf("  %d", em->pos_vec[i]);
    }
    Rprintf("\n");
  } else {
    Rprintf("pos_vec: NULL\n");
  }
  
  if (em->otu_vec) {
    Rprintf("otu_vec (safe=%d):", is_safe_ptr(em->otu_vec));
    for (int i = 0; i < em->nnz; i++) {
      Rprintf("  %d", em->otu_vec[i]);
    }
    Rprintf("\n");
  } else {
    Rprintf("otu_vec: NULL\n");
  }
  
  if (em->val_vec) {
    Rprintf("val_vec (safe=%d):", is_safe_ptr(em->val_vec));
    for (int i = 0; i < em->nnz; i++) {
      Rprintf("  %d", (int)em->val_vec[i]);
    }
    Rprintf("\n");
  } else {
    Rprintf("val_vec: NULL\n");
  }
  
  free_all();
  error("Halting execution after debug print.");
}


//=========================================================
// Ensure we don't overwrite values in other R objects.
//=========================================================

static void rw_vec (void **ptr, int n_bytes) {
  
  if (is_safe_ptr(*ptr)) return;
  
  void *new_ptr = safe_malloc(n_bytes);
  if (*ptr) memcpy(new_ptr, *ptr, n_bytes);
  *ptr = new_ptr;
}

static int* rw_sam_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->sam_vec), em->nnz * sizeof(int));
  return em->sam_vec;
}

static int* rw_pos_vec (ecomatrix_t *em) {
  int vec_len = em->n_samples + 1;
  rw_vec((void**)&(em->pos_vec), vec_len * sizeof(int));
  return em->pos_vec;
}

static int* rw_otu_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->otu_vec), em->nnz * sizeof(int));
  return em->otu_vec;
}

double* rw_val_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->val_vec), em->nnz * sizeof(double));
  return em->val_vec;
}

double* rw_clr_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->clr_vec), em->n_samples * sizeof(double));
  return em->clr_vec;
}



//=========================================================
// Accept double, integer, or logical values.
//=========================================================

static void assign_sparse_vals(ecomatrix_t *em, SEXP sexp_vals) {
  
  if (isReal(sexp_vals)) {
    em->val_vec = REAL(sexp_vals);
  }
  
  else if (isInteger(sexp_vals) || isLogical(sexp_vals)) {
    
    double *val_vec = rw_val_vec(em);
    
    int  nnz = em->nnz;
    int *tmp = INTEGER(sexp_vals);
    for (int i = 0; i < nnz; i++) {
      val_vec[i] = (double)(tmp[i]);
    }
  }
  
  else {
    error("Input must be numeric");
  }
  
}




static void assign_dense_vals(ecomatrix_t *em, SEXP sexp_vals, int margin) {
  
  int n_samples = em->n_samples;
  int n_otus    = em->n_otus;
  int mtx_len   = length(sexp_vals);
  
  double *dbl_vals = NULL;
  int    *int_vals = NULL;
  
  
  // Accept double, integer, or logical values.
  // --------------------------------------------
  
  if (isReal(sexp_vals)) {
    dbl_vals = REAL(sexp_vals);
  }
  else if (isInteger(sexp_vals) || isLogical(sexp_vals)) {
    int_vals = INTEGER(sexp_vals);
  }
  else {
    error("Input must be numeric");
  }
  
  
  // Count the number of non-zero entries.
  // --------------------------------------------
  
  int nnz = 0;
  
  if (dbl_vals != NULL) {
    for (int i = 0; i < mtx_len; i++) {
      if (dbl_vals[i]) nnz++;
    }
  }
  else {
    for (int i = 0; i < mtx_len; i++) {
      if (int_vals[i]) nnz++;
    }
  }
  em->nnz = nnz;
  
  
  // Allocate based on n_samples and nnz
  int    *pos_vec = rw_pos_vec(em);
  int    *otu_vec = rw_otu_vec(em);
  double *val_vec = rw_val_vec(em);
  if (!(pos_vec && otu_vec && val_vec)) return;
  
  pos_vec[n_samples] = nnz;
  
  
  
  // Four modes: int|dbl * margin1|margin2
  // --------------------------------------------
  
  int i = 0;
  
  if (margin == 1) {
    
    int n_rows = n_samples;
    int n_cols = n_otus;
    
    if (dbl_vals != NULL) { // Importing doubles
      
      for (int row = 0; row < n_rows; row++) {
        pos_vec[row] = i;
        for (int col = 0; col < n_cols; col++) {
          double v = dbl_vals[row + n_rows * col];
          if (v) {
            otu_vec[i] = col;
            val_vec[i] = v;
            i++;
          }
        }
      }
      
    }
    else { // Importing integers/logicals
      
      for (int row = 0; row < n_rows; row++) {
        pos_vec[row] = i;
        for (int col = 0; col < n_cols; col++) {
          int v = int_vals[row + n_rows * col];
          if (v) {
            otu_vec[i] = col;
            val_vec[i] = (double)v;
            i++;
          }
        }
      }
      
    }
    
  }
  else { // margin == 2
    
    int n_cols = n_samples;
    int n_rows = n_otus;
    
    if (dbl_vals != NULL) { // Importing doubles
    
      for (int col = 0; col < n_cols; col++) {
        pos_vec[col] = i;
        for (int row = 0; row < n_rows; row++) {
          double v = dbl_vals[row + n_rows * col];
          if (v) {
            otu_vec[i] = row;
            val_vec[i] = v;
            i++;
          }
        }
      }
      
    }
    else { // Importing integers/logicals
      
      for (int col = 0; col < n_cols; col++) {
        pos_vec[col] = i;
        for (int row = 0; row < n_rows; row++) {
          int v = int_vals[row + n_rows * col];
          if (v) {
            otu_vec[i] = row;
            val_vec[i] = (double)v;
            i++;
          }
        }
      }
      
    }
    
  }
  
  
}



//=========================================================
// Sorts an uncompressed sparse matrix's vectors.
//=========================================================

static void sort_triplet_q (int *sam_vec, int *otu_vec, double *val_vec, int lo, int hi) {
  
  int    sam, otu;
  double val;
  
  int pivot_sam = sam_vec[hi];
  int pivot_otu = otu_vec[hi];
  int i = lo;
  
  for (int j = lo; j <= hi - 1; j++) {
    
    sam = sam_vec[j];
    
    if (sam < pivot_sam || (sam == pivot_sam && otu_vec[j] < pivot_otu)) {
      sam = sam_vec[i]; sam_vec[i] = sam_vec[j]; sam_vec[j] = sam;
      otu = otu_vec[i]; otu_vec[i] = otu_vec[j]; otu_vec[j] = otu;
      val = val_vec[i]; val_vec[i] = val_vec[j]; val_vec[j] = val;
      i++;
    }
  }
  
  sam = sam_vec[i]; sam_vec[i] = sam_vec[hi]; sam_vec[hi] = sam;
  otu = otu_vec[i]; otu_vec[i] = otu_vec[hi]; otu_vec[hi] = otu;
  val = val_vec[i]; val_vec[i] = val_vec[hi]; val_vec[hi] = val;
  
  if (lo < i - 1) sort_triplet_q(sam_vec, otu_vec, val_vec, lo, i - 1);
  if (i + 1 < hi) sort_triplet_q(sam_vec, otu_vec, val_vec, i + 1, hi);
  
}

static void sort_triplet (ecomatrix_t *em) {
  
  int     nnz     = em->nnz;
  int    *sam_vec = em->sam_vec;
  int    *otu_vec = em->otu_vec;
  double *val_vec = em->val_vec;
  
  
  // Check if it's already sorted
  // --------------------------------------------
  
  if (nnz < 2) return;
  
  int sorted = 1;
  for (int i = 1; i < nnz; i++) {
    int s1 = sam_vec[i-1];
    int s2 = sam_vec[i];
    if (s1 == s2) {
      if (otu_vec[i-1] > otu_vec[i]) {
        sorted = 0;
        break;
      }
    }
    else if (s1 > s2) {
      sorted = 0;
      break;
    }
  }
  
  if (sorted) return;
  
  
  // Make sure we're not modifying original data.
  // --------------------------------------------
  sam_vec = rw_sam_vec(em);
  otu_vec = rw_otu_vec(em);
  val_vec = rw_val_vec(em);
  
  sort_triplet_q(sam_vec, otu_vec, val_vec, 0, nnz - 1);
}



//=========================================================
// Reorder sam/otu/val by sam/otu, then populate pos.
//=========================================================
static void compress_triplet (ecomatrix_t *em) {
  
  sort_triplet(em);
  
  int  n_samples = em->n_samples;
  int  nnz       = em->nnz;
  int *sam_vec   = em->sam_vec;
  int *pos_vec   = rw_pos_vec(em);
  
  int p = 0;
  for (int i = 0; i < n_samples; i++) {
    pos_vec[i] = p;
    while (p < nnz && sam_vec[p] == i) p++;
  }
  pos_vec[n_samples] = nnz;
  
  em->sam_vec = maybe_free_one(em->sam_vec);
}



//=========================================================
// Convert pos_vec to otu_vec.
//=========================================================
static void inflate_triplet_otus (ecomatrix_t *em) {
  
  int  n_otus  = em->n_otus;
  int *pos_vec = em->pos_vec;
  int *otu_vec = rw_otu_vec(em);
  
  for (int otu = 0; otu < n_otus; otu++) {
    int end = pos_vec[otu + 1];
    for (int i = pos_vec[otu]; i < end; i++) {
      otu_vec[i] = otu;
    }
  }
  
  em->pos_vec = maybe_free_one(em->pos_vec);
}



//=========================================================
// A matrix object from base R.
//=========================================================

static void parse_matrix (ecomatrix_t *em, SEXP sexp_matrix, int margin) {
  
  if (!isMatrix(sexp_matrix)) error("Input must be a matrix");
  
  if (margin == 1) {
    
    em->n_samples = nrows(sexp_matrix); // samples are in rows
    em->n_otus    = ncols(sexp_matrix); // OTUs are in columns
    
    SEXP sexp_dimnames = getAttrib(sexp_matrix, R_DimNamesSymbol);
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
    
  }
  else { // margin == 2
    
    em->n_samples = ncols(sexp_matrix); // samples are in columns
    em->n_otus    = nrows(sexp_matrix); // OTUs are in rows
    
    SEXP sexp_dimnames = getAttrib(sexp_matrix, R_DimNamesSymbol);
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
    
  }
  
  assign_dense_vals(em, sexp_matrix, margin);
}



//=========================================================
// A dense matrix from the Matrix package.
//=========================================================

static void parse_dgeMatrix (ecomatrix_t *em, SEXP sexp_dgeMatrix, int margin) {
  
  SEXP sexp_dge_x    = R_do_slot(sexp_dgeMatrix, install("x"));
  SEXP sexp_dim      = R_do_slot(sexp_dgeMatrix, install("Dim"));
  SEXP sexp_dimnames = R_do_slot(sexp_dgeMatrix, install("Dimnames"));
  int  n_rows        = INTEGER(sexp_dim)[0];
  int  n_cols        = INTEGER(sexp_dim)[1];
  
  
  if (margin == 1) {
    
    em->n_samples = n_rows; // samples are in rows
    em->n_otus    = n_cols; // OTUs are in columns
    
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
    
  }
  else { // margin == 2
    
    em->n_samples = n_cols; // samples are in columns
    em->n_otus    = n_rows; // OTUs are in rows
    
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
    
  }
  
  assign_dense_vals(em, sexp_dge_x, margin);
}



//=========================================================
// A simple_triplet_matrix from the slam package.
// Can NOT assume that i and j are sorted.
//=========================================================

static void parse_slam (ecomatrix_t *em, SEXP sexp_slam_mtx, int margin) {
  
  SEXP sexp_slam_i = get(sexp_slam_mtx, "i");
  SEXP sexp_slam_j = get(sexp_slam_mtx, "j");
  SEXP sexp_slam_v = get(sexp_slam_mtx, "v");
  int  n_rows      = asInteger(get(sexp_slam_mtx, "nrow"));
  int  n_cols      = asInteger(get(sexp_slam_mtx, "ncol"));
  int  nnz         = length(sexp_slam_v);
  
  
  // Import values as double precision
  em->nnz = nnz;
  assign_sparse_vals(em, sexp_slam_v);
  
  
  // The sample names will be the 1st or 2nd element of dimnames.
  SEXP sexp_dimnames = get(sexp_slam_mtx, "dimnames");
  
  
  // Accept samples in either rows or columns.
  //-------------------------------------------------------
  
  if (margin == 1) {
    
    em->n_samples = n_rows; // samples are in rows
    em->n_otus    = n_cols; // OTUs are in columns
    em->sam_vec   = INTEGER(sexp_slam_i);
    em->otu_vec   = INTEGER(sexp_slam_j);
    
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
  }
  
  else {
    
    em->n_samples = n_cols; // samples are in columns
    em->n_otus    = n_rows; // OTUs are in rows
    em->sam_vec   = INTEGER(sexp_slam_j);
    em->otu_vec   = INTEGER(sexp_slam_i);
    
    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
  }
  

  // debug_ecomatrix(em, "Before ingest");
  
  
  // Slam indices are 1-based. Convert to 0-based.
  //-------------------------------------------------------
  
  int *sam_vec = rw_sam_vec(em);
  int *otu_vec = rw_otu_vec(em);
  
  for (int i = 0; i < nnz; i++) {
    sam_vec[i]--;
    otu_vec[i]--;
  }
  
  
  // Compress and cleanup
  compress_triplet(em);
}




//=========================================================
// An sparse triplet matrix from the Matrix package.
// Can NOT assume that i and j are sorted.
//=========================================================

static void parse_dgTMatrix (ecomatrix_t *em, SEXP sexp_dgTMatrix, int margin) {
  
  SEXP sexp_dgt_i    = R_do_slot(sexp_dgTMatrix, install("i"));
  SEXP sexp_dgt_j    = R_do_slot(sexp_dgTMatrix, install("j"));
  SEXP sexp_dgt_x    = R_do_slot(sexp_dgTMatrix, install("x"));
  SEXP sexp_dim      = R_do_slot(sexp_dgTMatrix, install("Dim"));
  SEXP sexp_dimnames = R_do_slot(sexp_dgTMatrix, install("Dimnames"));
  int  n_rows        = INTEGER(sexp_dim)[0];
  int  n_cols        = INTEGER(sexp_dim)[1];
  int  nnz           = length(sexp_dgt_x);
  
  
  // Import values as double precision
  em->nnz = nnz;
  assign_sparse_vals(em, sexp_dgt_x);
  
  
  // Transpose if needed
  if (margin == 1) {
    
    em->n_samples         = n_rows; // samples are in rows
    em->n_otus            = n_cols; // OTUs are in columns
    em->sam_vec           = INTEGER(sexp_dgt_i);
    em->otu_vec           = INTEGER(sexp_dgt_j);
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
  }
  
  else { // margin == 2
    
    em->n_samples         = n_cols; // samples are in columns
    em->n_otus            = n_rows; // OTUs are in rows
    em->sam_vec           = INTEGER(sexp_dgt_j);
    em->otu_vec           = INTEGER(sexp_dgt_i);
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
  }
  
  
  compress_triplet(em);
}



//=========================================================
// A compressed sparse matrix from the Matrix package.
// CAN assume that i is sorted.
//=========================================================

static void parse_dgCMatrix (ecomatrix_t *em, SEXP sexp_dgCMatrix, int margin) {
  
  SEXP sexp_dgc_i    = R_do_slot(sexp_dgCMatrix, install("i"));
  SEXP sexp_dgc_p    = R_do_slot(sexp_dgCMatrix, install("p"));
  SEXP sexp_dgc_x    = R_do_slot(sexp_dgCMatrix, install("x"));
  SEXP sexp_dim      = R_do_slot(sexp_dgCMatrix, install("Dim"));
  SEXP sexp_dimnames = R_do_slot(sexp_dgCMatrix, install("Dimnames"));
  int  n_rows        = INTEGER(sexp_dim)[0];
  int  n_cols        = INTEGER(sexp_dim)[1];
  int  nnz           = length(sexp_dgc_x);
  
  
  // Import values as double precision
  em->nnz = nnz;
  assign_sparse_vals(em, sexp_dgc_x);
  
  
  // Transpose if needed
  if (margin == 1) {
    
    em->n_samples         = n_rows; // samples are in rows
    em->n_otus            = n_cols; // OTUs are in columns
    em->sam_vec           = INTEGER(sexp_dgc_i);
    em->pos_vec           = INTEGER(sexp_dgc_p);
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
    
    inflate_triplet_otus(em);
    compress_triplet(em);
  }
  
  else { // margin == 2
    
    em->n_samples         = n_cols; // samples are in columns
    em->n_otus            = n_rows; // OTUs are in rows
    em->pos_vec           = INTEGER(sexp_dgc_p);
    em->otu_vec           = INTEGER(sexp_dgc_i);
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
  }
  
}



//=========================================================
// Initialize a new ecomatrix_t struct.
//=========================================================
ecomatrix_t* new_ecomatrix(SEXP sexp_matrix, SEXP sexp_margin) {
  
  int margin = asInteger(sexp_margin);
  
  
  // function to run
  void (*parse_func)(ecomatrix_t*, SEXP, int) = NULL;
  
  if      (isMatrix(sexp_matrix))                          { parse_func = parse_matrix;    }
  else if (inherits(sexp_matrix, "simple_triplet_matrix")) { parse_func = parse_slam;      }
  else if (inherits(sexp_matrix, "dgCMatrix"))             { parse_func = parse_dgCMatrix; }
  else if (inherits(sexp_matrix, "dgTMatrix"))             { parse_func = parse_dgTMatrix; }
  else if (inherits(sexp_matrix, "dgeMatrix"))             { parse_func = parse_dgeMatrix; }
  else    { error("Unrecognized matrix format."); }
  
  
  ecomatrix_t *em       = (ecomatrix_t*) safe_malloc(sizeof(ecomatrix_t));
  em->n_samples         = 0;
  em->n_otus            = 0;
  em->nnz               = 0;
  em->sam_vec           = NULL;
  em->pos_vec           = NULL;
  em->otu_vec           = NULL;
  em->val_vec           = NULL;
  em->clr_vec           = NULL;
  em->sexp_sample_names = R_NilValue;
  
  parse_func(em, sexp_matrix, margin);
  
  return em;
}
