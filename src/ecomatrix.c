// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"

// An integer or logical value as a double, keeping NA as NA, as R does.
#define INT_TO_DBL(x) ((x) == NA_INTEGER ? NA_REAL : (double)(x))


//=========================================================
// Ensure we don't overwrite values in other R objects.
//=========================================================

static void rw_vec (void **ptr, size_t n_bytes) {

  if (is_safe_ptr(*ptr)) return;

  void *new_ptr = safe_malloc(n_bytes);
  if (*ptr) memcpy(new_ptr, *ptr, n_bytes);
  *ptr = new_ptr;
}

static int* rw_pos_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->pos_vec), ((size_t) em->n_samples + 1) * sizeof(int));
  return em->pos_vec;
}

static int* rw_otu_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->otu_vec), (size_t) em->nnz * sizeof(int));
  return em->otu_vec;
}

double* rw_val_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->val_vec), (size_t) em->nnz * sizeof(double));
  return em->val_vec;
}

double* rw_clr_vec (ecomatrix_t *em) {
  rw_vec((void**)&(em->clr_vec), (size_t) em->n_samples * sizeof(double));
  return em->clr_vec;
}


// Frees ecodive's allocations, then raises an R error.
static void import_error (const char *msg) {
  free_all();
  error("%s", msg);
}


// A sparse matrix's index vector must be an integer vector of
// the given length (that of its values).
static int* index_vec (SEXP sexp_idx, int nnz, const char *what) {
  if (TYPEOF(sexp_idx) != INTSXP || XLENGTH(sexp_idx) != nnz) {
    free_all();
    error("Sparse matrix is invalid: `%s` must be an integer vector as long as its values.", what);
  }
  return INTEGER(sexp_idx);
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
      val_vec[i] = INT_TO_DBL(tmp[i]);
    }
  }
  
  else {
    error("Input must be numeric"); // # nocov
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
            val_vec[i] = INT_TO_DBL(v);
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
            val_vec[i] = INT_TO_DBL(v);
            i++;
          }
        }
      }
      
    }
    
  }
  
  
}



//=========================================================
// Sparse input that isn't already sample-major (dgCMatrix
// with margin = 1, dgTMatrix, slam) is reordered so each
// sample's entries are contiguous with otus increasing.
//
// Stable counting sorts make this O(nnz + n_samples +
// n_otus). Since (sample, otu) pairs are unique, the order
// is fully determined. Duplicate pairs, if any, are kept
// as separate entries, adjacent and in input order.
//=========================================================


// Values may be stored as doubles, integers, or logicals.
static void sparse_val_ptrs (SEXP sexp_vals, double **dbl_in, int **int_in) {
  
  *dbl_in = NULL;
  *int_in = NULL;
  
  if      (isReal(sexp_vals))                            { *dbl_in = REAL(sexp_vals);    }
  else if (isInteger(sexp_vals) || isLogical(sexp_vals)) { *int_in = INTEGER(sexp_vals); }
  else    { error("Input must be numeric"); } // # nocov
}


// Every index must be in [base, base + n). Indices become
// memory offsets below, so out-of-range ones are fatal.
static void check_indices (int *idx_vec, int nnz, int n, int base) {
  
  for (int i = 0; i < nnz; i++) {
    int idx = idx_vec[i];
    if (idx < base || idx - base >= n) {
      free_all();
      error("Sparse matrix has an index out of range.");
    }
  }
}


// Allocate pos_vec and set pos_vec[s] to where sample s's
// entries start; pos_vec[n_samples] = nnz.
static int* count_samples (ecomatrix_t *em, int *sam_in, int base) {
  
  int  n_samples = em->n_samples;
  int  nnz       = em->nnz;
  int *pos_vec   = rw_pos_vec(em);
  
  memset(pos_vec, 0, ((size_t) n_samples + 1) * sizeof(int));
  
  for (int i = 0; i < nnz; i++)
    pos_vec[sam_in[i] - base + 1]++;
  
  for (int sam = 0; sam < n_samples; sam++)
    pos_vec[sam + 1] += pos_vec[sam];
  
  return pos_vec;
}


// Scattering with ptr_vec[b]++ as bucket b's write cursor
// leaves ptr_vec[b] at bucket b + 1's start. Shift back.
static void unshift_ptrs (int *ptr_vec, int n_buckets) {
  
  for (int b = n_buckets; b > 0; b--)
    ptr_vec[b] = ptr_vec[b - 1];
  
  ptr_vec[0] = 0;
}


// Entries grouped by otu (column-compressed: col_ptr[otu]
// is where otu's entries start) are scattered by sample.
// Visiting otus in increasing order and keeping each
// otu's entries in order makes the scatter stable.
static void scatter_by_sample (
    ecomatrix_t *em, int *sam_in, int base, int *col_ptr, 
    double *dbl_in, int *int_in ) {
  
  int n_samples = em->n_samples;
  int n_otus    = em->n_otus;
  
  int    *pos_vec = count_samples(em, sam_in, base);
  int    *otu_vec = rw_otu_vec(em);
  double *val_vec = rw_val_vec(em);
  
  for (int otu = 0; otu < n_otus; otu++) {
    int end = col_ptr[otu + 1];
    for (int i = col_ptr[otu]; i < end; i++) {
      int dst = pos_vec[sam_in[i] - base]++;
      otu_vec[dst] = otu;
      val_vec[dst] = dbl_in ? dbl_in[i] : INT_TO_DBL(int_in[i]);
    }
  }
  
  unshift_ptrs(pos_vec, n_samples);
}



//=========================================================
// Import (sample, otu, value) triplets given in any order.
// Indices are 0-based (Matrix) or 1-based (slam).
//=========================================================

static void import_triplet (ecomatrix_t *em, int *sam_in, int *otu_in, SEXP sexp_vals, int base) {
  
  int n_samples = em->n_samples;
  int n_otus    = em->n_otus;
  int nnz       = em->nnz;
  
  double *dbl_in;
  int    *int_in;
  sparse_val_ptrs(sexp_vals, &dbl_in, &int_in);
  
  check_indices(sam_in, nnz, n_samples, base);
  check_indices(otu_in, nnz, n_otus,    base);
  
  
  // Note how the input is already ordered.
  // --------------------------------------------
  
  int sam_sorted = 1; // by sample, then otu
  int otu_sorted = 1; // by otu
  
  for (int i = 1; i < nnz; i++) {
    int prev_sam = sam_in[i - 1], sam = sam_in[i];
    int prev_otu = otu_in[i - 1], otu = otu_in[i];
    if (prev_otu > otu) {
      otu_sorted = 0;
      if (prev_sam == sam) sam_sorted = 0;
    }
    if (prev_sam > sam) sam_sorted = 0;
  }
  
  
  // Already in (sample, otu) order.
  // --------------------------------------------
  
  if (sam_sorted) {
    
    count_samples(em, sam_in, base);
    
    if (base == 0) {
      em->otu_vec = otu_in;
    }
    else {
      int *otu_vec = rw_otu_vec(em);
      for (int i = 0; i < nnz; i++)
        otu_vec[i] = otu_in[i] - base;
    }
    
    assign_sparse_vals(em, sexp_vals);
    return;
  }
  
  
  // Where each otu's entries start, once grouped by otu.
  // --------------------------------------------
  
  int *col_ptr = (int*) safe_malloc(((size_t) n_otus + 1) * sizeof(int));
  memset(col_ptr, 0, ((size_t) n_otus + 1) * sizeof(int));
  
  for (int i = 0; i < nnz; i++)
    col_ptr[otu_in[i] - base + 1]++;
  
  for (int otu = 0; otu < n_otus; otu++)
    col_ptr[otu + 1] += col_ptr[otu];
  
  
  // Already grouped by otu, in order. Matrix and slam
  // build triplets this way (column-major).
  // --------------------------------------------
  
  if (otu_sorted) {
    scatter_by_sample(em, sam_in, base, col_ptr, dbl_in, int_in);
  }
  
  
  // Otherwise, a stable counting sort by otu moves the
  // (sample, value) pairs into column-compressed order.
  // --------------------------------------------
  
  else {
    
    int    *csc_sam = (int*)    safe_malloc((size_t) nnz * sizeof(int));
    double *csc_val = (double*) safe_malloc((size_t) nnz * sizeof(double));
    
    for (int i = 0; i < nnz; i++) {
      int dst = col_ptr[otu_in[i] - base]++;
      csc_sam[dst] = sam_in[i] - base;
      csc_val[dst] = dbl_in ? dbl_in[i] : INT_TO_DBL(int_in[i]);
    }
    unshift_ptrs(col_ptr, n_otus);
    
    scatter_by_sample(em, csc_sam, 0, col_ptr, csc_val, NULL);
    
    free_one(csc_sam);
    free_one(csc_val);
  }
  
  free_one(col_ptr);
}



//=========================================================
// Import a column-compressed matrix whose rows are samples.
//=========================================================

// Each entry must belong to exactly one column: col_ptr has
// n_cols + 1 elements, from 0 up to nnz, never decreasing.
static void check_col_ptr (SEXP sexp_ptr, int n_cols, int nnz) {

  if (TYPEOF(sexp_ptr) != INTSXP || XLENGTH(sexp_ptr) != (R_xlen_t) n_cols + 1)
    import_error("Sparse matrix has invalid column pointers.");

  int *col_ptr = INTEGER(sexp_ptr);
  int  bad_ptr = col_ptr[0] != 0 || col_ptr[n_cols] != nnz;
  for (int col = 0; col < n_cols && !bad_ptr; col++)
    if (col_ptr[col] > col_ptr[col + 1]) bad_ptr = 1;
  if (bad_ptr)
    import_error("Sparse matrix has invalid column pointers.");
}


static void import_csc_by_row (ecomatrix_t *em, int *row_in, int *col_ptr, SEXP sexp_vals) {

  int n_samples = em->n_samples;
  int nnz       = em->nnz;

  double *dbl_in;
  int    *int_in;
  sparse_val_ptrs(sexp_vals, &dbl_in, &int_in);

  check_indices(row_in, nnz, n_samples, 0);

  scatter_by_sample(em, row_in, 0, col_ptr, dbl_in, int_in);
}



//=========================================================
// Explicitly stored zeros are dropped, and repeated
// (sample, otu) entries summed, as the Matrix package
// defines them, so that sparse input gives the same
// results as dense input. Each sample's entries are in
// order of otu, so repeats are adjacent.
//=========================================================

static void compact_entries (ecomatrix_t *em) {

  int     n_samples = em->n_samples;
  int    *pos_vec   = em->pos_vec;
  int    *otu_vec   = em->otu_vec;
  double *val_vec   = em->val_vec;

  // Most inputs need nothing done: look before copying.
  int needed = 0;
  for (int sam = 0; sam < n_samples && !needed; sam++)
    for (int i = pos_vec[sam]; i < pos_vec[sam + 1]; i++)
      if (val_vec[i] == 0 || (i > pos_vec[sam] && otu_vec[i] == otu_vec[i - 1])) {
        needed = 1;
        break;
      }
  if (!needed) return;

  pos_vec = rw_pos_vec(em);
  otu_vec = rw_otu_vec(em);
  val_vec = rw_val_vec(em);

  int dst = 0;
  for (int sam = 0; sam < n_samples; sam++) {

    int start = dst;
    int end   = pos_vec[sam + 1];

    for (int i = pos_vec[sam]; i < end; i++) {
      if (dst > start && otu_vec[dst - 1] == otu_vec[i]) {
        val_vec[dst - 1] += val_vec[i];
      }
      else {
        otu_vec[dst] = otu_vec[i];
        val_vec[dst] = val_vec[i];
        dst++;
      }
    }

    // Drop zeros, including sums that cancel to zero.
    int keep = start;
    for (int i = start; i < dst; i++) {
      if (val_vec[i] != 0) {
        otu_vec[keep] = otu_vec[i];
        val_vec[keep] = val_vec[i];
        keep++;
      }
    }
    dst = keep;

    pos_vec[sam] = start;
  }

  pos_vec[n_samples] = dst;
  em->nnz            = dst;
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
  
  SEXP sexp_dge_x    = PROTECT(R_do_slot(sexp_dgeMatrix, install("x")));
  SEXP sexp_dim      = PROTECT(R_do_slot(sexp_dgeMatrix, install("Dim")));
  SEXP sexp_dimnames = PROTECT(R_do_slot(sexp_dgeMatrix, install("Dimnames")));
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
  UNPROTECT(3);
}



//=========================================================
// A simple_triplet_matrix from the slam package.
// Can NOT assume that i and j are sorted.
//=========================================================

static void parse_slam (ecomatrix_t *em, SEXP sexp_slam_mtx, int margin) {
  
  SEXP sexp_slam_i   = PROTECT(get(sexp_slam_mtx, "i"));
  SEXP sexp_slam_j   = PROTECT(get(sexp_slam_mtx, "j"));
  SEXP sexp_slam_v   = PROTECT(get(sexp_slam_mtx, "v"));
  SEXP sexp_dimnames = PROTECT(get(sexp_slam_mtx, "dimnames"));
  int  n_rows        = asInteger(get(sexp_slam_mtx, "nrow"));
  int  n_cols        = asInteger(get(sexp_slam_mtx, "ncol"));
  int  nnz           = length(sexp_slam_v);
  int *sam_in, *otu_in;

  if (n_rows == NA_INTEGER || n_cols == NA_INTEGER || n_rows < 0 || n_cols < 0)
    import_error("Sparse matrix is invalid: bad `nrow` or `ncol`.");

  int *i_in = index_vec(sexp_slam_i, nnz, "i");
  int *j_in = index_vec(sexp_slam_j, nnz, "j");

  em->nnz = nnz;


  // Accept samples in either rows or columns.
  //-------------------------------------------------------

  if (margin == 1) {

    em->n_samples = n_rows; // samples are in rows
    em->n_otus    = n_cols; // OTUs are in columns
    sam_in        = i_in;
    otu_in        = j_in;

    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = PROTECT(VECTOR_ELT(sexp_dimnames, 0));
  }

  else {

    em->n_samples = n_cols; // samples are in columns
    em->n_otus    = n_rows; // OTUs are in rows
    sam_in        = j_in;
    otu_in        = i_in;

    if (!isNull(sexp_dimnames))
      em->sexp_sample_names = PROTECT(VECTOR_ELT(sexp_dimnames, 1));
  }


  // Slam indices are 1-based.
  import_triplet(em, sam_in, otu_in, sexp_slam_v, 1);

  if (!isNull(sexp_dimnames)) UNPROTECT(1);
  UNPROTECT(4);
}




//=========================================================
// An sparse triplet matrix from the Matrix package.
// Can NOT assume that i and j are sorted.
//=========================================================

static void parse_dgTMatrix (ecomatrix_t *em, SEXP sexp_dgTMatrix, int margin) {
  
  SEXP sexp_dgt_i    = PROTECT(R_do_slot(sexp_dgTMatrix, install("i")));
  SEXP sexp_dgt_j    = PROTECT(R_do_slot(sexp_dgTMatrix, install("j")));
  SEXP sexp_dgt_x    = PROTECT(R_do_slot(sexp_dgTMatrix, install("x")));
  SEXP sexp_dim      = PROTECT(R_do_slot(sexp_dgTMatrix, install("Dim")));
  SEXP sexp_dimnames = PROTECT(R_do_slot(sexp_dgTMatrix, install("Dimnames")));
  int  n_rows        = INTEGER(sexp_dim)[0];
  int  n_cols        = INTEGER(sexp_dim)[1];
  int  nnz           = length(sexp_dgt_x);
  int *sam_in, *otu_in;

  int *i_in = index_vec(sexp_dgt_i, nnz, "i");
  int *j_in = index_vec(sexp_dgt_j, nnz, "j");

  em->nnz = nnz;


  // Transpose if needed
  if (margin == 1) {

    em->n_samples         = n_rows; // samples are in rows
    em->n_otus            = n_cols; // OTUs are in columns
    sam_in                = i_in;
    otu_in                = j_in;
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);
  }

  else { // margin == 2

    em->n_samples         = n_cols; // samples are in columns
    em->n_otus            = n_rows; // OTUs are in rows
    sam_in                = j_in;
    otu_in                = i_in;
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);
  }


  import_triplet(em, sam_in, otu_in, sexp_dgt_x, 0);
  UNPROTECT(5);
} 



//=========================================================
// A compressed sparse matrix from the Matrix package.
// CAN assume that i is sorted.
//=========================================================

static void parse_dgCMatrix (ecomatrix_t *em, SEXP sexp_dgCMatrix, int margin) {
  
  SEXP sexp_dgc_i    = PROTECT(R_do_slot(sexp_dgCMatrix, install("i")));
  SEXP sexp_dgc_p    = PROTECT(R_do_slot(sexp_dgCMatrix, install("p")));
  SEXP sexp_dgc_x    = PROTECT(R_do_slot(sexp_dgCMatrix, install("x")));
  SEXP sexp_dim      = PROTECT(R_do_slot(sexp_dgCMatrix, install("Dim")));
  SEXP sexp_dimnames = PROTECT(R_do_slot(sexp_dgCMatrix, install("Dimnames")));
  int  n_rows        = INTEGER(sexp_dim)[0];
  int  n_cols        = INTEGER(sexp_dim)[1];
  int  nnz           = length(sexp_dgc_x);

  // Pointers and indices become memory offsets: check them all.
  int *i_in = index_vec(sexp_dgc_i, nnz, "i");
  check_col_ptr(sexp_dgc_p, n_cols, nnz);

  em->nnz = nnz;


  // Transpose if needed
  if (margin == 1) {

    em->n_samples         = n_rows; // samples are in rows
    em->n_otus            = n_cols; // OTUs are in columns
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 0);

    import_csc_by_row(em, i_in, INTEGER(sexp_dgc_p), sexp_dgc_x);
  }

  else { // margin == 2

    check_indices(i_in, nnz, n_rows, 0);

    em->n_samples         = n_cols; // samples are in columns
    em->n_otus            = n_rows; // OTUs are in rows
    em->pos_vec           = INTEGER(sexp_dgc_p);
    em->otu_vec           = i_in;
    em->sexp_sample_names = VECTOR_ELT(sexp_dimnames, 1);

    assign_sparse_vals(em, sexp_dgc_x);
  }
  
  UNPROTECT(5);
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
  else    { error("Unrecognized matrix format."); } // # nocov
  
  
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

  // Dense input has neither zeros nor repeats among its entries.
  if (parse_func != parse_matrix && parse_func != parse_dgeMatrix)
    compact_entries(em);

  return em;
}

