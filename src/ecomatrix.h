// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#ifndef ECOMATRIX_H_INCLUDED
#define ECOMATRIX_H_INCLUDED

#include <R.h>
#include <Rinternals.h>

// ecomatrix data structure
typedef struct {
  int     n_samples;
  int     n_otus;
  int     nnz;
  int    *sam_vec;
  int    *pos_vec;
  int    *otu_vec;
  double *val_vec;
  double *clr_vec;
  SEXP    sexp_sample_names;
} ecomatrix_t;


// ecomatrix constructor
ecomatrix_t* new_ecomatrix(SEXP sexp_matrix, SEXP sexp_margin);

// avoid overwriting R object data
double* rw_val_vec(ecomatrix_t *em);
double* rw_clr_vec(ecomatrix_t *em);


#endif
