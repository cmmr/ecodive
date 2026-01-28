// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include "ecodive.h"


extern SEXP C_alpha_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pthreads(void);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 6},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  7},
  {"C_pthreads",  (DL_FUNC) &C_pthreads,  0},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    6},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_unifrac",   (DL_FUNC) &C_unifrac,   7},
  {NULL, NULL, 0}
};


// # nocov start
void R_init_ecodive(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
// # nocov end
