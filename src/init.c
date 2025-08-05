// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP C_alpha_div(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_faith(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gower(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pthreads(void);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 4},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  6},
  {"C_faith",     (DL_FUNC) &C_faith,     4},
  {"C_gower",     (DL_FUNC) &C_gower,     5},
  {"C_pthreads",  (DL_FUNC) &C_pthreads,  0},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    5},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_unifrac",   (DL_FUNC) &C_unifrac,   7},
  {NULL, NULL, 0}
};

void R_init_fastbiom(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
