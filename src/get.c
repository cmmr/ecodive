// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>

//======================================================
// Extract list elements by name.
//======================================================
SEXP get(SEXP sexp_list, const char *str) {
  
  SEXP names = getAttrib(sexp_list, R_NamesSymbol);
  
  for (int i = 0; i < length(sexp_list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      return VECTOR_ELT(sexp_list, i);
    }
  }
  
  return R_NilValue; // # nocov
}

//======================================================
// Set the value for a named list element.
//======================================================
void set(SEXP sexp_list, const char *str, SEXP sexp_value) {
  
  SEXP names = getAttrib(sexp_list, R_NamesSymbol);
  
  for (int i = 0; i < length(sexp_list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      SET_VECTOR_ELT(sexp_list, i, sexp_value);
      return;
    }
  }
  
  error("Element does not exist: %s\n", str); // # nocov
}
