// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"


//======================================================
// R interface. Returns if pthreads are available.
//======================================================
SEXP C_pthreads(void) {
#ifdef HAVE_PTHREAD
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}
