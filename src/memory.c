// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * All the active temporary memory allocations are tracked here.
 * If there is not enough free memory for a new allocation, then 
 * all the active allocations are freed and an R error is thrown 
 * with an appropriate error message.
 * 
 * init_n_ptrs() should be called at the beginning of .Call()-ed
 * C functions, and set to the number of safe_mallocs needed at
 * any one time during that function's run (including safe_mallocs 
 * by the functions it calls).
 * 
 * safe_malloc() should be called in place of malloc().
 * 
 * free_one() is optional, and can be called in place of free();
 * 
 * free_all() should be called at the end of .Call()-ed
 * C functions.
 */


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // free, malloc
#include "memory.h"

static void **ptr_vec;
static int    n_ptrs = 0;


void init_n_ptrs (int n) {
  
  if (n_ptrs || ptr_vec) {
    free_all(); // # nocov
  }
  
  n_ptrs  = n;
  ptr_vec = malloc(sizeof(void*) * n);
  if (!ptr_vec) {
    error("Insufficient memory."); // # nocov
  }
  
  for (int i = 0; i < n_ptrs; i++) {
    ptr_vec[i] = NULL;
  }
}


void* safe_malloc (size_t n_bytes) {
  
  // Find an open slot in ptr_vec
  int i = 0;
  for (; i < n_ptrs; i++) {
    if (!ptr_vec[i]) {
      break;
    }
  }
  
  // No slots available
  if (i == n_ptrs) {
    free_all(); // # nocov
    error("Insufficient n_ptrs"); // # nocov
  }
  
  // Allocate the requested memory
  ptr_vec[i] = malloc(n_bytes);
  if (!ptr_vec[i]) {
    free_all(); // # nocov
    error("Insufficient memory."); // # nocov
  }
  
  return ptr_vec[i];
}



void free_all () {
  
  if (!ptr_vec) return;
  
  for (int i = 0; i < n_ptrs; i++) {
    if (ptr_vec[i]) {
      free(ptr_vec[i]);
    }
  }
  
  free(ptr_vec);
  
  n_ptrs  = 0;
  ptr_vec = NULL;
}



int is_safe_ptr (void *ptr) {
  
  // Pointer already freed
  if (!ptr) return 0;
  
  // Find the matching slot in ptr_vec
  for (int i = 0; i < n_ptrs; i++) {
    if (ptr_vec[i] == ptr) {
      return 1;
    }
  }
  
  return 0; // # nocov
}



void* free_one (void *ptr) {
  
  // Pointer already freed
  if (!ptr) return NULL;
  
  // Find the matching slot in ptr_vec
  int i = 0;
  for (; i < n_ptrs; i++) {
    if (ptr_vec[i] == ptr) {
      break;
    }
  }
  
  // ptr is not tracked here
  if (i == n_ptrs) {
    free_all(); // # nocov
    error("free_one() cannot free that pointer"); // # nocov
  }
  
  // Free the given memory address
  free(ptr_vec[i]);
  ptr_vec[i] = NULL;
  
  return NULL;
}



void* maybe_free_one (void *ptr) {
  
  if (!ptr) return NULL;
  
  if (is_safe_ptr(ptr)) {
    return free_one(ptr);
  }
  
  return NULL; // # nocov
}


