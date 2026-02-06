// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"

//======================================================
// R interface. Returns if pthreads are available.
// Used by R code for `ecodive::n_cpus()` logic.
//======================================================
SEXP C_pthreads(void) {
  #ifdef HAVE_PTHREAD
    return ScalarLogical(1);
  #else
    return ScalarLogical(0);
  #endif
}

//======================================================
// run_parallel
// 
// Orchestrates the creation, execution, and cleanup
// of POSIX threads. Handles graceful degradation to
// single-threaded mode on failure.
//======================================================
void run_parallel(pthread_func_t func, int n_threads, int n_tasks) {
  
  // ---------------------------------------------------
  // Path A: Multithreading Attempt
  // ---------------------------------------------------
  #ifdef HAVE_PTHREAD
    // Only attempt threading if explicitly requested and compiled in.
    // The n_tasks check is a heuristic to avoid overhead for small workloads.
    if (n_threads > 1 && n_tasks >= 100) {
      
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      worker_t  *args = calloc(n_threads, sizeof(worker_t));
      
      // If out of memory, fall back to serial execution.
      if (!tids || !args) {
        free(tids); free(args);      // # nocov
        goto single_thread_fallback; // # nocov
      }
      
      // Initialize worker arguments.
      // Each thread gets its index (i) and total count (n).
      for (int i = 0; i < n_threads; i++) {
        args[i] = (worker_t){ .i = i, .n = n_threads };
      }
      
      // Create Threads.
      // We track 'created_count' to handle cases where the OS refuses
      // to create new threads (e.g., ulimit reached) halfway through the loop.
      int created_count = 0;
      for (int i = 0; i < n_threads; i++) {
        // pthread_create returns 0 on success. Non-zero is an error.
        if (pthread_create(&tids[i], NULL, func, &args[i]) != 0) {
            break; // Stop creating new threads  # nocov
        }
        created_count++;
      }
      
      // Only join threads that were actually created.
      for (int i = 0; i < created_count; i++) {
        pthread_join(tids[i], NULL);
      }
      
      // Clean up system heap memory explicitly.
      free(tids);
      free(args);
      
      
      // If we failed partially (created_count < n_threads), we fall through
      // to the fallback to ensure the job completes (by re-running it serially).
      if (created_count == n_threads) return;
    }
  #endif
  
  // ---------------------------------------------------
  // Path B: Single Thread Fallback
  // ---------------------------------------------------
  // This runs if:
  // 1. Pthreads are not supported on this OS.
  // 2. n_threads <= 1.
  // 3. n_tasks < 100.
  // 4. Memory allocation failed.
  // 5. Thread creation failed (partial or total).
  
  single_thread_fallback:
    
  // Arguments for processing the whole dataset
  worker_t args = { .i = 0, .n = 1 };
  func((void*)&args);
}
