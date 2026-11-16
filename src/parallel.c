// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

// Before R's headers, so that R's macros can't reach the
// system headers.
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h> // GetModuleHandleA, GetProcAddress
#else
#  include <dlfcn.h>   // dlopen, dlsym, dlclose, dlerror
#  include <unistd.h>  // getpid
#endif

#include "ecodive.h"

#include <limits.h> // INT_MAX
#include <stdio.h>  // sscanf

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
// BLAS thread count
//
// R's BLAS/LAPACK library may run threads of its own.
// C_blas_threads() reads that library's thread count,
// and optionally sets it, by calling the library's own
// functions. They are found by name among the
// libraries already loaded into R's process:
//
//   Unix:    The process's global symbols. R's BLAS is
//            a dependency of R itself, so it is there.
//   Windows: The exports of Rblas.dll and Rlapack.dll,
//            which users may replace with copies of an
//            optimized library. When both are copies,
//            each has a thread count of its own.
//
// Nothing is loaded that wasn't loaded before, and
// OpenMP's settings are left alone.
//
// Setting the count is optional, so it is only done
// where it is known to be safe:
//   - Only libraries whose builds honor the count (see
//     openblas_can_set()).
//   - Never above the highest count read, so that the
//     library doesn't have to start threads it didn't
//     already have. OpenBLAS, for one, doesn't check
//     that new threads were started. With no count
//     read, only 1.
//   - Only in the process that loaded ecodive. In a
//     child forked from it (e.g. by mclapply()),
//     OpenBLAS restarts its threads when its count is
//     set, and stops the process if that fails.
//======================================================

typedef void    (*sym_t)(void); // any function
typedef int     (*get_int_t)(void);
typedef void    (*set_int_t)(int);
typedef int64_t (*get_i64_t)(void);
typedef void    (*set_i64_t)(int64_t);
typedef char*   (*get_str_t)(void);

// Argument and return type of a library's functions.
#define BLAS_INT 0 // int
#define BLAS_I64 1 // int64_t (BLIS's dim_t)

// Which function names to look for, in order.
#define API_INT  0 // int get(void);   void set(int)
#define API_FLEX 1 // as API_INT, but see flexiblas_can_set()
#define API_OPEN 2 // as API_INT, but see openblas_can_set()
#define API_BLIS 3 // dim_t get(void); void set(dim_t)

static const struct {
  const char *get;
  const char *set;
  int         api;
} blas_apis[] = {

  // FlexiBLAS goes first: it also defines the names of
  // the libraries below, forwarding them to whichever
  // backend it has loaded.
  { "flexiblas_get_num_threads",  "flexiblas_set_num_threads",  API_FLEX },
  { "openblas_get_num_threads",   "openblas_set_num_threads",   API_OPEN },

  // MKL's lowercase mkl_get_max_threads() and
  // mkl_set_num_threads() are macros for these. As
  // symbols, the lowercase names can be Fortran
  // interfaces that take pointers.
  { "MKL_Get_Max_Threads",        "MKL_Set_Num_Threads",        API_INT  },
  { "bli_thread_get_num_threads", "bli_thread_set_num_threads", API_BLIS }
};

typedef struct {
  sym_t get;
  sym_t set;  // NULL when the count is only read.
  int   type;
  int   max;  // Highest count read; 0 if none yet.
} blas_lib_t;

typedef int (*get_name_t)(char *, size_t);

// FlexiBLAS's flexiblas_current_backend(), and the
// backend it named when FlexiBLAS was found.
static get_name_t flexiblas_backend_fn = NULL;
static char       flexiblas_backend[64];

// The libraries found, with each one's functions.
static blas_lib_t blas_libs[2];
static int        n_blas_libs = 0;

#ifndef _WIN32
  // The process that loaded ecodive.
  static pid_t blas_pid = 0;
#endif


//------------------------------------------------------
// Called when ecodive is loaded.
//------------------------------------------------------
void init_blas_threads(void) {
  #ifndef _WIN32
    blas_pid = getpid();
  #endif
}


//------------------------------------------------------
// Whether counts may be set from this process.
//------------------------------------------------------
static int blas_can_set_here(void) {
  #ifdef _WIN32
    return 1; // No fork() on Windows.
  #else
    return blas_pid != 0 && getpid() == blas_pid;
  #endif
}


//------------------------------------------------------
// Look up a function by name, or NULL if not found.
//------------------------------------------------------
static sym_t find_sym(void *handle, const char *name) {

  sym_t fn = NULL;

  #ifdef _WIN32
    FARPROC p = GetProcAddress((HMODULE) handle, name);
    if (p) fn = (sym_t) p;
  #else
    // POSIX guarantees that a function's address fits
    // in dlsym()'s void*; ISO C has no cast between the
    // two, so copy the bytes.
    void *p = dlsym(handle, name);
    if (p && sizeof(p) == sizeof(fn)) memcpy(&fn, &p, sizeof(fn));
  #endif

  return fn;
}


//------------------------------------------------------
// OpenBLAS's OpenMP builds before 0.3.34 take their
// count from OpenMP's:
//   - Up to 0.3.21, setting the count also calls
//     omp_set_num_threads(), changing how many threads
//     every later OpenMP region started from R's thread
//     uses.
//   - From 0.3.22 to 0.3.33, the next BLAS call resets
//     the count to omp_get_max_threads().
// Those builds' count is only read.
//------------------------------------------------------
static int openblas_can_set(void *handle) {

  // Very old versions lack it, so their build is unknown.
  sym_t parallel = find_sym(handle, "openblas_get_parallel");
  if (!parallel) return 0;

  // 0 = sequential (always one thread), 1 = pthreads,
  // 2 = OpenMP.
  switch (((get_int_t) parallel)()) {
    case 0: return 0;
    case 1: return 1;
  }

  // Starts with e.g. "OpenBLAS 0.3.34 " in recent versions.
  sym_t config = find_sym(handle, "openblas_get_config");
  const char *str = config ? ((get_str_t) config)() : NULL;

  int major, minor, patch;
  if (!str || sscanf(str, "OpenBLAS %d.%d.%d", &major, &minor, &patch) != 3)
    return 0;

  return major > 0 || minor > 3 || (minor == 3 && patch >= 34);
}


//------------------------------------------------------
// FlexiBLAS passes the count on to its backend. Only
// backends known to take it are set: OpenBLAS's builds
// other than OpenMP ones (whose version can't be
// checked here; see openblas_can_set()), BLIS, and MKL.
// They are known by FlexiBLAS's usual names for them;
// any other name, such as a library's path, is only
// read.
//------------------------------------------------------
static int flexiblas_name_ok(const char *name) {

  char up[sizeof(flexiblas_backend)];
  size_t i;

  for (i = 0; name[i] && i < sizeof(up) - 1; i++) {
    char c = name[i];
    if (c >= 'a' && c <= 'z') c = (char) (c - 'a' + 'A');
    if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-' || c == '_'))
      return 0;
    up[i] = c;
  }
  if (i == 0 || name[i]) return 0;
  up[i] = '\0';

  if (strstr(up, "OPENBLAS")) return !strstr(up, "OPENMP");
  return strstr(up, "BLIS") || strstr(up, "MKL");
}


// The name of FlexiBLAS's current backend, or "".
static void flexiblas_name(char *buf, size_t len) {
  buf[0] = '\0';
  if (flexiblas_backend_fn) {
    flexiblas_backend_fn(buf, len);
    buf[len - 1] = '\0';
  }
}


static int flexiblas_can_set(void *handle) {

  sym_t fn = find_sym(handle, "flexiblas_current_backend");
  if (!fn) return 0;

  flexiblas_backend_fn = (get_name_t) fn;
  flexiblas_name(flexiblas_backend, sizeof(flexiblas_backend));

  return flexiblas_name_ok(flexiblas_backend);
}


//------------------------------------------------------
// Whether FlexiBLAS, if it is the library, still uses
// the backend it did when found. It can be switched
// while running, e.g. with the flexiblas R package.
//------------------------------------------------------
static int flexiblas_unswitched(void) {

  if (!flexiblas_backend_fn) return 1;

  char now[sizeof(flexiblas_backend)];
  flexiblas_name(now, sizeof(now));

  return strcmp(now, flexiblas_backend) == 0;
}


//------------------------------------------------------
// BLIS's dim_t is 64 bits unless BLIS was configured
// with --int-size=32. Returns the BLAS_* type to call
// its functions with, or -1 if neither.
//------------------------------------------------------
static int blis_int_type(void *handle) {

  sym_t bits = find_sym(handle, "bli_info_get_int_type_size");
  if (!bits) return -1;

  // It returns 32 or 64 as a dim_t. Both fit in the
  // lower 32 bits, the part read as an int.
  switch (((get_int_t) bits)()) {
    case 32: return BLAS_INT;
    case 64: return BLAS_I64;
  }

  return -1;
}


//------------------------------------------------------
// Add the thread functions exported by `handle`'s
// library, if any, to blas_libs.
//------------------------------------------------------
static void add_blas_lib(void *handle) {

  int n_apis = (int) (sizeof(blas_apis) / sizeof(blas_apis[0]));

  for (int i = 0; i < n_apis; i++) {

    sym_t get = find_sym(handle, blas_apis[i].get);
    if (!get) continue;

    blas_lib_t lib = { get, find_sym(handle, blas_apis[i].set), BLAS_INT, 0 };

    if (blas_apis[i].api == API_FLEX && !flexiblas_can_set(handle))
      lib.set = NULL;

    if (blas_apis[i].api == API_OPEN && !openblas_can_set(handle))
      lib.set = NULL;

    if (blas_apis[i].api == API_BLIS) {
      lib.type = blis_int_type(handle);
      if (lib.type < 0) return;
    }

    // The same library can be reachable twice, e.g. when
    // Rlapack.dll forwards to Rblas.dll.
    for (int j = 0; j < n_blas_libs; j++)
      if (blas_libs[j].get == lib.get) return;

    if (n_blas_libs < 2) blas_libs[n_blas_libs++] = lib;

    return;
  }
}


//------------------------------------------------------
// Find R's BLAS/LAPACK libraries' thread functions, the
// first time it is called. On Windows, Rlapack.dll is
// only loaded when R first needs LAPACK, so it is
// looked for on each call until it is found.
//------------------------------------------------------
static void find_blas_libs(void) {

  #ifdef _WIN32

    static const char *dlls[]  = { "Rblas.dll", "Rlapack.dll" };
    static int         found[] = { 0, 0 };

    for (int i = 0; i < 2; i++) {
      if (found[i]) continue;
      HMODULE handle = GetModuleHandleA(dlls[i]);
      if (!handle) continue;
      found[i] = 1;
      add_blas_lib((void *) handle);
    }

  #else

    static int searched = 0;
    if (searched) return;
    searched = 1;

    // The libraries found stay loaded after dlclose():
    // R itself depends on them.
    void *handle = dlopen(NULL, RTLD_LAZY);
    if (handle) {
      add_blas_lib(handle);
      dlclose(handle);
    }

    dlerror(); // Don't leave failed lookups' messages.

  #endif
}


//------------------------------------------------------
// Read a library's count, or NA, and keep its highest.
//------------------------------------------------------
static int blas_get(blas_lib_t *lib) {

  int64_t n = (lib->type == BLAS_I64)
    ? ((get_i64_t) lib->get)()
    : ((get_int_t) lib->get)();

  // BLIS returns -1 when its count was never set.
  if (n < 1 || n > INT_MAX) return NA_INTEGER;

  if (n > lib->max) lib->max = (int) n;

  return (int) n;
}


//------------------------------------------------------
// Set a library's count to `n`, but never above the
// highest count read, and only to 1 if none was read.
//------------------------------------------------------
static void blas_set(const blas_lib_t *lib, int n) {

  if (!lib->set || !blas_can_set_here() || !flexiblas_unswitched()) return;

  if (lib->max == 0 && n != 1) return;
  if (lib->max >  0 && n > lib->max) n = lib->max;

  if (lib->type == BLAS_I64) { ((set_i64_t) lib->set)((int64_t) n); }
  else                       { ((set_int_t) lib->set)(n);           }
}


//======================================================
// R interface. Returns the BLAS thread count, or NA if
// it can't be read. If there is more than one library,
// returns the largest count. If `sexp_n` is an integer
// >= 1, then sets each library's count to it afterward,
// within blas_set()'s limits.
// Used by R code for `ecodive::blas_threads()`.
//======================================================
SEXP C_blas_threads(SEXP sexp_n) {

  find_blas_libs();

  int prev = NA_INTEGER;
  for (int i = 0; i < n_blas_libs; i++) {
    int n = blas_get(&blas_libs[i]);
    if (n != NA_INTEGER && (prev == NA_INTEGER || n > prev)) prev = n;
  }

  if (TYPEOF(sexp_n) == INTSXP && XLENGTH(sexp_n) == 1) {
    int n = INTEGER(sexp_n)[0];
    if (n != NA_INTEGER && n >= 1)
      for (int i = 0; i < n_blas_libs; i++)
        blas_set(&blas_libs[i], n);
  }

  return ScalarInteger(prev);
}


//======================================================
// R interface, for tests. Returns whether any BLAS
// library found has a count that C_blas_threads() can
// set from this process.
//======================================================
SEXP C_blas_settable(void) {

  find_blas_libs();

  for (int i = 0; i < n_blas_libs; i++)
    if (blas_libs[i].set && blas_can_set_here() && flexiblas_unswitched())
      return ScalarLogical(1);

  return ScalarLogical(0);
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
  {
    // Arguments for processing the whole dataset
    worker_t args = { .i = 0, .n = 1 };
    func((void*)&args);
  }
}
