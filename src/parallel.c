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
#  include <signal.h>  // pthread_sigmask, sigfillset
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
// The functions from here to the end marker only run
// when an optimized BLAS library is found. The coverage
// runner has R's reference BLAS, so they are left out
// of coverage. test-blas_threads.r tests them wherever
// such a library is present.
//------------------------------------------------------
// # nocov start


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
// # nocov end


//------------------------------------------------------
// Add the thread functions exported by `handle`'s
// library, if any, to blas_libs.
//------------------------------------------------------
static void add_blas_lib(void *handle) {

  int n_apis = (int) (sizeof(blas_apis) / sizeof(blas_apis[0]));

  for (int i = 0; i < n_apis; i++) {

    sym_t get = find_sym(handle, blas_apis[i].get);
    if (!get) continue;

    // # nocov start (only with an optimized BLAS; see above)
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
    // # nocov end
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
    int n = blas_get(&blas_libs[i]);                                   // # nocov
    if (n != NA_INTEGER && (prev == NA_INTEGER || n > prev)) prev = n; // # nocov
  }

  if (TYPEOF(sexp_n) == INTSXP && XLENGTH(sexp_n) == 1) {
    int n = INTEGER(sexp_n)[0];
    if (n != NA_INTEGER && n >= 1)
      for (int i = 0; i < n_blas_libs; i++)
        blas_set(&blas_libs[i], n); // # nocov
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

  int settable = 0;
  for (int i = 0; i < n_blas_libs; i++)
    if (blas_libs[i].set && blas_can_set_here() && flexiblas_unswitched()) // # nocov
      settable = 1;                                                      // # nocov

  return ScalarLogical(settable);
}

//======================================================
// run_parallel
//
// Runs func once for each share i in 0 .. n - 1 (n =
// n_threads), where each call does its own part (share i
// of n) of the job.
//
// Shares are claimed one at a time by whichever thread is
// free: the calling thread, and the threads of a pool that
// lives as long as the package is loaded, asleep between
// jobs. The calling thread waits only for shares that
// other threads have claimed, never for threads that
// haven't woken up yet. So on a machine whose cores are
// busy with other processes, a job takes at worst about as
// long as a serial run, rather than waiting for every
// thread to be scheduled (once for each of the hundreds of
// products in an iterative method), and threads aren't
// started and stopped for every job.
//
// So each share is done exactly once, by one thread, but
// which thread does it, and which shares run at the same
// time, vary. Workers may use their index i for scratch
// space of their own, but must not wait on each other (no
// barriers), and must not call R. Results are the same
// whatever the split, as long as each share's work depends
// only on (i, n).
//======================================================

#if defined(HAVE_PTHREAD) && (defined(__GNUC__) || defined(__clang__))
#define ECODIVE_POOL 1
#endif


#ifdef ECODIVE_POOL

static pthread_mutex_t pl_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  pl_wake  = PTHREAD_COND_INITIALIZER; // A new job, or stop.
static pthread_cond_t  pl_done  = PTHREAD_COND_INITIALIZER; // A job's last share is done.

static pthread_t      *pl_tids  = NULL; // The pool's threads.
static int             pl_size  = 0;
static int             pl_stop  = 0;

// The current job, set under pl_mutex.
static uint64_t        pl_gen   = 0;    // Job number.
static pthread_func_t  pl_func  = NULL;
static int             pl_n     = 0;    // Shares.

// Atomic: the job number (high 32 bits) and the next share
// to claim (low 32 bits), so that a thread still holding
// an old job can't claim a share of a newer one; and how
// many of the current job's shares aren't done yet.
static uint64_t        pl_claim = 0;
static int             pl_left  = 0;


// The calling thread's floating-point control state, which
// the pool's threads load before each job. POSIX threads
// inherit it when created, but Windows threads start with
// the system's default: x87 at 53-bit precision, where
// mingw's startup code sets R's thread to 64-bit. Math
// functions computed with x87 instructions (all of them on
// 32-bit Windows; exp(), log(), pow() in the msvcrt runtime
// of R < 4.2 on 64-bit) would then round differently on the
// pool's threads than with cpus = 1. The job's state is
// copied as a whole, so that changes R makes later (e.g. by
// another package) are followed too.
#if defined(__i386__) || defined(__x86_64__)
#define POOL_FPENV 1
typedef struct { unsigned short x87; unsigned int mxcsr; } fpenv_t;
static fpenv_t pl_fpenv;

static void fpenv_get (fpenv_t *e) {
  __asm__ __volatile__ ("fnstcw %0" : "=m" (e->x87));
  #if defined(__x86_64__) || defined(__SSE__)
    __asm__ __volatile__ ("stmxcsr %0" : "=m" (e->mxcsr));
  #else
    e->mxcsr = 0;
  #endif
}

static void fpenv_set (const fpenv_t *e) {
  __asm__ __volatile__ ("fldcw %0" : : "m" (e->x87));
  #if defined(__x86_64__) || defined(__SSE__)
    __asm__ __volatile__ ("ldmxcsr %0" : : "m" (e->mxcsr));
  #endif
}
#endif


static int claim_share (uint32_t gen, int n) {
  uint64_t v = __atomic_load_n(&pl_claim, __ATOMIC_ACQUIRE);
  for (;;) {
    if ((uint32_t) (v >> 32) != gen || (int) (uint32_t) v >= n) return -1;
    if (__atomic_compare_exchange_n(
        &pl_claim, &v, v + 1, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE ))
      return (int) (uint32_t) v;
  }
}


// Do shares of job `gen` until none are left to claim.
static void run_shares (uint32_t gen, pthread_func_t func, int n) {
  for (int i = claim_share(gen, n); i >= 0; i = claim_share(gen, n)) {
    worker_t w = { .i = i, .n = n };
    func(&w);
    if (__atomic_sub_fetch(&pl_left, 1, __ATOMIC_ACQ_REL) == 0) {
      pthread_mutex_lock(&pl_mutex);
      pthread_cond_broadcast(&pl_done);
      pthread_mutex_unlock(&pl_mutex);
    }
  }
}


static void *pool_thread (void *arg) {

  (void) arg;
  uint64_t seen = 0;

  pthread_mutex_lock(&pl_mutex);
  for (;;) {
    while (!pl_stop && pl_gen == seen)
      pthread_cond_wait(&pl_wake, &pl_mutex);
    if (pl_stop) break;
    seen = pl_gen;
    pthread_func_t func = pl_func;
    int            n    = pl_n;
    #ifdef POOL_FPENV
      fpenv_t env = pl_fpenv;
    #endif
    pthread_mutex_unlock(&pl_mutex);
    #ifdef POOL_FPENV
      fpenv_set(&env);
    #endif
    run_shares((uint32_t) seen, func, n);
    pthread_mutex_lock(&pl_mutex);
  }
  pthread_mutex_unlock(&pl_mutex);

  return NULL;
}


// A forked child (e.g., from parallel::mclapply()) has only
// the thread that forked, and copies of the pool's state:
// start over with an empty pool. Noticed by the process ID,
// as for BLAS above, rather than with pthread_atfork(),
// whose handler can't be removed if the DLL is unloaded.
#ifndef _WIN32
static pid_t pl_pid = 0; // The process the pool's threads belong to.

static void pool_check_fork (void) {
  if (pl_size > 0 && getpid() != pl_pid) {
    pthread_mutex_init(&pl_mutex, NULL);
    pthread_cond_init(&pl_wake, NULL);
    pthread_cond_init(&pl_done, NULL);
    pl_tids = NULL; // The parent's; not ours to free.
    pl_size = 0;
    pl_stop = 0;
  }
}
#endif


// Grow the pool to `size` threads, or as many as the OS
// allows. Main thread only, between jobs.
static void pool_grow (int size) {

  #ifndef _WIN32
    pool_check_fork();
  #endif

  if (size <= pl_size) return;

  pthread_t *tids = realloc(pl_tids, (size_t) size * sizeof(pthread_t));
  if (!tids) return; // # nocov
  pl_tids = tids;

  // The threads live as long as the package, so they block
  // every signal: process-wide signals such as SIGCHLD must
  // go to R's own thread, not run R's handlers on one of
  // these. (New threads inherit the creating thread's mask.)
  #ifndef _WIN32
    sigset_t all, old;
    sigfillset(&all);
    pthread_sigmask(SIG_SETMASK, &all, &old);
    pl_pid = getpid();
  #endif

  while (pl_size < size) {
    if (pthread_create(&pl_tids[pl_size], NULL, pool_thread, NULL) != 0)
      break; // # nocov: shares are done by the threads there are.
    pl_size++;
  }

  #ifndef _WIN32
    pthread_sigmask(SIG_SETMASK, &old, NULL);
  #endif
}


// Stop and join the pool's threads. Called from .onUnload(),
// before the DLL is unloaded: the threads run its code.
void stop_thread_pool (void) {

  #ifndef _WIN32
    pool_check_fork();
  #endif

  if (pl_size == 0) return;

  pthread_mutex_lock(&pl_mutex);
  pl_stop = 1;
  pthread_cond_broadcast(&pl_wake);
  pthread_mutex_unlock(&pl_mutex);

  for (int t = 0; t < pl_size; t++)
    pthread_join(pl_tids[t], NULL);

  free(pl_tids);
  pl_tids = NULL;
  pl_size = 0;
  pl_stop = 0;
}

#else

void stop_thread_pool (void) {}

#endif


void run_parallel(pthread_func_t func, int n_threads, int n_tasks) {

  #ifdef ECODIVE_POOL
    // Threads only pay off for enough work.
    if (n_threads > 1 && n_tasks >= 100) {

      pool_grow(n_threads - 1);

      if (pl_size > 0) {

        pthread_mutex_lock(&pl_mutex);
        uint32_t gen = (uint32_t) ++pl_gen;
        pl_func = func;
        pl_n    = n_threads;
        #ifdef POOL_FPENV
          fpenv_get(&pl_fpenv);
        #endif
        __atomic_store_n(&pl_left,  n_threads,              __ATOMIC_RELEASE);
        __atomic_store_n(&pl_claim, (uint64_t) gen << 32,   __ATOMIC_RELEASE);
        pthread_cond_broadcast(&pl_wake);
        pthread_mutex_unlock(&pl_mutex);

        run_shares(gen, func, n_threads);

        // Wait for shares that other threads claimed.
        pthread_mutex_lock(&pl_mutex);
        while (__atomic_load_n(&pl_left, __ATOMIC_ACQUIRE) > 0)
          pthread_cond_wait(&pl_done, &pl_mutex);
        pthread_mutex_unlock(&pl_mutex);

        return;
      }
    }
  #endif

  // ---------------------------------------------------
  // Serial: one thread, or not enough work, or no threads
  // could be started. The whole job as one share.
  // ---------------------------------------------------
  worker_t args = { .i = 0, .n = 1 };
  func((void*)&args);
}


//======================================================
// Whether the user has asked to interrupt (Esc, Ctrl+C).
// Unlike R_CheckUserInterrupt(), this returns instead of
// jumping out of the C code, so that the caller can free
// its malloc'd memory before raising an error. Main thread
// only.
//======================================================
static void check_interrupt_fn(void *dummy) {
  (void) dummy;
  R_CheckUserInterrupt();
}

int user_interrupted(void) {
  return R_ToplevelExec(check_interrupt_fn, NULL) == FALSE;
}


//======================================================
// R interface: stop the thread pool (from .onUnload()).
//======================================================
SEXP C_stop_thread_pool(void) {
  stop_thread_pool();
  return R_NilValue;
}
