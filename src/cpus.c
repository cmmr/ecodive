// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

//======================================================
// CPU limits
//
// What the operating system lets this process run on,
// for `ecodive::n_cpus()`. R code adds the limits that
// options, environment variables and job schedulers set
// (R/cpus.r).
//
//   logical   Logical cores online.
//   physical  Physical cores online.
//   affinity  Cores the process may run on: Linux's
//             sched_getaffinity() (taskset, numactl,
//             cgroup cpusets, Docker's --cpuset-cpus,
//             Slurm's and other schedulers' binding) or
//             Windows' process affinity mask (start
//             /affinity, job objects).
//   quota     CPU time the process may use, in cores:
//             Linux cgroup v1/v2 quotas (Docker's --cpus,
//             Kubernetes limits, systemd's CPUQuota) or a
//             Windows job object's hard CPU rate cap.
//
// Each is NA where it can't be read. The core counts and
// the quota are read once per process, and again in a
// forked child; the affinity is read on every call, since
// it can be changed while R runs (e.g. by taskset -p).
//======================================================

// Before any system header: sched_getaffinity() and the
// CPU_*_S() macros are GNU extensions.
#if defined(__linux__) && !defined(_GNU_SOURCE)
#  define _GNU_SOURCE
#endif

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#else
#  include <unistd.h> // getpid, sysconf
#endif

#ifdef __linux__
#  include <sched.h>  // sched_getaffinity, CPU_ALLOC
#  include <errno.h>  // errno, EINVAL
#endif

#ifdef __APPLE__
#  include <sys/types.h>
#  include <sys/sysctl.h> // sysctlbyname
#endif

#include <limits.h> // INT_MAX
#include <stdio.h>  // fopen, fgets, snprintf, sscanf
#include <stdlib.h> // getenv, strtol, strtod

#include "ecodive.h"


// Values read once per process.
static int    ci_logical  = 0;   // 0 = not yet read
static int    ci_physical = 0;   // set with ci_logical
static int    ci_quota    = 0;

#ifndef _WIN32
  static pid_t ci_pid   = 0;     // The process they were read in.
  static pid_t load_pid = 0;     // The process that loaded ecodive.
#endif


//------------------------------------------------------
// Called when ecodive is loaded.
//------------------------------------------------------
void init_cpu_info(void) {
  #ifndef _WIN32
    load_pid = getpid();
  #endif
}


//------------------------------------------------------
// Linux
//------------------------------------------------------
#ifdef __linux__

// Read a file's first line into buf. Returns 0 on failure.
static int read_line(const char *path, char *buf, size_t len) {
  FILE *f = fopen(path, "r");
  if (!f) return 0;
  char *ok = fgets(buf, (int) len, f);
  fclose(f);
  if (!ok) return 0;
  buf[strcspn(buf, "\r\n")] = '\0';
  return 1;
}


// Physical cores: online CPUs that are the first of
// their core's siblings.
static int linux_physical(void) {

  long n_conf = sysconf(_SC_NPROCESSORS_CONF);
  if (n_conf < 1 || n_conf > 1 << 20) return NA_INTEGER;

  char path[128], buf[512];
  int  n = 0;

  for (long i = 0; i < n_conf; i++) {
    snprintf(path, sizeof(path),
      "/sys/devices/system/cpu/cpu%ld/topology/core_cpus_list", i);
    if (!read_line(path, buf, sizeof(buf))) {
      snprintf(path, sizeof(path),
        "/sys/devices/system/cpu/cpu%ld/topology/thread_siblings_list", i);
      if (!read_line(path, buf, sizeof(buf))) continue; // offline
    }
    if (strtol(buf, NULL, 10) == i) n++;
  }

  return n > 0 ? n : NA_INTEGER;
}


// Cores the process may run on.
static int linux_affinity(void) {

  // Machines can have more CPUs than cpu_set_t's 1024.
  for (int n_cpus = 1024; n_cpus <= 1 << 20; n_cpus *= 4) {

    cpu_set_t *set  = CPU_ALLOC(n_cpus);
    if (!set) return NA_INTEGER; // # nocov
    size_t     size = CPU_ALLOC_SIZE(n_cpus);
    CPU_ZERO_S(size, set);

    if (sched_getaffinity(0, size, set) == 0) {
      int n = CPU_COUNT_S(size, set);
      CPU_FREE(set);
      return n > 0 ? n : NA_INTEGER;
    }

    int err = errno;
    CPU_FREE(set);
    if (err != EINVAL) return NA_INTEGER; // # nocov
  }

  return NA_INTEGER; // # nocov
}


// The directory of the process's cgroup in a hierarchy:
// from /proc/self/cgroup, the line whose controllers
// include `ctrl` (or the v2 line, "0::path", when ctrl is
// NULL); from /proc/self/mountinfo, where that hierarchy
// is mounted. The mount point goes in `mnt`.
static int cgroup_dir(const char *ctrl, char *dir, char *mnt, size_t len) {

  char  line[4096], path[4096] = "";
  FILE *f;

  // Its path in the hierarchy.
  if (!(f = fopen("/proc/self/cgroup", "r"))) return 0;
  while (fgets(line, sizeof(line), f)) {
    line[strcspn(line, "\r\n")] = '\0';
    char *c1 = strchr(line, ':');
    char *c2 = c1 ? strchr(c1 + 1, ':') : NULL;
    if (!c2) continue;
    *c2 = '\0';
    const char *ctrls = c1 + 1;
    int hit = 0;
    if (!ctrl) {
      hit = (*ctrls == '\0');
    }
    else {
      char tmp[256];
      snprintf(tmp, sizeof(tmp), "%s", ctrls);
      for (char *t = strtok(tmp, ","); t; t = strtok(NULL, ","))
        if (strcmp(t, ctrl) == 0) hit = 1;
    }
    if (hit) { snprintf(path, sizeof(path), "%s", c2 + 1); break; }
  }
  fclose(f);
  if (!path[0]) return 0;

  // Where the hierarchy is mounted: mountinfo's fields are
  // "id parent dev root mountpoint opts ... - type src superopts".
  int found = 0;
  char root[4096] = "";
  if (!(f = fopen("/proc/self/mountinfo", "r"))) return 0;
  while (!found && fgets(line, sizeof(line), f)) {
    line[strcspn(line, "\r\n")] = '\0';
    char *sep = strstr(line, " - ");
    if (!sep) continue;
    *sep = '\0';
    char type[64], src[1024], super[1024] = "";
    if (sscanf(sep + 3, "%63s %1023s %1023s", type, src, super) < 2) continue;
    if (ctrl) {
      if (strcmp(type, "cgroup") != 0) continue;
      int hit = 0;
      for (char *t = strtok(super, ","); t; t = strtok(NULL, ","))
        if (strcmp(t, ctrl) == 0) hit = 1;
      if (!hit) continue;
    }
    else if (strcmp(type, "cgroup2") != 0) continue;
    char r[4096], m[4096];
    if (sscanf(line, "%*s %*s %*s %4095s %4095s", r, m) != 2) continue;
    snprintf(root, sizeof(root), "%s", r);
    snprintf(mnt,  len,          "%s", m);
    found = 1;
  }
  fclose(f);
  if (!found) return 0;

  // The path relative to the mount's root.
  const char *rel  = path;
  size_t      rlen = strlen(root);
  if (strcmp(root, "/") != 0 && strncmp(path, root, rlen) == 0 &&
      (path[rlen] == '/' || path[rlen] == '\0'))
    rel = path + rlen;
  if (strcmp(rel, "/") == 0) rel = "";

  if ((size_t) snprintf(dir, len, "%s%s", mnt, rel) >= len) return 0;
  return 1;
}


// Remove dir's last component. Returns 0 once dir is
// the mount point, or above it.
static int cgroup_up(char *dir, const char *mnt) {
  if (strlen(dir) <= strlen(mnt)) return 0;
  char *slash = strrchr(dir, '/');
  if (!slash) return 0;
  *slash = '\0';
  return strlen(dir) >= strlen(mnt);
}


// The smallest CPU quota, in cores, set on the process's
// cgroup or any of its parents, or < 0 if none.
static double linux_quota(void) {

  char   dir[4096], mnt[4096], file[4200], buf[256];
  double best = -1;

  // cgroup v2: cpu.max is "max period" or "quota period".
  if (cgroup_dir(NULL, dir, mnt, sizeof(dir))) {
    do {
      snprintf(file, sizeof(file), "%s/cpu.max", dir);
      if (!read_line(file, buf, sizeof(buf)) || strncmp(buf, "max", 3) == 0)
        continue;
      char  *end;
      double quota  = strtod(buf, &end);
      double period = strtod(end, NULL);
      if (quota > 0 && period > 0 && (best < 0 || quota / period < best))
        best = quota / period;
    } while (cgroup_up(dir, mnt));
  }

  // cgroup v1: cpu.cfs_quota_us is -1 when unlimited.
  if (cgroup_dir("cpu", dir, mnt, sizeof(dir))) {
    do {
      snprintf(file, sizeof(file), "%s/cpu.cfs_quota_us", dir);
      if (!read_line(file, buf, sizeof(buf))) continue;
      double quota = strtod(buf, NULL);
      snprintf(file, sizeof(file), "%s/cpu.cfs_period_us", dir);
      if (!read_line(file, buf, sizeof(buf))) continue;
      double period = strtod(buf, NULL);
      if (quota > 0 && period > 0 && (best < 0 || quota / period < best))
        best = quota / period;
    } while (cgroup_up(dir, mnt));
  }

  return best;
}

#endif // __linux__


//------------------------------------------------------
// Windows
//------------------------------------------------------
#ifdef _WIN32

static int popcount(DWORD_PTR x) {
  int n = 0;
  for (; x; x &= x - 1) n++;
  return n;
}


static int win_physical(void) {

  DWORD len = 0;
  GetLogicalProcessorInformationEx(RelationProcessorCore, NULL, &len);
  if (len == 0) return NA_INTEGER; // # nocov

  char *buf = malloc(len);
  if (!buf) return NA_INTEGER; // # nocov

  int n = 0;
  if (GetLogicalProcessorInformationEx(
      RelationProcessorCore, (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) buf, &len )) {
    for (DWORD off = 0; off < len; ) {
      PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX p =
        (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) (buf + off);
      if (p->Size == 0) break; // # nocov
      if (p->Relationship == RelationProcessorCore) n++;
      off += p->Size;
    }
  }
  free(buf);

  return n > 0 ? n : NA_INTEGER;
}


// The process's affinity mask covers one processor group
// (up to 64 cores). On a machine with several groups, a
// mask that is the whole group is no limit: Windows 11
// and Server 2022 run a process's threads on every group.
static int win_affinity(void) {

  DWORD_PTR proc, sys;
  if (!GetProcessAffinityMask(GetCurrentProcess(), &proc, &sys) || !proc)
    return NA_INTEGER; // # nocov: threads in several groups already

  if (proc == sys && GetActiveProcessorGroupCount() > 1)
    return NA_INTEGER; // # nocov

  return popcount(proc);
}


// A job object's hard cap on CPU time, as used by Windows
// containers' --cpus. Declared here rather than taken
// from <winnt.h>, which only has it for newer targets.
typedef struct {
  DWORD ControlFlags;
  union {
    DWORD CpuRate;
    DWORD Weight;
    struct { WORD MinRate; WORD MaxRate; } mm;
  } u;
} cpu_rate_t;

#define CPU_RATE_INFO       15     // JobObjectCpuRateControlInformation
#define CPU_RATE_ENABLE     0x1
#define CPU_RATE_WEIGHT     0x2
#define CPU_RATE_HARD_CAP   0x4
#define CPU_RATE_MIN_MAX    0x10

static double win_quota(int logical) {

  cpu_rate_t r;
  memset(&r, 0, sizeof(r));

  // NULL: the job the process is in, if any.
  if (!QueryInformationJobObject(NULL, (JOBOBJECTINFOCLASS) CPU_RATE_INFO, &r, sizeof(r), NULL))
    return -1;

  // # nocov start (only inside a job object with a CPU cap)
  if (!(r.ControlFlags & CPU_RATE_ENABLE) || (r.ControlFlags & CPU_RATE_WEIGHT))
    return -1;

  // In hundredths of a percent of all the cores.
  double rate = -1;
  if      (r.ControlFlags & CPU_RATE_MIN_MAX)  rate = r.u.mm.MaxRate;
  else if (r.ControlFlags & CPU_RATE_HARD_CAP) rate = r.u.CpuRate;

  return rate > 0 ? rate / 10000.0 * logical : -1;
  // # nocov end
}

#endif // _WIN32


//------------------------------------------------------
// Read the values that are kept for the process.
//------------------------------------------------------
static void read_cpu_info(void) {

  #ifndef _WIN32
    if (ci_logical > 0 && ci_pid == getpid()) return;
    ci_pid = getpid();
  #else
    if (ci_logical > 0) return;
  #endif

  int    logical  = NA_INTEGER;
  int    physical = NA_INTEGER;
  double quota    = -1;

  #if defined(_WIN32)

    DWORD n = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    if (n > 0 && n < INT_MAX) logical = (int) n;
    physical = win_physical();
    if (logical != NA_INTEGER) quota = win_quota(logical);

  #elif defined(__APPLE__)

    int    v;
    size_t len = sizeof(v);
    if (sysctlbyname("hw.logicalcpu", &v, &len, NULL, 0) == 0 && v > 0) logical = v;
    len = sizeof(v);
    if (sysctlbyname("hw.physicalcpu", &v, &len, NULL, 0) == 0 && v > 0) physical = v;

  #else

    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n > 0 && n < INT_MAX) logical = (int) n;
    #ifdef __linux__
      physical = linux_physical();
      quota    = linux_quota();
    #endif

  #endif

  if (logical == NA_INTEGER) logical = 1; // # nocov

  // A quota is rounded to whole cores, at least 1, and
  // ignored when it is more than the machine has.
  int q = NA_INTEGER;
  if (quota > 0 && quota <= logical) {
    q = (int) floor(quota + 0.5);
    if (q < 1) q = 1;
  }

  ci_logical  = logical;
  ci_physical = physical;
  ci_quota    = q;
}


//------------------------------------------------------
// An environment variable's value, or NA_STRING. On
// Windows, read from the process's environment block,
// which Sys.setenv() updates: the C runtime's getenv()
// takes microseconds there.
//------------------------------------------------------
static SEXP env_value(const char *name) {

  #ifdef _WIN32

    char  buf[256];
    DWORD len = GetEnvironmentVariableA(name, buf, sizeof(buf));
    if (len == 0)          return NA_STRING;
    if (len < sizeof(buf)) return mkChar(buf);

    // # nocov start (values of 256 bytes or more)
    char *big = malloc(len);
    if (!big) return NA_STRING;
    DWORD got = GetEnvironmentVariableA(name, big, len);
    SEXP  res = (got > 0 && got < len) ? mkChar(big) : NA_STRING;
    free(big);
    return res;
    // # nocov end

  #else

    const char *val = getenv(name);
    return val ? mkChar(val) : NA_STRING;

  #endif
}


//======================================================
// R interface. Returns list(
//   c(logical, physical, affinity, quota, forked, pthreads),
//   the values of the environment variables named in
//   `sexp_env` (NA when unset)
// ). `forked` is 1 in a process forked from the one that
// loaded ecodive, e.g. by parallel::mclapply().
// Used by R code for `ecodive::n_cpus()`.
//======================================================
SEXP C_cpu_info(SEXP sexp_env) {

  read_cpu_info();

  int affinity = NA_INTEGER;
  #if defined(__linux__)
    affinity = linux_affinity();
  #elif defined(_WIN32)
    affinity = win_affinity();
  #endif

  int forked = 0;
  #ifndef _WIN32
    forked = load_pid != 0 && getpid() != load_pid;
  #endif

  int pthreads = 0;
  #ifdef HAVE_PTHREAD
    pthreads = 1;
  #endif

  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SEXP info   = allocVector(INTSXP, 6);
  SET_VECTOR_ELT(result, 0, info);

  int *p = INTEGER(info);
  p[0] = ci_logical;
  p[1] = ci_physical;
  p[2] = affinity;
  p[3] = ci_quota;
  p[4] = forked;
  p[5] = pthreads;

  // Environment variables, much faster than Sys.getenv().
  R_xlen_t n   = XLENGTH(sexp_env);
  SEXP     env = allocVector(STRSXP, n);
  SET_VECTOR_ELT(result, 1, env);
  for (R_xlen_t i = 0; i < n; i++)
    SET_STRING_ELT(env, i, env_value(CHAR(STRING_ELT(sexp_env, i))));
  setAttrib(env, R_NamesSymbol, sexp_env);

  UNPROTECT(1);
  return result;
}
