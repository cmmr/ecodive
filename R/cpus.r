# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Number of CPU Cores
#'
#' How many threads ecodive's functions use by default: the CPU cores this
#' R process may use, within any limits set on it.
#'
#' @details
#' `n_cpus()` is the default for every function's `cpus` argument. It
#' follows the checks of `parallelly::availableCores()`, without needing the
#' `parallelly` package, and takes about as long as reading a few
#' environment variables. It is the smallest of the values below that
#' apply.
#'
#' **Resource limits.** A `cpus` given explicitly is capped at the smallest
#' of these alone:
#'
#' * The logical cores online (physical cores with
#'   `options(parallelly.availableCores.logical = FALSE)`).
#' * The cores the process may run on: its CPU affinity, as set by
#'   `taskset`, `numactl`, cgroup cpusets (Docker's `--cpuset-cpus`), job
#'   schedulers that bind tasks to cores, or Windows' `start /affinity` and
#'   job objects.
#' * A CPU quota, in whole cores: Linux cgroup v1 and v2 quotas (Docker's
#'   `--cpus`, Kubernetes CPU limits, systemd's `CPUQuota`) or a Windows job
#'   object's hard CPU rate cap.
#' * The cores a job scheduler allocated on this node: Slurm
#'   (`SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, ...), PBS/Torque
#'   (`PBS_NUM_PPN`, `NCPUS`), LSF (`LSB_DJOB_NUMPROC`), SGE (the
#'   `PE_HOSTFILE`, `NSLOTS`), Fujitsu PJM (`PJM_VNODE_CORE`,
#'   `PJM_PROC_BY_NODE`) and HTCondor (`_CONDOR_MACHINE_AD`,
#'   `_CONDOR_JOB_AD`).
#' * 2 under `R CMD check --as-cran` (`_R_CHECK_LIMIT_CORES_`), and 4 on
#'   Bioconductor's build machines.
#' * `options(ecodive.cpus = n)`, or the environment variable
#'   `R_ECODIVE_CPUS`.
#'
#' **Defaults.** These lower `n_cpus()`, but not the cap on an explicit
#' `cpus`:
#'
#' * The `mc.cores` option, or the `MC_CORES` environment variable.
#' * `BIOCPARALLEL_WORKER_NUMBER`.
#' * 1 in a worker process: one forked by, e.g., `parallel::mclapply()`,
#'   future's multicore plan or BiocParallel's `MulticoreParam` after
#'   ecodive was loaded, or a worker of a socket cluster (`parallel`,
#'   `parallelly`, future's multisession plan, `snow`), `mirai` or
#'   `clustermq`. Those workers already run side by side, so more threads
#'   in each would compete for the same cores. Give `cpus` to use more.
#' * These options of `parallelly` (`parallelly.availableCores.*`, also
#'   read as `future.availableCores.*`, or the environment variables
#'   `R_PARALLELLY_AVAILABLECORES_*`), whether or not `parallelly` is
#'   installed: `methods` and `methods.excludes` (which of parallelly's
#'   methods above apply), `fallback`, `custom`, `min`, `system`, `logical`,
#'   `threadsPerWorker`, `fraction`, `omit` and `max`.
#'
#' Settings that are missing, not numbers or less than 1 are ignored.
#' Results are the same for any number of threads.
#'
#' @return   A scalar integer, at least `1`. Always `1` where ecodive was
#'           built without `pthread` support.
#'
#' @export
#'
#' @examples
#'     n_cpus()
#'
n_cpus <- function () {
  cpu_limits()[[2L]]
}


# The largest `cpus` that functions accept: the resource limits alone.
max_cpus <- function () {
  cpu_limits()[[1L]]
}



# The environment variables cpu_limits() reads, all in one .Call().
PLY_SETTINGS <- c(
  'methods', 'methods.excludes', 'fallback', 'min', 'system', 'logical',
  'threadsPerWorker', 'fraction', 'omit', 'max', 'custom' )

PLY_ENV <- toupper(gsub('.', '_', PLY_SETTINGS, fixed = TRUE))

CPU_ENV <- c(
  '_R_CHECK_LIMIT_CORES_', 'R_ECODIVE_CPUS', 'MC_CORES',
  'BIOCPARALLEL_WORKER_NUMBER', 'IS_BIOC_BUILD_MACHINE', 'BBS_HOME',
  'SLURM_CPUS_ON_NODE', 'SLURM_STEP_ID', 'SLURM_STEPID',
  'SLURM_CPUS_PER_TASK', 'SLURM_STEP_TASKS_PER_NODE', 'SLURM_STEP_NODELIST',
  'SLURMD_NODENAME', 'SLURM_JOB_NUM_NODES', 'SLURM_NNODES',
  'SLURM_TASKS_PER_NODE',
  'PBS_NUM_PPN', 'NCPUS', 'LSB_DJOB_NUMPROC',
  'PJM_VNODE_CORE', 'PJM_PROC_BY_NODE',
  'PE_HOSTFILE', 'SGE_JOB_SPOOL_DIR', 'NSLOTS',
  '_CONDOR_MACHINE_AD', '_CONDOR_JOB_AD',
  paste0('R_PARALLELLY_AVAILABLECORES_', PLY_ENV),
  paste0('R_FUTURE_AVAILABLECORES_',     PLY_ENV) )

# parallelly's methods, in its default order.
PLY_METHODS <- c(
  'system', '/proc/self/status', 'cgroups.cpuset', 'cgroups.cpuquota',
  'cgroups2.cpuset.cpus', 'cgroups2.cpuset.cpus.effective',
  'cgroups2.cpu.max', 'nproc', 'mc.cores', 'BiocParallel',
  '_R_CHECK_LIMIT_CORES_', 'Bioconductor', 'LSF', 'PJM', 'PBS', 'SGE',
  'Slurm', 'HTCondor', 'fallback', 'custom' )

# parallelly's methods that read the process's CPU affinity or quota.
PLY_AFFINITY <- c(
  '/proc/self/status', 'nproc', 'cgroups.cpuset', 'cgroups2.cpuset.cpus',
  'cgroups2.cpuset.cpus.effective' )
PLY_QUOTA <- c('cgroups.cpuquota', 'cgroups2.cpu.max')

# Command lines of cluster workers: parallel's and parallelly's socket
# workers (.workRSOCK, .slaveRSOCK), snow's (RSOCKnode.R, RMPInode.R,
# RMPISNOW), mirai's daemons and clustermq's workers.
WORKER_ARGS <- 'workRSOCK|slaveRSOCK|RSOCKnode|RMPInode|RMPISNOW|mirai::daemon|clustermq:::worker'



#------------------------------------------------------------------------------
# c(cap, default): the largest `cpus` accepted, and n_cpus().
#
# `info` and `worker` replace what the process reports, for tests: info is
# c(logical, physical, affinity, quota, forked, pthreads), as C_cpu_info()
# returns.
#------------------------------------------------------------------------------
cpu_limits <- function (info = NULL, worker = NULL) {

  res <- .Call(C_cpu_info, CPU_ENV)
  env <- res[[2L]]
  if (is.null(info)) info <- res[[1L]]

  if (!info[[6L]]) return (c(1L, 1L))

  if (is.null(worker)) worker <- info[[5L]] == 1L || cluster_worker()

  # Usually no option or environment variable is set: then only the
  # machine's limits apply.
  opts <- names(.Options)
  opts <- opts[startsWith(opts, 'parallelly.availableCores.') |
               startsWith(opts, 'future.availableCores.') |
               opts == 'mc.cores' | opts == 'ecodive.cpus' ]
  if (!length(opts) && all(is.na(env))) {
    cap <- min(info[c(1L, 3L, 4L)], na.rm = TRUE) # logical, affinity, quota
    return (c(cap, if (worker) 1L else cap))
  }

  logical <- info[[1L]]
  ply     <- ply_settings(env, opts)

  methods <- if (is.null(ply[['methods']])) PLY_METHODS else ply[['methods']]
  if (!is.null(ply[['methods.excludes']]))
    methods <- methods[!methods %in% ply[['methods.excludes']]]


  # Resource limits ===========================================================

  hard <- c(system = logical)

  if (!is.null(ply[['system']])) {
    hard[['system']] <- ply[['system']]
  }
  else if (isFALSE(ply[['logical']]) && !is.na(info[[2L]])) {
    hard[['system']] <- info[[2L]]
  }

  if (!'system' %in% methods)            hard <- hard[0]
  if (any(PLY_AFFINITY %in% methods))    hard[['affinity']] <- info[[3L]]
  if (any(PLY_QUOTA    %in% methods))    hard[['quota']]    <- info[[4L]]

  if ('_R_CHECK_LIMIT_CORES_' %in% methods) {
    chk <- tolower(env[['_R_CHECK_LIMIT_CORES_']])
    if (!is.na(chk) && nzchar(chk) && chk != 'false')
      hard[['_R_CHECK_LIMIT_CORES_']] <- 2L
  }

  if ('Bioconductor' %in% methods) {
    if (isTRUE(as.logical(env[['IS_BIOC_BUILD_MACHINE']])) ||
        isTRUE(as.logical(env[['BBS_HOME']])) )
      hard[['Bioconductor']] <- 4L
  }

  if ('LSF'      %in% methods) hard[['LSF']]      <- env_int(env, 'LSB_DJOB_NUMPROC')
  if ('PJM'      %in% methods) hard[['PJM']]      <- env_int(env, c('PJM_VNODE_CORE', 'PJM_PROC_BY_NODE'))
  if ('PBS'      %in% methods) hard[['PBS']]      <- env_int(env, c('PBS_NUM_PPN', 'NCPUS'))
  if ('SGE'      %in% methods) hard[['SGE']]      <- cores_sge(env)
  if ('Slurm'    %in% methods) hard[['Slurm']]    <- cores_slurm(env)
  if ('HTCondor' %in% methods) hard[['HTCondor']] <- cores_htcondor(env)

  ecodive <- getOption('ecodive.cpus')
  if (is.null(ecodive)) ecodive <- env[['R_ECODIVE_CPUS']]
  hard[['ecodive']] <- as_count(ecodive)


  # Defaults ==================================================================

  soft <- integer(0)

  if ('mc.cores' %in% methods || 'mc.cores+1' %in% methods) {
    mc <- getOption('mc.cores')
    if (is.null(mc)) mc <- env[['MC_CORES']]
    mc <- as_count(mc, zero = TRUE)
    if ('mc.cores' %in% methods)
      soft[['mc.cores']] <- if (isTRUE(mc == 0L)) 1L else mc
    if ('mc.cores+1' %in% methods)
      soft[['mc.cores+1']] <- mc + 1L
  }

  if ('BiocParallel' %in% methods)
    soft[['BiocParallel']] <- env_int(env, 'BIOCPARALLEL_WORKER_NUMBER')

  if ('custom' %in% methods && is.function(ply[['custom']]))
    soft[['custom']] <- as_count(tryCatch(ply[['custom']](), error = function (e) NA))

  # Any other method names an option or environment variable.
  for (m in methods[!methods %in% c(PLY_METHODS, 'mc.cores+1')]) {
    v <- getOption(m)
    if (is.null(v)) v <- Sys.getenv(m, NA)
    soft[[m]] <- as_count(v)
  }


  # Combined, as parallelly::availableCores() does. =========================

  hard <- hard[!is.na(hard) & hard >= 1]
  soft <- soft[!is.na(soft) & soft >= 1]
  all  <- c(hard, soft)

  cap  <- if (length(hard)) min(hard) else logical

  # The fallback replaces the machine's limits when nothing else is set and
  # the affinity and quota don't limit them.
  dflt <- if (length(all)) min(all) else 1L
  fb   <- if ('fallback' %in% methods) as_count(ply[['fallback']]) else NA
  if (!is.na(fb)) {
    machine <- c('system', 'affinity', 'quota')
    others  <- setdiff(names(all), machine)
    limited <- any(all[intersect(c('affinity', 'quota'), names(all))] < all['system'])
    if (!length(others) && !isTRUE(limited)) dflt <- fb
  }

  if (!is.null(ply[['min']]) && isTRUE(dflt < ply[['min']]))
    dflt <- min(as.integer(ply[['min']]), logical)
  if (!is.null(ply[['threadsPerWorker']]) && isTRUE(ply[['threadsPerWorker']] > 1))
    dflt <- dflt %/% as.integer(ply[['threadsPerWorker']])
  if (!is.null(ply[['fraction']]) && isTRUE(ply[['fraction']] > 0 && ply[['fraction']] < 1))
    dflt <- floor(ply[['fraction']] * dflt)
  if (!is.null(ply[['omit']]) && isTRUE(ply[['omit']] > 0))
    dflt <- dflt - as.integer(ply[['omit']])
  if (!is.null(ply[['max']]) && isTRUE(ply[['max']] >= 1))
    dflt <- min(dflt, ply[['max']])

  if (worker) dflt <- 1L

  cap  <- as.integer(max(1, cap))
  dflt <- as.integer(max(1, min(dflt, cap)))

  c(cap, dflt)
}



#------------------------------------------------------------------------------
# The availableCores settings of parallelly that are set, as a named list:
# from options parallelly.availableCores.X or future.availableCores.X, else
# environment variables R_PARALLELLY_AVAILABLECORES_X or
# R_FUTURE_AVAILABLECORES_X.
#------------------------------------------------------------------------------
ply_settings <- function (env, opts) {

  vars <- env[startsWith(names(env), 'R_PARALLELLY_') | startsWith(names(env), 'R_FUTURE_')]
  if (!any(startsWith(opts, 'parallelly.') | startsWith(opts, 'future.')) &&
      all(is.na(vars)) )
    return (list())

  res <- list()
  for (i in seq_along(PLY_SETTINGS)) {

    x <- PLY_SETTINGS[[i]]
    v <- getOption(paste0('parallelly.availableCores.', x))
    if (is.null(v)) v <- getOption(paste0('future.availableCores.', x))

    if (is.null(v)) {
      v <- vars[[paste0('R_PARALLELLY_AVAILABLECORES_', PLY_ENV[[i]])]]
      if (is.na(v)) v <- vars[[paste0('R_FUTURE_AVAILABLECORES_', PLY_ENV[[i]])]]
      if (is.na(v) || !nzchar(v <- trimws(v))) next
      if (x %in% c('methods', 'methods.excludes'))
        v <- trimws(strsplit(v, ',', fixed = TRUE)[[1]])
    }

    v <- switch(
      EXPR             = x,
      methods          = ,
      methods.excludes = as.character(v),
      custom           = if (is.function(v)) v else NULL,
      logical          = as.logical(v)[1],
      suppressWarnings(as.numeric(v)[1]) )

    if (is.function(v) || (length(v) && !anyNA(v))) res[[x]] <- v
  }

  res
}


# Whether this process is a cluster worker, by its command line. Fixed for
# the process, so kept.
cluster_worker <- function () {
  if (is.null(CACHE$cluster_worker))
    CACHE$cluster_worker <- any(grepl(WORKER_ARGS, commandArgs()))
  CACHE$cluster_worker
}


# A count from an option or environment variable: a whole number, or NA.
as_count <- function (x, zero = FALSE) {
  if (length(x) != 1L || is.na(x)) return (NA_integer_)
  # as.numeric() allows the spaces around a number.
  if (is.character(x)) x <- suppressWarnings(as.numeric(x))
  if (is.na(x) || !is.finite(x) || x < (if (zero) 0 else 1) || x > .Machine$integer.max)
    return (NA_integer_)
  as.integer(x)
}


# The first of the environment variables `nms` that is a count.
env_int <- function (env, nms) {
  for (nm in nms) {
    v <- env[[nm]]
    if (!is.na(v)) {
      v <- suppressWarnings(as.integer(trimws(v)))
      if (!is.na(v)) return (v)
    }
  }
  NA_integer_
}



#------------------------------------------------------------------------------
# Slurm: the CPUs of this task on this node. As parallelly's
# availableCoresSlurm(): within a job step (srun), SLURM_CPUS_PER_TASK, or
# the node's CPUs shared among the step's tasks on it; in the batch script,
# SLURM_CPUS_ON_NODE.
#------------------------------------------------------------------------------
cores_slurm <- function (env) {

  ncpus <- env_int(env, 'SLURM_CPUS_ON_NODE')

  step <- suppressWarnings(as.numeric(env[['SLURM_STEP_ID']]))
  if (is.na(step)) step <- suppressWarnings(as.numeric(env[['SLURM_STEPID']]))
  # The batch and extern steps' IDs.
  if (!is.na(step) && step >= 4294967280) step <- NA

  if (is.na(step)) {
    n <- ncpus
  }
  else {
    n <- env_int(env, 'SLURM_CPUS_PER_TASK')
    if (is.na(n) && !is.na(ncpus)) {
      ntasks <- slurm_step_ntasks(env)
      if (!is.na(ntasks) && ntasks > 0L) n <- max(1L, ncpus %/% ntasks)
    }
  }

  if (is.na(n)) n <- env_int(env, 'SLURM_CPUS_PER_TASK')

  if (is.na(n)) {
    nnodes <- env_int(env, c('SLURM_JOB_NUM_NODES', 'SLURM_NNODES'))
    if (is.na(nnodes) || nnodes == 1L) {
      n <- ncpus
    }
    else {
      counts <- env[['SLURM_TASKS_PER_NODE']]
      if (!is.na(counts)) {
        counts <- slurm_counts(counts)
        n <- if (length(counts) && !anyNA(counts)) counts[[1]] else NA_integer_
      }
    }
  }

  n
}


# Tasks of the current job step on this node.
slurm_step_ntasks <- function (env) {

  counts <- env[['SLURM_STEP_TASKS_PER_NODE']]
  if (is.na(counts)) return (NA_integer_)
  counts <- slurm_counts(counts)
  if (!length(counts) || anyNA(counts)) return (NA_integer_)
  if (length(unique(counts)) == 1L) return (counts[[1]])

  nodelist <- env[['SLURM_STEP_NODELIST']]
  nodename <- env[['SLURMD_NODENAME']]
  if (is.na(nodelist) || is.na(nodename)) return (NA_integer_)
  nodes <- slurm_nodes(nodelist)
  if (length(nodes) != length(counts)) return (NA_integer_)
  counts[match(nodename, nodes)]
}


# Slurm's per-node counts: "2(x3),1" is c(2, 2, 2, 1).
slurm_counts <- function (x) {
  x <- gsub('[[:space:]]', '', strsplit(x, ',', fixed = TRUE)[[1]])
  x <- x[nzchar(x)]
  unlist(lapply(x, function (count) {
    m <- regmatches(count, regexec('^([0-9]+)\\(x([0-9]+)\\)$', count))[[1]]
    if (length(m)) rep(as.integer(m[[2]]), as.integer(m[[3]]))
    else           suppressWarnings(as.integer(count))
  }))
}


# Slurm's host lists: "a[01-03,7],b" is c("a01", "a02", "a03", "a7", "b").
slurm_nodes <- function (x) {

  # Split at the commas outside brackets.
  chars <- strsplit(x, '')[[1]]
  depth <- cumsum((chars == '[') - (chars == ']'))
  cut   <- which(chars == ',' & depth == 0L)
  parts <- substring(x, c(1L, cut + 1L), c(cut - 1L, nchar(x)))
  parts <- parts[nzchar(parts)]

  unlist(lapply(parts, function (p) {
    m <- regmatches(p, regexec('^([^[]*)\\[([^]]*)\\](.*)$', p))[[1]]
    if (!length(m)) return (p)
    ids <- unlist(lapply(strsplit(m[[3]], ',', fixed = TRUE)[[1]], function (r) {
      ab <- strsplit(r, '-', fixed = TRUE)[[1]]
      if (length(ab) == 1L) return (ab)
      sprintf('%0*d', nchar(ab[[1]]), seq(as.integer(ab[[1]]), as.integer(ab[[2]])))
    }))
    paste0(m[[2]], ids, m[[4]])
  }))
}



#------------------------------------------------------------------------------
# SGE: this host's slots in the parallel environment's host file, else
# NSLOTS. The file is read once.
#------------------------------------------------------------------------------
cores_sge <- function (env) {

  path  <- env[['PE_HOSTFILE']]
  spool <- env[['SGE_JOB_SPOOL_DIR']]
  if (is.na(path) && is.na(spool)) return (env_int(env, 'NSLOTS'))

  key <- paste(path, spool, sep = '\r')
  if (!identical(CACHE$sge_key, key)) {
    CACHE$sge_key <- key
    CACHE$sge     <- sge_slots(env, path, spool)
  }

  if (is.na(CACHE$sge)) env_int(env, 'NSLOTS') else CACHE$sge
}


sge_slots <- function (env, path, spool) {

  if (is.na(path) || !file.exists(path)) {
    if (is.na(spool)) return (NA_integer_)
    path <- file.path(c(spool, dirname(spool)), 'pe_hostfile')
    path <- path[file.exists(path)][1]
    if (is.na(path)) return (NA_integer_)
  }

  lines <- tryCatch(readLines(path, warn = FALSE), error = function (e) NULL)
  if (!length(lines)) return (NA_integer_)
  fields <- strsplit(trimws(lines), '[[:space:]]+')
  nodes  <- vapply(fields, `[`, '', 1L)
  counts <- suppressWarnings(as.integer(vapply(fields, `[`, '', 2L)))
  if (anyNA(counts)) return (NA_integer_)

  host <- Sys.getenv('HOSTNAME', NA)
  if (is.na(host)) host <- Sys.info()[['nodename']]
  short <- function (x) sub('[.].*', '', x)
  here  <- short(nodes) == short(host)
  if (!any(here)) return (NA_integer_)
  sum(counts[here])
}



#------------------------------------------------------------------------------
# HTCondor: Cpus from the machine ad, else CpusProvisioned or Cpus from the
# job ad. The files are read once.
#------------------------------------------------------------------------------
cores_htcondor <- function (env) {

  machine <- env[['_CONDOR_MACHINE_AD']]
  job     <- env[['_CONDOR_JOB_AD']]
  if (is.na(machine) && is.na(job)) return (NA_integer_)

  key <- paste(machine, job, sep = '\r')
  if (!identical(CACHE$condor_key, key)) {
    n <- classad_count(machine, 'Cpus')
    if (is.na(n)) n <- classad_count(job, 'CpusProvisioned')
    if (is.na(n)) n <- classad_count(job, 'Cpus')
    CACHE$condor_key <- key
    CACHE$condor     <- n
  }

  CACHE$condor
}


# A ClassAd file's `attr = value`, rounded to a whole number, or NA.
classad_count <- function (path, attr) {

  if (is.na(path) || !nzchar(path) || !file.exists(path)) return (NA_integer_)

  lines <- tryCatch(readLines(path, warn = FALSE), error = function (e) NULL)
  lines <- lines[!grepl('^\\s*(#|//)', lines)]
  pat   <- sprintf('^\\s*%s\\s*=\\s*(.*)$', attr)
  hit   <- grep(pat, lines, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) return (NA_integer_)

  val <- trimws(sub(pat, '\\1', hit[[1]], ignore.case = TRUE))
  val <- gsub('^["\']|["\']$', '', val)
  val <- suppressWarnings(as.numeric(val))
  if (is.na(val) || !is.finite(val)) return (NA_integer_)
  val <- as.integer(floor(val + 0.5))
  if (val <= 0L) NA_integer_ else val
}
