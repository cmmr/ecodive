#test_that("n_cpus", {

  # The probe. ====

  res <- .Call(C_cpu_info, c('PATH', 'ECODIVE_NO_SUCH_VAR'))
  expect_true(is.list(res) && length(res) == 2L)
  expect_true(is.integer(res[[1]]) && length(res[[1]]) == 6L)
  expect_true(res[[1]][[1]] >= 1L)
  expect_identical(names(res[[2]]), c('PATH', 'ECODIVE_NO_SUCH_VAR'))
  expect_identical(res[[2]][[2]], NA_character_)
  expect_identical(res[[1]][[5]], 0L) # not forked

  # Sys.setenv() is seen.
  Sys.setenv(ECODIVE_TEST_VAR = 'abc')
  expect_identical(.Call(C_cpu_info, 'ECODIVE_TEST_VAR')[[2]][[1]], 'abc')
  Sys.unsetenv('ECODIVE_TEST_VAR')
  expect_identical(.Call(C_cpu_info, 'ECODIVE_TEST_VAR')[[2]][[1]], NA_character_)

  # This process.
  n <- n_cpus()
  m <- max_cpus()
  expect_true(is.integer(n) && length(n) == 1L && n >= 1L)
  expect_true(is.integer(m) && length(m) == 1L && n <= m)
  expect_false(cluster_worker())



  # cpu_limits(), with the machine replaced. ====

  # 8 logical cores, 4 physical; affinity and quota unknown.
  hw <- c(8L, 4L, NA, NA, 0L, 1L)

  # c(cap, default) with only `env` and `opts` set.
  ply_opts <- function ()
    grep('^(parallelly|future)[.]availableCores[.]', names(.Options), value = TRUE)

  case <- function (env = list(), opts = list(), info = hw, worker = FALSE) {

    vars <- unique(c(CPU_ENV, names(env)))
    old  <- Sys.getenv(vars, NA, names = TRUE)
    nms  <- unique(c('mc.cores', 'ecodive.cpus', ply_opts(), names(opts)))
    oldo <- options(stats::setNames(rep(list(NULL), length(nms)), nms))

    on.exit({
      options(oldo)
      Sys.unsetenv(vars)
      set <- old[!is.na(old)]
      if (length(set)) do.call(Sys.setenv, as.list(set))
    })

    Sys.unsetenv(vars)
    if (length(env))  do.call(Sys.setenv, env)
    if (length(opts)) options(opts)

    cpu_limits(info = info, worker = worker)
  }

  # The machine.
  expect_identical(case(), c(8L, 8L))
  expect_identical(case(info = c(8L, 4L, 3L, NA, 0L, 1L)), c(3L, 3L))
  expect_identical(case(info = c(8L, 4L, 3L,  2L, 0L, 1L)), c(2L, 2L))
  expect_identical(case(info = c(8L, 4L, NA, NA, 0L, 0L)), c(1L, 1L)) # no pthreads

  # Workers use one thread by default.
  expect_identical(case(worker = TRUE), c(8L, 1L))
  expect_identical(case(env = list(R_ECODIVE_CPUS = '4'), worker = TRUE), c(4L, 1L))
  expect_identical(case(info = c(8L, 4L, NA, NA, 1L, 1L), worker = NULL), c(8L, 1L))
  expect_true(grepl(WORKER_ARGS, 'tryCatch(parallel:::.workRSOCK,error=function(e)parallel:::.slaveRSOCK)()'))
  expect_true(grepl(WORKER_ARGS, "mirai::daemon('tcp://localhost:5555')"))
  expect_false(grepl(WORKER_ARGS, '--file=analysis.R'))

  # R CMD check and Bioconductor.
  expect_identical(case(list(`_R_CHECK_LIMIT_CORES_` = 'TRUE')),  c(2L, 2L))
  expect_identical(case(list(`_R_CHECK_LIMIT_CORES_` = 'warn')),  c(2L, 2L))
  expect_identical(case(list(`_R_CHECK_LIMIT_CORES_` = 'false')), c(8L, 8L))
  expect_identical(case(list(IS_BIOC_BUILD_MACHINE = 'true')),    c(4L, 4L))

  # ecodive's own limit.
  expect_identical(case(list(R_ECODIVE_CPUS = '3')),            c(3L, 3L))
  expect_identical(case(list(R_ECODIVE_CPUS = ' 3 ')),          c(3L, 3L))
  expect_identical(case(opts = list(ecodive.cpus = 5)),         c(5L, 5L))
  expect_identical(case(list(R_ECODIVE_CPUS = '3'), list(ecodive.cpus = 5)), c(5L, 5L))
  expect_identical(case(list(R_ECODIVE_CPUS = 'abc')),          c(8L, 8L))
  expect_identical(case(list(R_ECODIVE_CPUS = '0')),            c(8L, 8L))
  expect_identical(case(opts = list(ecodive.cpus = 100)),       c(8L, 8L))

  # mc.cores and BiocParallel lower the default only.
  expect_identical(case(opts = list(mc.cores = 2)),             c(8L, 2L))
  expect_identical(case(opts = list(mc.cores = 0)),             c(8L, 1L))
  expect_identical(case(list(MC_CORES = '3')),                  c(8L, 3L))
  expect_identical(case(list(MC_CORES = '3'), list(mc.cores = 5)), c(8L, 5L))
  expect_identical(case(list(BIOCPARALLEL_WORKER_NUMBER = '3')), c(8L, 3L))

  # Schedulers.
  expect_identical(case(list(LSB_DJOB_NUMPROC = '5')),          c(5L, 5L))
  expect_identical(case(list(PBS_NUM_PPN = '6')),               c(6L, 6L))
  expect_identical(case(list(NCPUS = '7')),                     c(7L, 7L))
  expect_identical(case(list(PBS_NUM_PPN = '6', NCPUS = '7')),  c(6L, 6L))
  expect_identical(case(list(PJM_VNODE_CORE = '3')),            c(3L, 3L))
  expect_identical(case(list(PJM_PROC_BY_NODE = '2')),          c(2L, 2L))
  expect_identical(case(list(NSLOTS = '4')),                    c(4L, 4L))

  # Slurm: in the batch script; in a job step; in a step whose tasks
  # share the node; in the batch or extern step; across nodes.
  expect_identical(case(list(SLURM_CPUS_ON_NODE = '6')), c(6L, 6L))
  expect_identical(case(list(
    SLURM_CPUS_ON_NODE = '6', SLURM_STEP_ID = '0', SLURM_CPUS_PER_TASK = '2' )), c(2L, 2L))
  expect_identical(case(list(
    SLURM_CPUS_ON_NODE = '6', SLURM_STEP_ID = '0',
    SLURM_STEP_TASKS_PER_NODE = '3(x2)' )), c(2L, 2L))
  expect_identical(case(list(
    SLURM_CPUS_ON_NODE = '8', SLURM_STEP_ID = '1',
    SLURM_STEP_TASKS_PER_NODE = '2,4', SLURM_STEP_NODELIST = 'n[1-2]',
    SLURMD_NODENAME = 'n2' )), c(2L, 2L))
  expect_identical(case(list(
    SLURM_CPUS_ON_NODE = '6', SLURM_STEP_ID = '4294967294',
    SLURM_CPUS_PER_TASK = '2' )), c(6L, 6L))
  expect_identical(case(list(
    SLURM_JOB_NUM_NODES = '2', SLURM_TASKS_PER_NODE = '4(x2)' )), c(4L, 4L))

  expect_identical(slurm_counts('2(x3),1'), c(2L, 2L, 2L, 1L))
  expect_identical(slurm_counts(' 5 '), 5L)
  expect_identical(slurm_nodes('a[01-03,7],b'), c('a01', 'a02', 'a03', 'a7', 'b'))
  expect_identical(slurm_nodes('node1,node2'), c('node1', 'node2'))

  # SGE's host file.
  host <- Sys.getenv('HOSTNAME', Sys.info()[['nodename']])
  pe   <- tempfile()
  writeLines(c(paste(host, '3 all.q@x UNDEFINED'), 'other.host 5 all.q@y UNDEFINED'), pe)
  expect_identical(case(list(PE_HOSTFILE = pe, NSLOTS = '8')), c(3L, 3L))
  pe2 <- tempfile()
  writeLines('other.host 5 all.q@y UNDEFINED', pe2)
  expect_identical(case(list(PE_HOSTFILE = pe2, NSLOTS = '6')), c(6L, 6L))
  unlink(c(pe, pe2))

  # HTCondor's ads.
  machine <- tempfile()
  job     <- tempfile()
  writeLines(c('# Cpus = 1', 'Cpus = 3'), machine)
  writeLines(c('CpusProvisioned = 2', 'Cpus = "5"'), job)
  expect_identical(case(list(`_CONDOR_MACHINE_AD` = machine)), c(3L, 3L))
  expect_identical(case(list(`_CONDOR_JOB_AD` = job)),         c(2L, 2L))
  job2 <- tempfile()
  writeLines('Cpus = 4.6', job2)
  expect_identical(case(list(`_CONDOR_JOB_AD` = job2)),        c(5L, 5L))
  unlink(c(machine, job, job2))

  # parallelly's settings, as options or environment variables.
  expect_identical(case(opts = list(parallelly.availableCores.fallback = 1)), c(8L, 1L))
  expect_identical(case(list(R_PARALLELLY_AVAILABLECORES_FALLBACK = '2')),   c(8L, 2L))
  expect_identical(case(list(R_FUTURE_AVAILABLECORES_FALLBACK = '2')),       c(8L, 2L))
  expect_identical(case(opts = list(future.availableCores.fallback = 3)),    c(8L, 3L))
  # The fallback gives way to other settings, and to a limiting affinity.
  expect_identical(case(opts = list(parallelly.availableCores.fallback = 1, mc.cores = 4)), c(8L, 4L))
  expect_identical(case(
    opts = list(parallelly.availableCores.fallback = 1),
    info = c(8L, 4L, 3L, NA, 0L, 1L) ), c(3L, 3L))

  expect_identical(case(opts = list(
    parallelly.availableCores.methods.excludes = 'mc.cores', mc.cores = 2 )), c(8L, 8L))
  expect_identical(case(
    list(R_PARALLELLY_AVAILABLECORES_METHODS_EXCLUDES = 'Slurm, mc.cores', SLURM_CPUS_ON_NODE = '2'),
    list(mc.cores = 3) ), c(8L, 8L))
  expect_identical(case(
    opts = list(parallelly.availableCores.methods = c('system', 'mc.cores'), mc.cores = 2),
    info = c(8L, 4L, 3L, NA, 0L, 1L) ), c(8L, 2L))
  expect_identical(case(
    list(MY_CPUS = '3'),
    list(parallelly.availableCores.methods = c('system', 'MY_CPUS')) ), c(8L, 3L))

  expect_identical(case(opts = list(parallelly.availableCores.system = 4)),         c(4L, 4L))
  expect_identical(case(opts = list(parallelly.availableCores.logical = FALSE)),    c(4L, 4L))
  expect_identical(case(list(R_PARALLELLY_AVAILABLECORES_LOGICAL = 'FALSE')),      c(4L, 4L))
  expect_identical(case(opts = list(parallelly.availableCores.threadsPerWorker = 2)), c(8L, 4L))
  expect_identical(case(opts = list(parallelly.availableCores.fraction = 0.5)),     c(8L, 4L))
  expect_identical(case(opts = list(parallelly.availableCores.omit = 3)),           c(8L, 5L))
  expect_identical(case(opts = list(parallelly.availableCores.omit = 30)),          c(8L, 1L))
  expect_identical(case(opts = list(parallelly.availableCores.max = 6)),            c(8L, 6L))
  expect_identical(case(opts = list(parallelly.availableCores.min = 3, mc.cores = 1)), c(8L, 3L))
  expect_identical(case(opts = list(parallelly.availableCores.custom = function () 3L)), c(8L, 3L))
  expect_identical(case(opts = list(parallelly.availableCores.custom = function () stop('x'))), c(8L, 8L))
  expect_identical(case(list(R_PARALLELLY_AVAILABLECORES_OMIT = 'abc')),           c(8L, 8L))



  # validate_cpus(). ====

  f <- function (cpus = n_cpus()) { validate_cpus(); cpus }
  g <- function (cpus = n_cpus()) f(cpus = cpus)

  old <- options(mc.cores = 1)
  expect_identical(f(), 1L)
  expect_identical(g(), 1L)
  expect_identical(f(2), min(2L, max_cpus()))
  options(old)

  old <- options(ecodive.cpus = 1)
  expect_identical(f(4), 1L)
  expect_identical(f(), 1L)
  options(old)

  expect_identical(f(.Machine$integer.max), max_cpus())
  expect_error(f(0), '`cpus`')

  env <- new.env()
  env$cpus <- 1e6
  validate_cpus(env)
  expect_identical(env$cpus, max_cpus())



  # Worker processes. ====

  if (requireNamespace('parallel', quietly = TRUE)) {

    # A forked child.
    if (.Platform$OS.type == 'unix') {
      .Call(C_stop_thread_pool)
      job <- parallel::mcparallel(c(n_cpus(), max_cpus()))
      expect_identical(parallel::mccollect(job)[[1]], c(1L, max_cpus()))
    }

    # A socket cluster's worker.
    cl <- parallel::makeCluster(1L)
    res <- tryCatch(
      expr = {
        # The worker's own .libPaths(), not a serialized copy of ours.
        parallel::clusterCall(cl, eval, bquote(.libPaths(.(.libPaths()))))
        parallel::clusterEvalQ(cl, c(ecodive::n_cpus(), ecodive:::max_cpus()))[[1]]
      },
      finally = parallel::stopCluster(cl) )
    expect_identical(res[[1]], 1L)
    expect_true(res[[2]] >= 1L)
  }

#})
