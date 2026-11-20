# Number of CPU Cores

How many threads ecodive's functions use by default: the CPU cores this
R process may use, within any limits set on it.

## Usage

``` r
n_cpus()
```

## Value

A scalar integer, at least `1`. Always `1` where ecodive was built
without `pthread` support.

## Details

`n_cpus()` is the default for every function's `cpus` argument. It
follows the checks of `parallelly::availableCores()`, without needing
the `parallelly` package, and takes about as long as reading a few
environment variables. It is the smallest of the values below that
apply.

**Resource limits.** A `cpus` given explicitly is capped at the smallest
of these alone:

- The logical cores online (physical cores with
  `options(parallelly.availableCores.logical = FALSE)`).

- The cores the process may run on: its CPU affinity, as set by
  `taskset`, `numactl`, cgroup cpusets (Docker's `--cpuset-cpus`), job
  schedulers that bind tasks to cores, or Windows' `start /affinity` and
  job objects.

- A CPU quota, in whole cores: Linux cgroup v1 and v2 quotas (Docker's
  `--cpus`, Kubernetes CPU limits, systemd's `CPUQuota`) or a Windows
  job object's hard CPU rate cap.

- The cores a job scheduler allocated on this node: Slurm
  (`SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, ...), PBS/Torque
  (`PBS_NUM_PPN`, `NCPUS`), LSF (`LSB_DJOB_NUMPROC`), SGE (the
  `PE_HOSTFILE`, `NSLOTS`), Fujitsu PJM (`PJM_VNODE_CORE`,
  `PJM_PROC_BY_NODE`) and HTCondor (`_CONDOR_MACHINE_AD`,
  `_CONDOR_JOB_AD`).

- 2 under `R CMD check --as-cran` (`_R_CHECK_LIMIT_CORES_`), and 4 on
  Bioconductor's build machines.

- `options(ecodive.cpus = n)`, or the environment variable
  `R_ECODIVE_CPUS`.

**Defaults.** These lower `n_cpus()`, but not the cap on an explicit
`cpus`:

- The `mc.cores` option, or the `MC_CORES` environment variable.

- `BIOCPARALLEL_WORKER_NUMBER`.

- 1 in a worker process: one forked by, e.g.,
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  future's multicore plan or BiocParallel's `MulticoreParam` after
  ecodive was loaded, or a worker of a socket cluster (`parallel`,
  `parallelly`, future's multisession plan, `snow`), `mirai` or
  `clustermq`. Those workers already run side by side, so more threads
  in each would compete for the same cores. Give `cpus` to use more.

- These options of `parallelly` (`parallelly.availableCores.*`, also
  read as `future.availableCores.*`, or the environment variables
  `R_PARALLELLY_AVAILABLECORES_*`), whether or not `parallelly` is
  installed: `methods` and `methods.excludes` (which of parallelly's
  methods above apply), `fallback`, `custom`, `min`, `system`,
  `logical`, `threadsPerWorker`, `fraction`, `omit` and `max`.

Settings that are missing, not numbers or less than 1 are ignored.
Results are the same for any number of threads.

## Examples

``` r
    n_cpus()
#> [1] 4
```
