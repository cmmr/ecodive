# Writes tests/tinytests/expected/sourcetracker.rds: inputs and reference
# results for test-sourcetracker.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     python -I data-raw/expected/sourcetracker.py <sourcetracker2 dir>
#     SOURCETRACKER_V1=<SourceTracker.r> Rscript data-raw/expected/sourcetracker.r
#
# (The v1 part takes about 3 minutes on 4 cores.)
#
# SourceTracker2 2.0.1-dev (BSD-3): collapse_source_data(),
# ConditionalProbability and gibbs_sampler()'s draw schedule on fixed
# inputs, and gibbs() on simulated mixtures with 40 seeds per setting (see
# the Python script). Every real number in its outputs is a hexadecimal
# float, read back exactly.
#
# SourceTracker 1.0 (GPL; github.com/danknights/sourcetracker,
# src/SourceTracker.r): run for its outputs only, on the same simulated
# tables, for a loose cross-check. Give the path of SourceTracker.r in
# SOURCETRACKER_V1 (it is source()d, and uses R's own random numbers);
# without it, the v1 part is left out.

dir <- file.path('data-raw', 'expected', 'sourcetracker')

csv   <- function (f, ...) read.csv(file.path(dir, f), ...)
words <- function (x) strsplit(x, ' ', fixed = TRUE)[[1]]
hex   <- function (x) as.numeric(words(x))
mat   <- function (f) {
  m <- as.matrix(csv(f, header = FALSE))
  dimnames(m) <- NULL
  storage.mode(m) <- 'double'
  m
}


#--------------------------------------------------------
# Deterministic parts.
#--------------------------------------------------------
col <- csv('collapse-result.csv', check.names = FALSE)
collapse <- list(
  table    = mat('collapse-table.csv'),
  env      = csv('collapse-env.csv', colClasses = 'character')$env,
  expected = `dimnames<-`(as.matrix(col[, -1]), list(col$env, names(col)[-1])) )
storage.mode(collapse$expected) <- 'double'

cp <- csv('cp.csv', colClasses = 'character')
cond_prob <- list(
  source = mat('cp-source.csv'),
  cases  = lapply(seq_len(nrow(cp)), function (i) {
    r <- cp[i, ]
    list(
      alpha1 = as.numeric(r$alpha1), alpha2 = as.numeric(r$alpha2),
      beta = as.numeric(r$beta), n = as.numeric(r$n), xi = as.integer(r$xi),
      m_xiV = as.numeric(r$m_xiV), m_V = as.numeric(r$m_V),
      n_vnoti = as.integer(words(r$n_vnoti)), expected = hex(r$result) )
  }) )

sch <- csv('schedule.csv', colClasses = 'character')
schedule <- lapply(seq_len(nrow(sch)), function (i) list(
  burnin = as.integer(sch$burnin[[i]]),
  draws  = as.integer(sch$draws_per_restart[[i]]),
  delay  = as.integer(sch$delay[[i]]),
  passes = as.integer(words(sch$passes[[i]])) ))


#--------------------------------------------------------
# Simulated mixtures: SourceTracker2's proportions over
# seeds, as the mean and sd between seeds per sink, the
# mean of its reported sd (mps), and seed 1's feature
# tables (sums over the draws).
#--------------------------------------------------------
envs   <- c('soil', 'water', 'skin')
sinks  <- mat('sim-sinks.csv')
srcs   <- mat('sim-sources.csv')
truth  <- mat('sim-truth.csv')
F      <- ncol(sinks)
feats  <- sprintf('f%d', seq_len(F) - 1L)
snames <- sprintf('sink%d', seq_len(nrow(sinks)))
dimnames(sinks) <- list(snames, feats)
dimnames(srcs)  <- list(envs, feats)
dimnames(truth) <- list(snames, c(envs, 'Unknown'))

runs <- list(
  default = list(n_restarts = 10L, n_draws = 1L, n_burnin = 100L, thin = 1L, loo = FALSE),
  thinned = list(n_restarts = 4L,  n_draws = 3L, n_burnin = 20L,  thin = 5L, loo = FALSE),
  loo     = list(n_restarts = 10L, n_draws = 1L, n_burnin = 100L, thin = 1L, loo = TRUE) )

for (run in names(runs)) {
  r    <- csv(sprintf('sim-%s.csv', run), colClasses = 'character')
  rows <- unique(r$sink)
  cols <- c(envs, 'Unknown')
  P    <- t(sapply(r$proportions, hex, USE.NAMES = FALSE))
  S    <- t(sapply(r$sd, hex, USE.NAMES = FALSE))
  agg  <- function (M, f) {
    out <- t(sapply(rows, function (s) apply(M[r$sink == s, , drop = FALSE], 2, f)))
    dimnames(out) <- list(rows, cols)
    out
  }
  runs[[run]]$n_seeds <- length(unique(r$seed))
  runs[[run]]$mean    <- agg(P, mean)
  runs[[run]]$sd      <- agg(P, stats::sd)
  runs[[run]]$mps     <- agg(S, mean)

  fas <- csv(sprintf('sim-%s-fas.csv', run), check.names = FALSE, colClasses = c('character', 'character', rep('numeric', F)))
  runs[[run]]$fas <- lapply(setNames(nm = rows), function (s) {
    m <- as.matrix(fas[fas$sink == s, -(1:2)])
    dimnames(m) <- list(fas$source[fas$sink == s], feats)
    m
  })
}


#--------------------------------------------------------
# SourceTracker 1 (output only): 20 seeds of R's own
# random numbers, with SourceTracker2's defaults (delay 1).
#--------------------------------------------------------
v1 <- NULL
v1_path <- Sys.getenv('SOURCETRACKER_V1')
if (nzchar(v1_path)) {
  st1 <- new.env()
  sys.source(v1_path, envir = st1)
  # Each seed seeds R's generator itself, so forked workers give the same
  # results as a serial run.
  props <- parallel::mclapply(1:20, mc.cores = 4L, function (seed) {
    set.seed(seed)
    fit <- st1$sourcetracker(srcs, envs, rarefaction_depth = 2000)
    res <- st1$predict.sourcetracker(
      fit, sinks, burnin = 100, nrestarts = 10, ndraws.per.restart = 1,
      delay = 1, alpha1 = 1e-3, alpha2 = 0.1, beta = 10,
      rarefaction_depth = 400, verbosity = 0 )
    res$proportions[, c(envs, 'Unknown')]
  })
  v1 <- list(
    n_seeds = length(props),
    mean    = Reduce(`+`, props) / length(props),
    sd      = sqrt(Reduce(`+`, lapply(props, function (p) p^2)) / length(props) -
                     (Reduce(`+`, props) / length(props))^2) )
}


saveRDS(
  list(
    collapse  = collapse,
    cond_prob = cond_prob,
    schedule  = schedule,
    sim       = list(sources = srcs, sinks = sinks, truth = truth, runs = runs, v1 = v1) ),
  file.path('tests', 'tinytests', 'expected', 'sourcetracker.rds'),
  version = 2 )
