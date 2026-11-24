# Writes tests/tinytests/expected/decontam.rds: inputs and reference
# results for test-decontam.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     DECONTAM_LIB=<library with decontam> Rscript data-raw/expected/decontam.r
#
# Two stages. With decontam installed (set DECONTAM_LIB to the library
# holding it), the tables and sample data are simulated and written to
# data-raw/expected/decontam/, and decontam's isContaminant() and
# isNotContaminant() (Artistic-2.0; Bioconductor) are run on them, their
# results written there too. Every real number in those files is a
# hexadecimal float, read back exactly. Then, with or without decontam, the
# files are folded into the .rds.
#
# The committed files come from decontam 1.32.0 (Bioconductor 3.23),
# installed from its source tarball into R 4.3.3; 1.22.0 (Bioconductor
# 3.18, the newest release for R 4.3) has the same R/decontam.R, byte for
# byte, and gives the same files.
#
# Only simulated tables are used: decontam's own example data
# (inst/extdata) are not copied into ecodive.

dir  <- file.path('data-raw', 'expected', 'decontam')
hexs <- function (x) ifelse(is.na(x), 'NA', sprintf('%a', x))
hex  <- function (x) suppressWarnings(as.numeric(x))
csv  <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

lib <- Sys.getenv('DECONTAM_LIB')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))

if (requireNamespace('decontam', quietly = TRUE)) {

  cat('decontam', format(packageVersion('decontam')), '\n')
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(format(packageVersion('decontam')), file.path(dir, 'version.txt'))


  #--------------------------------------------------------
  # Tables (samples x features, as decontam takes them) and
  # sample data.
  #--------------------------------------------------------
  sim <- function (seed, n, n_neg, D, n_contam, mu = 50, size = 0.5, low = FALSE) {
    set.seed(seed)
    neg  <- rep(c(FALSE, TRUE), c(n - n_neg, n_neg))
    conc <- signif(exp(runif(n, 0, 4)), 6)
    conc[neg] <- signif(runif(n_neg, 0.2, 1), 6)
    X <- matrix(rnbinom(n * D, mu = mu, size = size), n, D)
    # Contaminants: a fixed amount of DNA, diluted by the sample's own.
    X[, seq_len(n_contam)] <- rpois(n * n_contam, 200 / conc)
    # True features: mostly absent from the controls.
    keep <- (n_contam + 1):D
    X[neg, keep] <- X[neg, keep] * rbinom(n_neg * length(keep), 1, if (low) 0.6 else 0.1)
    dimnames(X) <- list(sprintf('S%02d', seq_len(n)), sprintf('F%02d', seq_len(D)))
    batch <- rep(c('run1', 'run2'), length.out = n)
    list(X = X, conc = conc, neg = neg, batch = batch)
  }

  tabs <- list(
    main = sim(1, 40, 8, 60, 8),
    low  = sim(2, 24, 8, 40, 15, mu = 5, size = 0.3, low = TRUE) )

  # Three batches, one of them without controls, and features present in
  # only one batch (so that 'product' and 'fisher' meet NA scores).
  b3 <- sim(3, 30, 6, 30, 5)
  b3$batch <- rep(c('A', 'B', 'C'), each = 10)
  b3$neg   <- rep(c(FALSE, TRUE, FALSE, TRUE, FALSE), c(8, 2, 6, 4, 10))
  b3$X[11:30, 26:30] <- 0
  b3$X[c(1, 5), 29] <- c(3, 4)
  tabs$b3 <- b3

  # Edge cases: an empty sample, a feature in one sample, one in none, one
  # in every sample.
  ed <- sim(4, 16, 4, 12, 2)
  ed$X[, 10] <- 0; ed$X[7, 10] <- 5
  ed$X[, 11] <- 0
  ed$X[, 12] <- ed$X[, 12] + 1
  ed$X[3, ] <- 0
  ed$batch <- rep(c('x', 'y'), each = 8)
  tabs$edge <- ed

  for (nm in names(tabs)) {
    t <- tabs[[nm]]
    write.csv(t$X, file.path(dir, paste0(nm, '-counts.csv')))
    write.csv(
      data.frame(sample = rownames(t$X), conc = hexs(t$conc), neg = t$neg, batch = t$batch),
      file.path(dir, paste0(nm, '-samples.csv')), row.names = FALSE )
  }


  #--------------------------------------------------------
  # isContaminant(): every method, batch combination and
  # threshold form; with and without batches; normalized
  # and not.
  #--------------------------------------------------------
  run <- function (nm, method, combine, batched, thr, normalize = TRUE, use = c('conc', 'neg')) {
    t <- tabs[[nm]]
    X <- t$X
    if (!normalize) X <- X / rowSums(X) * 1000   # equal depth
    X[!is.finite(X)] <- 0
    r <- suppressWarnings(suppressMessages(decontam::isContaminant(
      X,
      conc = if ('conc' %in% use) t$conc,
      neg  = if ('neg' %in% use) t$neg,
      method = method, batch = if (batched) t$batch,
      batch.combine = combine, threshold = thr, normalize = normalize, detailed = TRUE )))
    data.frame(
      case = nm, method = method, batch_combine = combine, batched = batched,
      threshold = paste(thr, collapse = ' '), normalize = normalize,
      use = paste(use, collapse = ' '),
      feature = rownames(r), freq = hexs(r$freq), prev = r$prev,
      p_freq = hexs(as.numeric(r$p.freq)), p_prev = hexs(as.numeric(r$p.prev)),
      p = hexs(as.numeric(r$p)), contaminant = r$contaminant )
  }

  out <- list()
  for (nm in c('main', 'b3', 'edge'))
    for (m in c('frequency', 'prevalence', 'combined', 'minimum', 'either', 'both'))
      for (bc in c('minimum', 'product', 'fisher'))
        for (batched in c(FALSE, TRUE)) {
          thr <- if (m %in% c('either', 'both')) c(0.1, 0.5) else 0.1
          out[[length(out) + 1]] <- run(nm, m, bc, batched, thr)
        }
  # One threshold for 'either' and 'both'; another threshold; one input only.
  out[[length(out) + 1]] <- run('main', 'either', 'minimum', FALSE, 0.3)
  out[[length(out) + 1]] <- run('main', 'both',   'minimum', TRUE,  0.3)
  out[[length(out) + 1]] <- run('main', 'prevalence', 'minimum', FALSE, 0.5)
  out[[length(out) + 1]] <- run('main', 'auto', 'minimum', FALSE, 0.2, use = 'conc')
  out[[length(out) + 1]] <- run('main', 'auto', 'fisher',  TRUE,  0.2, use = 'neg')
  out[[length(out) + 1]] <- run('low',  'prevalence', 'minimum', FALSE, 0.5)
  # Not normalized (equal-depth input).
  for (m in c('frequency', 'prevalence', 'combined'))
    out[[length(out) + 1]] <- run('main', m, 'minimum', TRUE, 0.1, normalize = FALSE)

  write.csv(do.call(rbind, out), file.path(dir, 'contaminant.csv'), row.names = FALSE)


  #--------------------------------------------------------
  # isNotContaminant().
  #--------------------------------------------------------
  out <- list()
  for (nm in c('main', 'low', 'edge'))
    for (thr in c(0.5, 0.05)) {
      t <- tabs[[nm]]
      r <- suppressWarnings(decontam::isNotContaminant(t$X, neg = t$neg, threshold = thr, detailed = TRUE))
      out[[length(out) + 1]] <- data.frame(
        case = nm, threshold = thr, feature = rownames(r), freq = hexs(r$freq),
        prev = r$prev, p_prev = hexs(r$p.prev), p = hexs(r$p),
        not_contaminant = r$not.contaminant )
    }
  write.csv(do.call(rbind, out), file.path(dir, 'not_contaminant.csv'), row.names = FALSE)
}


#--------------------------------------------------------
# Fold the files into the .rds.
#--------------------------------------------------------
cases <- c('main', 'low', 'b3', 'edge')
tables <- lapply(stats::setNames(cases, cases), function (nm) {
  X <- read.csv(file.path(dir, paste0(nm, '-counts.csv')), row.names = 1, check.names = FALSE)
  X <- as.matrix(X)
  storage.mode(X) <- 'integer'
  s <- csv(paste0(nm, '-samples.csv'))
  list(
    X     = X,
    conc  = stats::setNames(hex(s$conc), s$sample),
    neg   = stats::setNames(as.logical(s$neg), s$sample),
    batch = stats::setNames(s$batch, s$sample) )
})

num <- function (d) {
  for (k in c('freq', 'p_freq', 'p_prev', 'p'))
    if (k %in% names(d)) d[[k]] <- hex(d[[k]])
  for (k in c('prev')) d[[k]] <- as.integer(d[[k]])
  for (k in c('contaminant', 'not_contaminant', 'batched', 'normalize'))
    if (k %in% names(d)) d[[k]] <- as.logical(d[[k]])
  d
}

saveRDS(
  list(
    version         = readLines(file.path(dir, 'version.txt')),
    tables          = tables,
    contaminant     = num(csv('contaminant.csv')),
    not_contaminant = num(csv('not_contaminant.csv')) ),
  file.path('tests', 'tinytests', 'expected', 'decontam.rds'),
  version = 2 )
