# Writes tests/tinytests/expected/fastancom.rds: inputs and reference
# results for test-fastancom.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/fastancom.r
#
# Two stages. With fastANCOM 0.0.4 installed (GitHub only, GPL-3: run for
# its outputs only; set FASTANCOM_LIB to the library holding it), the
# tables are simulated and written to data-raw/expected/fastancom/, and
# fastANCOM() is run on them, its results written there too. Every real
# number in those files is a hexadecimal float, read back exactly. Then,
# with or without fastANCOM, the files are folded into the .rds.
#
# Every case runs fastANCOM() with zero_cut = 1 (no feature filter) and its
# default pseudo = 0.5 unless given. Its standard errors average the
# reference set's residuals over ref.rate * K rather than over the set's
# size (see plans/fastancom.md), so every case's ref.rate * K is a whole
# number or below 1, where the two agree. With covariates, fastANCOM tests
# on n - 2 degrees of freedom; the covariate case also stores its pairwise
# statistics (`joint$statistic`), from which the tests count W on n - p.
#
# decisions.csv holds, for 600 random two-group tables, fastANCOM's W, its
# effect = FALSE `detected` (the original ANCOM's cutoff) and its REJECT,
# for the decision rules alone.

dir  <- file.path('data-raw', 'expected', 'fastancom')
hexs <- function (x) sprintf('%a', x)
hex  <- function (x) as.numeric(x)
csv  <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

lib <- Sys.getenv('FASTANCOM_LIB')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))

if (requireNamespace('fastANCOM', quietly = TRUE)) {

  stopifnot(packageVersion('fastANCOM') == '0.0.4')
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  fa <- function (...) suppressMessages(suppressWarnings(
    fastANCOM::fastANCOM(..., zero_cut = 1) ))

  write_case <- function (name, Y, meta, fit, stat = FALSE) {
    write.csv(
      data.frame(sample = rownames(Y), apply(Y, 2, hexs), check.names = FALSE),
      file.path(dir, paste0(name, '_counts.csv')), row.names = FALSE )
    m <- meta
    for (v in names(m)) if (is.numeric(m[[v]])) m[[v]] <- hexs(m[[v]])
    write.csv(m, file.path(dir, paste0(name, '_meta.csv')), row.names = FALSE)
    r <- fit$results
    f <- r$final
    write.csv(
      data.frame(
        feature = colnames(Y),
        beta    = hexs(r$marginal$beta),
        log2FC  = hexs(f$log2FC),
        SD      = hexs(f$log2FC.SD),
        pval    = hexs(f$log2FC.pval),
        qval    = hexs(f$log2FC.qval),
        W       = as.integer(r$joint$W),
        REJECT  = f$REJECT ),
      file.path(dir, paste0(name, '_final.csv')), row.names = FALSE )
    if (stat)
      write.csv(
        apply(r$joint$statistic, 2, hexs),
        file.path(dir, paste0(name, '_stat.csv')), row.names = FALSE )
  }

  args <- list()


  #--------------------------------------------------------
  # The manual's example: 60 samples, 100 features, the
  # first 10 scaled by 2^b in group 1.
  #--------------------------------------------------------
  set.seed(123)
  Y <- matrix(rpois(100 * 60, 10), 60)
  b <- rnorm(10, 4)
  Y[31:60, 1:10] <- t(t(Y[31:60, 1:10]) * 2^b)
  dimnames(Y) <- list(sprintf('sample%d', 1:60), sprintf('microbe%d', 1:100))
  grp <- rep(c(0, 1), each = 30)
  write_case('manual', Y, data.frame(group = grp), fa(Y = Y, x = grp))
  args$manual <- list(pseudo = 0.5, sig = 0.05, detect = 0.7, ref = 0.05)

  # A numeric variable of interest, on the same table.
  set.seed(7)
  age <- round(rnorm(60, 40, 10), 1)
  write_case('numeric', Y, data.frame(age = age), fa(Y = Y, x = age))
  args$numeric <- list(pseudo = 0.5, sig = 0.05, detect = 0.7, ref = 0.05)


  #--------------------------------------------------------
  # Two groups and a numeric covariate.
  #--------------------------------------------------------
  set.seed(11)
  n <- 30; K <- 40
  Y <- matrix(rnbinom(n * K, mu = 30, size = 2), n)
  g <- rep(c('ctrl', 'trt'), each = 15)
  z <- round(rnorm(n), 3)
  Y[g == 'trt', 1:5] <- round(Y[g == 'trt', 1:5] * 3)
  Y <- round(Y * exp(0.3 * z))
  dimnames(Y) <- list(sprintf('S%02d', 1:n), sprintf('F%02d', 1:K))
  write_case(
    'covariate', Y, data.frame(grp = g, z = z),
    fa(Y = Y, x = factor(g), Z = z), stat = TRUE )
  args$covariate <- list(pseudo = 0.5, sig = 0.05, detect = 0.7, ref = 0.05)


  #--------------------------------------------------------
  # Zeros, other settings, few features.
  #--------------------------------------------------------
  set.seed(5)
  n <- 16; K <- 20
  Y <- matrix(rpois(n * K, 8), n)
  Y[sample(length(Y), 40)] <- 0
  g <- rep(c('a', 'b'), each = 8)
  Y[g == 'b', 1:3] <- Y[g == 'b', 1:3] * 4
  dimnames(Y) <- list(sprintf('S%02d', 1:n), sprintf('F%02d', 1:K))
  write_case(
    'zeros', Y, data.frame(grp = g),
    fa(Y = Y, x = factor(g), sig = 0.1, detect.rate = 0.6, ref.rate = 0.1) )
  args$zeros <- list(pseudo = 0.5, sig = 0.1, detect = 0.6, ref = 0.1)

  write_case('pseudo1', Y, data.frame(grp = g), fa(Y = Y, x = factor(g), pseudo = 1))
  args$pseudo1 <- list(pseudo = 1, sig = 0.05, detect = 0.7, ref = 0.05)

  Yp <- Y + 1
  write_case('nopseudo', Yp, data.frame(grp = g), fa(Y = Yp, x = factor(g), pseudo = 0))
  args$nopseudo <- list(pseudo = 0, sig = 0.05, detect = 0.7, ref = 0.05)

  P <- Y / rowSums(Y)
  write_case('props', P, data.frame(grp = g), fa(Y = P, x = factor(g)))
  args$props <- list(pseudo = 0.5, sig = 0.05, detect = 0.7, ref = 0.05)

  set.seed(9)
  n <- 12; K <- 8
  Y <- matrix(rpois(n * K, 30), n)
  g <- rep(c('a', 'b'), each = 6)
  Y[g == 'b', 1:2] <- Y[g == 'b', 1:2] * 3
  dimnames(Y) <- list(sprintf('S%02d', 1:n), sprintf('F%02d', 1:K))
  write_case('small', Y, data.frame(grp = g), fa(Y = Y, x = factor(g), detect.rate = 0.5))
  args$small <- list(pseudo = 0.5, sig = 0.05, detect = 0.5, ref = 0.05)

  # ref.rate * K = 28.999999999999996: 29 reference features (R's `:`).
  set.seed(21)
  n <- 20; K <- 100
  Y <- matrix(rpois(n * K, 25), n)
  g <- rep(c('a', 'b'), each = 10)
  Y[g == 'b', 1:15] <- Y[g == 'b', 1:15] * 2
  dimnames(Y) <- list(sprintf('S%02d', 1:n), sprintf('F%03d', 1:K))
  write_case('ref29', Y, data.frame(grp = g), fa(Y = Y, x = factor(g), ref.rate = 0.29))
  args$ref29 <- list(pseudo = 0.5, sig = 0.05, detect = 0.7, ref = 0.29)

  # Ties in W around the reference set, and W on the detect.rate line.
  set.seed(31)
  n <- 10; K <- 30
  Y <- matrix(rpois(n * K, 40), n)
  g <- rep(c('a', 'b'), each = 5)
  Y[g == 'b', 1:9] <- Y[g == 'b', 1:9] * 3
  dimnames(Y) <- list(sprintf('S%02d', 1:n), sprintf('F%02d', 1:K))
  write_case('ties', Y, data.frame(grp = g), fa(Y = Y, x = factor(g), detect.rate = 0.6, ref.rate = 0.1))
  args$ties <- list(pseudo = 0.5, sig = 0.05, detect = 0.6, ref = 0.1)

  saveRDS(args, file.path(dir, 'args.rds'), version = 2)


  #--------------------------------------------------------
  # The decision rules: 600 random tables.
  #--------------------------------------------------------
  dec <- vector('list', 600)
  for (s in seq_along(dec)) {
    set.seed(s)
    K  <- sample(c(4:60, 100, 127, 177), 1)
    n  <- sample(6:16, 1)
    a  <- round(runif(1, 0.05, 0.95), 2)
    Y  <- matrix(rpois(n * K, 30), n)
    nd <- sample(0:max(1, K %/% 3), 1)
    g  <- rep(0:1, length.out = n)
    if (nd) Y[g == 1, 1:nd] <- Y[g == 1, 1:nd] * sample(c(1.5, 2, 3), 1)
    f0 <- fa(Y = Y, x = g, detect.rate = a, effect = FALSE)$results$final
    f1 <- fa(Y = Y, x = g, detect.rate = a)$results$final
    stopifnot(all(f0$W == f1$Reject.number))
    dec[[s]] <- data.frame(
      K = K, detect_rate = a,
      W        = paste(f0$W, collapse = ' '),
      detected = paste(as.integer(f0$detected), collapse = ''),
      reject   = paste(as.integer(f1$REJECT), collapse = '') )
  }
  write.csv(do.call(rbind, dec), file.path(dir, 'decisions.csv'), row.names = FALSE)
}


#--------------------------------------------------------
# Fold the files into the .rds.
#--------------------------------------------------------
args  <- readRDS(file.path(dir, 'args.rds'))
cases <- list()
for (name in names(args)) {
  cnt <- csv(paste0(name, '_counts.csv'))
  Y   <- sapply(cnt[-1], hex)
  dimnames(Y) <- list(cnt$sample, names(cnt)[-1])
  meta <- read.csv(file.path(dir, paste0(name, '_meta.csv')), colClasses = 'character')
  for (v in names(meta)) if (grepl('^-?0x', meta[[v]][[1]])) meta[[v]] <- hex(meta[[v]])
  rownames(meta) <- rownames(Y)
  fin <- csv(paste0(name, '_final.csv'))
  out <- list(
    counts = Y, meta = meta, args = args[[name]],
    beta   = hex(fin$beta),    log2FC = hex(fin$log2FC), SD = hex(fin$SD),
    pval   = hex(fin$pval),    qval   = hex(fin$qval),
    W      = as.integer(fin$W), REJECT = as.logical(fin$REJECT) )
  sf <- file.path(dir, paste0(name, '_stat.csv'))
  if (file.exists(sf)) out$stat <- unname(as.matrix(sapply(csv(basename(sf)), hex)))
  cases[[name]] <- out
}

dec <- read.csv(file.path(dir, 'decisions.csv'), colClasses = 'character')
decisions <- lapply(seq_len(nrow(dec)), function (i) list(
  K           = as.integer(dec$K[[i]]),
  detect_rate = as.numeric(dec$detect_rate[[i]]),
  W           = as.integer(strsplit(dec$W[[i]], ' ')[[1]]),
  detected    = strsplit(dec$detected[[i]], '')[[1]] == '1',
  reject      = strsplit(dec$reject[[i]], '')[[1]] == '1' ))

saveRDS(
  list(cases = cases, decisions = decisions),
  file.path('tests', 'tinytests', 'expected', 'fastancom.rds'), version = 2 )
