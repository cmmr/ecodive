# Writes tests/tinytests/expected/linda.rds: inputs and reference results
# for test-linda.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/linda.r
#
# Two stages. With MicrobiomeStat installed (set MICROBIOMESTAT_LIB to the
# library holding it), the tables and metadata are simulated and written to
# data-raw/expected/linda/, and MicrobiomeStat 1.4's linda() (CRAN, GPL-3:
# run for its outputs only) and modeest 2.5.0's shorth() and meanshift()
# (GPL-3, likewise) are run on them, their results written there too. Every
# real number in those files is a hexadecimal float, read back exactly.
# Then, with or without MicrobiomeStat, the files are folded into the .rds.
#
# Each case runs linda() with adaptive = FALSE and an explicit
# zero.handling (MicrobiomeStat 1.4's adaptive mode reports its choice but
# adds the pseudocount whatever it chose), pseudo.cnt = 0.5, no filters,
# and p.adj.method = 'BH'. The adaptive cases also record the choice that
# linda(adaptive = TRUE) reports ("Imputation approach is used." or
# "Pseudo-count approach is used."), and their results are those of the
# explicit run with that choice.

dir <- file.path('data-raw', 'expected', 'linda')
hexs <- function (x) sprintf('%a', x)
hex  <- function (x) as.numeric(x)
csv  <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

lib <- Sys.getenv('MICROBIOMESTAT_LIB')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))

if (requireNamespace('MicrobiomeStat', quietly = TRUE)) {

  stopifnot(packageVersion('MicrobiomeStat') == '1.4')
  stopifnot(packageVersion('modeest') == '2.5.0')
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)


  #--------------------------------------------------------
  # Tables (features x samples, as LinDA takes them) and
  # metadata.
  #--------------------------------------------------------
  sim <- function (seed, n, D, mu = 20, size = 0.7, depth_shift = 0, sparse = 0) {
    set.seed(seed)
    g   <- rep(c('ctrl', 'trt'), length.out = n)
    lib <- exp(rnorm(n, 0, 0.5) + depth_shift * (g == 'trt'))
    Y   <- matrix(rnbinom(n * D, mu = mu, size = size), D, n) *
             rep(1, D) %o% lib
    Y   <- round(Y)
    # Differentially abundant features (10%).
    da  <- seq_len(max(1, D %/% 10))
    Y[da, g == 'trt'] <- round(Y[da, g == 'trt'] * 4)
    if (sparse > 0) Y[matrix(runif(n * D) < sparse, D, n)] <- 0
    dimnames(Y) <- list(sprintf('F%02d', seq_len(D)), sprintf('S%02d', seq_len(n)))
    meta <- data.frame(
      g   = g,
      x   = round(rnorm(n, 50, 10), 2),
      c1  = round(runif(n, 18, 70)),
      c2  = sample(c('north', 'east', 'south'), n, TRUE),
      flg = sample(c(TRUE, FALSE), n, TRUE),
      row.names = colnames(Y), stringsAsFactors = FALSE )
    list(Y = Y, meta = meta)
  }

  tabs <- list(
    two    = sim(1, 40, 80),
    conf   = sim(2, 36, 60, depth_shift = 0),
    depth  = sim(3, 40, 70, depth_shift = 1.2),
    sparse = sim(4, 30, 90, mu = 5, size = 0.3, sparse = 0.4),
    small  = sim(5, 12, 2) )
  # Proportions: closed, with their zeros.
  tabs$prop <- tabs$two
  tabs$prop$Y <- t(t(tabs$two$Y) / colSums(tabs$two$Y))
  # Rare features: a single non-zero count each, many alike, so that the
  # mode's starting interval has ties.
  tabs$rare <- sim(6, 24, 50)
  tabs$rare$Y[31:50, ] <- 0
  for (f in 31:50) tabs$rare$Y[f, 1 + (f %% 3)] <- 1

  for (nm in names(tabs)) {
    Y <- tabs[[nm]]$Y
    write.csv(matrix(hexs(Y), nrow(Y)), file.path(dir, sprintf('table-%s.csv', nm)), row.names = FALSE)
    m <- tabs[[nm]]$meta
    m$x <- hexs(m$x)
    write.csv(m, file.path(dir, sprintf('meta-%s.csv', nm)), row.names = FALSE)
  }


  #--------------------------------------------------------
  # The cases.
  #--------------------------------------------------------
  cases <- read.table(header = TRUE, stringsAsFactors = FALSE, text = '
    case          table   formula           zero         winsor
    two_pc        two     ~g                pseudo-count TRUE
    two_imp       two     ~g                imputation   TRUE
    two_nowin     two     ~g                pseudo-count FALSE
    num_pc        two     ~x                pseudo-count TRUE
    num_imp       two     ~x                imputation   TRUE
    conf_pc       conf    ~g+c1+c2          pseudo-count TRUE
    conf_imp      conf    ~g+c1+c2          imputation   TRUE
    inter         conf    ~g*x              pseudo-count TRUE
    flag          conf    ~flg+c1           imputation   FALSE
    depth_adapt   depth   ~g                adaptive     TRUE
    depth_adapt2  depth   ~g+x              adaptive     TRUE
    two_adapt     two     ~g+x              adaptive     TRUE
    conf_adapt    two     ~g+c1+c2          adaptive     TRUE
    sparse_adapt  sparse  ~g                adaptive     TRUE
    sparse_pc     sparse  ~g+x              pseudo-count TRUE
    sparse_imp    sparse  ~g+x              imputation   TRUE
    prop          prop    ~g+x              proportion   TRUE
    prop_nowin    prop    ~g                proportion   FALSE
    rare          rare    ~g                pseudo-count TRUE
    small         small   ~g+x              pseudo-count TRUE
  ')

  results <- NULL
  extra   <- NULL
  for (i in seq_len(nrow(cases))) {

    cs  <- cases[i, ]
    tb  <- tabs[[cs$table]]
    typ <- if (cs$zero == 'proportion') 'proportion' else 'count'

    chosen <- cs$zero
    if (cs$zero == 'adaptive') {
      out <- capture.output(suppressWarnings(MicrobiomeStat::linda(
        tb$Y, tb$meta, cs$formula, feature.dat.type = 'count',
        is.winsor = cs$winsor, adaptive = TRUE, verbose = TRUE )))
      chosen <- if (any(grepl('Imputation approach is used', out))) 'imputation'
                else if (any(grepl('Pseudo-count approach is used', out))) 'pseudo-count'
                else stop('No choice reported.')
    }

    r <- suppressWarnings(MicrobiomeStat::linda(
      tb$Y, tb$meta, cs$formula, feature.dat.type = typ,
      is.winsor = cs$winsor, outlier.pct = 0.03, adaptive = FALSE,
      zero.handling = if (typ == 'count') chosen else 'pseudo-count',
      pseudo.cnt = 0.5, p.adj.method = 'BH', alpha = 0.05, verbose = FALSE ))

    for (v in r$variables) {
      o <- r$output[[v]]
      results <- rbind(results, data.frame(
        case = cs$case, covariate = v, feature = rownames(o),
        log2_fc = hexs(o$log2FoldChange), se = hexs(o$lfcSE),
        t = hexs(o$stat), p_value = hexs(o$pvalue), p_adjusted = hexs(o$padj),
        reject = as.integer(o$reject), dof = o$df,
        biased = hexs(r$wald$beta[v, rownames(o)]),
        stringsAsFactors = FALSE ))
    }

    # Covariances between coefficients (LinDA: only with two or more
    # variables), and each feature's residual standard deviation.
    cv <- ''
    if (length(r$variables) > 1L) {
      v1 <- r$variables[[1]]; v2 <- r$variables[[2]]
      cv <- paste(hexs(r$covariance[[v1]][[v2]]), collapse = ' ')
    }
    extra <- rbind(extra, data.frame(
      case    = cs$case, chosen = chosen,
      bias    = paste(hexs(r$bias), collapse = ' '),
      cov12   = cv,
      sig     = paste(hexs(r$wald$sig), collapse = ' '),
      stringsAsFactors = FALSE ))
  }

  cases$chosen <- extra$chosen
  write.csv(cases,   file.path(dir, 'cases.csv'),   row.names = FALSE)
  write.csv(results, file.path(dir, 'results.csv'), row.names = FALSE)
  write.csv(extra[, -2], file.path(dir, 'extra.csv'), row.names = FALSE)


  #--------------------------------------------------------
  # modeest's shorth() and meanshift() (Gaussian kernel,
  # bw.nrd0, from the shorth), the mode LinDA uses, on
  # vectors with and without ties.
  #--------------------------------------------------------
  set.seed(7)
  vecs <- list()
  for (k in 1:40) {
    n <- sample(c(2:12, 30, 100, 500), 1)
    v <- switch(
      1 + k %% 4,
      rnorm(n),
      c(rnorm(n, 0, 0.3), rnorm(n %/% 4 + 1, 3)),
      round(rnorm(n), 1),
      rexp(n) * 10^runif(1, -4, 3) )
    vecs[[k]] <- v
  }
  vecs[[41]] <- c(-2, -1, 0, 1, 2)
  vecs[[42]] <- rep(0.25, 7)
  modes <- data.frame(
    values  = vapply(vecs, function (v) paste(hexs(v), collapse = ' '), ''),
    shorth  = vapply(vecs, function (v) hexs(suppressWarnings(modeest::shorth(v))), ''),
    mode    = vapply(vecs, function (v) hexs(suppressWarnings(
                as.numeric(modeest::meanshift(v, bw = stats::bw.nrd0(v))))), ''),
    stringsAsFactors = FALSE )
  write.csv(modes, file.path(dir, 'modes.csv'), row.names = FALSE)
}



#--------------------------------------------------------
# Fold the text files into the .rds.
#--------------------------------------------------------
words <- function (x) if (nzchar(x)) hex(strsplit(x, ' ', fixed = TRUE)[[1]]) else numeric(0)

cases   <- csv('cases.csv')
results <- csv('results.csv')
extra   <- csv('extra.csv')
modes   <- csv('modes.csv')

tables <- list()
for (nm in unique(cases$table)) {
  m <- as.matrix(read.csv(file.path(dir, sprintf('table-%s.csv', nm)), colClasses = 'character'))
  Y <- matrix(hex(m), nrow(m), ncol(m))
  D <- nrow(Y); n <- ncol(Y)
  x <- t(Y)   # samples x features, as ecodive takes them
  dimnames(x) <- list(sprintf('S%02d', seq_len(n)), sprintf('F%02d', seq_len(D)))
  meta <- read.csv(file.path(dir, sprintf('meta-%s.csv', nm)), stringsAsFactors = FALSE)
  meta$x   <- hex(meta$x)
  meta$flg <- as.logical(meta$flg)
  rownames(meta) <- rownames(x)
  tables[[nm]] <- list(counts = x, meta = meta)
}

ref <- lapply(setNames(nm = cases$case), function (cs) {
  cc <- cases[cases$case == cs, ]
  r  <- results[results$case == cs, ]
  e  <- extra[extra$case == cs, ]
  list(
    table   = cc$table,
    formula = cc$formula,
    zero    = cc$zero,
    chosen  = sub('-', '', cc$chosen),
    winsor  = as.logical(cc$winsor),
    result  = data.frame(
      feature    = r$feature,
      covariate  = r$covariate,
      log2_fc    = hex(r$log2_fc),
      se         = hex(r$se),
      t          = hex(r$t),
      dof        = as.numeric(r$dof),
      p_value    = hex(r$p_value),
      p_adjusted = hex(r$p_adjusted),
      reject     = r$reject == '1',
      biased     = hex(r$biased),
      stringsAsFactors = FALSE ),
    bias  = words(e$bias),
    cov12 = words(e$cov12),
    sig   = words(e$sig) )
})

modes <- lapply(seq_len(nrow(modes)), function (i) list(
  values = words(modes$values[[i]]),
  shorth = hex(modes$shorth[[i]]),
  mode   = hex(modes$mode[[i]]) ))

saveRDS(
  list(tables = tables, cases = ref, modes = modes),
  file.path('tests', 'tinytests', 'expected', 'linda.rds'),
  version = 2 )
