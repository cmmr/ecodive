# Writes tests/tinytests/expected/ancombc.rds: inputs and reference results
# for test-ancombc.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     python data-raw/expected/ancombc.py   # scikit-bio's outputs
#     Rscript data-raw/expected/ancombc.r   # folds them in, and runs ANCOMBC
#
# scikit-bio 0.7.4 (BSD-3, derived from ANCOMBC): ancombc() and ancombc2()
# on each case of data-raw/expected/ancombc/cases.csv, with their biases,
# sampling fractions and fitted quantities. Every number in those files is a
# hexadecimal float, read back exactly.
#
# ANCOMBC (Bioconductor, Artistic-2.0): ancombc() and ancombc2() on the same
# tables, where their settings can match: ancombc() always adds 1, so only
# the cases with `added` = 1; ancombc2() with prv_cut = 0, lib_cut = 0 and
# no sensitivity analysis. ANCOMBC 2.14.0 (the newest release) needs
# R >= 4.5; this was run with 2.8.1, the newest that installs on R 4.3.3,
# whose estimation code (R/ancombc_bias_correct.R: .iter_mle(), .bias_em())
# is the same as 2.14.0's. Set ANCOMBC_LIB to the library holding it.

dir <- file.path('data-raw', 'expected', 'ancombc')

hex <- function (x) as.numeric(x)
csv <- function (f) read.csv(file.path(dir, f), colClasses = 'character')

tabs    <- csv('tables.csv')
cases   <- csv('cases.csv')
results <- csv('results.csv')
bias    <- csv('bias.csv')
fit     <- csv('fit.csv')
vcv     <- csv('vcov.csv')
dofs    <- csv('dof.csv')
sfs     <- csv('sf.csv')

read_table <- function (nm) {
  m <- as.matrix(read.csv(
    file.path(dir, sprintf('table-%s.csv', nm)),
    header = FALSE, colClasses = 'character' ))
  m <- matrix(hex(m), nrow(m), ncol(m))
  if (all(m == round(m))) storage.mode(m) <- 'integer'
  dimnames(m) <- list(sprintf('S%02d', seq_len(nrow(m))), sprintf('F%02d', seq_len(ncol(m))))
  m
}

read_meta <- function (nm) {
  tr <- tabs[tabs$name == nm, ]
  m  <- read.csv(file.path(dir, sprintf('meta-%s.csv', nm)), colClasses = 'character')
  ty <- strsplit(strsplit(tr$types, ' ', fixed = TRUE)[[1]], ':', fixed = TRUE)
  lv <- if (nzchar(tr$levels)) strsplit(strsplit(tr$levels, ' ', fixed = TRUE)[[1]], '=', fixed = TRUE) else list()
  for (t in ty) {
    v <- t[[1]]
    m[[v]] <- switch(
      t[[2]],
      str  = m[[v]],
      num  = hex(m[[v]]),
      bool = as.logical(m[[v]]),
      cat  = factor(m[[v]], levels = strsplit(
        Filter(function (l) l[[1]] == v, lv)[[1]][[2]], '|', fixed = TRUE)[[1]]) )
  }
  rownames(m) <- sprintf('S%02d', seq_len(nrow(m)))
  m
}

# patsy's column names, as model.matrix() names them.
r_name <- function (x) {
  x <- gsub('[T.True]',  'TRUE',  x, fixed = TRUE)
  x <- gsub('[T.False]', 'FALSE', x, fixed = TRUE)
  x <- gsub('([A-Za-z_][A-Za-z0-9_.]*)\\[T\\.([^]]*)\\]', '\\1\\2', x)
  x <- gsub('([A-Za-z_][A-Za-z0-9_.]*)\\[([^]]*)\\]',     '\\1\\2', x)
  x[x == 'Intercept'] <- '(Intercept)'
  x
}

# A long table (feature, covariate, value) as a features x covariates matrix,
# covariates in patsy's order.
as_mat <- function (f, cov, v, D) {
  lev <- unique(cov)
  m   <- matrix(NA_real_, D, length(lev), dimnames = list(NULL, r_name(lev)))
  m[cbind(f, match(cov, lev))] <- v
  m
}


#--------------------------------------------------------
# scikit-bio.
#--------------------------------------------------------
tables <- lapply(setNames(tabs$name, tabs$name), function (nm)
  list(counts = read_table(nm), data = read_meta(nm)) )

ref <- lapply(setNames(cases$case, cases$case), function (cs) {

  cr <- cases[cases$case == cs, ]
  x  <- tables[[cr$table]]$counts
  D  <- ncol(x)
  r  <- results[results$case == cs, ]
  f  <- as.integer(r$feature)

  out <- list(
    table        = cr$table,
    method       = cr$method,
    formula      = cr$formula,
    groups       = if (nzchar(cr$grouping)) cr$grouping,
    pseudocount  = if (cr$method == 'ancombc') hex(cr$added) else hex(cr$pseudocount),
    var_quantile = hex(cr$var_quantile),
    p_adjust     = cr$p_adjust,
    alpha        = hex(cr$alpha),
    max_iters    = as.integer(cr$max_iter),
    tol          = hex(cr$tol),
    log_fc       = as_mat(f, r$covariate, hex(r$lfc), D),
    se           = as_mat(f, r$covariate, hex(r$se),  D),
    w            = as_mat(f, r$covariate, hex(r$w),   D),
    p_value      = as_mat(f, r$covariate, hex(r$p),   D),
    p_adjusted   = as_mat(f, r$covariate, hex(r$q),   D),
    reject       = as_mat(f, r$covariate, r$signif == '1', D) == 1 )

  b <- bias[bias$case == cs, ]
  out$bias <- cbind(
    delta_em  = hex(b$delta_em),
    delta_wls = hex(b$delta_wls),
    var_delta = hex(b$var_delta) )
  rownames(out$bias) <- colnames(out$log_fc)
  out$n_iters <- setNames(as.integer(b$n_iter), colnames(out$log_fc))

  s <- sfs[sfs$case == cs, ]
  if (nrow(s)) out$sample_fractions <- hex(s$theta)[order(as.integer(s$sample))]

  if (!is.null(out$groups)) {
    fr <- fit[fit$case == cs, ]
    P  <- ncol(out$log_fc)
    out$coef <- matrix(hex(fr$beta), D, P, byrow = TRUE, dimnames = list(NULL, colnames(out$log_fc)))
    out$var  <- matrix(hex(fr$var),  D, P, byrow = TRUE, dimnames = list(NULL, colnames(out$log_fc)))

    vr <- vcv[vcv$case == cs, ]
    k  <- max(as.integer(vr$i)) + 1L
    out$vcov <- array(NA_real_, c(D, k, k))
    out$vcov[cbind(as.integer(vr$feature) + 1L, as.integer(vr$i) + 1L, as.integer(vr$j) + 1L)] <- hex(vr$value)

    dr <- dofs[dofs$case == cs, ]
    dr <- dr[order(as.integer(dr$feature)), ]
    out$dof       <- if (cr$method == 'ancombc2') hex(dr$dof)
    out$estimable <- dr$estimable == '1'
  }

  out
})


#--------------------------------------------------------
# ANCOMBC.
#--------------------------------------------------------
lib <- Sys.getenv('ANCOMBC_LIB')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ANCOMBC))
ancombc_version <- as.character(packageVersion('ANCOMBC'))

ancombc_ref <- list()

for (cs in cases$case) {

  cr <- cases[cases$case == cs, ]
  if (as.integer(cr$max_iter) != 100L || hex(cr$tol) != 1e-5) next
  if (grepl('0 +', cr$formula, fixed = TRUE)) next   # ANCOMBC prepends `~`

  x    <- tables[[cr$table]]$counts
  meta <- tables[[cr$table]]$data
  for (v in names(meta)) if (is.logical(meta[[v]])) meta[[v]] <- factor(meta[[v]], levels = c(FALSE, TRUE))
  meta$sample_id <- rownames(meta)   # ANCOMBC drops a one-column data frame
  grp  <- if (nzchar(cr$grouping)) cr$grouping

  if (cr$method == 'ancombc') {
    if (hex(cr$added) != 1) next
    res <- suppressWarnings(suppressMessages(ANCOMBC::ancombc(
      data = x, taxa_are_rows = FALSE, meta_data = meta, formula = cr$formula,
      p_adj_method = cr$p_adjust, prv_cut = 0, lib_cut = 0, group = grp,
      struc_zero = FALSE, tol = 1e-5, max_iter = 100, conserve = FALSE,
      alpha = hex(cr$alpha), global = FALSE, n_cl = 1, verbose = FALSE )))
    rr  <- res$res
    get <- function (what) {
      m <- as.matrix(rr[[what]][, -1, drop = FALSE])
      rownames(m) <- rr[[what]]$taxon
      m
    }
    out <- list(
      log_fc = get('lfc'), se = get('se'), w = get('W'), p_value = get('p_val'),
      p_adjusted = get('q_val'), reject = get('diff_abn'),
      bias = as.matrix(res$delta_em), sample_fractions = res$samp_frac )
  }
  else {
    set.seed(1)   # ANCOMBC's %dorng% (nothing random is used)
    res <- suppressWarnings(suppressMessages(ANCOMBC::ancombc2(
      data = x, taxa_are_rows = FALSE, meta_data = meta,
      fix_formula = cr$formula, p_adj_method = cr$p_adjust,
      pseudo = hex(cr$pseudocount), pseudo_sens = FALSE, prv_cut = 0,
      lib_cut = 0, s0_perc = hex(cr$var_quantile), group = grp,
      struc_zero = FALSE, alpha = hex(cr$alpha), n_cl = 1, verbose = FALSE,
      global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
      em_control = list(tol = 1e-5, max_iter = 100) )))
    rr  <- res$res
    get <- function (prefix) {
      cols <- grep(paste0('^', prefix, '_'), names(rr), value = TRUE)
      m    <- as.matrix(rr[cols])
      colnames(m) <- sub(paste0('^', prefix, '_'), '', cols)
      rownames(m) <- rr$taxon
      m
    }
    out <- list(
      log_fc = get('lfc'), se = get('se'), w = get('W'), p_value = get('p'),
      p_adjusted = get('q'), reject = get('diff'),
      sample_fractions = res$samp_frac )
  }

  ancombc_ref[[cs]] <- out
}


#--------------------------------------------------------
# Save. Only the ANCOMBC cases whose bias EM converges for every covariate
# within max_iters (100) are compared by the tests: elsewhere its estimates
# depend on the path of nloptr's minimisations rather than SciPy's (see
# plans/ancombc.md for the sizes, and for the cases that differ by design:
# var_quantile = 0, interactions with zeros as missing, and the badly scaled
# `days` covariate).
#--------------------------------------------------------
converged <- c('doc1', 'doc2', 'doc3_pc0.1', 'doc3_pc0.5', 'doc3_pc1.0', 'sim_bc2_pc', 'prop_bc2')
for (cs in names(ancombc_ref)) {
  r <- ancombc_ref[[cs]]
  cat(sprintf(
    '%-15s ANCOMBC vs scikit-bio: max |log_fc diff| %.1e, max rel se diff %.1e, reject differs %d\n',
    cs, max(abs(r$log_fc[, colnames(ref[[cs]]$log_fc)] - ref[[cs]]$log_fc), na.rm = TRUE),
    max(abs(r$se[, colnames(ref[[cs]]$se)] / ref[[cs]]$se - 1), na.rm = TRUE),
    sum(r$reject[, colnames(ref[[cs]]$reject)] != ref[[cs]]$reject, na.rm = TRUE) ))
}
ancombc_ref <- ancombc_ref[converged]

saveRDS(
  list(
    tables  = tables,
    skbio   = ref,
    ancombc = ancombc_ref,
    versions = c(skbio = '0.7.4', ANCOMBC = ancombc_version) ),
  file.path('tests', 'tinytests', 'expected', 'ancombc.rds'),
  version = 2, compress = 'xz' )
