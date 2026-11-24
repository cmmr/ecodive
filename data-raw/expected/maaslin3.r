# Writes tests/tinytests/expected/maaslin3.rds: inputs and reference
# results for test-maaslin3.r, so that the tests need no reference package.
#
# Run from the package root, with maaslin3 1.4.0 (Bioconductor 3.23, MIT;
# installed from its source tarball on R 4.3.3 by relaxing its `Depends: R
# (>= 4.4)`, with collapse 2.1.9 and data.table 1.18 from CRAN), lme4
# 2.1.0, lmerTest 3.2-1 and multcomp 1.4-25 (GPL, run for their outputs
# only) installed:
#
#     Rscript data-raw/expected/maaslin3.r
#
# Every case runs maaslin3() with its filters off (min_variance = -1, so
# that features with constant non-zero values are kept; the tables are
# made so that maaslin3 keeps every feature, which is checked), and no
# plots. maaslin3's median comparison draws 10,000 sets of coefficients
# per covariate with rnorm(), in its covariates' order, after the seed set
# here: the tests regenerate those standard normals with the same
# set.seed() and give them to maaslin3_details(sims = ), so the median
# comparisons can be compared exactly.
#
# The mixed cases fit lme4's lmer() (abundance, through lmerTest) and
# glmer() (prevalence). maaslin3 tries lme4's optimizers in turn and keeps
# the first fit without a warning; for the logistic model nloptwrap often
# stops short of the optimum on flat (singular) fits, so each feature's
# glmer() model is also refitted here with nloptwrap, bobyqa and
# Nelder_Mead, and the fit with the lowest criterion (from glmer(...,
# devFunOnly = TRUE)) is stored, with its theta, coefficients and standard
# errors.

library(lme4)

stopifnot(
  packageVersion('maaslin3') == '1.4.0',
  packageVersion('lme4') == '2.1.0' )

out <- file.path('tests', 'tinytests', 'expected', 'maaslin3.rds')

cols <- c(
  'feature', 'name', 'coef', 'stderr', 'null_hypothesis', 'pval_individual',
  'qval_individual', 'pval_joint', 'qval_joint', 'error', 'N_not_zero' )


#--------------------------------------------------------
# Tables (samples x features) and metadata.
#--------------------------------------------------------
sim <- function (seed, n, D, prev = 0.65, mu = 6) {
  set.seed(seed)
  x <- matrix(rpois(n * D, mu) * rbinom(n * D, 1, prev), n, D)
  # At least two distinct non-zero values per feature, so that maaslin3's
  # variance filter keeps every feature.
  for (f in seq_len(D)) {
    nz <- which(x[, f] > 0)
    if (length(unique(x[nz, f])) < 2) x[sample(n, 2), f] <- c(3, 7)
  }
  dimnames(x) <- list(sprintf('s%02d', seq_len(n)), sprintf('f%02d', seq_len(D)))
  x
}

# No feature present in as many samples of one half as of the other: the
# prevalence coefficient of a two-group comparison would be 0, and the
# sign of its rounding error (which warn_prevalence's flag compares)
# arbitrary.
unbalance <- function (x) {
  h <- seq_len(nrow(x) / 2)
  for (f in seq_len(ncol(x)))
    if (sum(x[h, f] > 0) == sum(x[-h, f] > 0)) x[which(x[h, f] == 0)[1], f] <- 4
  x
}

tables <- list()

x <- sim(1, 40, 25)
x[21:40, 1:3] <- x[21:40, 1:3] * 5
x[21:40, 4] <- 0
x[c(25, 31, 38), 4] <- c(2, 5, 9)                    # rarer in group b
x <- unbalance(x)
tables$two <- list(x = x, meta = data.frame(
  grp = rep(c('a', 'b'), each = 20), row.names = rownames(x) ))

set.seed(2)
x <- sim(3, 45, 30)
x[1:15, 1:4] <- x[1:15, 1:4] * 6
x[, 5] <- c(rep(0, 15), rpois(30, 5) + 1)           # absent at the reference level
x[, 6] <- c(rep(0, 39), c(1, 2, 3, 4, 5, 6))         # rare
x[, 7] <- c(rpois(15, 3) + 1, rep(0, 30))            # only at the reference level
tables$three <- list(x = x, meta = data.frame(
  trt = factor(rep(c('ctl', 'lo', 'hi'), each = 15), levels = c('ctl', 'lo', 'hi')),
  age = round(rnorm(45, 50, 10)),
  sex = rep(c('F', 'M'), length.out = 45),
  row.names = rownames(x) ))

set.seed(4)
x <- sim(5, 48, 20, prev = 0.6, mu = 5)
x[25:48, 1:3] <- x[25:48, 1:3] * 4
x <- unbalance(x)
tables$mixed <- list(x = x, meta = data.frame(
  grp = rep(c('a', 'b'), each = 24), age = rnorm(48),
  id  = rep(sprintf('p%02d', 1:12), 4), row.names = rownames(x) ))

set.seed(6)
x <- sim(7, 45, 20, prev = 0.7)
x[1:15, 1:3] <- x[1:15, 1:3] * 5
tables$small <- list(x = x, meta = data.frame(
  trt = factor(rep(c('ctl', 'lo', 'hi'), each = 15), levels = c('ctl', 'lo', 'hi')),
  age = rnorm(45), id = rep(sprintf('p%02d', 1:15), 3), row.names = rownames(x) ))

# A table whose every feature has a positive CLR value in some sample
# (maaslin3's min_abundance = 0 filter drops the others with CLR).
x <- sim(8, 40, 20, prev = 0.8, mu = 8)
x[21:40, 1:3] <- x[21:40, 1:3] * 4
x <- unbalance(x)
tables$clr <- list(x = x, meta = data.frame(
  grp = rep(c('a', 'b'), each = 20), age = rnorm(40), row.names = rownames(x) ))


#--------------------------------------------------------
# Cases: a table, a formula, maaslin3's arguments and
# the matching ecodive arguments.
#--------------------------------------------------------
cases <- list(
  two = list(table = 'two', formula = '~ grp',
    ref = list(correction = 'BH'), eco = list(p_adjust = 'BH', alpha = 0.1)),
  three = list(table = 'three', formula = '~ trt + age + sex',
    ref = list(correction = 'holm', max_significance = 0.05), eco = list(p_adjust = 'holm', alpha = 0.05)),
  clr = list(table = 'clr', formula = '~ grp + age',
    ref = list(normalization = 'CLR', transform = 'NONE', warn_prevalence = FALSE),
    eco = list(norm = 'rclr', transform = 'none', warn_prevalence = FALSE, p_adjust = 'BH', alpha = 0.1)),
  none_log = list(table = 'three', formula = '~ trt + age',
    ref = list(normalization = 'NONE', transform = 'LOG'),
    eco = list(norm = 'none', transform = 'log', p_adjust = 'BH', alpha = 0.1)),
  tss_none = list(table = 'two', formula = '~ grp',
    ref = list(transform = 'NONE', warn_prevalence = FALSE),
    eco = list(transform = 'none', warn_prevalence = FALSE, p_adjust = 'BH', alpha = 0.1)),
  median_both = list(table = 'three', formula = '~ trt + age',
    ref = list(median_comparison_prevalence = TRUE, median_comparison_abundance_threshold = 0.3,
               median_comparison_prevalence_threshold = 0.2, subtract_median = TRUE),
    eco = list(median_comparison_prevalence = TRUE, median_comparison_abundance_threshold = 0.3,
               median_comparison_prevalence_threshold = 0.2, subtract_median = TRUE, p_adjust = 'BH', alpha = 0.1)),
  no_augment = list(table = 'three', formula = '~ trt + age',
    ref = list(augment = FALSE, warn_prevalence = FALSE, median_comparison_abundance = FALSE, standardize = FALSE),
    eco = list(augment = FALSE, warn_prevalence = FALSE, median_comparison_abundance = FALSE,
               standardize = FALSE, p_adjust = 'BH', alpha = 0.1)),
  mixed = list(table = 'mixed', formula = '~ grp + age + (1 | id)',
    ref = list(), eco = list(p_adjust = 'BH', alpha = 0.1)),
  small_re = list(table = 'small', formula = '~ trt + age + (1 | id)',
    ref = list(small_random_effects = TRUE), eco = list(small_random_effects = TRUE, p_adjust = 'BH', alpha = 0.1)),
  small_groups = list(table = 'small', formula = '~ trt + age + (1 | id)',
    ref = list(), eco = list(p_adjust = 'BH', alpha = 0.1)) )


#--------------------------------------------------------
# glmer() refitted with three optimizers: the best fit
# of each feature's augmented logistic model.
#--------------------------------------------------------
best_glmer <- function (x, meta, formula) {
  f   <- stats::as.formula(paste('expr', formula))
  fx  <- reformulas::nobars(f)
  n   <- nrow(x)
  # maaslin3 standardizes numeric metadata and makes others factors.
  for (v in names(meta)) {
    if (is.numeric(meta[[v]])) meta[[v]] <- as.vector(scale(meta[[v]]))
  }
  p <- ncol(stats::model.matrix(fx[-2L], meta)) - 1
  a <- p / (2 * n)
  w <- c(rep(1, n), rep(a, 2 * n))
  res <- lapply(colnames(x), function (ft) {
    d   <- data.frame(expr = as.numeric(x[, ft] > 0), meta)
    if (length(unique(d$expr)) < 2) return (NULL)
    aug <- rbind(d, transform(d, expr = 1), transform(d, expr = 0))
    dfun <- suppressWarnings(suppressMessages(
      glmer(f, data = aug, family = binomial, weights = w, devFunOnly = TRUE) ))
    fits <- lapply(c('nloptwrap', 'bobyqa', 'Nelder_Mead'), function (o) tryCatch(
      suppressWarnings(suppressMessages(glmer(
        f, data = aug, family = binomial, weights = w,
        control = glmerControl(optimizer = o) ))),
      error = function (e) NULL ))
    fits <- Filter(Negate(is.null), fits)
    dev  <- vapply(fits, function (m) dfun(c(getME(m, 'theta'), fixef(m))), 0)
    m    <- fits[[which.min(dev)]]
    co   <- stats::coef(summary(m))
    data.frame(
      feature = ft, name = rownames(co)[-1], coef = co[-1, 1], stderr = co[-1, 2],
      theta = getME(m, 'theta'), deviance = min(dev), row.names = NULL )
  })
  do.call(rbind, res)
}


#--------------------------------------------------------
# Run.
#--------------------------------------------------------
seed <- 20261124L
for (nm in names(cases)) {
  cs <- cases[[nm]]
  tb <- tables[[cs$table]]
  set.seed(seed)
  fit <- suppressWarnings(suppressMessages(do.call(maaslin3::maaslin3, c(list(
    input_data = tb$x, input_metadata = tb$meta, output = tempfile(),
    formula = stats::as.formula(cs$formula), min_variance = -1,
    plot_summary_plot = FALSE, plot_associations = FALSE, verbosity = 'ERROR'), cs$ref))))
  ra <- fit$fit_data_abundance$results
  rp <- fit$fit_data_prevalence$results
  stopifnot(setequal(ra$feature, colnames(tb$x)), setequal(rp$feature, colnames(tb$x)))
  cases[[nm]]$seed       <- seed
  cases[[nm]]$abundance  <- as.data.frame(ra[, cols])
  cases[[nm]]$prevalence <- as.data.frame(rp[, cols])
  if (grepl('|', cs$formula, fixed = TRUE) && !isTRUE(cs$ref$small_random_effects))
    cases[[nm]]$glmer <- best_glmer(tb$x, tb$meta, cs$formula)
}

versions <- c(
  maaslin3 = as.character(packageVersion('maaslin3')),
  lme4     = as.character(packageVersion('lme4')),
  lmerTest = as.character(packageVersion('lmerTest')),
  multcomp = as.character(packageVersion('multcomp')),
  R        = paste(R.version$major, R.version$minor, sep = '.') )

saveRDS(list(tables = tables, cases = cases, versions = versions), out, version = 2)
