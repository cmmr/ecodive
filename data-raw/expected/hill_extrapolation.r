# Writes tests/tinytests/expected/hill_extrapolation.rds, and the same
# values as text under data-raw/expected/hill_extrapolation/: reference
# values for test-hill_extrapolation.r, so that the tests need no reference
# package.
#
# Run from the package root, with ecodive and iNEXT installed:
#
#     Rscript data-raw/expected/hill_extrapolation.r
#
# iNEXT is GPL: only its outputs are used, never its source.
# Generated with R 4.3.3 and iNEXT 3.0.2 (versions stored in `versions`).
#
# Tables (`tables`, each samples x features, whole numbers):
#   spider   iNEXT's spider data (Chao et al. 2014's worked example):
#            Girdled (n = 168) and Logged (n = 252).
#   ex       ecodive's ex_counts (4 samples).
#   edge     small samples: f2 = 0; f1 = 0; all singletons; one feature;
#            one singleton and no doubletons; c(1, 1, 2); c(2, 2, 3).
#   deep     one simulated sample of 20,000 observations (lognormal
#            abundances over 600 features; R's RNG, seed 64).
#
# Fields, for each table name (`curves`, `asy`, `info`), from
# iNEXT::iNEXT(rows, q = c(0, 1, 2), datatype = 'abundance', nboot = 0):
#   curves   size_based (sample, q, depth, method, qD, SC).
#   cov      coverage_based (sample, q, SC, depth, method, qD).
#   asy      AsyEst (sample, q, observed, estimate).
#   info     DataInfo (sample, n, s_obs, SC, f1 .. f10).
#   grids    The depths iNEXT chose for spider's Girdled (n = 168) with
#            other `endpoint` and `knots`: a list of (endpoint, knots, m).
#   estd     iNEXT::estimateD() on spider and ex: base 'size' and
#            'coverage', default and given levels (sample, q, base, level,
#            depth, method, SC, qD).
#   boot     The bootstrap, with iNEXT's own draws: for spider's Girdled,
#            set.seed(7), then iNEXT(..., size = boot$sizes, nboot = 200),
#            whose draws are rmultinom(200, n, p) with p the Appendix G
#            assemblage, most abundant first (settled from the outputs):
#            the first 200 (draws1) give the size- and coverage-based
#            standard errors, the next 100 (draws2) the asymptotes' (iNEXT
#            always uses 100 replicates for these, whatever nboot). Then
#            set.seed(5) and estimateD(..., nboot = 60) at size 100 and at
#            coverage 0.9, each with draws rmultinom(60, n, p) (draws_e).
#            Standard errors are (UCL - estimate) / qnorm(0.975), where the
#            limits aren't clipped.
#   se2000   iNEXT's standard errors with nboot = 2000 (set.seed(11)), for
#            spider's Girdled and ex's Stool: size_based at each sample's
#            half size, size and double size; and for the asymptotes, whose
#            bootstrap iNEXT always runs with 100 replicates, the root mean
#            square of 20 runs' standard errors (set.seed(11 + k)).

library(ecodive)
library(iNEXT)
stopifnot(packageVersion('iNEXT') == '3.0.2')

dir <- file.path('data-raw', 'expected', 'hill_extrapolation')
dir.create(dir, showWarnings = FALSE)


# Tables -------------------------------------------------------------------

rows_to_matrix <- function (lst) {
  k <- max(lengths(lst))
  m <- t(vapply(lst, function (v) c(v, rep(0, k - length(v))), numeric(k)))
  rownames(m) <- names(lst)
  m
}

data(spider, package = 'iNEXT')
set.seed(64)
deep <- as.vector(stats::rmultinom(1, 20000, local({ p <- stats::rlnorm(600, 0, 2); p / sum(p) })))

tables <- list(
  spider = rows_to_matrix(spider),
  ex     = ex_counts,
  edge   = rows_to_matrix(list(
    f2_zero    = c(1, 1, 1, 3, 5),
    f1_zero    = c(2, 2, 3, 5),
    singletons = c(1, 1, 1, 1, 1),
    one        = c(10),
    single_one = c(1, 9),
    small_a    = c(1, 1, 2),
    small_b    = c(2, 2, 3) )),
  deep   = matrix(deep, 1, dimnames = list('deep', NULL)) )

as_list <- function (m) {
  lst <- lapply(seq_len(nrow(m)), function (i) { v <- m[i, ]; v[v > 0] })
  names(lst) <- rownames(m)
  lst
}


# Curves, asymptotes and data info ------------------------------------------

curves <- cov <- asy <- info <- list()
for (nm in names(tables)) {
  o  <- iNEXT(as_list(tables[[nm]]), q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
  sb <- o$iNextEst$size_based
  cb <- o$iNextEst$coverage_based
  curves[[nm]] <- data.frame(
    sample = as.character(sb$Assemblage), q = sb$Order.q, depth = sb$m,
    method = tolower(as.character(sb$Method)), qD = sb$qD, SC = sb$SC,
    stringsAsFactors = FALSE )
  cov[[nm]] <- data.frame(
    sample = as.character(cb$Assemblage), q = cb$Order.q, SC = cb$SC, depth = cb$m,
    method = tolower(as.character(cb$Method)), qD = cb$qD, stringsAsFactors = FALSE )
  utils::write.csv(cov[[nm]], file.path(dir, paste0('cov_', nm, '.csv')), row.names = FALSE)
  a <- o$AsyEst
  asy[[nm]] <- data.frame(
    sample = as.character(a$Assemblage),
    q = match(a$Diversity, c('Species richness', 'Shannon diversity', 'Simpson diversity')) - 1,
    observed = a$Observed, estimate = a$Estimator, stringsAsFactors = FALSE )
  di <- DataInfo(as_list(tables[[nm]]), datatype = 'abundance')
  info[[nm]] <- data.frame(
    sample = as.character(di$Assemblage), n = di$n, s_obs = di$S.obs, SC = di$SC,
    di[, paste0('f', 1:10)], stringsAsFactors = FALSE )
  utils::write.csv(curves[[nm]], file.path(dir, paste0('curves_', nm, '.csv')), row.names = FALSE)
  utils::write.csv(asy[[nm]],    file.path(dir, paste0('asy_', nm, '.csv')),    row.names = FALSE)
  utils::write.csv(info[[nm]],   file.path(dir, paste0('info_', nm, '.csv')),   row.names = FALSE)
}


# Depth grids ----------------------------------------------------------------

girdled <- list(Girdled = spider$Girdled)
grid_cases <- list(
  c(336, 40), c(336, 41), c(336, 77), c(500, 21), c(309, 39), c(200, 10),
  c(100, 15), c(168, 12), c(169, 4), c(1000, 5) )
grids <- lapply(grid_cases, function (g) {
  m <- iNEXT(girdled, q = 0, datatype = 'abundance', endpoint = g[[1]], knots = g[[2]], nboot = 0)$iNextEst$size_based$m
  list(endpoint = g[[1]], knots = g[[2]], m = m)
})
writeLines(
  vapply(grids, function (g) paste(g$endpoint, g$knots, paste(g$m, collapse = ' ')), ''),
  file.path(dir, 'grids.txt') )


# estimateD -------------------------------------------------------------------

estd <- list()
for (nm in c('spider', 'ex')) {
  lst <- as_list(tables[[nm]])
  for (base in c('size', 'coverage')) {
    levels <- if (base == 'size') list(NULL, c(50, 300)) else list(NULL, c(0.8, 0.95, 0.99, 1))
    for (lv in levels) {
      e <- estimateD(lst, q = c(0, 1, 2), datatype = 'abundance', base = base, level = lv, nboot = 0)
      estd[[length(estd) + 1L]] <- data.frame(
        table = nm, sample = as.character(e$Assemblage), q = e$Order.q, base = base,
        given = !is.null(lv), depth = e$m, method = tolower(as.character(e$Method)),
        SC = e$SC, qD = e$qD, stringsAsFactors = FALSE )
    }
  }
}
estd <- do.call(rbind, estd)
utils::write.csv(estd, file.path(dir, 'estimateD.csv'), row.names = FALSE)


# The bootstrap, with iNEXT's own draws -----------------------------------------

x  <- spider$Girdled
n  <- sum(x)
f1 <- sum(x == 1); f2 <- sum(x == 2)
C  <- 1 - f1 / n * (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
f0 <- ceiling((n - 1) / n * f1^2 / (2 * f2))
lam <- (1 - C) / sum(x / n * (1 - x / n)^n)
p   <- sort(c(x / n * (1 - lam * (1 - x / n)^n), rep((1 - C) / f0, f0)), decreasing = TRUE)

z     <- stats::qnorm(0.975)
sizes <- c(1, 10, 50, 167, 168, 169, 300)
set.seed(7)
o <- iNEXT(girdled, q = c(0, 1, 2), datatype = 'abundance', size = sizes, nboot = 200)
set.seed(7)
draws1 <- stats::rmultinom(200, n, p)
draws2 <- stats::rmultinom(100, n, p)
sb <- o$iNextEst$size_based
cb <- o$iNextEst$coverage_based
boot <- list(
  sizes  = sort(unique(sb$m)),
  p      = p,
  draws1 = draws1,
  draws2 = draws2,
  size   = data.frame(
    q = sb$Order.q, depth = sb$m, se = (sb$qD.UCL - sb$qD) / z,
    sc_se = (sb$SC.UCL - sb$SC) / z, sc_clipped = sb$SC.UCL >= 1 | sb$SC.LCL <= 0 ),
  coverage = data.frame(q = cb$Order.q, depth = cb$m, SC = cb$SC, se = (cb$qD.UCL - cb$qD) / z),
  asy_se = o$AsyEst$s.e. )

set.seed(5)
e_size <- estimateD(girdled, q = c(0, 1, 2), datatype = 'abundance', base = 'size', level = 100, nboot = 60)
set.seed(5)
e_cov <- estimateD(girdled, q = c(0, 1, 2), datatype = 'abundance', base = 'coverage', level = 0.9, nboot = 60)
set.seed(5)
boot$draws_e <- stats::rmultinom(60, n, p)
boot$estd_size_se <- (e_size$qD.UCL - e_size$qD) / z
boot$estd_cov_se  <- (e_cov$qD.UCL - e_cov$qD) / z
storage.mode(boot$draws1) <- storage.mode(boot$draws2) <- storage.mode(boot$draws_e) <- 'integer'
utils::write.csv(boot$size,     file.path(dir, 'boot_size.csv'),     row.names = FALSE)
utils::write.csv(boot$coverage, file.path(dir, 'boot_coverage.csv'), row.names = FALSE)


# Standard errors from 2,000 replicates ---------------------------------------

se2000 <- list()
for (s in list(list('spider', 'Girdled'), list('ex', 'Stool'))) {
  v  <- tables[[s[[1]]]][s[[2]], ]
  v  <- v[v > 0]
  nn <- sum(v)
  set.seed(11)
  o  <- iNEXT(stats::setNames(list(v), s[[2]]), q = c(0, 1, 2), datatype = 'abundance',
              size = c(floor(nn / 2), nn, 2 * nn), nboot = 2000)
  sb <- o$iNextEst$size_based
  sb <- sb[sb$m %in% c(floor(nn / 2), nn, 2 * nn), ]
  asy_se <- sapply(1:20, function (k) {
    set.seed(11 + k)
    iNEXT(stats::setNames(list(v), s[[2]]), q = c(0, 1, 2), datatype = 'abundance', size = nn, nboot = 2)$AsyEst$s.e.
  })
  se2000[[s[[2]]]] <- list(
    asy  = data.frame(q = 0:2, se = sqrt(rowMeans(asy_se^2))),
    size = data.frame(q = sb$Order.q, depth = sb$m, se = (sb$qD.UCL - sb$qD) / z,
                      clipped = sb$qD.LCL <= 0) )
}


versions <- c(R = as.character(getRversion()), iNEXT = as.character(packageVersion('iNEXT')))
writeLines(paste(names(versions), versions), file.path(dir, 'versions.txt'))

saveRDS(
  list(
    tables = tables, curves = curves, cov = cov, asy = asy, info = info, grids = grids,
    estd = estd, boot = boot, se2000 = se2000, versions = versions ),
  file.path('tests', 'tinytests', 'expected', 'hill_extrapolation.rds'),
  version = 2 )
