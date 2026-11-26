#test_that("mantel_correlogram", {

  ex <- readRDS('expected/mantel_correlogram.rds')
  x  <- ex$x
  xy <- ex$xy
  n  <- attr(x, 'Size')



  # Same classes, statistics and p-values as vegan::mantel.correlog(),
  # given the same permutations: stored results, so at home only. ====

  if (at_home()) {
    for (cs in ex$cases) {

      a    <- cs$args
      info <- paste(deparse(a), collapse = '')

      res <- mantel_correlogram(
        x, xy,
        breaks       = a$break.pts,
        n_classes    = a$n.class,
        cutoff       = if (is.null(a$cutoff)) TRUE else a$cutoff,
        method       = if (is.null(a$r.type)) 'pearson' else a$r.type,
        p_adjust     = if (is.null(a$mult)) 'holm' else a$mult,
        progressive  = if (is.null(a$progressive)) TRUE else a$progressive,
        permutations = ex$perms,
        cpus         = 1 )

      # vegan lowers its first break by .Machine$double.eps.
      expect_equal(attr(res, 'breaks')[-1], cs$breaks[-1], tolerance = 1e-14, info = info)
      expect_equal(attr(res, 'breaks')[[1]], cs$breaks[[1]], tolerance = 1e-12, info = info)
      expect_equal(res$class_mid, cs$class_mid, tolerance = 1e-12, info = info)
      expect_identical(2L * res$n_pairs, as.integer(cs$n_dist), info = info)
      expect_identical(is.na(res$r), is.na(cs$r), info = info)
      expect_identical(sum(!is.na(res$r)), cs$n_tests, info = info)
      expect_equal(res$r, cs$r, tolerance = 1e-12, info = info)
      expect_identical(res$p_value, cs$p_value, info = info)
      # Sidak's 1 - (1 - p)^k goes through a power function, whose last bits
      # differ between R builds and C libraries (R 4.6 and musl differ from
      # the R 4.3 / glibc that made the fixture); p.adjust()'s methods only
      # multiply and take minima, so they match exactly everywhere.
      if (identical(a$mult, 'sidak'))
        expect_equal(res$p_adjusted, cs$p_adjusted, tolerance = 1e-14, info = info)
      else
        expect_identical(res$p_adjusted, cs$p_adjusted, info = info)
    }
  }



  # Result format. ====

  res <- mantel_correlogram(x, xy, permutations = 99, cpus = 1)

  expect_inherits(res, c('ecodive_mantel_correlogram', 'data.frame'))
  expect_identical(class(res), c('ecodive_mantel_correlogram', 'data.frame'))
  expect_identical(names(res), c(
    'class_from', 'class_to', 'class_mid', 'n_pairs', 'r', 'p_value', 'p_adjusted' ))
  expect_identical(nrow(res), 13L)   # Sturges: ceiling(log2(2415) + 1)
  expect_identical(attr(res, 'breaks'), seq(min(dist(xy)), max(dist(xy)), length.out = 14))
  expect_identical(res$class_from, attr(res, 'breaks')[-14])
  expect_identical(res$class_to,   attr(res, 'breaks')[-1])
  expect_identical(sum(res$n_pairs), (n * (n - 1L)) %/% 2L)
  expect_identical(attr(res, 'method'), 'pearson')
  expect_identical(attr(res, 'alternative'), 'sign')
  expect_identical(attr(res, 'p_adjust'), 'holm')
  expect_true(attr(res, 'progressive'))
  expect_identical(attr(res, 'n_perm'), 99L)
  expect_identical(attr(res, 'permutation'), 'free')
  expect_identical(dim(attr(res, 'r_perm')), c(13L, 99L))
  expect_identical(is.na(attr(res, 'r_perm')[, 1]), is.na(res$r))
  expect_true(all(res$p_value >= 0.01 & res$p_value <= 1, na.rm = TRUE))
  expect_true(all(res$p_adjusted >= res$p_value, na.rm = TRUE))
  expect_identical(is.na(res$p_value), is.na(res$r))
  expect_identical(attr(res, 'display')$title, 'Mantel correlogram')

  out <- capture.output(print(res))
  expect_true(any(grepl('Mantel correlogram', out, fixed = TRUE)))
  expect_true(any(grepl('A data.frame: 13 ', out, fixed = TRUE)))



  # Titles: of both inputs when both have them, else of dm. ====

  set.seed(1)
  counts <- matrix(rpois(20 * 8, 5), 20, 8, dimnames = list(paste0('S', 1:20), NULL))
  geo    <- matrix(runif(40), 20, 2, dimnames = list(paste0('S', 1:20), c('lat', 'lon')))
  meta   <- as.data.frame(geo)
  dm     <- bray(counts)

  expect_identical(
    attr(mantel_correlogram(dm, euclidean(geo), permutations = 9), 'display')$title,
    'Mantel correlogram of Bray-Curtis Dissimilarity and Euclidean Distance' )
  expect_identical(
    attr(mantel_correlogram(dm, geo, permutations = 9), 'display')$title,
    'Mantel correlogram of Bray-Curtis Dissimilarity' )

  titles <- function (x, ...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control('enable')
    plot(x, ...)
    txt <- unlist(lapply(grDevices::recordPlot()[[1]], function (op) {
      a <- unlist(op[[2]], recursive = TRUE)
      a[vapply(a, is.character, logical(1))]
    }))
    gsub('\n', ' ', unique(unlist(txt)))
  }
  r2 <- mantel_correlogram(dm, geo, permutations = 9)
  expect_true('Mantel correlogram of Bray-Curtis Dissimilarity' %in% titles(r2))
  expect_true('Mantel r (Pearson)' %in% titles(r2))
  expect_true('x' %in% titles(r2, main = 'x'))
  r3 <- r2; attr(r3, 'display') <- NULL
  expect_true('Mantel correlogram' %in% titles(r3))



  # Each class is a Mantel test of dm against the class's 0/1 model matrix
  # (0 in the class), with the classes permuted: mantel() of the model
  # matrix against dm, given the same permutations. ====

  check_classes <- function (dm, y, method, perms, info, ...) {
    res <- mantel_correlogram(dm, y, method = method, permutations = perms, cpus = 1, ...)
    yd  <- as.matrix(if (inherits(y, 'dist')) y else dist(y))
    br  <- attr(res, 'breaks')
    rp  <- attr(res, 'r_perm')
    for (k in which(!is.na(res$r))) {
      inside <- (yd > br[[k]] & yd <= br[[k + 1]]) | (k == 1 & yd == br[[1]])
      m      <- suppressWarnings(mantel(
        as.dist(1 - inside), dm, method = method, permutations = perms,
        alternative = if (res$r[[k]] > 0) 'greater' else 'less', cpus = 1 ))
      expect_equal(res$r[[k]], m$statistic[['r']], tolerance = 1e-12, info = paste(info, k))
      expect_equal(rp[k, ], attr(m, 'r_perm'), tolerance = 1e-12, info = paste(info, k))
      expect_identical(res$p_value[[k]], m$p_value, info = paste(info, k))
    }
  }

  set.seed(2)
  perms <- t(replicate(49, sample(n)))
  for (method in c('pearson', 'spearman', 'kendall'))
    check_classes(x, xy, method, perms, method, n_classes = 6, cutoff = FALSE)

  # Tied distances in both.
  set.seed(3)
  ti <- matrix(sample(0:3, 30 * 2, replace = TRUE), 30, 2)
  td <- dist(matrix(sample(0:2, 30 * 3, replace = TRUE), 30, 3))
  tp <- t(replicate(39, sample(30)))
  for (method in c('pearson', 'spearman', 'kendall'))
    check_classes(td, ti, method, tp, paste('ties', method), breaks = c(0, 1, 2, 3, 5))



  # Arguments. ====

  # alternative: the p-values of the permuted statistics.
  for (alt in c('greater', 'less', 'two.sided')) {
    ra <- mantel_correlogram(x, xy, alternative = alt, permutations = 99, cpus = 1)
    expect_identical(ra$r, res$r, info = alt)
    k  <- which(!is.na(ra$r))
    expect_identical(
      ra$p_value[k],
      vapply(k, function (i) ecodive:::perm_p(ra$r[[i]], attr(ra, 'r_perm')[i, ], alt), 0),
      info = alt )
  }
  k <- which(!is.na(res$r))
  expect_identical(
    res$p_value[k],
    vapply(k, function (i) ecodive:::perm_p(
      res$r[[i]], attr(res, 'r_perm')[i, ], if (res$r[[i]] > 0) 'greater' else 'less'), 0) )

  # p_adjust and progressive.
  ra <- mantel_correlogram(x, xy, p_adjust = 'none', permutations = 99, cpus = 1)
  expect_identical(ra$p_adjusted, ra$p_value)
  ra <- mantel_correlogram(x, xy, p_adjust = 'BH', progressive = FALSE, permutations = 99, cpus = 1)
  expect_identical(ra$p_adjusted[k], p.adjust(ra$p_value[k], 'BH'))
  ra <- mantel_correlogram(x, xy, p_adjust = 'bonf', permutations = 99, cpus = 1)
  expect_identical(attr(ra, 'p_adjust'), 'bonferroni')
  expect_identical(ra$p_adjusted[k], vapply(seq_along(k), function (j)
    if (j == 1) ra$p_value[k[[1]]] else min(1, j * ra$p_value[k[[j]]]), 0))
  ra <- mantel_correlogram(x, xy, p_adjust = 'sidak', progressive = FALSE, permutations = 99, cpus = 1)
  expect_identical(ra$p_adjusted[k], 1 - (1 - ra$p_value[k]) ^ length(k))

  # cutoff = FALSE tests every class with pairs.
  ra <- mantel_correlogram(x, xy, cutoff = FALSE, permutations = 9, cpus = 1)
  expect_true(all(!is.na(ra$r)))
  expect_identical(ra$r[k], res$r[k])

  # n_classes and breaks.
  ra <- mantel_correlogram(x, xy, n_classes = 4, permutations = 9, cpus = 1)
  expect_identical(nrow(ra), 4L)
  rb <- mantel_correlogram(x, xy, breaks = attr(ra, 'breaks'), permutations = 9, cpus = 1)
  expect_identical(rb$r, ra$r)
  expect_identical(rb$n_pairs, ra$n_pairs)

  # Distances outside the breaks are in no class; a class with no pairs is
  # not tested.
  d  <- as.vector(dist(xy))
  u  <- sort(unique(d[d > 2]))[1:2]
  rb <- mantel_correlogram(x, xy, breaks = c(1, u[[1]], mean(u), 4, 6), cutoff = FALSE, permutations = 9, cpus = 1)
  expect_identical(rb$n_pairs, c(
    sum(d >= 1 & d <= u[[1]]), 0L, sum(d > mean(u) & d <= 4), sum(d > 4 & d <= 6) ))
  expect_identical(rb$n_pairs[[2]], 0L)
  expect_true(is.na(rb$r[[2]]) && is.na(rb$p_value[[2]]) && is.na(rb$p_adjusted[[2]]))
  expect_true(all(!is.na(rb$r[-2])))
  expect_identical(rb$p_adjusted[[3]], p.adjust(rb$p_value[c(1, 3)], 'holm')[[2]])

  # The first class holds its lower boundary, even for large distances.
  rk <- mantel_correlogram(x, xy * 1000, permutations = 9, cpus = 1)
  expect_identical(sum(rk$n_pairs), (n * (n - 1L)) %/% 2L)
  expect_identical(rk$n_pairs, res$n_pairs)
  expect_equal(rk$r, res$r, tolerance = 1e-12)

  # Permutations: none, within strata.
  r0 <- mantel_correlogram(x, xy, permutations = 0)
  expect_identical(r0$r, res$r)
  expect_true(all(is.na(r0$p_value)) && all(is.na(r0$p_adjusted)))
  expect_identical(dim(attr(r0, 'r_perm')), c(13L, 0L))
  rs <- mantel_correlogram(x, xy, strata = rep(1:2, n / 2), permutations = 19, method = 'kendall')
  expect_identical(attr(rs, 'permutation'), 'within strata')
  expect_identical(attr(mantel_correlogram(x, xy, permutations = ex$perms[1:9, ]), 'permutation'), 'supplied matrix')



  # The forms of y, and matching samples by name. ====

  ref <- mantel_correlogram(dm, dist(geo), permutations = 19, cpus = 1)
  same <- function (a, b) expect_identical(no_display(a), no_display(b))
  same(mantel_correlogram(dm, geo, permutations = 19, cpus = 1), ref)
  same(mantel_correlogram(dm, meta, permutations = 19, cpus = 1), ref)
  same(mantel_correlogram(dm, as.matrix(dist(geo)), permutations = 19, cpus = 1), ref)
  same(mantel_correlogram(dm, c('lat', 'lon'), meta, permutations = 19, cpus = 1), ref)
  same(mantel_correlogram(dm, geo[20:1, ], permutations = 19, cpus = 1), ref)
  same(mantel_correlogram(dm, c('lat', 'lon'), meta[20:1, ], permutations = 19, cpus = 1), ref)

  ref1 <- mantel_correlogram(dm, dist(geo[, 'lat']), permutations = 19, cpus = 1)
  same(mantel_correlogram(dm, geo[, 'lat'], permutations = 19, cpus = 1), ref1)
  same(mantel_correlogram(dm, 'lat', meta, permutations = 19, cpus = 1), ref1)
  same(mantel_correlogram(dm, meta$lat, meta, permutations = 19, cpus = 1), ref1)

  # Unnamed coordinates are in dm's order.
  same(mantel_correlogram(dm, unname(geo), permutations = 19, cpus = 1), ref)



  # Reproducible, and identical for any number of threads. ====

  a <- mantel_correlogram(x, xy, permutations = 199, seed = 3, cpus = 1)
  expect_identical(a, mantel_correlogram(x, xy, permutations = 199, seed = 3, cpus = 1))
  expect_false(identical(attr(a, 'r_perm'), attr(mantel_correlogram(x, xy, permutations = 199, seed = 4, cpus = 1), 'r_perm')))
  for (cpus in if (at_home()) c(2, 3, 4, 8) else 2)
    for (method in c('pearson', 'kendall'))
      expect_identical(
        mantel_correlogram(x, xy, method = method, permutations = 199, seed = 3, cpus = cpus),
        mantel_correlogram(x, xy, method = method, permutations = 199, seed = 3, cpus = 1),
        info = paste(method, cpus) )



  # Edge cases. ====

  # Three samples.
  d3 <- dist(c(a = 0, b = 1, c = 3))
  r3 <- mantel_correlogram(d3, c(a = 0, b = 1, c = 5), n_classes = 2, cutoff = FALSE, permutations = 5)
  expect_identical(r3$n_pairs, c(1L, 2L))
  expect_true(all(is.finite(r3$r)))

  # Every pair in one class: not tested.
  r1 <- mantel_correlogram(x, xy, n_classes = 1, permutations = 9)
  expect_identical(r1$n_pairs, as.integer(n * (n - 1) / 2))
  expect_true(is.na(r1$r))

  # Constant dm.
  expect_warning(
    rc <- mantel_correlogram(dist(rep(1, 10)) + 1, matrix(runif(20), 10), permutations = 9),
    'same distance between every pair' )
  expect_true(all(is.na(rc$r)))



  # Plots. ====

  drawn <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  grDevices::pdf(NULL)
  expect_true(drawn(res))
  expect_true(drawn(res, alpha = 0.01, main = 'x', col = 'blue'))
  expect_true(drawn(r0))
  expect_error(plot(res, alpha = 2), '`alpha`')
  expect_error(plot(rc), 'no tested classes')
  grDevices::dev.off()



  # Invalid arguments. ====

  expect_error(mantel_correlogram(x), '`y` is required')
  expect_error(mantel_correlogram(x, xy, breaks = c(2, 1)), '`breaks` must be')
  expect_error(mantel_correlogram(x, xy, breaks = 1), '`breaks` must be')
  expect_error(mantel_correlogram(x, xy, breaks = c(0, 1), n_classes = 3), 'not both')
  expect_error(mantel_correlogram(x, xy, n_classes = 0), '`n_classes` must be')
  expect_error(mantel_correlogram(x, xy, n_classes = 2.5), '`n_classes` must be')
  expect_error(mantel_correlogram(x, xy, cutoff = NA), '`cutoff` must be TRUE or FALSE')
  expect_error(mantel_correlogram(x, xy, progressive = 'yes'), '`progressive` must be TRUE or FALSE')
  expect_error(mantel_correlogram(x, xy, p_adjust = 'nope'), '`p_adjust` must be one of')
  expect_error(mantel_correlogram(x, xy, alternative = 'nope'), "'arg' should be one of")
  expect_error(mantel_correlogram(x, xy, method = 'nope'), "'arg' should be one of")
  expect_error(mantel_correlogram(x, dist(rep(1, n))), 'no distance classes')
  expect_error(mantel_correlogram(x, xy[-1, ]), 'samples')
  expect_error(mantel_correlogram(dm, c('lat', 'nope'), meta), "not in `data`: 'nope'")
  expect_error(mantel_correlogram(dm, c('lat', 'lon')), '`y`')
  expect_error(mantel_correlogram(dm, 'lat'), 'names a column, but `data` is not given')
  bad <- meta; bad$lat[[2]] <- NA
  expect_error(mantel_correlogram(dm, c('lat', 'lon'), bad), 'must not contain NA')
  bad <- meta; bad$lat <- letters[1:20]
  expect_error(mantel_correlogram(dm, c('lat', 'lon'), bad), 'must be numeric')
  expect_error(mantel_correlogram(dm, geo, strata = rep(1:2, 10), permutations = t(replicate(3, 1:20))), 'not both')
  expect_error(mantel_correlogram(dist(1:2), 1:2), 'at least 3 samples')

#})
