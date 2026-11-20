#test_that("enterotype", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  props <- function (m) m / rowSums(m)

  # Two partitions are the same, whatever the clusters' numbers.
  same_partition <- function (a, b) {
    tab <- table(a, b)
    all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
  }

  # PAM through enterotype()'s internals, on a dist.
  pam_of <- function (d, k, cpus = 1L) {
    dv <- dist_vector(d)
    pam_fit(dv$dm, dv$n, k, cpus)
  }



  # Stored results from cluster::pam(), cluster::silhouette() and
  # clusterSim::index.G1() (data-raw/expected/enterotype.r). ====

  if (at_home()) {

    ref <- readRDS(file.path('expected', 'enterotype.rds'))

    inputs <- list(
      dirichlet = list(d = jensen(ref$dirichlet_counts), x = ref$dirichlet_counts),
      ex_counts = list(d = jensen(ecodive::ex_counts),   x = ecodive::ex_counts),
      varespec  = list(d = jensen(ref$varespec_counts),  x = ref$varespec_counts),
      ties      = list(d = dist(ref$ties_x, 'manhattan'), x = NULL) )

    for (nm in names(inputs)) {

      d  <- inputs[[nm]]$d
      ks <- vapply(ref[[nm]]$fits, `[[`, 0L, 'k')

      # The whole procedure: the CH index measures spread between the
      # samples' proportions. (For `ties`, index.G1() was given the points,
      # so the CH index is checked below, from the points.)
      et <- if (!is.null(inputs[[nm]]$x)) enterotype(inputs[[nm]]$x, n_types = ks)
      sel <- if (!is.null(et)) attr(et, 'selection')

      for (i in seq_along(ks)) {

        f    <- ref[[nm]]$fits[[i]]
        k    <- f$k
        info <- paste(nm, k)
        fit  <- pam_of(d, k)

        expect_identical(fit$medoids, sort(as.integer(f$id_med)), info = info)
        expect_true(same_partition(fit$clustering, f$clustering), info = info)
        expect_equal(fit$objective, f$objective, tolerance = 1e-12, info = info)

        # Types are numbered by their medoids' positions.
        expect_identical(fit$clustering[fit$medoids], seq_len(k), info = info)

        dv  <- dist_vector(d)
        sil <- .Call(C_silhouette, dv$dm, dv$n, fit$clustering, k, 1L)[[1]]
        expect_equal(sil, f$silhouette, tolerance = 1e-12, info = info)

        if (is.null(et)) {
          ch <- .Call(C_calinski, dv$dm, dv$n, fit$clustering, k, 1L, ref$ties_x * 1)
          expect_equal(ch, f$ch_medoids, tolerance = 1e-12, info = info)
          next
        }

        expect_equal(sel$ch[[i]], f$ch_medoids, tolerance = 1e-12, info = info)
        expect_equal(sel$silhouette[[i]], mean(f$silhouette), tolerance = 1e-12, info = info)
        expect_equal(sel$objective[[i]], f$objective[[2]], tolerance = 1e-12, info = info)
      }

      if (is.null(et)) next

      best <- which.max(vapply(ref[[nm]]$fits, `[[`, 0, 'ch_medoids'))
      expect_identical(which(sel$selected), best, info = nm)
      expect_equal(et$silhouette, ref[[nm]]$fits[[best]]$silhouette, tolerance = 1e-12, info = nm)
      expect_true(same_partition(et$type, ref[[nm]]$fits[[best]]$clustering), info = nm)
    }

    # The centroid form: index.G1() on the data (Calinski & Harabasz's
    # index) for Euclidean distances.
    et <- enterotype(dist(iris[, 1:4]), n_types = 2:5, ch = 'centroids')
    expect_equal(
      attr(et, 'selection')$ch,
      vapply(ref$euclid$fits, `[[`, 0, 'ch_centroids'),
      tolerance = 1e-12 )
  }



  # Hand-made cases. ====

  # A line of points with an obvious split.
  x  <- matrix(c(1, 2, 3, 4, 20, 21, 22), ncol = 1, dimnames = list(letters[1:7], NULL))
  et <- enterotype(dist(x), n_types = 2:3)
  expect_identical(attr(et, 'selection')$n_types[attr(et, 'selection')$selected], 2L)
  expect_identical(as.integer(et$type), rep(1:2, c(4, 3)))
  expect_identical(attr(et, 'types')$medoid, c('b', 'f'))
  expect_identical(et$sample, letters[1:7])
  expect_identical(rownames(et), letters[1:7])
  # Medoids b (or c) and f: distances 1, 0, 1, 2 and 1, 0, 1.
  expect_equal(unname(attr(et, 'objective')[['swap']]), 6 / 7)

  # The centroid form is Calinski & Harabasz's index on Euclidean data.
  set.seed(1)
  pts <- rbind(matrix(rnorm(30), 15), matrix(rnorm(30, 3), 15), matrix(rnorm(30, c(0, 5)), 15, byrow = TRUE))
  et  <- enterotype(dist(pts), n_types = 2:4, ch = 'centroids')
  for (i in 1:3) {
    k   <- i + 1L
    cl  <- pam_of(dist(pts), k)$clustering
    cen <- colMeans(pts)
    W   <- sum(vapply(seq_len(k), function (c) {
      y <- pts[cl == c, , drop = FALSE]
      sum(sweep(y, 2, colMeans(y))^2) }, 0))
    B   <- sum(vapply(seq_len(k), function (c) sum(cl == c) * sum((colMeans(pts[cl == c, , drop = FALSE]) - cen)^2), 0))
    expect_equal(attr(et, 'selection')$ch[[i]], (B / (k - 1)) / (W / (nrow(pts) - k)), tolerance = 1e-12)
  }

  # The medoid form with no abundances measures spread by the distances,
  # so for Euclidean distances it matches the same form on the points.
  for (i in 1:3) {
    k  <- i + 1L
    cl <- pam_of(dist(pts), k)$clustering
    D  <- as.matrix(dist(pts))
    M  <- which.min(rowSums(D))
    md <- vapply(seq_len(k), function (c) { j <- which(cl == c); j[which.min(rowSums(D[j, j, drop = FALSE]))] }, 1L)
    Tt <- sum(D[, M]^2)
    W  <- sum(D[cbind(seq_along(cl), md[cl])]^2)
    expect_equal(
      attr(enterotype(dist(pts), n_types = k), 'selection')$ch,
      ((Tt - W) / (k - 1)) / (W / (nrow(pts) - k)), tolerance = 1e-12 )
  }

  # Silhouettes against a direct base-R implementation.
  d  <- dist(pts)
  et <- enterotype(d, n_types = 3)
  D  <- as.matrix(d)
  cl <- as.integer(et$type)
  s  <- vapply(seq_along(cl), function (i) {
    own <- cl == cl[i]
    if (sum(own) == 1) return (0)
    a <- sum(D[i, own]) / (sum(own) - 1)
    b <- min(vapply(setdiff(unique(cl), cl[i]), function (c) mean(D[i, cl == c]), 0))
    (b - a) / max(a, b)
  }, 0)
  expect_equal(et$silhouette, s, tolerance = 1e-12)
  expect_equal(attr(et, 'types')$silhouette, as.vector(tapply(s, cl, mean)), tolerance = 1e-12)

  # PAM finds the best medoids on small cases (checked by brute force).
  set.seed(3)
  y <- matrix(runif(24), 12)
  D <- as.matrix(dist(y))
  for (k in 2:3) {
    cmb  <- utils::combn(12, k)
    cost <- apply(cmb, 2, function (m) sum(apply(D[, m, drop = FALSE], 1, min)))
    expect_equal(pam_of(dist(y), k)$objective[[2]] * 12, min(cost), tolerance = 1e-12)
  }



  # The result. ====

  set.seed(45)
  base   <- rep(0.3, 12)
  types  <- list(replace(base, 1, 12), replace(base, 2, 12), replace(base, 3:4, 6))
  grp    <- rep(1:3, c(40, 35, 30))
  counts <- t(vapply(grp, function (g) {
    a <- rgamma(12, types[[g]] * 3)
    as.integer(round(a / sum(a) * 2000)) + 1L
  }, integer(12)))
  dimnames(counts) <- list(sprintf('S%03d', seq_along(grp)), sprintf('G%02d', 1:12))

  et <- enterotype(counts, n_types = 2:5, cpus = 1L)

  expect_inherits(et, c('ecodive_enterotype', 'ecodive_types', 'data.frame'))
  expect_identical(names(et), c('sample', 'type', 'silhouette'))
  expect_identical(et$sample, rownames(counts))
  expect_identical(rownames(et), rownames(counts))
  expect_true(is.factor(et$type))
  expect_identical(levels(et$type), c('1', '2', '3'))
  expect_true(same_partition(et$type, grp))
  expect_identical(attr(et, 'method'), 'enterotype')
  expect_identical(attr(et, 'ch'), 'medoids')
  expect_identical(attr(et, 'display')$title, 'Enterotypes of Jensen-Shannon Distance')
  expect_identical(plot_title(et), 'Enterotypes of Jensen-Shannon Distance')
  expect_identical(names(attr(et, 'objective')), c('build', 'swap'))

  ty <- attr(et, 'types')
  expect_identical(names(ty), c('type', 'n_samples', 'top_feature', 'medoid', 'silhouette'))
  expect_identical(ty$n_samples, as.vector(table(et$type)))
  expect_identical(levels(ty$type), levels(et$type))
  expect_true(all(ty$medoid %in% rownames(counts)))
  expect_identical(et$type[match(ty$medoid, et$sample)], ty$type)

  pr <- attr(et, 'profiles')
  expect_identical(dim(pr), c(3L, 12L))
  expect_identical(dimnames(pr), list(levels(et$type), colnames(counts)))
  p <- props(counts)
  expect_equal(pr[2, ], colMeans(p[et$type == '2', ]))
  expect_equal(unname(rowSums(pr)), rep(1, 3))
  expect_identical(ty$top_feature, colnames(counts)[apply(pr, 1, which.max)])
  expect_true(setequal(ty$top_feature, c('G01', 'G02', 'G03')) || setequal(ty$top_feature, c('G01', 'G02', 'G04')))

  sel <- attr(et, 'selection')
  expect_identical(names(sel), c('n_types', 'ch', 'silhouette', 'objective', 'selected'))
  expect_identical(sel$n_types, 2:5)
  expect_identical(sum(sel$selected), 1L)
  expect_identical(which(sel$selected), which.max(sel$ch))
  expect_identical(sel$n_types[sel$selected], 3L)
  expect_equal(sel$silhouette[sel$selected], mean(et$silhouette))

  # The types as groups in a test.
  pm <- permanova(jensen(counts), et$type, permutations = 9, cpus = 1L)
  expect_true(is.finite(pm$statistic))

  # A dist gives the same types, and profiles only with `counts`.
  d   <- jensen(counts)
  et2 <- enterotype(d, n_types = 2:5)
  expect_identical(et2$type, et$type)
  expect_null(attr(et2, 'profiles'))
  expect_true(all(is.na(attr(et2, 'types')$top_feature)))
  expect_identical(attr(et2, 'display')$title, 'Enterotypes of Jensen-Shannon Distance')
  et3 <- enterotype(d, n_types = 2:5, counts = counts)
  expect_identical(attr(et3, 'profiles'), attr(et, 'profiles'))
  expect_identical(attr(et3, 'selection'), attr(et, 'selection'))

  # `counts` is matched to the distances by name.
  rev_counts <- counts[rev(seq_len(nrow(counts))), ]
  expect_identical(attr(enterotype(d, n_types = 3, counts = rev_counts), 'profiles'), pr)

  # Another distance's title, and an untitled one.
  expect_identical(attr(enterotype(bray(counts), n_types = 3), 'display')$title, 'Enterotypes of Bray-Curtis Dissimilarity')
  expect_identical(attr(enterotype(dist(pts), n_types = 3), 'display')$title, 'Enterotypes')

  # margin = 2.
  et_t <- enterotype(t(counts), n_types = 2:5, margin = 2L)
  expect_identical(no_display(et_t), no_display(et))

  # The default n_types is 2:10, cut to 2:(n - 1).
  expect_identical(attr(enterotype(counts[1:6, ]), 'selection')$n_types, 2:5)
  expect_identical(attr(enterotype(counts[1:30, ]), 'selection')$n_types, 2:10)

  # A single value fixes the number of types.
  et4 <- enterotype(counts, n_types = 4)
  expect_identical(nlevels(et4$type), 4L)
  expect_identical(attr(et4, 'selection')$selected, TRUE)

  # ch = 'centroids'.
  et5 <- enterotype(counts, n_types = 2:5, ch = 'centroids')
  expect_identical(attr(et5, 'ch'), 'centroids')
  expect_identical(et5$type, et$type)
  expect_equal(
    attr(et5, 'selection')$ch[[2]],
    permanova(d, et$type, permutations = 0, cpus = 1L)$statistic[[1]], tolerance = 1e-12 )



  # Invariance. ====

  # Feature order.
  perm_f <- c(7, 3, 12, 1, 5, 9, 2, 11, 4, 8, 10, 6)
  etf <- enterotype(counts[, perm_f], n_types = 2:5)
  expect_identical(etf$type, et$type)
  expect_equal(etf$silhouette, et$silhouette, tolerance = 1e-12)
  expect_equal(attr(etf, 'profiles')[, colnames(counts)], pr, tolerance = 1e-12)
  expect_equal(attr(etf, 'selection'), attr(et, 'selection'), tolerance = 1e-12)

  # Sample order: assignments permute with the samples.
  set.seed(2)
  perm_s <- sample(nrow(counts))
  ets <- enterotype(counts[perm_s, ], n_types = 2:5)
  expect_true(same_partition(ets$type, et$type[perm_s]))
  expect_equal(ets$silhouette, et$silhouette[perm_s], tolerance = 1e-12)
  expect_equal(attr(ets, 'selection')$ch, attr(et, 'selection')$ch, tolerance = 1e-10)



  # Threads: identical results for any number. ====

  big <- counts[rep(seq_len(nrow(counts)), 2), ]
  big <- big + matrix(rpois(length(big), 3), nrow(big))
  rownames(big) <- sprintf('R%03d', seq_len(nrow(big)))
  r1 <- enterotype(big, n_types = 2:6, cpus = 1L)
  for (cpus in 2:max_cpus)
    expect_identical(enterotype(big, n_types = 2:6, cpus = cpus), r1, info = cpus)
  for (cpus in 2:max_cpus)
    expect_identical(enterotype(jensen(big), n_types = 3:4, ch = 'centroids', cpus = cpus),
                     enterotype(jensen(big), n_types = 3:4, ch = 'centroids', cpus = 1L), info = cpus)



  # Edge cases. ====

  # n = 3.
  et3s <- enterotype(counts[c(1, 50, 100), ])
  expect_identical(attr(et3s, 'selection')$n_types, 2L)
  expect_identical(nrow(et3s), 3L)

  # Duplicate samples (zero distances).
  dup <- counts[c(1, 1, 1, 50, 50, 100, 2), ]
  rownames(dup) <- paste0('D', 1:7)
  etd <- enterotype(dup, n_types = 2:4)
  expect_identical(as.integer(etd$type)[1:3], rep(as.integer(etd$type)[[1]], 3))
  expect_true(all(is.finite(etd$silhouette)))
  expect_true(all(table(etd$type) > 0))
  expect_true(all(is.finite(attr(etd, 'selection')$ch) | attr(etd, 'selection')$ch == Inf))

  # Two medoids that are duplicates each keep their own type.
  same <- counts[c(1, 1, 1, 1, 50, 100), ]
  rownames(same) <- paste0('Q', 1:6)
  ets <- enterotype(same, n_types = 4)
  expect_identical(nlevels(ets$type), 4L)
  expect_true(all(table(ets$type) > 0))

  # A feature absent everywhere changes nothing.
  z <- cbind(counts, G13 = 0L)
  etz <- enterotype(z, n_types = 2:5)
  expect_identical(etz$type, et$type)
  expect_equal(attr(etz, 'selection'), attr(et, 'selection'), tolerance = 1e-12)
  expect_identical(attr(etz, 'profiles')[, 'G13'], c(`1` = 0, `2` = 0, `3` = 0))

  # Errors.
  expect_error(enterotype(counts, n_types = 1), '`n_types` must be between 2 and 104')
  expect_error(enterotype(counts, n_types = 105), '`n_types` must be between 2 and 104')
  expect_error(enterotype(counts, n_types = 2.5), 'whole numbers')
  expect_error(enterotype(counts, n_types = NA), 'whole numbers')
  expect_error(enterotype(counts[1:2, ]), 'at least 3 samples')
  zero <- counts; zero[5, ] <- 0L
  expect_error(enterotype(zero), 'positive total.*S005')
  na <- counts; na[5, 1] <- NA
  expect_error(enterotype(na), 'NA')
  neg <- counts; neg[5, 1] <- -1L
  expect_error(enterotype(neg), 'negative')
  expect_error(enterotype(counts, counts = counts), 'only with distances')
  expect_error(enterotype(d, counts = counts[1:10, ]), 'same 105 samples')
  other <- counts; rownames(other)[1] <- 'nobody'
  expect_error(enterotype(d, counts = other), 'lacks samples.*S001')
  expect_error(enterotype(counts, ch = 'nope'))
  expect_error(enterotype(counts, margin = 3), '`margin`')
  expect_error(enterotype(counts, cpus = 0), '`cpus`')



  # print() and plot(). ====

  out <- capture.output(print(et))
  expect_identical(out[[1]], '# Enterotypes of Jensen-Shannon Distance')
  expect_identical(out[[2]], '# 105 samples in 3 types, chosen by the Calinski-Harabasz index from 2 to 5')
  expect_true(any(out == '# Types:'))
  expect_true(any(out == '# Selection:'))
  expect_true(any(grepl('^# Samples: a data.frame: 105 . 3$', out)))
  expect_true(any(grepl('^# Attributes: types, profiles, selection, objective, ch, method$', out)))
  expect_identical(capture.output(r <- print(et))[1:2], out[1:2])
  expect_identical(r, et)
  out4 <- capture.output(print(et4))
  expect_identical(out4[[2]], '# 105 samples in 4 types')
  expect_false(any(out4 == '# Selection:'))
  expect_identical(print_range(c(2, 4, 7)), '2, 4 and 7')

  grDevices::pdf(NULL)
  expect_identical(plot(et), et)
  expect_identical(plot(et, 'selection'), et)
  expect_identical(plot(et, n_features = 2, main = 'Mine'), et)
  expect_identical(plot(et, n_features = Inf), et)
  expect_identical(plot(et4, 'selection'), et4)
  expect_error(plot(et2), 'no profiles')
  expect_identical(plot(et2, 'selection'), et2)
  expect_error(plot(et, 'scores'))
  grDevices::dev.off()

#})
