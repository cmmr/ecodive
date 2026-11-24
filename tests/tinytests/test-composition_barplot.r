#test_that("composition_barplot", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  # Draws on a null device; returns the visible-or-not result too.
  cb <- function (...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    v <- withVisible(composition_barplot(...))
    v$value$visible <- v$visible
    v$value
  }

  set.seed(1)
  n <- 24; p <- 30
  x <- matrix(stats::rpois(n * p, rep(exp(seq(3, -1, length.out = p)), each = n)), n, p,
              dimnames = list(paste0('S', 1:n), paste0('F', 1:p)))
  x[3, ] <- x[3, ] + 1L   # some variety
  prop  <- x / rowSums(x)
  means <- colMeans(prop)
  pc    <- min(x[x > 0]) / 2



  # Proportions and the features drawn. ====

  r <- cb(x, n_features = 6)
  expect_false(r$visible)
  expect_identical(names(r)[1:4], c('proportions', 'col', 'order', 'x'))
  P <- r$proportions
  top <- order(-means, seq_len(p))[1:6]
  expect_identical(rownames(P), c(colnames(x)[top], 'Other'))
  expect_identical(colnames(P), rownames(x)[r$order])
  expect_equal(unname(colSums(P)), rep(1, n), tolerance = 1e-15)
  expect_equal(P[1:6, ], t(prop[r$order, top]), tolerance = 1e-15)
  expect_equal(unname(P['Other', ]), unname(1 - colSums(P[1:6, ])), tolerance = 1e-15)
  expect_equal(unname(P['Other', ]), unname(rowSums(prop[r$order, -top])), tolerance = 1e-14)

  # The default: 12 features.
  expect_identical(nrow(cb(x)$proportions), 13L)

  # Ties go to the feature that comes first.
  tie <- matrix(c(5, 5, 1, 5, 5, 1, 5, 5, 1), 3, 3, byrow = TRUE,
                dimnames = list(c('a', 'b', 'c'), c('X', 'Y', 'Z')))
  expect_identical(rownames(cb(tie, n_features = 1)$proportions), c('X', 'Other'))
  tie2 <- tie[, c(2, 1, 3)]
  expect_identical(rownames(cb(tie2, n_features = 1)$proportions), c('Y', 'Other'))

  # All features drawn: no "Other". Features that are zero everywhere are not drawn.
  z <- cbind(x[, 1:4], Empty = 0)
  rz <- cb(z, n_features = Inf)
  expect_identical(rownames(rz$proportions), colnames(x)[order(-colMeans(x[, 1:4] / rowSums(x[, 1:4])), 1:4)])
  expect_identical(rownames(cb(z, n_features = 4)$proportions), rownames(rz$proportions))

  # Proportions are drawn as they are (e.g. a table of proportions).
  expect_equal(cb(prop, n_features = 6, clust = 'none')$proportions, cb(x, n_features = 6, clust = 'none')$proportions, tolerance = 1e-15)

  # Samples of zeros are left out, with a warning; `order` is in the input's positions.
  x0 <- x
  x0[c(2, 7), ] <- 0L
  expect_warning(r0 <- cb(x0, n_features = 6, clust = 'none'), "all zero: 'S2', 'S7'")
  expect_identical(r0$order, setdiff(1:n, c(2L, 7L)))
  expect_identical(colnames(r0$proportions), rownames(x)[-c(2, 7)])



  # Feature groups. ====

  fg <- rep(c('Firm', 'Bact', 'Prot', 'Actino', 'Verr', 'Fuso'), each = 5)
  fg[c(4, 12)] <- NA
  ga  <- ifelse(is.na(fg), 'Unassigned', fg)
  gm  <- tapply(means, ga, sum)[unique(ga)]
  sel <- names(gm)[order(-gm, seq_along(gm))][1:3]
  r   <- cb(x, feature_groups = fg, n_feature_groups = 3, n_features = 2)
  P   <- r$proportions
  want <- unlist(lapply(sel, function (g) {
    mem <- which(ga == g)
    mem <- mem[order(-means[mem], mem)]
    c(colnames(x)[mem[1:2]], if (length(mem) > 2) paste('Other', g))
  }))
  expect_identical(rownames(P), c(want, 'Other'))
  for (g in sel) {
    mem  <- which(ga == g)
    mem  <- mem[order(-means[mem], mem)]
    rest <- mem[-(1:2)]
    if (length(rest)) expect_equal(unname(P[paste('Other', g), ]), unname(rowSums(prop[r$order, rest, drop = FALSE])), tolerance = 1e-14)
  }
  expect_equal(unname(P['Other', ]), unname(rowSums(prop[r$order, !ga %in% sel])), tolerance = 1e-14)
  expect_equal(unname(colSums(P)), rep(1, n), tolerance = 1e-15)

  # NA is "Unassigned", a group like any other.
  fu <- rep(NA, p)
  fu[1:5] <- 'G'
  ru <- cb(x, feature_groups = fu, n_features = 1)
  expect_true(all(c('Other Unassigned', 'Other G') %in% rownames(ru$proportions)))

  # Matched by name, in any order; by position otherwise.
  named <- stats::setNames(fg, colnames(x))
  perm  <- sample(p)
  expect_identical(cb(x, feature_groups = named[perm], n_feature_groups = 3, n_features = 2)$proportions, P)
  expect_identical(cb(x, feature_groups = unname(fg), n_feature_groups = 3, n_features = 2)$proportions, P)

  # A group with no more members than n_features has no "Other <group>".
  rs <- cb(x, feature_groups = fg, n_feature_groups = Inf, n_features = 5)
  expect_false(any(grepl('^Other ', rownames(rs$proportions))))
  expect_false('Other' %in% rownames(rs$proportions))



  # microshades 1.13's choice of members, and their shades. ====

  if (at_home()) {
    ref <- readRDS(file.path('expected', 'composition_barplot.rds'))
    for (nm in names(ref)) {
      cs <- ref[[nm]]
      rr <- cb(cs$counts, feature_groups = cs$phylum, n_feature_groups = cs$n_feature_groups,
               n_features = cs$n_features, clust = 'none')
      seg <- rownames(rr$proportions)
      # The groups, in order (microshades' are the largest summed proportions).
      grp <- vapply(seg, function (s) {
        if (s == 'Other') return (NA_character_)
        if (startsWith(s, 'Other ')) return (substring(s, 7))
        ph <- cs$phylum[match(s, colnames(cs$counts))]
        if (is.na(ph)) 'Unassigned' else ph
      }, '')
      expect_identical(unique(grp[!is.na(grp)]), cs$selected, info = nm)
      # The members, from the darkest shade to the lightest, then "Other".
      ms <- ifelse(cs$cdf$subgroup == 'Other', paste('Other', cs$cdf$group), cs$cdf$subgroup)
      expect_identical(seg[seg != 'Other'], ms, info = nm)
      # Shades darken with abundance, and the group's "Other" is lightest:
      # the same order of shade as microshades'.
      lum <- function (cl) as.vector(grDevices::convertColor(
        t(grDevices::col2rgb(cl)) / 255, from = 'sRGB', to = 'Lab')[, 1])
      for (g in cs$selected) {
        i <- which(grp == g)
        expect_identical(order(lum(rr$col[i])), order(cs$cdf$shade[cs$cdf$group == g]), info = paste(nm, g))
      }
    }
  }



  # Order. ====

  d <- aitchison(x, pseudocount = pc)
  expect_identical(cb(x)$order, order(pcoa(d, n_axes = 1L)[, 1]))
  expect_identical(cb(x)$order, ecodive:::dist_clust_order(d, 'pcoa'))
  expect_identical(cb(x, clust = 'average')$order, stats::hclust(stats::as.dist(as.matrix(d)), 'average')$order)
  expect_identical(cb(x, clust = 'ward.D2', pseudocount = 1)$order,
                   stats::hclust(stats::as.dist(as.matrix(aitchison(x, pseudocount = 1))), 'ward.D2')$order)
  expect_identical(cb(x, clust = 'none')$order, 1:n)
  expect_identical(cb(x, clust = FALSE)$order, 1:n)
  expect_identical(cb(x, clust = NULL)$order, 1:n)

  # Left-out samples: ordered by the distances between the others.
  keep <- setdiff(1:n, c(2L, 7L))
  d0   <- aitchison(x0[keep, ], pseudocount = pc)
  expect_identical(suppressWarnings(cb(x0)$order), keep[order(pcoa(d0, n_axes = 1L)[, 1])])

  # Groups: sorted by group, the order kept within each (a stable sort).
  grp <- rep(c('b', 'a', 'c'), length.out = n)
  base <- cb(x)$order
  rg <- cb(x, groups = grp)
  expect_identical(rg$order, base[order(grp[base])])
  meta <- data.frame(row.names = rev(rownames(x)), site = rev(grp))
  expect_identical(cb(x, groups = 'site', data = meta)$order, rg$order)
  expect_identical(cb(x, groups = factor(grp, c('c', 'b', 'a')))$order,
                   base[order(factor(grp, c('c', 'b', 'a'))[base])])
  # Numbers are categories too; NA is a group of its own, last.
  gn <- rep(c(10, 2, NA), length.out = n)
  expect_identical(cb(x, groups = gn)$order, base[order(gn[base], na.last = TRUE)])

  # Bars 1 wide, 0.1 apart, a bar's width more between groups.
  gaps <- diff(rg$x)
  expect_equal(rg$x[[1]], 0.6)
  expect_equal(sort(unique(round(gaps, 10))), c(1.1, 2.1))
  expect_identical(sum(abs(gaps - 2.1) < 1e-10), 2L)
  # Above 60 samples, no space between bars (but between groups).
  big <- matrix(stats::rpois(200 * 8, 5), 200, 8)
  rb  <- cb(big, groups = rep(1:2, each = 100), clust = 'none', heatmap = TRUE)
  expect_equal(diff(rb$x), c(rep(1, 99), 2, rep(1, 99)))



  # Colours. ====

  r <- cb(x, n_features = 4)
  expect_identical(names(r$col), rownames(r$proportions))
  expect_identical(unname(r$col), c(grDevices::hcl.colors(4, 'Dark 3'), 'grey80'))

  r <- cb(x, feature_groups = fg, n_feature_groups = 3, n_features = 2)
  hue <- seq(15, 375, length.out = 4)
  # Unassigned (third) has only 2 members, so no "Other Unassigned": its
  # shades are the same as with one.
  expect_identical(sel[[3]], 'Unassigned')
  expect_identical(
    unname(r$col),
    c(unlist(lapply(1:3, function (g) grDevices::hcl(hue[[g]], 55, seq(35, 80, length.out = 3))[if (g < 3) 1:3 else 1:2])), 'grey80') )

  r2 <- cb(x, feature_groups = fg, n_feature_groups = 3, n_features = 2,
           col = c(Other = 'black', F1 = 'red'))
  expect_identical(r2$col[['Other']], 'black')
  expect_identical(r2$col[['F1']], 'red')
  expect_identical(r2$col[setdiff(names(r2$col), c('Other', 'F1'))], r$col[setdiff(names(r$col), c('Other', 'F1'))])
  r3 <- cb(x, n_features = 2, col = c('red', 'green', 'blue'))
  expect_identical(unname(r3$col), c('red', 'green', 'blue'))
  expect_error(cb(x, n_features = 2, col = c('red', 'green')), 'one colour per segment \\(3\\)')
  expect_error(cb(x, n_features = 2, col = c(Nope = 'red')), "not drawn: 'Nope'")



  # The heatmap: centred clr of the drawn features. ====

  r <- cb(x, n_features = 5, heatmap = TRUE)
  lx  <- log(x + pc)
  clr <- (lx - rowMeans(lx))[r$order, rownames(r$proportions)[1:5]]
  clr <- t(clr) - colMeans(clr)
  expect_equal(r$clr, clr, tolerance = 1e-12)
  expect_identical(dimnames(r$clr), list(rownames(r$proportions)[1:5], colnames(r$proportions)))
  r <- cb(x, feature_groups = fg, n_feature_groups = 3, n_features = 2, heatmap = TRUE, pseudocount = 1)
  expect_identical(rownames(r$clr), rownames(r$proportions)[!grepl('^Other', rownames(r$proportions))])
  # Without zeros, no pseudocount.
  xp  <- x + 1L
  r   <- cb(xp, n_features = 3, heatmap = TRUE, clust = 'none')
  lx  <- log(xp)
  clr <- (lx - rowMeans(lx))[, rownames(r$proportions)[1:3]]
  expect_equal(r$clr, t(clr) - colMeans(clr), tolerance = 1e-12)

  # The pseudocount's default warns of nothing.
  expect_silent(cb(x, heatmap = TRUE))



  # The same for any cpus, sparse and dense input, either margin. ====

  args <- list(feature_groups = fg, n_feature_groups = 3, n_features = 2, heatmap = TRUE)
  ref  <- do.call(cb, c(list(x, cpus = 1L), args))
  for (cpus in 2:max_cpus)
    expect_identical(do.call(cb, c(list(x, cpus = cpus), args)), ref, info = cpus)
  expect_identical(do.call(cb, c(list(t(x), margin = 2L), args)), ref)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(x, sparse = TRUE)
    expect_identical(do.call(cb, c(list(sp), args)), ref)
    expect_identical(do.call(cb, c(list(methods::as(sp, 'TsparseMatrix')), args)), ref)
    expect_identical(do.call(cb, c(list(Matrix::t(sp), margin = 2L), args)), ref)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(do.call(cb, c(list(slam::as.simple_triplet_matrix(x)), args)), ref)
  # Many samples, with threads.
  bx <- matrix(stats::rpois(500 * 40, 3), 500, 40)
  rb1 <- cb(bx, heatmap = TRUE, groups = rep(1:4, 125), cpus = 1L)
  for (cpus in 2:max_cpus)
    expect_identical(cb(bx, heatmap = TRUE, groups = rep(1:4, 125), cpus = cpus), rb1, info = cpus)



  # Drawing. ====

  grDevices::pdf(NULL)
  op   <- graphics::par(c('mfrow', 'mai', 'cex'))
  seed <- .Random.seed
  v <- withVisible(composition_barplot(x, heatmap = TRUE, groups = grp, legend = NULL))
  expect_false(v$visible)
  expect_identical(graphics::par(c('mfrow', 'mai', 'cex')), op)
  composition_barplot(x, legend = 'bottom', feature_groups = fg, heatmap = TRUE, main = 'Mine', sub = 'Sub')
  composition_barplot(big)
  composition_barplot(x[1:2, ], heatmap = TRUE)
  composition_barplot(x[1, , drop = FALSE])
  composition_barplot(matrix(5, 4, 3))
  composition_barplot(x[, 1, drop = FALSE], heatmap = TRUE)
  composition_barplot(ecodive::ex_counts, groups = c('a', 'a', 'b', 'b'), legend = 'bottom')
  expect_identical(graphics::par(c('mfrow', 'mai', 'cex')), op)
  expect_identical(.Random.seed, seed)
  # Inside the user's own panels.
  graphics::par(mfrow = c(1, 2))
  composition_barplot(x)
  expect_identical(graphics::par('mfrow'), c(1L, 2L))
  graphics::par(op)
  grDevices::dev.off()



  # Errors. ====

  expect_error(cb(x, n_features = 0), '`n_features` must be at least 1')
  expect_error(cb(x, n_feature_groups = 0, feature_groups = fg), '`n_feature_groups` must be at least 1')
  expect_error(cb(x, feature_groups = fg[-1]), '`feature_groups` has 29 values, but `counts` has 30 features')
  expect_error(cb(x, feature_groups = list(1)), '`feature_groups` must be a vector')
  expect_error(cb(x, clust = 'nope'), '`clust` must be')
  expect_error(cb(x * 0L), 'no sample with counts above 0')
  expect_error(cb(-x), 'negative')
  expect_error(cb(x, legend = 'left'), '`legend` must be')
  expect_error(cb(x, heatmap = NA), '`heatmap` must be TRUE or FALSE')
  expect_error(cb(x, groups = 1:3), '`groups` has 3 values, but `counts` has 24 samples')
  expect_error(cb(x, pseudocount = -1), '`pseudocount`')

#})
