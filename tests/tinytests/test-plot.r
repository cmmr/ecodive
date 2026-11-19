#test_that("plot", {

  # The plot() methods draw on a null device. They return their argument
  # invisibly, take their labels from the result, and leave R's random
  # number generator alone.

  grDevices::pdf(NULL)

  cpus <- 1L
  dm   <- bray(iris[, 1:4], cpus = cpus)
  cnt  <- round(iris[, 1:4] * 10)
  two  <- droplevels(subset(iris, Species != 'setosa'))

  # Returns its argument invisibly.
  invisibly <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }



  # Ordinations. ====

  ord <- pcoa(dm, cpus = cpus)
  expect_true(invisibly(ord))
  expect_true(invisibly(ord, 'Species', iris))
  expect_true(invisibly(ord, iris$Species, axes = c(1, 3)))
  expect_true(invisibly(ord, 'Sepal.Length', iris))           # a gradient
  expect_true(invisibly(ord, col = iris$Species))              # as for a matrix
  expect_true(invisibly(ord, 'Species', iris, legend = NULL, main = 'x', asp = NA))
  expect_true(invisibly(ord, 'Species', iris, col = c(setosa = 'red', versicolor = 'blue', virginica = 'green')))

  # Every method, with its overlays.
  ords <- list(
    pcoa      = pcoa(dm, correction = 'cailliez', cpus = cpus),
    mds       = mds(dm, cpus = cpus),
    nmds      = nmds(dm, cpus = cpus),
    isomap    = isomap(dm, n_neighbors = 50, cpus = cpus),
    tsne      = tsne(dm, cpus = cpus),
    umap      = umap(dm, cpus = cpus),
    dbrda     = dbrda(dm, 'Species', iris, cpus = cpus),
    capscale  = capscale(dm, data = iris, formula = ~ Species + Sepal.Width, cpus = cpus),
    rda       = rda(cnt, data = iris, formula = ~ Species + Sepal.Width, norm = 'hellinger', cpus = cpus),
    rda_u     = rda(cnt, norm = 'hellinger', cpus = cpus),
    cca       = cca(cnt, 'Species', iris, cpus = cpus),
    ca        = ca(cnt, cpus = cpus),
    pca       = pca(iris[, 1:4], scale = TRUE, cpus = cpus),
    rpca      = rpca(cnt, cpus = cpus),
    rpca_phy  = rpca_phylo(ecodive::ex_counts, ecodive::ex_tree, n_axes = 2, cpus = cpus),
    lda       = lda(iris[, 1:4], 'Species', iris),
    plsda     = plsda(cnt, 'Species', iris),
    oplsda    = oplsda(round(two[, 1:4] * 10), 'Species', two),
    procr     = procrustes(pcoa(dm, n_axes = 2, cpus = cpus), pca(iris[, 1:4], n_axes = 2, cpus = cpus)) )

  for (nm in names(ords)) {
    o <- ords[[nm]]
    expect_true(invisibly(o), info = nm)
    if (nrow(o) == nrow(iris))
      expect_true(invisibly(o, 'Species', iris), info = nm)
    expect_true(invisibly(o, loadings = Inf, biplot = FALSE, centroids = FALSE, target = FALSE), info = nm)
  }

  # One axis: a strip chart, by group.
  ord1 <- lda(two[, 1:4], 'Species', two)
  expect_identical(ncol(ord1), 1L)
  expect_true(invisibly(ord1))
  expect_true(invisibly(ord1, groups = NULL, legend = NULL))
  expect_true(invisibly(pca(iris[, 1:4], n_axes = 1, cpus = cpus)))

  # Titles and labels, from the attributes.
  expect_identical(ord_method_name('pcoa'), 'Principal Coordinates Analysis')
  expect_identical(ord_method_name('oplsda'), 'Orthogonal Partial Least Squares Discriminant Analysis')
  expect_identical(ord_method_name('procrustes'), 'Procrustes Rotation')
  expect_identical(ord_method_name('other'), 'other')

  expect_identical(ord_axis_labels(ord)[1:2], c('PCoA1 (89.0%)', 'PCoA2 (9.3%)'))
  expect_identical(ord_axis_labels(ords$nmds), c('NMDS1', 'NMDS2'))

  expect_null(ord_subtitle(ord))
  expect_true(grepl('^Cailliez correction, constant [0-9.]+$', ord_subtitle(ords$pcoa)))
  expect_true(grepl('^Stress = [0-9.]+$', ord_subtitle(ords$nmds)))
  expect_identical(ord_subtitle(ords$isomap), '50 nearest neighbours')
  expect_true(grepl('^KL divergence = [0-9.]+, perplexity = 30$', ord_subtitle(ords$tsne)))
  expect_identical(ord_subtitle(ords$umap), '15 nearest neighbours')
  expect_identical(
    ord_subtitle(ords$dbrda),
    sprintf('R-squared = %s, adj. R-squared = %s',
            plot_num(attr(ords$dbrda, 'r_squared')), plot_num(attr(ords$dbrda, 'adj_r_squared'))) )
  expect_null(ord_subtitle(ords$rda_u))   # no constraints
  expect_true(grepl('^R2X = [0-9.]+, R2Y = [0-9.]+$', ord_subtitle(ords$oplsda)))
  expect_true(grepl('^Procrustes sum of squares = [0-9.e+-]+$', ord_subtitle(ords$procr)))

  # Groups: by column name (which titles the legend), or by vector, matched
  # by name when both are named.
  g <- ord_plot_groups(ord, 'Species', iris)
  expect_identical(g$f, iris$Species)
  expect_identical(g$title, 'Species')
  g <- ord_plot_groups(ord, as.character(iris$Species), NULL)
  expect_identical(levels(g$f), levels(iris$Species))
  expect_null(g$title)
  g <- ord_plot_groups(ords$plsda, NULL, NULL)   # the stored groups
  expect_identical(g$f, iris$Species)

  named <- ord
  rownames(named) <- paste0('s', 1:150)
  v <- setNames(as.character(iris$Species), paste0('s', 1:150))[150:1]
  expect_identical(as.character(ord_plot_groups(named, v, NULL)$f), as.character(iris$Species))

  g <- ord_plot_groups(ord, c(NA, as.character(iris$Species[-1])), NULL)
  expect_identical(nlevels(g$f), 4L)   # NA is a group of its own

  g <- ord_plot_groups(ord, 'Sepal.Length', iris)
  expect_null(g$f)
  expect_identical(g$values, iris$Sepal.Length)

  # Colours: one per level, the user's (by name or in order) or the palette.
  g    <- ord_plot_groups(ord, 'Species', iris)
  look <- ord_plot_look(g, NULL, NULL)
  expect_identical(look$pal, plot_palette(3))
  expect_identical(look$col, plot_palette(3)[as.integer(iris$Species)])
  expect_identical(look$legend, levels(iris$Species))
  expect_identical(look$pch, 16)
  look <- ord_plot_look(g, c(virginica = 'c', setosa = 'a', versicolor = 'b'), 1)
  expect_identical(look$pal, c('a', 'b', 'c'))
  expect_identical(look$pch, 1)
  look <- ord_plot_look(g, c('x', 'y'), NULL)
  expect_identical(look$pal, c('x', 'y', 'x'))
  look <- ord_plot_look(list(f = NULL, values = NULL), NULL, NULL)
  expect_identical(look$pch, 1)

  # A gradient: low to high, with a legend of round values in the range.
  g    <- ord_plot_groups(ord, 'Sepal.Length', iris)
  look <- ord_plot_look(g, c('white', 'black'), NULL)
  expect_identical(look$col[which.min(iris$Sepal.Length)], '#FFFFFF')
  expect_identical(look$col[which.max(iris$Sepal.Length)], '#000000')
  at <- as.numeric(look$legend)
  expect_true(all(at >= min(iris$Sepal.Length) & at <= max(iris$Sepal.Length)))

  # Centroids of the plotted groups are labelled with the group.
  g    <- ord_plot_groups(ords$dbrda, 'Species', iris)
  look <- ord_plot_look(g, NULL, NULL)
  lc   <- ord_centroid_labels(rownames(attr(ords$dbrda, 'centroids')), g, look)
  expect_identical(lc$labels, c('setosa', 'versicolor', 'virginica'))
  expect_identical(lc$col, plot_palette(3))
  lc <- ord_centroid_labels(c('Speciessetosa', 'Other'), list(f = NULL), NULL)
  expect_identical(lc$labels, c('Speciessetosa', 'Other'))
  g  <- ord_plot_groups(ords$plsda, NULL, NULL)
  lc <- ord_centroid_labels(rownames(attr(ords$plsda, 'centroids')), g, ord_plot_look(g, NULL, NULL))
  expect_identical(lc$labels, levels(iris$Species))

  # Arrows: constraints of factor levels are left out with their centroids,
  # and the longest loadings are kept.
  bp <- ord_plot_biplot(attr(ords$capscale, 'biplot'), 1:2, attr(ords$capscale, 'centroids'))
  expect_identical(rownames(bp), 'Sepal.Width')
  bp <- ord_plot_biplot(attr(ords$capscale, 'biplot'), 1:2, NULL)
  expect_identical(nrow(bp), 3L)
  ld <- ord_plot_loadings(attr(ords$pca, 'loadings'), 1:2, 2)
  len <- rowSums(attr(ords$pca, 'loadings')[, 1:2] ^ 2)
  expect_identical(rownames(ld), names(sort(len, decreasing = TRUE))[1:2])

  m <- plot_arrow_mul(cbind(c(1, -2), c(0.5, 1)), c(4, 1))
  expect_equal(m, 0.8 * min(4 / 2, 1 / 1))

  # Errors.
  expect_error(plot(ord, axes = c(1, 1)), 'two different axis numbers from 1 to 3')
  expect_error(plot(ord, axes = c(1, 4)), 'two different axis numbers from 1 to 3')
  expect_error(plot(ord, axes = 1), 'two different axis numbers')
  expect_error(plot(ord, legend = 'middle'), '`legend` must be one of')
  expect_error(plot(ord, loadings = -1), '`loadings` must be a number')
  expect_error(plot(ord, biplot = NA), '`biplot` must be TRUE or FALSE')
  expect_error(plot(ord, iris$Species[1:10]), 'has 10 values, but `x` has 150 samples')
  expect_error(plot(ord, 'Species'), 'names a column, but `data` is not given')
  expect_error(plot(ord, 'Nope', iris), 'not one')
  expect_error(plot(ord, 'Species', iris[1:10, ]), 'has 10 rows, but `x` has 150 samples')



  # Permutation tests. ====

  res <- list(
    permanova = permanova(dm, 'Species', iris, permutations = 19, cpus = cpus),
    anosim    = anosim(dm, 'Species', iris, permutations = 19, cpus = cpus),
    mrpp      = mrpp(dm, 'Species', iris, permutations = 19, cpus = cpus),
    mantel    = mantel(dm, euclidean(iris[, 1:2]), permutations = 19, cpus = cpus),
    partial   = mantel_partial(dm, 'Sepal.Length', given = 'Petal.Width', data = iris, permutations = 19, cpus = cpus),
    permdisp  = suppressWarnings(permdisp(dm, 'Species', iris, permutations = 19, cpus = cpus)),
    protest   = procrustes_test(ords$procr, pca(iris[, 1:4], n_axes = 2, cpus = cpus), permutations = 19, cpus = cpus) )

  for (nm in names(res)) {
    expect_true(invisibly(res[[nm]]), info = nm)
    expect_true(invisibly(res[[nm]], type = 'null'), info = nm)
    expect_true(invisibly(res[[nm]], type = 'null', main = 'x', breaks = 5, col = 'white'), info = nm)
  }
  expect_true(invisibly(res$mrpp, 'groups', col = 'white'))
  expect_true(invisibly(res$permdisp, 'groups'))
  expect_true(invisibly(res$protest, 'fit', groups = 'Species', data = iris))

  expect_identical(test_subtitle(res$permanova), sprintf(
    'F = %s, p = 0.05, z = %s', plot_num(res$permanova$statistic), plot_num(res$permanova$z) ))
  expect_identical(
    test_subtitle(mrpp(dm, 'Species', iris, permutations = 0, cpus = cpus)),
    sprintf('delta = %s', plot_num(res$mrpp$statistic)) )

  expect_error(plot(res$permanova, 'groups'), '`type` must be "null"')
  expect_error(plot(res$mrpp, 'fit'), '`type` must be "groups" or "null"')
  expect_error(plot(res$protest, 'groups'), '`type` must be "fit" or "null"')
  expect_error(
    plot(permanova(dm, 'Species', iris, permutations = 0, cpus = cpus)),
    'no permutations to plot' )
  # Without permutations, MRPP and PERMDISP still have their groups.
  expect_true(invisibly(mrpp(dm, 'Species', iris, permutations = 0, cpus = cpus)))



  # envfit(). ====

  fit <- envfit(ord, iris, permutations = 19, cpus = cpus)
  expect_true(invisibly(fit))
  plot(ord, 'Species', iris)
  expect_true(invisibly(fit, add = TRUE))
  expect_true(invisibly(fit, add = TRUE, p_max = 0.01))  # nothing left
  expect_true(invisibly(fit, axes = c(1, 3)))
  expect_true(invisibly(envfit(ord, iris[, 1:4], permutations = 0, cpus = cpus)))  # no factors
  expect_true(invisibly(envfit(ord, iris[, 5, drop = FALSE], permutations = 0, cpus = cpus)))  # no vectors
  expect_error(plot(fit, add = NA), '`add` must be TRUE or FALSE')
  expect_error(plot(fit, axes = c(1, 4)), 'two different axis numbers from 1 to 3')
  expect_error(plot(envfit(pcoa(dm, n_axes = 1, cpus = cpus), iris, permutations = 0, cpus = cpus)), 'only one axis')



  # simper(), dirmult_ttest(), ancom(). ====

  sp <- simper(cnt, 'Species', iris, permutations = 19, cpus = cpus)
  expect_true(invisibly(sp))
  expect_true(invisibly(sp, n_taxa = 2, pairs = 'setosa-virginica', sub = 'x'))
  expect_true(invisibly(simper(cnt, 'Species', iris, permutations = 0, cpus = cpus)))
  expect_error(plot(sp, pairs = 'a-b'), '`pairs` must be among')
  expect_error(plot(sp, n_taxa = NA), '`n_taxa` must be a number')

  dt <- dirmult_ttest(cnt, 'Species', iris, cpus = cpus)
  expect_true(invisibly(dt))
  expect_true(invisibly(dt, labels = 0, main = 'x'))
  expect_true(invisibly(dirmult_ttest(cnt, 'Species', iris, contrast = c('setosa', 'virginica'), cpus = cpus)))

  set.seed(1)
  x <- matrix(rpois(200, 20), 10, dimnames = list(NULL, paste0('F', 1:20)))
  x[1:5, 1:3] <- x[1:5, 1:3] * 8
  an <- ancom(x, rep(c('a', 'b'), each = 5), cpus = cpus)
  expect_true(invisibly(an))
  expect_true(invisibly(an, n_features = 5, sub = 'x'))
  an0 <- an
  attr(an0, 'cutoff') <- NA_real_   # nothing could be rejected
  an0$reject <- FALSE
  expect_true(invisibly(an0))



  # bioenv(), forward_select(), varpart(). ====

  be <- bioenv(dm, iris[, 1:4], cpus = cpus)
  expect_true(invisibly(be))
  expect_true(invisibly(bioenv(dm, iris[, 1:4], given = 'Sepal.Width', cpus = cpus)))

  fs <- forward_select(
    iris[, 3:4], ~ Species + Sepal.Length + Sepal.Width, iris,
    norm = 'hellinger', permutations = 19, cpus = cpus )
  expect_true(invisibly(fs))
  expect_identical(forward_select_how(fs), 'RDA (hellinger)')
  fs0 <- fs[0, ]
  attributes(fs0)[c('global', 'next_best', 'stop', 'method', 'norm')] <- attributes(fs)[c('global', 'next_best', 'stop', 'method', 'norm')]
  expect_true(invisibly(fs0))   # nothing selected

  sets <- list(
    a = ~ Species, b = ~ Sepal.Length, c = ~ Sepal.Width, d = ~ Petal.Width )
  for (k in 2:4) {
    vp <- do.call(varpart, c(list(iris[, 3:4]), sets[seq_len(k)], list(data = iris, norm = 'hellinger', cpus = cpus)))
    expect_true(invisibly(vp), info = k)
  }
  expect_true(invisibly(vp, digits = 3, col = 'grey', main = 'x'))
  expect_error(plot(vp, digits = -1), '`digits` must be a number')

  # Each fraction's label sits in its sets' shapes, and no others.
  inside <- function (e, p) {
    r  <- e[['angle']] * pi / 180
    dx <- p[[1]] - e[['x']]
    dy <- p[[2]] - e[['y']]
    u  <-  dx * cos(r) + dy * sin(r)
    v  <- -dx * sin(r) + dy * cos(r)
    (u / e[['a']]) ^ 2 + (v / e[['b']]) ^ 2 < 1
  }
  for (k in 2:4) {
    lay <- varpart_layout(k)
    reg <- varpart_regions(k)
    expect_identical(sort(names(lay$labels)), sort(vapply(reg, paste, '', collapse = ',')), info = k)
    for (r in reg) {
      p   <- lay$labels[[paste(r, collapse = ',')]]
      hit <- vapply(seq_len(k), function (i) inside(lay$shapes[i, ], p), TRUE)
      expect_identical(which(hit), as.integer(r), info = paste(k, paste(r, collapse = ',')))
    }
    # The names lie outside every shape.
    for (i in seq_len(k)) {
      hit <- vapply(seq_len(k), function (j) inside(lay$shapes[j, ], lay$names[i, ]), TRUE)
      expect_false(any(hit), info = paste(k, i))
    }
  }
  expect_null(varpart_layout(5))



  # Helpers. ====

  expect_identical(plot_num(30), '30')
  expect_identical(plot_num(0.0378123), '0.0378')
  expect_identical(plot_num(NA), 'NA')
  expect_identical(plot_p(0.001), '0.001')
  expect_identical(plot_cap('spearman'), 'Spearman')
  expect_identical(plot_count(TRUE, 'n'), Inf)
  expect_identical(plot_count(FALSE, 'n'), 0)
  expect_identical(plot_count(2.7, 'n'), 2)
  expect_null(plot_legend_pos(NULL))
  expect_null(plot_legend_pos(FALSE))
  expect_identical(plot_legend_pos('topleft'), 'topleft')

  # Defaults give way to the user's arguments.
  f <- function (...) list(...)
  expect_identical(plot_do(f, list(main = 'u'), list(main = 'd', xlab = 'x')), list(main = 'u', xlab = 'x'))
  expect_identical(plot_do(f, list(main = NULL), list(main = 'd')), list(main = NULL))
  expect_identical(plot_arg(list(sub = 's'), 'sub', 'd'), 's')
  expect_identical(plot_arg(list(), 'sub', 'd'), 'd')

  grDevices::dev.off()

  # A title too wide for the plot is wrapped.
  grDevices::pdf(NULL, width = 3, height = 3)
  plot.new()
  long <- 'Permutational multivariate analysis of variance (PERMANOVA)'
  expect_true(grepl('\n', plot_fit_title(long), fixed = TRUE))
  expect_identical(gsub('\n', ' ', plot_fit_title(long), fixed = TRUE), long)
  expect_identical(plot_fit_title('PCoA'), 'PCoA')
  expect_identical(plot_fit_title(NULL), NULL)
  grDevices::dev.off()

  grDevices::pdf(NULL, width = 12, height = 7)
  plot.new()
  expect_identical(plot_fit_title(long), long)
  grDevices::dev.off()



  # R's random number generator is untouched. ====

  grDevices::pdf(NULL)
  set.seed(42)
  seed <- .Random.seed
  plot(ord, 'Species', iris)
  plot(ords$lda)
  plot(ord1)
  plot(res$permanova)
  plot(sp)
  plot(vp)
  expect_identical(.Random.seed, seed)
  grDevices::dev.off()

#})
