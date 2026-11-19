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

  # One axis: a strip chart, by group, with the same colours, symbols and
  # legends as two.
  ord1 <- lda(two[, 1:4], 'Species', two)
  expect_identical(ncol(ord1), 1L)
  expect_true(invisibly(ord1))
  expect_true(invisibly(ord1, groups = NULL, legend = NULL))
  expect_silent(plot(ord1, pch = c(1, 2)))
  pc1 <- pca(iris[, 1:4], n_axes = 1, cpus = cpus)
  expect_true(invisibly(pc1))
  expect_silent(plot(pc1, iris$Sepal.Length))
  expect_silent(plot(pc1, col = iris$Species))
  expect_silent(plot(pc1, 'Species', iris, xaxt = 'n'))

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

  # Unused levels keep their colours, and are left out of the legend.
  f    <- factor(rep(c('a', 'c'), 75), levels = c('a', 'b', 'c'))
  look <- ord_plot_look(ord_plot_groups(ord, f, NULL), c('red', 'green', 'blue'), NULL)
  expect_identical(unique(look$col), c('red', 'blue'))
  expect_identical(look$legend, c('a', 'c'))
  expect_identical(look$legend_col, c('red', 'blue'))
  expect_identical(look$used, c(TRUE, FALSE, TRUE))

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

  # A constant one has one value in its legend; one of NAs, none.
  look <- ord_plot_look(ord_plot_groups(ord, rep(3, 150), NULL), NULL, NULL)
  expect_identical(look$legend, '3')
  look <- ord_plot_look(ord_plot_groups(ord, rep(NA_real_, 150), NULL), NULL, NULL)
  expect_null(look$legend)
  expect_identical(look$col, rep(PLOT_MUTED, 150))
  expect_silent(plot(ord, rep(NA_real_, 150)))

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

  # The x axis spans the whole histogram, and the observed statistic.
  plot(res$mantel)
  usr <- graphics::par('usr')
  brk <- graphics::hist(attr(res$mantel, 'r_perm'), plot = FALSE)$breaks
  expect_true(usr[[1]] <= min(brk) && usr[[2]] >= max(brk, res$mantel$statistic))
  plot(res$mantel, breaks = 3)
  usr <- graphics::par('usr')
  brk <- graphics::hist(attr(res$mantel, 'r_perm'), breaks = 3, plot = FALSE)$breaks
  expect_true(usr[[1]] <= min(brk) && usr[[2]] >= max(brk))
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
  expect_silent(plot(sp, col = 'grey40'))
  expect_error(plot(sp, pairs = 'a-b'), '`pairs` must be among')
  expect_error(plot(sp, n_taxa = NA), '`n_taxa` must be a number')
  expect_error(plot(sp, n_taxa = 0), '`n_taxa` must be at least 1')

  # The rows go to dotchart() in the order it sorts them to, last pair
  # first, so that it moves none of them (R < 3.6.2 didn't move `pch`), each
  # pair's largest last (at the top).
  rows <- simper_plot_rows(sp$pair, unique(sp$pair), 2)
  expect_identical(sort.list(as.numeric(rows$groups), decreasing = TRUE), seq_along(rows$groups))
  expect_identical(levels(rows$groups), unique(sp$pair))
  for (p in unique(sp$pair)) {
    a <- sp$average[rows$idx[rows$groups == p]]
    expect_identical(a, sort(sp$average[sp$pair == p], decreasing = TRUE)[2:1], info = p)
  }

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
  expect_silent(plot(an, col = 'red'))
  expect_error(plot(an, n_features = 0), '`n_features` must be at least 1')
  an0 <- an
  attr(an0, 'cutoff') <- NA_real_   # nothing could be rejected
  an0$reject <- FALSE
  expect_true(invisibly(an0))



  # ancombc() and its post-hoc tests. ====

  ax <- rbind(
    c(1, 10, 16, 0, 0, 12, 0, 7), c(0, 10, 45, 2, 1, 25, 1, 9),
    c(4,  8, 22, 1, 5, 14, 1, 5), c(1, 11, 22, 1, 4, 14, 1, 3),
    c(4,  7, 36, 0, 0, 19, 0, 9), c(2, 30, 26, 0, 0, 12, 1, 0),
    c(1, 36, 10, 1, 5, 15, 2, 9), c(1, 25, 11, 2, 2,  8, 0, 5),
    c(1, 36, 23, 0, 2, 11, 1, 2), c(3, 60, 33, 3, 4, 14, 1, 6),
    c(3, 57, 30, 3, 1,  3, 0, 4), c(0, 65, 23, 0, 2,  2, 0, 8),
    c(2, 46, 19, 0, 1,  2, 0, 1), c(0, 41, 26, 0, 2,  0, 0, 3),
    c(3, 46, 20, 0, 0,  3, 0, 0) )
  dimnames(ax) <- list(paste0('S', 1:15), paste0('F', 1:8))
  am <- data.frame(
    status = rep(c('mild', 'moderate', 'severe'), each = 5),
    age    = c(25, 48, 23, 35, 52, 51, 18, 29, 40, 29, 44, 39, 26, 37, 46),
    row.names = rownames(ax) )
  fit  <- ancombc2(ax, ~ status + age, am, groups = 'status', cpus = cpus)
  fit1 <- ancombc(ax + 1, ~ status + age, am, cpus = cpus)

  # Volcano plots, of the first covariate after the intercept by default.
  expect_true(invisibly(fit))
  expect_true(invisibly(fit1, 'age', labels = 0, main = 'x'))
  expect_true(invisibly(fit, '(Intercept)'))
  expect_identical(ancombc_plot_name(fit), 'ANCOM-BC2')
  expect_identical(ancombc_plot_name(fit1), 'ANCOM-BC')
  expect_identical(ancombc_plot_groups(fit), 'status')
  expect_identical(ancombc_plot_groups(fit1), 'the groups')
  expect_error(plot(fit, 'nope'), '`covariate` must be one of')
  expect_error(plot(fit, c('age', 'statussevere')), '`covariate` must be one of')

  # Heatmaps of the comparisons; dot charts of the global and trend tests.
  ap <- ancombc_pairwise(fit, cpus = cpus)
  ad <- ancombc_dunnett(fit, n_bootstraps = 20L, cpus = cpus)
  ag <- ancombc_global(fit, cpus = cpus)
  at <- ancombc_trend(fit, n_bootstraps = 20L, cpus = cpus)
  for (r in list(ap, ad, ag, at)) {
    expect_true(invisibly(r), info = class(r)[[1]])
    expect_true(invisibly(r, n_features = 3, main = 'x', sub = 'y'), info = class(r)[[1]])
    expect_error(plot(r, n_features = 0), '`n_features` must be at least 1', info = class(r)[[1]])
  }
  expect_true(invisibly(ap, digits = 0, col = grDevices::gray.colors(5)))
  expect_error(plot(ap, digits = -1), '`digits` must be a number')

  # The heatmap restores the margins it widened.
  mai <- graphics::par('mai')
  plot(ap)
  expect_equal(graphics::par('mai'), mai)

  # Features with no estimate, and a feature that's all NA, are drawn blank.
  ap_na <- ap
  ap_na$log_fc[ap_na$feature == 'F1'] <- NaN
  ap_na$p_adjusted[ap_na$feature == 'F1'] <- NA
  expect_silent(plot(ap_na))
  ap_na$log_fc <- NaN
  expect_silent(plot(ap_na))



  # bioenv(), forward_select(), varpart(). ====

  be <- bioenv(dm, iris[, 1:4], cpus = cpus)
  expect_true(invisibly(be))
  expect_true(invisibly(be, xaxt = 'n'))
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
  expect_true(plot_x_axis(list()))
  expect_false(plot_x_axis(list(axes = FALSE)))
  expect_false(plot_x_axis(list(xaxt = 'n')))
  expect_identical(plot_dot_args(list(col = 'r', sub = 's', cex = 1)), list(cex = 1, color = 'r'))
  expect_identical(plot_dot_args(list(col = 'r', color = 'b')), list(color = 'b'))

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
