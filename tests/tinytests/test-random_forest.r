#test_that("random_forest", {

# Inputs and reference results, written by data-raw/expected/random_forest.r
# (scikit-learn 1.9.1, q2-sample-classifier 2026.10.0, randomForest, ranger).
ref <- readRDS(file.path('expected', 'random_forest.rds'))
rf_details <- ecodive:::random_forest_details

max_cpus <- if (at_home()) 4L else 2L

rel_ok <- function (a, b, tol = 1e-12) {
  a <- as.vector(a); b <- as.vector(b)
  length(a) == length(b) && identical(is.na(a), is.na(b)) &&
    all(a[!is.na(a)] == b[!is.na(b)] |
        abs(a - b)[!is.na(a)] <= tol * pmax(abs(a), abs(b), 1)[!is.na(a)])
}

# Arguments of a scikit-learn case, in ecodive's names.
case_args <- function (cs, cpus) {
  dat  <- ref[[cs$data]]
  x    <- dat$x
  if (cs$sparse) x <- slam_or_matrix(x)
  args <- list(counts = x, y = dat$y, norm = 'none', n_trees = ncol(cs$inbag),
               n_folds = 0, numpy_random_state = cs$random_state,
               sparse = cs$sparse, cpus = cpus)
  pr <- cs$params
  if (!is.null(pr$max_features))
    args$max_features <- if (identical(pr$max_features, 'log2')) 'log2' else as.numeric(pr$max_features)
  for (a in c('criterion', 'class_weight', 'max_samples', 'min_samples_leaf',
              'min_samples_split', 'max_depth'))
    if (!is.null(pr[[a]])) args[[a]] <- pr[[a]]
  args
}

# A sparse table without Matrix or slam installed: a simple_triplet_matrix.
slam_or_matrix <- function (x) {
  nz <- which(x != 0, arr.ind = TRUE)
  structure(
    list(i = nz[, 1], j = nz[, 2], v = x[nz], nrow = nrow(x), ncol = ncol(x),
         dimnames = dimnames(x)),
    class = 'simple_triplet_matrix' )
}



# The random draws: numpy's legacy RandomState and our_rand_r(). ====

draws <- function (kind, seed, n, count) .Call(ecodive:::C_rf_rng_test, kind, seed, n, count)
expect_identical(draws(0L, 0, 0, 1000L), ref$rng$random_sample)
for (n in names(ref$rng$randint))
  expect_identical(draws(1L, 1, as.numeric(n), 1000L), ref$rng$randint[[n]], info = n)
expect_identical(draws(3L, 7, 0, 1000L), ref$rng$tree_seeds)
expect_identical(draws(2L, 12345, 0, 1000L), ref$rng$our_rand_r)
for (z in ref$rng$np_sum)
  expect_identical(.Call(ecodive:::C_rf_rng_test, 4L, 0, z$x, 0L), z$sum)



# Identical trees to scikit-learn's, given its draws. ====

if (at_home()) {
  for (nm in names(ref$cases)) {
    cs  <- ref$cases[[nm]]
    fit <- suppressWarnings(do.call(rf_details, case_args(cs, max_cpus)))
    f   <- fit$forests[[1]]
    tr  <- cs$trees
    expect_identical(f$inbag, cs$inbag, info = nm)
    expect_identical(f$offsets, tr$offsets, info = nm)
    expect_identical(f$left, tr$left, info = nm)
    expect_identical(f$right, tr$right, info = nm)
    expect_identical(f$feature, tr$feature, info = nm)
    expect_identical(f$n_node_samples, tr$n_node_samples, info = nm)
    expect_true(rel_ok(f$threshold, tr$threshold), info = nm)
    expect_true(rel_ok(f$value, tr$value), info = nm)
    expect_true(rel_ok(f$impurity, tr$impurity), info = nm)
    expect_true(rel_ok(f$weighted_n_node_samples, tr$weighted_n_node_samples), info = nm)

    feats <- colnames(ref[[cs$data]]$x)
    imp   <- fit$importance
    expect_true(rel_ok(imp$impurity[match(feats, imp$feature)], cs$importance), info = nm)

    # OOB predictions (where any tree left the sample out) and score.
    ok <- !is.na(fit$oob_value[1, ])
    expect_true(rel_ok(fit$oob_value[, ok], cs$oob[, ok]), info = nm)
    if (all(ok)) {
      sc <- fit$performance[[if (cs$data == 'reg') 'oob_r_squared' else 'oob_accuracy']]
      expect_true(rel_ok(sc, cs$oob_score, 1e-10), info = nm)
    }

    # Predictions of the training table, through predict().
    x  <- ref[[cs$data]]$x
    pr <- if (cs$data == 'reg') predict(fit, x, cpus = 1) else t(predict(fit, x, 'prob', cpus = 1))
    expect_true(rel_ok(unname(pr), cs$pred), info = nm)
  }
}



# q2-sample-classifier's nested cross-validation, given its folds. ====

if (at_home()) {
  q <- ref$q2_cls
  fit <- rf_details(slam_or_matrix(ref$three$x), ref$three$y, norm = 'none',
                    n_trees = 20, folds = q$folds, numpy_random_state = 5,
                    sparse = TRUE, cpus = max_cpus)
  pr <- fit$predictions
  expect_identical(as.character(pr$predicted), q$predictions)
  expect_true(rel_ok(as.matrix(pr[, paste0('prob_', q$classes)]), q$probabilities))
  expect_true(rel_ok(fit$importance$impurity_cv_mean[match(names(q$importances), fit$importance$feature)],
                     q$importances))
  acc <- vapply(1:5, function (f) mean(pr$predicted[pr$fold == f] == pr$observed[pr$fold == f]), 0)
  expect_true(rel_ok(acc, q$scores))
  expect_true(rel_ok(fit$performance[['accuracy_fold_mean']], mean(q$scores)))
  for (cl in c(q$classes, 'micro', 'macro')) {
    r <- fit$roc[fit$roc$class == cl, ]
    expect_true(rel_ok(r$fpr, q$roc[[cl]]$fpr), info = cl)
    expect_true(rel_ok(r$tpr, q$roc[[cl]]$tpr), info = cl)
  }
  expect_true(rel_ok(fit$performance[['auc_micro']], q$roc$micro$auc))
  expect_true(rel_ok(fit$performance[['auc_macro']], q$roc$macro$auc))

  q <- ref$q2_reg
  fit <- rf_details(slam_or_matrix(ref$reg$x), ref$reg$y, norm = 'none',
                    n_trees = 20, folds = q$folds, numpy_random_state = 5,
                    sparse = TRUE, cpus = max_cpus)
  expect_true(rel_ok(fit$predictions$predicted, q$predictions))
  expect_true(rel_ok(fit$importance$impurity_cv_mean[match(names(q$importances), fit$importance$feature)],
                     q$importances))
  mse <- vapply(1:5, function (f) with(fit$predictions[fit$predictions$fold == f, ],
                                       mean((predicted - observed)^2)), 0)
  expect_true(rel_ok(mse, q$scores))
}



# roc_curve() and auc() on scores with ties. ====

rb <- ref$roc$bin
rc <- ecodive:::rf_roc_curve(rb$truth, rb$score)
expect_identical(rc$fpr, rb$fpr)
expect_identical(rc$tpr, rb$tpr)
expect_identical(rc$threshold, rb$thresholds)
expect_true(rel_ok(ecodive:::rf_auc(rc$fpr, rc$tpr), rb$auc))
# Binary AUC is the Mann-Whitney statistic, ties counted 1/2.
pos <- rb$score[rb$truth == 1]; neg <- rb$score[rb$truth == 0]
expect_true(rel_ok(rb$auc, mean(outer(pos, neg, '>') + 0.5 * outer(pos, neg, '=='))))



# Folds: scikit-learn's StratifiedKFold allocation; KFold's sizes. ====

y3 <- ref$three$y
f3 <- ecodive:::rf_folds(y3, TRUE, 5L, 0L)
expect_identical(sort(unique(f3)), 1:5)
for (cl in levels(y3)) {
  cnt <- tabulate(f3[y3 == cl], 5)
  expect_true(max(cnt) - min(cnt) <= 1, info = cl)
}
# The allocation is sklearn's: the same counts per fold and class.
for (cl in levels(y3))
  expect_identical(sort(as.vector(table(f3, y3)[, cl])),
                   sort(as.vector(table(ref$folds$three, y3)[, cl])), info = cl)
fr <- ecodive:::rf_folds(ref$reg$y, FALSE, 7L, 0L)
expect_identical(as.vector(table(fr)), c(8L, rep(7L, 6)))
expect_identical(ecodive:::rf_folds(ref$reg$y, FALSE, 7L, 0L), fr)
expect_false(identical(ecodive:::rf_folds(ref$reg$y, FALSE, 7L, 1L), fr))



# The same forest for any number of threads, sparse or dense. ====

xb <- ref$bin$x; yb <- ref$bin$y
base <- random_forest(xb, yb, n_trees = 40, importance = 'both',
                      max_features = c('sqrt', 0.3), min_samples_leaf = 1:2, seed = 3, cpus = 1)
for (cp in unique(c(2L, max_cpus, if (at_home()) c(3L, 8L))))
  expect_identical(
    random_forest(xb, yb, n_trees = 40, importance = 'both', max_features = c('sqrt', 0.3),
                  min_samples_leaf = 1:2, seed = 3, cpus = cp), base, info = cp )
expect_identical(
  random_forest(slam_or_matrix(xb), yb, n_trees = 40, importance = 'both',
                max_features = c('sqrt', 0.3), min_samples_leaf = 1:2, seed = 3, cpus = max_cpus),
  base )
expect_false(identical(random_forest(xb, yb, n_trees = 40, seed = 4, cpus = 1)$importance,
                       random_forest(xb, yb, n_trees = 40, seed = 3, cpus = 1)$importance))

# Regression: sparse and dense input sum in different orders, so tied splits
# (common in small nodes) can go either way, as in scikit-learn; the
# forests are equally good.
xr <- ref$reg$x; yr <- ref$reg$y
rr  <- random_forest(xr, yr, n_trees = 30, norm = 'none', seed = 1, cpus = max_cpus)
rrs <- random_forest(slam_or_matrix(xr), yr, n_trees = 30, norm = 'none', seed = 1, cpus = max_cpus)
expect_true(abs(rr$performance[['r_squared']] - rrs$performance[['r_squared']]) < 0.1)
expect_true(stats::cor(rr$predictions$predicted, rrs$predictions$predicted) > 0.95)
expect_identical(random_forest(xr, yr, n_trees = 30, norm = 'none', seed = 1, cpus = 1), rr)

# Each tree has its own stream: tree t is the same in a larger forest.
a <- rf_details(xb, yb, n_trees = 30, n_folds = 0, seed = 9, cpus = 1)$forests[[1]]
b <- rf_details(xb, yb, n_trees = 50, n_folds = 0, seed = 9, cpus = max_cpus)$forests[[1]]
expect_identical(a$inbag, b$inbag[, 1:30])
expect_identical(a$left, b$left[seq_along(a$left)])
expect_identical(a$threshold, b$threshold[seq_along(a$threshold)])



# Like scikit-learn's forests over many seeds (distributional). ====

if (at_home()) {
  acc <- imp <- NULL
  for (s in 0:19) {
    f <- random_forest(xb, yb, n_trees = 500, n_folds = 0, norm = 'none', seed = s, cpus = max_cpus)
    acc <- c(acc, f$performance[['oob_accuracy']])
    imp <- cbind(imp, f$importance$impurity[match(colnames(xb), f$importance$feature)])
  }
  sk  <- ref$dist$sklearn_oob_accuracy
  ski <- ref$dist$sklearn_importance
  expect_true(abs(mean(acc) - mean(sk)) < 3 * sqrt(stats::var(acc) / 20 + stats::var(sk) / 20))
  se <- sqrt(apply(imp, 1, stats::var) / 20 + apply(ski, 1, stats::var) / 20)
  expect_true(all(abs(rowMeans(imp) - rowMeans(ski)) <= 3 * pmax(se, 1e-12)))
  # randomForest's and ranger's OOB error (500 trees, 20 seeds) is similar.
  rf_err <- c(ref$dist$randomForest_oob_error, ref$dist$ranger_oob_error)
  expect_true(abs(mean(1 - acc) - mean(rf_err)) < 0.03)
  # Permutation importances rank like randomForest's and ranger's.
  pm <- random_forest(xb, yb, n_trees = 2000, n_folds = 0, norm = 'none',
                      importance = 'permutation', seed = 1, cpus = max_cpus)$importance
  pmv <- pm$permutation[match(colnames(xb), pm$feature)]
  expect_true(stats::cor(pmv, ref$dist$randomForest_perm, method = 'spearman') > 0.8)
  expect_true(stats::cor(pmv, ref$dist$ranger_perm, method = 'spearman') > 0.8)
}

# A planted signal: two informative features among 50 rank first.
sig <- local({
  n <- 80; p <- 50
  set.seed(42)
  x <- matrix(stats::rpois(n * p, 30), n, p,
              dimnames = list(paste0('S', 1:n), sprintf('F%02d', 1:p)))
  y <- factor(rep(c('a', 'b'), each = n / 2))
  x[y == 'b', 7]  <- x[y == 'b', 7] + 60
  x[y == 'b', 23] <- x[y == 'b', 23] + 60
  list(x = x, y = y)
})
fs <- random_forest(sig$x, sig$y, n_trees = 200, norm = 'none', importance = 'both',
                    n_folds = 4, cpus = max_cpus)
expect_identical(sort(fs$importance$feature[1:2]), c('F07', 'F23'))
expect_identical(sort(fs$importance$feature[order(-fs$importance$impurity)][1:2]), c('F07', 'F23'))
expect_true(fs$performance[['accuracy']] > 0.9)
expect_true(fs$performance[['auc']] > 0.95)



# The result. ====

res <- random_forest(ref$three$x, ref$three$y, n_trees = 30, cpus = max_cpus)
expect_identical(class(res), 'ecodive_random_forest')
expect_identical(attr(res, 'method'), 'random_forest')
expect_identical(attr(res, 'type'), 'classification')
expect_identical(attr(res, 'n_trees'), 30L)
expect_identical(attr(res, 'max_features'), 5L)      # floor(sqrt(30))
expect_identical(attr(res, 'min_samples_leaf'), 1L)
expect_identical(attr(res, 'criterion'), 'gini')
expect_identical(attr(res, 'n_folds'), 5L)
expect_identical(names(res), c('importance', 'predictions', 'performance', 'confusion', 'roc', 'forest'))
expect_identical(names(res$importance), c('feature', 'impurity', 'impurity_cv_mean', 'impurity_cv_sd'))
expect_true(abs(sum(res$importance$impurity) - 1) < 1e-12)
expect_false(is.unsorted(-res$importance$impurity))
expect_identical(names(res$predictions),
  c('sample', 'observed', 'predicted', 'fold', 'oob_predicted', 'prob_A', 'prob_B', 'prob_C'))
expect_identical(res$predictions$sample, rownames(ref$three$x))
expect_true(all(abs(rowSums(res$predictions[, 6:8]) - 1) < 1e-12))
expect_identical(names(res$performance), c(
  'accuracy', 'accuracy_fold_mean', 'accuracy_fold_sd', 'balanced_accuracy',
  'baseline_accuracy', 'accuracy_ratio', 'auc_macro', 'auc_micro', 'oob_accuracy' ))
expect_equal(res$performance[['baseline_accuracy']], 20 / 45)
expect_equal(res$performance[['accuracy']],
             mean(res$predictions$predicted == res$predictions$observed))
expect_identical(sum(res$confusion), 45L)
expect_identical(dimnames(res$confusion), list(observed = c('A', 'B', 'C'), predicted = c('A', 'B', 'C')))
expect_identical(as.vector(res$confusion),
                 as.vector(table(res$predictions$observed, res$predictions$predicted)))
expect_identical(unique(res$roc$class), c('A', 'B', 'C', 'micro', 'macro'))
# balanced accuracy is the mean recall.
cm <- res$confusion
expect_equal(res$performance[['balanced_accuracy']], mean(diag(cm) / rowSums(cm)))

# Regression.
rg <- random_forest(xr, yr, n_trees = 30, importance = 'permutation', cpus = max_cpus)
expect_identical(attr(rg, 'type'), 'regression')
expect_identical(attr(rg, 'max_features'), 30L)      # all features
expect_identical(attr(rg, 'criterion'), 'squared_error')
expect_identical(names(rg), c('importance', 'predictions', 'performance', 'forest'))
expect_identical(names(rg$importance), c('feature', 'impurity', 'impurity_cv_mean',
                                         'impurity_cv_sd', 'permutation', 'permutation_sd'))
expect_false(is.unsorted(-rg$importance$permutation))
expect_identical(names(rg$performance), c('mse', 'r_squared', 'mse_fold_mean', 'mse_fold_sd',
                                          'slope', 'intercept', 'r', 'oob_r_squared'))
pp <- rg$predictions
expect_equal(rg$performance[['mse']], mean((pp$predicted - pp$observed)^2))
expect_equal(rg$performance[['slope']], unname(stats::coef(stats::lm(predicted ~ observed, pp))[[2]]))
expect_equal(rg$performance[['r']], stats::cor(pp$observed, pp$predicted))
expect_true(rg$performance[['r_squared']] > 0.3)

# Binary: one ROC curve, of the second class.
rb2 <- random_forest(xb, yb, n_trees = 30, cpus = max_cpus)
expect_identical(unique(rb2$roc$class), 'B')
expect_true('auc' %in% names(rb2$performance))
expect_equal(rb2$performance[['auc']],
             ecodive:::rf_auc(rb2$roc$fpr, rb2$roc$tpr))

# Without CV: OOB results only.
nc <- random_forest(xb, yb, n_trees = 30, n_folds = 0, cpus = max_cpus)
expect_true(all(is.na(nc$predictions$predicted)) && all(is.na(nc$predictions$fold)))
expect_true(is.na(nc$performance[['accuracy']]) && !is.na(nc$performance[['oob_accuracy']]))
expect_true(all(is.na(nc$importance$impurity_cv_mean)))
expect_identical(sum(nc$confusion), sum(!is.na(nc$predictions$oob_predicted)))
expect_identical(attr(random_forest(xb, yb, n_trees = 5, n_folds = 1, cpus = 1), 'n_folds'), 0L)

# A grid: chosen by OOB error inside each training set.
gr <- random_forest(xb, yb, n_trees = 30, max_features = c(1, 'sqrt', 0.5),
                    min_samples_leaf = c(1, 3), n_folds = 3, cpus = max_cpus)
expect_identical(nrow(gr$tuning), 6L * 4L)
expect_identical(unique(gr$tuning$training_set), c('full', 'fold 1', 'fold 2', 'fold 3'))
for (s in unique(gr$tuning$training_set)) {
  t1 <- gr$tuning[gr$tuning$training_set == s, ]
  expect_identical(sum(t1$chosen), 1L, info = s)
  expect_identical(which(t1$chosen), which.min(t1$oob_error), info = s)
}
full <- gr$tuning[gr$tuning$training_set == 'full' & gr$tuning$chosen, ]
expect_identical(attr(gr, 'max_features'), full$max_features)
expect_identical(attr(gr, 'min_samples_leaf'), full$min_samples_leaf)

# keep_forest = FALSE; data columns; folds given.
expect_null(random_forest(xb, yb, n_trees = 5, keep_forest = FALSE, cpus = 1)$forest)
meta <- data.frame(row.names = rownames(xb), cls = yb, fold = rep(1:3, 20))
fd <- random_forest(xb, 'cls', meta, n_trees = 10, folds = 'fold', cpus = 1)
expect_identical(fd$predictions$fold, rep(1:3, 20))
expect_identical(attr(fd, 'n_folds'), 3L)
fd2 <- random_forest(xb[60:1, ], 'cls', meta, n_trees = 10, folds = 'fold', cpus = 1)
expect_identical(fd2$predictions$sample, rownames(xb)[60:1])
expect_identical(random_forest(t(xb), yb, margin = 2, n_trees = 10, cpus = 1),
                 random_forest(xb, yb, n_trees = 10, cpus = 1))

# class_weight, criterion and max_samples are accepted (parity above).
expect_identical(attr(random_forest(xb, yb, n_trees = 5, class_weight = 'balanced', cpus = 1), 'class_weight'), 'balanced')
expect_identical(attr(random_forest(xb, yb, n_trees = 5, class_weight = c(B = 2), cpus = 1), 'class_weight'), c(B = 2))
expect_identical(attr(random_forest(xb, yb, n_trees = 5, criterion = 'entropy', cpus = 1), 'criterion'), 'entropy')
expect_silent(random_forest(xb, yb, n_trees = 5, max_samples = 0.5, max_depth = 2, cpus = 1))
expect_silent(random_forest(xb, yb, n_trees = 25, norm = 'clr', pseudocount = 1, cpus = 1))
expect_silent(random_forest(xr, yr, n_trees = 25, criterion = 'mse', cpus = 1))



# predict(). ====

# The same as a walk of the stored nodes in base R.
walk <- function (fo, x) {
  tr <- fo$trees
  k  <- nrow(tr$value)
  out <- matrix(0, nrow(x), k)
  x32 <- x   # float32 rounding matters only at a threshold; none are this close here
  for (t in seq_len(length(tr$offsets) - 1)) {
    o <- tr$offsets[[t]]
    for (i in seq_len(nrow(x))) {
      node <- 0L
      while (tr$left[[o + node + 1]] != -1L) {
        f <- tr$feature[[o + node + 1]] + 1L
        node <- if (x32[i, f] <= tr$threshold[[o + node + 1]]) tr$left[[o + node + 1]] else tr$right[[o + node + 1]]
      }
      out[i, ] <- out[i, ] + tr$value[, o + node + 1]
    }
  }
  out / (length(tr$offsets) - 1)
}
f3 <- random_forest(ref$three$x, ref$three$y, n_trees = 20, norm = 'none', n_folds = 0, cpus = 1)
newx <- ref$three$x[c(3, 9, 30), ]
p3 <- predict(f3, newx, 'prob', cpus = 1)
expect_true(rel_ok(p3, walk(f3$forest, newx)))
expect_identical(dimnames(p3), list(rownames(newx), c('A', 'B', 'C')))
pc <- predict(f3, newx, cpus = 1)
expect_identical(levels(pc), c('A', 'B', 'C'))
expect_identical(as.vector(unclass(pc)), max.col(p3, ties.method = 'first'))
expect_identical(names(pc), rownames(newx))
# Matched by name: shuffled and extra columns give the same predictions.
sh <- cbind(newx[, 30:1], extra = 5)
expect_identical(predict(f3, sh, 'prob', cpus = 1), p3)
expect_identical(predict(f3, slam_or_matrix(newx), 'prob', cpus = 1), p3)
expect_identical(predict(f3, t(newx), 'prob', margin = 2, cpus = 1), p3)
# Missing features count as 0.
miss <- newx[, -(1:2)]
z <- newx; z[, 1:2] <- 0
expect_warning(pm <- predict(f3, miss, 'prob', cpus = 1), '2 of the forest.s 30 features are not')
expect_identical(pm, predict(f3, z, 'prob', cpus = 1))
expect_error(predict(f3, matrix(1, 2, 2, dimnames = list(NULL, c('x', 'y')))), 'none of the forest')
expect_error(predict(random_forest(xb, yb, n_trees = 5, keep_forest = FALSE, cpus = 1), xb), 'keep_forest')
# Normalized as the forest was: percent over the forest's features.
fp <- random_forest(ref$three$x, ref$three$y, n_trees = 20, n_folds = 0, cpus = 1)
nx <- newx / rowSums(newx)
expect_true(rel_ok(predict(fp, cbind(newx, extra = 100), 'prob', cpus = 1), walk(fp$forest, nx)))
# Regression returns numbers.
pr <- predict(rr, xr[1:4, ], cpus = 1)
expect_true(is.numeric(pr) && identical(names(pr), rownames(xr)[1:4]))
expect_true(rel_ok(pr, as.vector(walk(rr$forest, xr[1:4, ]))))



# Edge cases and errors. ====

y1 <- yb; y1[c(2, 5)] <- NA
expect_warning(r1 <- random_forest(xb, y1, n_trees = 10, cpus = 1), 'Left out 2 samples')
expect_identical(nrow(r1$predictions), 58L)
expect_identical(r1$predictions$sample, rownames(xb)[-c(2, 5)])

yf <- factor(as.character(yb), levels = c('A', 'B', 'Z'))
expect_warning(r2 <- random_forest(xb, yf, n_trees = 10, cpus = 1), "unused `y` level: 'Z'")
expect_identical(levels(r2$predictions$observed), c('A', 'B'))

expect_error(random_forest(xb, rep('A', 60), cpus = 1), 'at least two classes')
expect_error(random_forest(xb, cpus = 1), '`y` is required')
expect_error(random_forest(xb, yb[-1], cpus = 1), 'y')
ys <- yb; ys[1:57] <- 'A'; ys[58:60] <- 'B'
expect_warning(random_forest(xb, ys, n_trees = 5, cpus = 1), 'smallest class .* 3 samples')
ys2 <- rep(c('A', 'B', 'C'), c(2, 2, 2))
expect_error(random_forest(xb[1:6, ], ys2, n_trees = 5, n_folds = 3, cpus = 1), 'larger than every class')
expect_error(random_forest(xb, yb, n_folds = 61, cpus = 1), 'at most the number of samples')
expect_error(random_forest(xb, yb, max_features = 41, cpus = 1), 'at most the number of features')
expect_error(random_forest(xb, yb, max_features = 'cube', cpus = 1), 'max_features')
expect_error(random_forest(xb, yb, max_features = 2.5, cpus = 1), 'fraction below 1')
expect_error(random_forest(xb, yb, min_samples_leaf = 0, cpus = 1), 'min_samples_leaf')
expect_error(random_forest(xb, yb, n_trees = 0, cpus = 1), 'n_trees')
expect_error(random_forest(xb, yb, criterion = 'squared_error', cpus = 1), "'gini' or 'entropy'")
expect_error(random_forest(xr, yr, criterion = 'gini', cpus = 1), "'squared_error'")
expect_error(random_forest(xr, yr, class_weight = 'balanced', cpus = 1), 'classification only')
expect_error(random_forest(xb, yb, class_weight = c(Q = 2), cpus = 1), "lacks: 'Q'")
expect_error(random_forest(xb, yb, class_weight = 'equal', cpus = 1), 'class_weight')
expect_error(random_forest(xb, yb, max_samples = 1.5, cpus = 1), 'max_samples')
expect_error(random_forest(xb, yb, importance = 'gain', cpus = 1))
expect_error(random_forest(xb, yb, folds = rep(1, 60), cpus = 1), 'at least two folds')
expect_error(random_forest(xb, yb, folds = factor(rep(1:2, 30), levels = 1:3), cpus = 1), 'must have samples')
expect_error(random_forest(xb, yb, folds = ifelse(yb == 'A', 1, 2), cpus = 1), 'only one class')
expect_error(random_forest(xb, yb, keep_forest = NA, cpus = 1), 'keep_forest')
bad <- xb; bad[1, 1] <- NA
expect_error(random_forest(bad, yb, cpus = 1), 'NA')
expect_error(random_forest(xb, list(1), cpus = 1), 'y')

# p = 1; n = 2 without CV; a constant table (single-node trees).
p1 <- random_forest(xb[, 1, drop = FALSE], yb, n_trees = 10, cpus = 1)
expect_identical(p1$importance$impurity, 1)
expect_identical(attr(p1, 'max_features'), 1L)
two <- c(which(yb == 'A')[[1]], which(yb == 'B')[[1]])
expect_warning(n2 <- random_forest(xb[two, ], yb[two], n_trees = 1, n_folds = 0, seed = 2, cpus = 1),
               'never out of bag')
expect_identical(nrow(n2$predictions), 2L)
const <- matrix(1, 20, 5, dimnames = list(NULL, paste0('F', 1:5)))
cf <- random_forest(const, rep(c('a', 'b'), 10), n_trees = 10, norm = 'none', cpus = 1)
expect_identical(cf$importance$impurity, rep(0, 5))
f0 <- rf_details(const, rep(c('a', 'b'), 10), n_trees = 10, norm = 'none', n_folds = 0, cpus = 1)
expect_identical(f0$forests[[1]]$offsets, 0:10)
# Logical y; samples without names.
lg <- random_forest(unname(xb), as.vector(yb == 'B'), n_trees = 10, n_folds = 2, cpus = 1)
expect_identical(levels(lg$predictions$observed), c('FALSE', 'TRUE'))
expect_identical(lg$predictions$sample, as.character(1:60))



# Printing and plots. ====

out <- capture.output(print(res))
expect_true(any(grepl('Random forest classification: 45 samples x 30 features', out)))
expect_true(any(grepl('5-fold CV accuracy', out)))
expect_true(any(grepl('Top 10 features by impurity importance', out)))
expect_identical(withVisible(print(res))$visible, FALSE)
out <- capture.output(print(rg))
expect_true(any(grepl('Random forest regression', out)) && any(grepl('permutation importance', out)))
expect_true(any(grepl('Out-of-bag accuracy', capture.output(print(nc)))))

drawn <- function (x, ...) {
  v <- withVisible(plot(x, ...))
  identical(v$value, x) && !v$visible
}
had_seed <- exists('.Random.seed', globalenv())
if (had_seed) seed0 <- get('.Random.seed', globalenv())
grDevices::pdf(NULL)
par0 <- graphics::par(no.readonly = TRUE)
expect_true(drawn(res))
expect_true(drawn(res, 'roc'))
expect_true(drawn(res, 'confusion'))
expect_true(drawn(rb2, 'roc'))
expect_true(drawn(nc, 'roc'))
expect_true(drawn(rg))
expect_true(drawn(rg, importance = 'impurity', n_features = 3, main = 'x'))
expect_true(drawn(rg, 'predictions'))
expect_true(drawn(random_forest(xr, yr, n_trees = 10, n_folds = 0, cpus = 1), 'predictions'))
expect_true(drawn(gr, 'tuning'))
par1 <- graphics::par(no.readonly = TRUE)
drop <- c('usr', 'xaxp', 'yaxp')
expect_identical(par1[setdiff(names(par1), drop)], par0[setdiff(names(par0), drop)])
grDevices::dev.off()
expect_identical(exists('.Random.seed', globalenv()), had_seed)
if (had_seed) expect_identical(get('.Random.seed', globalenv()), seed0)
expect_error(plot(rg, 'roc'), 'classification forests')
expect_error(plot(res, 'predictions'), 'regression forests')
expect_error(plot(res, 'tuning'), 'grid')
expect_error(plot(res, importance = 'permutation'), 'no permutation importances')
expect_error(plot(res, n_features = 0), 'at least 1')

#})
