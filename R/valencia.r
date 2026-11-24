# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# valencia() ports VALENCIA's nearest-centroid classification (Valencia.py,
# commit 8559d45 of github.com/ravel-lab/VALENCIA; see src/valencia.c), and
# ships its reference centroids (data/valencia_centroids.rda, built by
# data-raw/valencia_centroids.r). The reference scores' means and standard
# deviations and the sub-CST colours in VALENCIA_SUBCST below are taken from
# Valencia.py's plot. VALENCIA is under this licence:
#
# MIT License
#
# Copyright (c) 2023 ravel-lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


#' VALENCIA community state types
#'
#' Assigns vaginal microbiome samples to community state types (CSTs) with
#' VALENCIA (France et al. 2020): each sample goes to the most similar of 13
#' reference centroids (sub-CSTs), by Yue and Clayton's (2005) similarity of
#' relative abundances.
#'
#' @param counts  A count table (or relative abundances) whose features are
#'        named as VALENCIA names its taxa (see "Taxon names"). Samples are
#'        rows unless `margin = 2`. A named vector is one sample.
#'
#' @param centroids  The reference centroids: a matrix of relative
#'        abundances, one row per sub-CST (named in the row names) and one
#'        column per taxon (named in the column names). Default:
#'        [valencia_centroids], VALENCIA's 13 published centroids. A data
#'        frame read from VALENCIA's centroids file, whose `sub_CST` column
#'        names the rows, is also accepted.
#'
#' @param totals  Each sample's total read count, by which its counts are
#'        divided, as VALENCIA's `read_count` column: a vector with one value
#'        per sample (matched by name, or else by position), normally at
#'        least the sample's sum (it warns otherwise). Default `NULL` uses
#'        each sample's sum. Give
#'        `totals` when reads were left out of `counts` (such as unclassified
#'        ones) and VALENCIA would have counted them.
#'
#' @param margin  If your samples are in the table's rows, set to `1L`
#'        (the default). If they are in columns, set to `2L`.
#'
#' @param cpus  How many parallel processing threads to use. Results are
#'        the same for any number. Default: all available cores.
#'
#' @details
#' For each sample, the relative abundances \eqn{p} are its counts divided
#' by its total. Its similarity to each centroid \eqn{q} is Yue and
#' Clayton's \eqn{\theta}, over the union of the sample's and the centroid's
#' taxa (a taxon missing from either is 0 there):
#' \deqn{\theta = \frac{\sum p q}{\sum (p - q)^2 + \sum p q}}
#' It is 1 when \eqn{p = q}, and 0 when they share no taxa. Taxa that the
#' centroids don't name still count, through \eqn{\sum p^2}, and lower
#' every similarity. The sample's sub-CST (`type`) is the centroid of
#' largest similarity (the first in the centroids' order, on ties), and its
#' `score` is that similarity. Sub-CSTs are then grouped into CSTs (`cst`):
#' I-A and I-B into I, III-A and III-B into III, and IV-C0 to IV-C4 into
#' IV-C.
#'
#' **The community state types** (France et al. 2020):
#' * **I** (I-A, I-B): dominated by *Lactobacillus crispatus*; almost
#'   entirely in I-A, less so in I-B.
#' * **II**: dominated by *L. gasseri*.
#' * **III** (III-A, III-B): dominated by *L. iners*.
#' * **IV-A**: high to moderate BVAB1 and *Gardnerella vaginalis*.
#' * **IV-B**: high to moderate *G. vaginalis* and *Atopobium vaginae*.
#' * **IV-C**: few *Lactobacillus*, *G. vaginalis* or BVAB1: IV-C0, an even
#'   community with *Prevotella*; IV-C1 to IV-C4, dominated by
#'   *Streptococcus*, *Enterococcus*, *Bifidobacterium* or *Staphylococcus*.
#' * **V**: dominated by *L. jensenii*.
#'
#' **Scores.** VALENCIA's authors compare each sub-CST's scores with those of
#' their reference samples (the means and standard deviations in the
#' `types` attribute, and in `plot(x, "scores")`): a sample scoring far
#' below its sub-CST's reference fits it poorly.
#'
#' @section Taxon names:
#' Features are matched to the centroids' taxa by their exact names, so
#' name them as VALENCIA does (see `colnames(valencia_centroids)`):
#' * species of *Lactobacillus*, *Gardnerella*, *Prevotella*, *Atopobium* and
#'   *Sneathia* as `Genus_species`, such as `"Lactobacillus_crispatus"`,
#'   `"Gardnerella_vaginalis"`;
#' * other taxa summarized to genus, or a higher rank, with a prefix:
#'   `"g_Streptococcus"`, `"f_Enterobacteriaceae"`, `"o_Micrococcales"`;
#' * a few names of their own: `"BVAB1"`, `"BVAB2/BVAB3"`, `"BVAB_TM7"`.
#'
#' Sum features that share a name first (duplicated names are an error).
#' `valencia()` warns when `counts` has none of the taxa that dominate the
#' centroids (*L. crispatus*, *L. iners*, *L. gasseri*, *L. jensenii*,
#' BVAB1, *G. vaginalis*, and the genera *Streptococcus*, *Enterococcus*,
#' *Bifidobacterium* and *Staphylococcus*), since a misnamed one leaves
#' samples in the wrong type; the `matched` column shows how much of each
#' sample was matched. *Prevotella* species are often named differently:
#' VALENCIA's authors suggest merging them all into `"g_Prevotella"`, with
#' centroids merged the same way (`centroids`). VALENCIA's
#' `convert_qiime.py`, which renames QIIME 2 taxonomies, has no equivalent
#' here.
#'
#' @section Differences from Valencia.py:
#' * Valencia.py divides counts by a `read_count` column, which can exceed
#'   the sum of the taxa; `valencia()` divides by each sample's sum unless
#'   `totals` gives the read counts.
#' * A sample of total 0 has `NaN` scores in both, and a type of `NA` here,
#'   with a warning.
#' * A sample that shares no taxa with any centroid scores 0 against all of
#'   them. Valencia.py then assigns it to the first, I-A, with score 0;
#'   `valencia()` does the same, with a warning.
#' * The similarities are computed over each sample's non-zero taxa only,
#'   as \eqn{\sum pq / (\sum p^2 + \sum q^2 - \sum pq)}, which equals
#'   Valencia.py's form to rounding.
#'
#' @return A data frame of class
#'   `c("ecodive_valencia", "ecodive_types", "data.frame")`, one row per
#'   sample (see [ecodive_types]), with columns:
#'   * `sample`: the sample's name.
#'   * `type`: its sub-CST, a factor whose levels are the centroids' names
#'     (`"I-A"` to `"V"`), in the centroids' order.
#'   * `cst`: its CST, a factor (levels `"I"`, `"II"`, `"III"`, `"IV-A"`,
#'     `"IV-B"`, `"IV-C"`, `"V"`).
#'   * `score`: its similarity to its sub-CST's centroid, from 0 to 1.
#'   * `matched`: the proportion of its reads in taxa that the centroids
#'     name.
#'
#'   Its attributes:
#'   * `types`: a data frame, one row per sub-CST, with `type`, `n_samples`,
#'     `top_feature` (the centroid's most abundant taxon), `cst`,
#'     `score_mean` (the mean score of the samples assigned to it), and
#'     `reference_score_mean` and `reference_score_sd`, VALENCIA's mean and
#'     standard deviation of the scores of its reference samples (`NA` for
#'     centroids of other names).
#'   * `profiles`: the centroids, a sub-CSTs \eqn{\times} taxa matrix.
#'   * `scores`: a samples \eqn{\times} sub-CSTs matrix of every
#'     similarity.
#'   * `method`: `"valencia"`.
#'
#'   There is no `selection`: the number of types is fixed by the
#'   centroids. `plot()` draws the centroids' profiles or the scores (see
#'   [plot.ecodive_types()]); its title is "VALENCIA Community State Types".
#'
#' @section Same names in other packages:
#' None known.
#'
#' @references
#' France, M. T., Ma, B., Gajer, P., Brown, S., Humphrys, M. S., Holm, J. B.,
#' Waetjen, L. E., Brotman, R. M., & Ravel, J. (2020). VALENCIA: a nearest
#' centroid classification method for vaginal microbial communities based on
#' composition. *Microbiome*, 8, 166. \doi{10.1186/s40168-020-00934-6}
#'
#' Yue, J. C., & Clayton, M. K. (2005). A similarity measure based on species
#' proportions. *Communications in Statistics - Theory and Methods*, 34,
#' 2123-2131. \doi{10.1080/STA-200066418}
#'
#' Ravel, J., Gajer, P., Abdo, Z., et al. (2011). Vaginal microbiome of
#' reproductive-age women. *Proceedings of the National Academy of
#' Sciences*, 108, 4680-4687. \doi{10.1073/pnas.1002611107}
#'
#' @family community types
#' @seealso [valencia_centroids], [ecodive_types], [plot.ecodive_types()]
#' @export
#' @examples
#'     # Read counts of four samples, with VALENCIA's taxon names.
#'     counts <- rbind(
#'       S1 = c(Lactobacillus_crispatus = 9500, Lactobacillus_iners = 300, g_Streptococcus = 200),
#'       S2 = c(Lactobacillus_crispatus = 400, Lactobacillus_iners = 7000, g_Streptococcus = 100),
#'       S3 = c(Lactobacillus_crispatus = 0, Lactobacillus_iners = 50, g_Streptococcus = 9000),
#'       S4 = c(Lactobacillus_crispatus = 3000, Lactobacillus_iners = 3000, g_Streptococcus = 3000) )
#'     counts <- cbind(counts, Gardnerella_vaginalis = c(0, 500, 300, 2000))
#'
#'     # It warns that some of the centroids' key taxa are missing.
#'     vl <- suppressWarnings(valencia(counts))
#'     vl
#'     attr(vl, 'scores')
#'
#'     plot(vl, 'scores')
valencia <- function (
    counts, centroids = valencia_centroids, totals = NULL, margin = 1L,
    cpus = n_cpus() ) {

  if (missing(centroids))
    centroids <- getExportedValue('ecodive', 'valencia_centroids')

  # A named vector is one sample.
  if (is.atomic(counts) && is.null(dim(counts)) && !is.null(names(counts))) {
    counts <- matrix(counts, nrow = 1L, dimnames = list(NULL, names(counts)))
    margin <- 1L
  }

  validate_counts()
  validate_margin()
  validate_cpus()

  if (is.logical(counts) && is.matrix(counts))
    storage.mode(counts) <- 'double'

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  remove('vals')

  smp     <- count_samples(counts, margin)
  n       <- smp$n
  samples <- if (is.null(smp$labels)) as.character(seq_len(n)) else smp$labels

  features <- count_samples(counts, 3L - margin)$labels
  if (is.null(features) || anyNA(features) || any(!nzchar(features)))
    stop(
      "`counts` must name its features, as VALENCIA names its taxa ",
      "(see ?valencia).", call. = FALSE )
  if (anyDuplicated(features))
    stop(
      '`counts` has duplicated feature names: ',
      name_list(unique(features[duplicated(features)])),
      '. Sum the features that share a name first.', call. = FALSE )

  cents <- valencia_centroid_matrix(centroids)
  sub   <- rownames(cents)
  K     <- nrow(cents)

  tot <- sample_totals(counts, margin)
  if (!is.null(totals)) totals <- valencia_totals(totals, tot, smp$labels, n)


  #--------------------------------------------------
  # Similarities, sub-CSTs and scores (src/valencia.c).
  #--------------------------------------------------
  col <- match(features, colnames(cents), nomatch = 0L)
  res <- .Call(C_valencia, counts, margin, col, cents, totals, cpus)

  scores <- res[[1]]
  best   <- res[[2]]
  dimnames(scores) <- list(samples, sub)

  score <- rep(NaN, n)
  ok    <- which(!is.na(best))
  score[ok] <- scores[cbind(ok, best[ok])]

  cst_of <- valencia_cst(sub)
  type   <- factor(sub[best], levels = sub)
  cst    <- factor(cst_of[best], levels = unique(cst_of))


  #--------------------------------------------------
  # Warnings.
  #--------------------------------------------------
  empty <- is.na(res[[3]])
  if (any(empty))
    warning(
      sum(empty), ' sample', if (sum(empty) == 1) ' has' else 's have',
      ' no reads, so ', if (sum(empty) == 1) 'its' else 'their',
      ' scores are NaN and type NA: ', name_list(samples[empty]), '.', call. = FALSE )

  none <- !empty & !is.na(best) & score == 0
  if (any(none))
    warning(
      sum(none), ' sample', if (sum(none) == 1) ' shares' else 's share',
      ' no taxa with the centroids, so ',
      if (sum(none) == 1) 'it scores' else 'they score',
      ' 0 against every one and, as in VALENCIA, ',
      if (sum(none) == 1) 'is' else 'are', ' assigned to the first, ', sub[[1]],
      ': ', name_list(samples[none]), '. Check the taxon names (see ?valencia).',
      call. = FALSE )

  key  <- valencia_key_taxa(cents)
  miss <- setdiff(key, features)
  if (length(miss))
    warning(
      '`counts` has no feature named ', name_list(miss, max = Inf),
      if (length(miss) == 1) ', which dominates' else ', which dominate',
      ' some of the centroids. Check the taxon names (see ?valencia).', call. = FALSE )


  #--------------------------------------------------
  # The result.
  #--------------------------------------------------
  ref <- VALENCIA_SUBCST[match(sub, VALENCIA_SUBCST$id), , drop = FALSE]

  types <- data.frame(
    stringsAsFactors = FALSE,
    type        = factor(sub, levels = sub),
    n_samples   = tabulate(best[ok], K),
    top_feature = types_top_feature(cents, K),
    cst         = factor(cst_of, levels = levels(cst)),
    score_mean  = vapply(seq_len(K), function (k) {
      s <- score[which(best == k)]
      if (length(s)) mean(s) else NA_real_ }, 0),
    reference_score_mean = ref$mean,
    reference_score_sd   = ref$sd )
  rownames(types) <- NULL

  df <- data.frame(
    stringsAsFactors = FALSE,
    sample  = samples,
    type    = type,
    cst     = cst,
    score   = score,
    matched = res[[3]] )

  new_ecodive_types(df, 'valencia', types, cents, NULL, scores = scores)
}



#' VALENCIA's reference centroids
#'
#' The 13 reference centroids of VALENCIA (France et al. 2020), against which
#' [valencia()] classifies vaginal microbiome samples: the mean relative
#' abundances of 199 taxa in the reference samples of each sub-CST.
#'
#' @format A numeric matrix of 13 sub-CSTs (rows: `"I-A"`, `"I-B"`, `"II"`,
#'   `"III-A"`, `"III-B"`, `"IV-A"`, `"IV-B"`, `"IV-C0"` to `"IV-C4"`, `"V"`)
#'   \eqn{\times} 199 taxa (columns, named as VALENCIA names them; see
#'   [valencia()]). Each row sums to 1, to rounding.
#'
#' @details
#' The centroids were computed from 13,160 taxonomic profiles of 1,975
#' North American women, from 16S rRNA gene amplicon sequencing; only these
#' means are included, no individual data. A revised set of centroids (for
#' instance with all *Prevotella* merged into `"g_Prevotella"`, and the
#' samples' taxa merged the same way) can be given to [valencia()] as its
#' `centroids`.
#'
#' @source `CST_centroids_012920.csv` from VALENCIA, commit 8559d45 of
#'   <https://github.com/ravel-lab/VALENCIA>, built by
#'   `data-raw/valencia_centroids.r` in ecodive's source repository.
#'   Released under the MIT licence:
#'
#'   Copyright (c) 2023 ravel-lab
#'
#'   Permission is hereby granted, free of charge, to any person obtaining a
#'   copy of this software and associated documentation files (the
#'   "Software"), to deal in the Software without restriction, including
#'   without limitation the rights to use, copy, modify, merge, publish,
#'   distribute, sublicense, and/or sell copies of the Software, and to
#'   permit persons to whom the Software is furnished to do so, subject to
#'   the following conditions:
#'
#'   The above copyright notice and this permission notice shall be included
#'   in all copies or substantial portions of the Software.
#'
#'   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#'   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#'   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
#'   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#'   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#'   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
#'   USE OR OTHER DEALINGS IN THE SOFTWARE.
#'
#' @references
#' France, M. T., Ma, B., Gajer, P., Brown, S., Humphrys, M. S., Holm, J. B.,
#' Waetjen, L. E., Brotman, R. M., & Ravel, J. (2020). VALENCIA: a nearest
#' centroid classification method for vaginal microbial communities based on
#' composition. *Microbiome*, 8, 166. \doi{10.1186/s40168-020-00934-6}
#'
#' @seealso [valencia()]
#' @examples
#'     dim(valencia_centroids)
#'     round(valencia_centroids[, 1:5], 3)
"valencia_centroids"



# VALENCIA's sub-CSTs (MIT; from Valencia.py's plot): the mean and standard
# deviation of the reference samples' scores against their own sub-CST
# (`ref_ave`, `ref_std`), and the colours of `CST_color_scheme`.
VALENCIA_SUBCST <- data.frame(
  stringsAsFactors = FALSE,
  id   = c('I-A', 'I-B', 'II', 'III-A', 'III-B', 'IV-A', 'IV-B',
           'IV-C0', 'IV-C1', 'IV-C2', 'IV-C3', 'IV-C4', 'V'),
  mean = c(0.995678, 0.858443, 0.811734, 0.972983, 0.810466, 0.718100,
           0.659592, 0.321432, 0.745991, 0.695203, 0.758735, 0.681502, 0.734973),
  sd   = c(0.007615, 0.151048, 0.161987, 0.058698, 0.133990, 0.150829,
           0.162577, 0.137558, 0.231891, 0.246322, 0.213653, 0.186462, 0.140276),
  colour = c('#ff6868', '#ffd4da', '#b4ff68', '#ffbc6b', '#e4a67b', '#c1adec',
             '#91a8ed', '#989898', '#ffc0cb', '#a8e5e5', '#9acc9a', '#800080',
             '#ffff71') )


# The higher-order CST of each sub-CST name, as Valencia.py maps them; other
# names are kept.
valencia_cst <- function (sub) {
  map <- c(
    'I-A' = 'I', 'I-B' = 'I', 'III-A' = 'III', 'III-B' = 'III',
    'IV-C0' = 'IV-C', 'IV-C1' = 'IV-C', 'IV-C2' = 'IV-C', 'IV-C3' = 'IV-C',
    'IV-C4' = 'IV-C' )
  ifelse(sub %in% names(map), map[sub], sub)
}


# The colour of each sub-CST name: VALENCIA's, or else from plot_palette().
valencia_colours <- function (sub) {
  col   <- VALENCIA_SUBCST$colour[match(sub, VALENCIA_SUBCST$id)]
  other <- which(is.na(col))
  if (length(other)) col[other] <- plot_palette(length(other))
  col
}


# The taxa that dominate the centroids: each centroid's most abundant taxon,
# where it makes up at least a quarter (for VALENCIA's centroids, all but
# the even IV-C0).
valencia_key_taxa <- function (cents) {
  top <- apply(cents, 1L, which.max)
  big <- cents[cbind(seq_len(nrow(cents)), top)] >= 0.25
  unique(colnames(cents)[top[big]])
}


# `centroids` as a double matrix with named, unique rows and columns.
valencia_centroid_matrix <- function (centroids) {

  if (is.data.frame(centroids)) {
    if ('sub_CST' %in% names(centroids)) {
      rn <- as.character(centroids[['sub_CST']])
      centroids <- centroids[, names(centroids) != 'sub_CST', drop = FALSE]
      rownames(centroids) <- NULL
      centroids <- as.matrix(centroids)
      rownames(centroids) <- rn
    }
    else {
      centroids <- as.matrix(centroids)
    }
  }

  if (!is.matrix(centroids) || !(is.numeric(centroids) || is.logical(centroids)) ||
      nrow(centroids) == 0 || ncol(centroids) == 0)
    stop('`centroids` must be a numeric matrix: one row per sub-CST, one column per taxon.', call. = FALSE)

  rn <- rownames(centroids)
  cn <- colnames(centroids)
  if (is.null(rn) || anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn))
    stop('`centroids` must have unique row names: the sub-CSTs.', call. = FALSE)
  if (is.null(cn) || anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn))
    stop('`centroids` must have unique column names: the taxa.', call. = FALSE)

  storage.mode(centroids) <- 'double'
  if (!all(is.finite(centroids)) || any(centroids < 0))
    stop('`centroids` must hold finite, non-negative values.', call. = FALSE)
  if (any(rowSums(centroids) <= 0))
    stop(
      '`centroids` has centroids of zeros: ',
      name_list(rn[rowSums(centroids) <= 0]), '.', call. = FALSE )

  centroids
}


# `totals`, one per sample, in the samples' order: matched by name when both
# are named, else by position; each at least the sample's sum `tot`.
valencia_totals <- function (totals, tot, labels, n) {

  if (!is.numeric(totals) || !is.null(dim(totals)) || anyNA(totals) || any(!is.finite(totals)))
    stop('`totals` must be a numeric vector of finite values, one per sample.', call. = FALSE)
  if (length(totals) != n)
    stop('`totals` must have one value per sample (', n, ').', call. = FALSE)

  if (!is.null(names(totals)) && !is.null(labels)) {
    i <- match_rows(labels, n, names(totals), length(totals), ref = '`counts`', what = '`totals`', unit = 'value')
    totals <- totals[i]
  }
  totals <- as.vector(as.double(totals))

  nm <- function (i) name_list(if (is.null(labels)) which(i) else labels[i])

  zero <- totals <= 0 & tot > 0
  if (any(totals < 0) || any(zero))
    stop('`totals` must be positive for samples with counts: ', nm(totals < 0 | zero), '.', call. = FALSE)

  # Below the sum, the proportions add up to more than 1; allowed (as in
  # VALENCIA), with tolerance for rounded relative abundances.
  low <- totals < tot * (1 - 1e-6)
  if (any(low))
    warning(
      '`totals` is smaller than the sum of counts for ', sum(low),
      ' sample', if (sum(low) == 1) '' else 's', ': ', nm(low), '.', call. = FALSE )

  totals
}


# valencia()'s default `centroids` names the lazy-loaded dataset, which the
# body fetches with getExportedValue() when `centroids` is missing.
utils::globalVariables('valencia_centroids')
