# Expected results for tests/tinytests/test-valencia.r, from VALENCIA's own
# Valencia.py. Run from the package root, in three steps:
#
#     Rscript data-raw/expected/valencia.r inputs <all_samples_csv>
#     python3 data-raw/expected/valencia.py <Valencia.py> <CST_centroids_012920.csv>
#     Rscript data-raw/expected/valencia.r
#
# The first step writes the inputs, in Valencia.py's format (sampleID,
# read_count, then taxa), to data-raw/expected/valencia/: <all_samples_csv> is
# Publication_materials/Data_and_metadata/all_samples_taxonomic_composition_data.csv
# from VALENCIA (commit 8559d45), from which 50 samples are taken. The second
# runs Valencia.py (VALENCIA commit 8559d45; Python 3, pandas 2.2.3, numpy
# 1.26.4) on each input and writes its outputs there, with every similarity
# as a hexadecimal float. The third folds inputs and outputs into
# tests/tinytests/expected/valencia.rds.
#
# VALENCIA, its centroids and its publication data are MIT-licensed
# (Copyright (c) 2023 ravel-lab; see data-raw/expected/valencia/LICENSE).
#
# Cases:
#   synthetic:   195 samples drawn as Dirichlet-multinomial counts around
#                each of the 13 centroids (15 each), with two taxa that the
#                centroids lack added to some, and some read counts above
#                the taxa's sum (unclassified reads); plus an all-zero
#                sample, a sample of unmatched taxa only (scores 0 everywhere:
#                a tie, which goes to I-A), and mixtures of two centroids.
#                Taxa are in random order, and centroid taxa absent from
#                every sample are left out.
#   ties:        30 of those samples against centroids in which I-B is a
#                copy of I-A and III-A a copy of III-B, so that those pairs
#                tie exactly (the first wins).
#   publication: 50 samples of VALENCIA's publication data, about 4 per
#                sub-CST as VALENCIA assigned them.

dir <- file.path('data-raw', 'expected', 'valencia')
args <- commandArgs(TRUE)

load(file.path('data', 'valencia_centroids.rda'))
cents <- valencia_centroids
sub   <- rownames(cents)


write_input <- function (counts, read_count, file) {
  df <- data.frame(
    sampleID = rownames(counts), read_count = read_count,
    as.data.frame(counts, optional = TRUE), check.names = FALSE )
  utils::write.csv(df, file.path(dir, file), row.names = FALSE, quote = FALSE)
}

read_input <- function (file) {
  df <- utils::read.csv(file.path(dir, file), check.names = FALSE, stringsAsFactors = FALSE)
  m  <- as.matrix(df[, -(1:2), drop = FALSE])
  storage.mode(m) <- 'integer'
  rownames(m) <- df$sampleID
  list(counts = m, totals = stats::setNames(df$read_count, df$sampleID))
}



if (length(args) && args[[1]] == 'inputs') {

  dir.create(dir, showWarnings = FALSE)

  rdirichlet <- function (alpha) {
    g <- stats::rgamma(length(alpha), alpha)
    g / sum(g)
  }

  #--------------------------------------------------
  # synthetic
  #--------------------------------------------------
  set.seed(52)
  taxa <- c(colnames(cents), 'Unknown_taxon_A', 'g_NotInCentroids')
  rows <- list()
  for (k in seq_along(sub)) for (r in 1:15) {
    p     <- rdirichlet(cents[k, ] * stats::runif(1, 20, 200))
    depth <- sample(2000:20000, 1)
    x     <- c(stats::rmultinom(1, depth, p), 0L, 0L)
    if (stats::runif(1) < 0.3) x[length(x) - 1:0] <- sample(0:300, 2)
    rows[[sprintf('%s_%02d', sub[k], r)]] <- x
  }
  rows[['zero']]      <- integer(length(taxa))
  rows[['unmatched']] <- c(integer(ncol(cents)), 120L, 30L)
  mix <- function (a, b, w) round(c(cents[a, ] * w + cents[b, ] * (1 - w), 0, 0) * 10000)
  rows[['mix_IA_IIIA']]  <- mix('I-A', 'III-A', 0.5)
  rows[['mix_IVB_IIIB']] <- mix('IV-B', 'III-B', 0.4)
  rows[['mix_V_II']]     <- mix('V', 'II', 0.55)

  counts <- do.call(rbind, rows)
  colnames(counts) <- taxa
  storage.mode(counts) <- 'integer'
  counts <- counts[, colSums(counts) > 0]
  counts <- counts[, sample(ncol(counts))]

  rs    <- rowSums(counts)
  extra <- ifelse(stats::runif(nrow(counts)) < 0.25 & rs > 0, sample(1:2000, nrow(counts), TRUE), 0L)
  write_input(counts, rs + extra, 'synthetic-input.csv')

  #--------------------------------------------------
  # ties: centroids with two duplicated pairs.
  #--------------------------------------------------
  tc <- cents
  tc['I-B', ]   <- tc['I-A', ]
  tc['III-A', ] <- tc['III-B', ]
  df <- data.frame(sub_CST = rownames(tc), tc, check.names = FALSE)
  utils::write.csv(df, file.path(dir, 'ties-centroids.csv'), row.names = FALSE, quote = FALSE)

  pick <- grep('^(I-A|I-B|III-A|III-B)_', rownames(counts))[c(1:5, 16:20, 31:35, 46:50)]
  pick <- c(pick, match(c('unmatched', 'mix_IA_IIIA', 'mix_IVB_IIIB', 'mix_V_II'), rownames(counts)))
  pick <- c(pick, grep('^(II|V)_', rownames(counts))[1:6])
  tcnt <- counts[pick, ]
  tcnt <- tcnt[, colSums(tcnt) > 0]
  write_input(tcnt, rowSums(tcnt), 'ties-input.csv')

  #--------------------------------------------------
  # publication: 50 samples of VALENCIA's data.
  #--------------------------------------------------
  pub  <- utils::read.csv(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
  # Its 199 taxa: those of the centroids, but for g_Fastidiosipila, which the
  # centroids call BVAB2/BVAB3 (so it is a taxon they lack, here).
  taxa <- setdiff(names(pub)[-(1:7)], paste0(sub, '_sim'))
  stopifnot(length(taxa) == 199)
  set.seed(7)
  per  <- c(rep(4L, 11), 3L, 3L)  # 50 samples
  pick <- unlist(lapply(seq_along(sub), function (k) {
    i <- which(pub$Val_subCST == sub[k])
    i[sample.int(length(i), per[k])]
  }))
  pc <- as.matrix(pub[pick, taxa])
  storage.mode(pc) <- 'integer'
  rownames(pc) <- sprintf('P%05d', pick)
  stopifnot(all(rowSums(pc) == pub$total_reads[pick]))
  pc <- pc[, colSums(pc) > 0]
  write_input(pc, pub$total_reads[pick], 'publication-input.csv')

  # VALENCIA's published assignments and (rounded) similarities, to check
  # Valencia.py's outputs against in the third step.
  pubout <- pub[pick, c('Val_subCST', paste0(sub, '_sim'))]
  pubout <- data.frame(sampleID = rownames(pc), pubout, check.names = FALSE)
  utils::write.csv(pubout, file.path(dir, 'publication-published.csv'), row.names = FALSE, quote = FALSE)

  quit(save = 'no')
}



#--------------------------------------------------
# Fold the outputs in.
#--------------------------------------------------

read_output <- function (file) {
  df <- utils::read.csv(file.path(dir, file), check.names = FALSE, stringsAsFactors = FALSE,
                        colClasses = 'character', na.strings = '')
  sc <- sapply(paste0(sub, '_sim'), function (k) as.numeric(df[[k]]))
  dimnames(sc) <- list(df$sampleID, sub)
  list(
    scores = sc,
    subCST = df$subCST,
    CST    = df$CST,
    score  = as.numeric(df$score) )
}

ref <- list()
for (case in c('synthetic', 'ties', 'publication')) {
  ref[[case]] <- c(read_input(paste0(case, '-input.csv')), read_output(paste0(case, '-output.csv')))
  stopifnot(identical(rownames(ref[[case]]$counts), rownames(ref[[case]]$scores)))
}

tc <- utils::read.csv(file.path(dir, 'ties-centroids.csv'), check.names = FALSE)
ref$ties$centroids <- as.matrix(tc[, -1])
rownames(ref$ties$centroids) <- tc$sub_CST

# Valencia.py agrees with VALENCIA's published assignments (50 of 50). Its
# similarities differ from the published ones (6 decimals) by up to 0.07:
# those were computed with g_Fastidiosipila named BVAB2/BVAB3, as in the
# centroids (renamed so, ecodive's valencia() gives the published values to
# 6.5e-8).
pubout <- utils::read.csv(file.path(dir, 'publication-published.csv'), check.names = FALSE)
stopifnot(identical(pubout$Val_subCST, ref$publication$subCST))

saveRDS(ref, file.path('tests', 'tinytests', 'expected', 'valencia.rds'), version = 2)
