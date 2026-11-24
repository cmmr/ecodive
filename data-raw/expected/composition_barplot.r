# Expected results for tests/tinytests/test-composition_barplot.r, from
# microshades 1.13 (github.com/KarstensLab/microshades, MIT, Copyright (c)
# 2021 Erin Dahl): which members of each selected group create_color_dfs()
# keeps, in which order of shade, and that the group's "Other" takes the
# lightest shade. Run from the package root:
#
#     Rscript data-raw/expected/composition_barplot.r <microshades R directory>
#
# with <microshades R directory> holding microshades' R/color_mapping_functions.R
# and R/palettes.R (version 1.13, from the GitHub repository), and with dplyr,
# tidyr, forcats and rlang installed (run with R 4.3.3, dplyr 1.2.1, tidyr
# 1.3.2, forcats 1.0.1, rlang 1.3.0). Only those two files are sourced: prep_mdf()'s
# phyloseq step is replaced by building its melted data frame (one row per
# sample and genus, Abundance = the genus's proportion in the sample)
# directly. Nothing of microshades is copied into ecodive; this script only
# runs it.
#
# microshades leaves the choice of groups to the user (`selected_groups`):
# the groups given to it here are the `n_feature_groups` with the largest
# summed proportion, computed below in base R, which composition_barplot()
# must choose too (it ranks by the mean, which ranks the same).
#
# Cases (all proportions continuous, so that no two features or groups tie):
#   main:     30 samples x 65 genera in 7 phyla plus 5 genera of unknown
#             phylum ("Unassigned"), with zeros; 4 groups, 4 members each.
#   two:      the same table, 3 groups of 2 members each.
#   small:    20 samples x 14 genera in 3 phyla, one phylum with only 2
#             genera (so it has no "Other"); 3 groups, 4 members each.

args <- commandArgs(TRUE)
stopifnot(length(args) == 1)
ms_dir <- args[[1]]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(rlang)
})
source(file.path(ms_dir, 'palettes.R'))
source(file.path(ms_dir, 'color_mapping_functions.R'))


set.seed(47)

make_table <- function (n, phyla, sizes, n_unassigned) {
  p      <- sum(sizes) + n_unassigned
  scale  <- exp(stats::rnorm(p, sd = 1.5))
  counts <- matrix(0, n, p)
  for (j in seq_len(p))
    counts[, j] <- stats::rpois(n, scale[[j]] * stats::rgamma(n, 2, 2) * 40)
  dimnames(counts) <- list(paste0('S', seq_len(n)), paste0('Genus', seq_len(p)))
  phylum <- c(rep(phyla, sizes), rep(NA, n_unassigned))
  perm   <- sample(p)
  list(counts = counts[, perm], phylum = phylum[perm])
}

# microshades' melted data frame, and the result of create_color_dfs().
run_microshades <- function (counts, phylum, n_groups, n_members) {

  prop <- counts / rowSums(counts)
  ph   <- ifelse(is.na(phylum), 'Unassigned', phylum)

  mdf <- data.frame(
    Sample    = rep(rownames(prop), times = ncol(prop)),
    Genus     = rep(colnames(prop), each  = nrow(prop)),
    Phylum    = rep(ph,             each  = nrow(prop)),
    Abundance = as.vector(prop),
    stringsAsFactors = FALSE )
  rownames(mdf) <- paste0('r', seq_len(nrow(mdf)))

  sums     <- tapply(mdf$Abundance, mdf$Phylum, sum)
  selected <- names(sort(sums, decreasing = TRUE))[seq_len(n_groups)]

  res <- create_color_dfs(
    mdf, selected_groups = selected, top_n_subgroups = n_members,
    group_level = 'Phylum', subgroup_level = 'Genus' )
  cdf <- as.data.frame(res$cdf, stringsAsFactors = FALSE)

  # Each colour's shade, 1 for the darkest of its palette.
  pal   <- lapply(microshades_palettes, rev)
  shade <- vapply(cdf$hex, function (h) {
    for (q in pal) if (h %in% q) return (match(h, q))
    NA_integer_
  }, 1L)

  out <- data.frame(
    group    = as.character(cdf$Top_Phylum),
    subgroup = as.character(cdf$Top_Genus),
    shade    = unname(shade),
    stringsAsFactors = FALSE )
  out <- out[out$group %in% selected, , drop = FALSE]
  out <- out[order(match(out$group, selected), out$shade), , drop = FALSE]
  # With fewer than 4 members kept, microshades lists the group's "Other"
  # once per shade it has left (e.g. shades 3, 4 and 5 for 2 members); the
  # first of them follows the members.
  out <- out[!duplicated(out[c('group', 'subgroup')]), , drop = FALSE]
  rownames(out) <- NULL

  list(selected = selected, cdf = out)
}


cases <- list()

t1 <- make_table(30, paste0('Phylum', LETTERS[1:7]), c(15, 12, 10, 9, 8, 4, 2), 5)
cases$main <- c(t1, list(n_feature_groups = 4L, n_features = 4L),
                run_microshades(t1$counts, t1$phylum, 4L, 4L))
cases$two  <- c(t1, list(n_feature_groups = 3L, n_features = 2L),
                run_microshades(t1$counts, t1$phylum, 3L, 2L))

t3 <- make_table(20, c('Firmicutes', 'Bacteroidota', 'Proteobacteria'), c(7, 5, 2), 0)
cases$small <- c(t3, list(n_feature_groups = 3L, n_features = 4L),
                 run_microshades(t3$counts, t3$phylum, 3L, 4L))

for (nm in names(cases)) {
  cat(nm, ':', paste(cases[[nm]]$selected, collapse = ', '), '\n')
  print(cases[[nm]]$cdf)
}

saveRDS(
  cases, file.path('tests', 'tinytests', 'expected', 'composition_barplot.rds'),
  version = 2 )
