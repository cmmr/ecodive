# Builds data/valencia_centroids.rda: VALENCIA's 13 reference centroids.
# Run from the package root:
#
#     Rscript data-raw/valencia_centroids.r
#
# Source: CST_centroids_012920.csv from VALENCIA (France et al. 2020), at
# commit 8559d45 of github.com/ravel-lab/VALENCIA, fetched through
# jsDelivr's mirror of GitHub. (VALENCIA's README calls the file
# CST_profiles_012920.csv; in the repository it is CST_centroids_012920.csv.)
# The file's MD5 is checked, so the dataset is exactly the published one.
#
# The centroids are the mean relative abundances of the samples assigned to
# each sub-CST among 13,160 vaginal microbiome profiles of 1,975 women; only
# these means are shipped, no individual data. VALENCIA is released under the
# MIT licence:
#
#   MIT License
#
#   Copyright (c) 2023 ravel-lab
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to permit
#   persons to whom the Software is furnished to do so, subject to the
#   following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
#   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
#   USE OR OTHER DEALINGS IN THE SOFTWARE.

url <- paste0(
  'https://cdn.jsdelivr.net/gh/ravel-lab/VALENCIA@8559d45/',
  'CST_centroids_012920.csv' )
md5 <- '40029ff8c590f919a8506ad3bff4269f'

csv <- tempfile(fileext = '.csv')
utils::download.file(url, csv, mode = 'wb', quiet = TRUE)
stopifnot(unname(tools::md5sum(csv)) == md5)

# check.names = FALSE keeps the taxon names as they are ('BVAB2/BVAB3',
# 'g_Escherichia.Shigella', ...): samples are matched to them exactly.
df <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(names(df)[[1]] == 'sub_CST')

valencia_centroids <- as.matrix(df[, -1])
dimnames(valencia_centroids) <- list(df$sub_CST, names(df)[-1])

stopifnot(
  identical(rownames(valencia_centroids), c(
    'I-A', 'I-B', 'II', 'III-A', 'III-B', 'IV-A', 'IV-B',
    'IV-C0', 'IV-C1', 'IV-C2', 'IV-C3', 'IV-C4', 'V' )),
  identical(dim(valencia_centroids), c(13L, 199L)),
  !anyDuplicated(colnames(valencia_centroids)),
  is.double(valencia_centroids),
  all(is.finite(valencia_centroids) & valencia_centroids >= 0) )

save(valencia_centroids, file = 'data/valencia_centroids.rda', version = 2, compress = 'xz')
