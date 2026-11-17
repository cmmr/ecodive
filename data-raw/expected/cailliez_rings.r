# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

# "Ring" distances on which an earlier version of pcoa()'s Cailliez
# iteration settled on a root that wasn't the largest, for
# test-pcoa_correct.r (through data-raw/expected/pcoa_correct.r). g
# points on a circle; the distance between two
# depends only on their circular separation k = 1, ..., g/2, through a
# profile h[k]. The profile is hill-climbed (deterministically) to put the
# eigen-mode that sets the Cailliez constant far up B's spectrum, where
# Ritz vectors from B's smallest eigenvalues miss it. Run with:
#
#     Rscript data-raw/expected/cailliez_rings.r
#
# Writes each profile, as hexadecimal floats, to
# data-raw/expected/pcoa_correct/cailliez_rings.csv (committed), which
# data-raw/expected/pcoa_correct.r reads.

# For a circulant distance matrix, B and B2 share the Fourier modes: per
# mode, B's eigenvalue b, B2's b2, and the larger root of c^2/2 + 2 b2 c + b.
mode_stats <- function (h, g) {
  dk   <- c(0, h[pmin(1:(g - 1), g - 1:(g - 1))])
  b    <- -0.5 * Re(stats::fft(dk^2))[-1]
  b2   <- -0.5 * Re(stats::fft(dk))[-1]
  disc <- 4 * b2^2 - 2 * b
  list(b = b, cf = ifelse(disc >= 0, -2 * b2 + sqrt(pmax(disc, 0)), -Inf))
}

# How many of B's eigenvalues lie below the Cailliez mode's.
objective <- function (h, g) {
  ms    <- mode_stats(h, g)
  fstar <- which.max(ms$cf)
  sum(ms$b < ms$b[fstar] - 1e-9 * max(abs(ms$b))) +
    1e-3 * (ms$cf[fstar] - max(ms$cf[-fstar]))
}

profile <- function (g) {
  set.seed(11)
  best_h <- 1 + stats::runif(g / 2, 0, 1)
  best   <- objective(best_h, g)
  for (it in 1:20000) {
    h    <- best_h
    j    <- sample(g / 2, sample(1:3, 1))
    h[j] <- pmax(0.001, h[j] * exp(stats::rnorm(length(j), 0, 0.3)))
    v    <- objective(h, g)
    if (v >= best) { best <- v; best_h <- h }
  }
  best_h
}

rows <- do.call(rbind, lapply(c(120L, 200L), function (g)
  data.frame(g = g, k = seq_len(g / 2), h = sprintf('%a', profile(g))) ))

write.csv(rows, file.path('data-raw', 'expected', 'pcoa_correct', 'cailliez_rings.csv'), row.names = FALSE)
