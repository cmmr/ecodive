# Machine learning with log-ratio transforms

Machine learning with log-ratio transforms

## Machine learning

Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
and compositional geometry; don't expect them to improve prediction.
With random forests predicting 65 metadata variables in four datasets,
Yerke et al. (2024) found that alr, clr and ilr (PhILR) "performed
generally slightly worse or only as well as compositionally naïve
transformations", and relative abundances did best, by a small margin.
Compare with `norm = 'percent'` (see
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md))
before choosing a log-ratio predictor table.
