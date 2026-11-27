// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

// https://pdodds.w3.uvm.edu/research/papers/others/everything/cha2007a.pdf

// The partition of beta diversity (the *_part() kernels below) follows
// abdiv 0.2.0's R/beta_components.R, under the MIT License:
//
//   Copyright (c) 2020 Kyle Bittinger
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// The multiple-site values (C_beta_multi()) are written from Baselga
// (2010, 2012, 2017), and checked against betapart's outputs only.

#include "ecodive.h"
#include <Rmath.h> // lgammafn

#define BDIV_BHATTACHARYYA  1
#define BDIV_BRAY           2
#define BDIV_CANBERRA       3
#define BDIV_CHEBYSHEV      4
#define BDIV_CLARK          6
#define BDIV_DIVERGENCE     7
#define BDIV_EUCLIDEAN      8
#define BDIV_GOWER          9
#define BDIV_HAMMING       10
#define BDIV_HORN          11
#define BDIV_JACCARD       12
#define BDIV_JSD           13
#define BDIV_LORENTZIAN    14
#define BDIV_MANHATTAN     15
#define BDIV_MINKOWSKI     16
#define BDIV_MORISITA      17
#define BDIV_MOTYKA        18
#define BDIV_OCHIAI        19
#define BDIV_SOERGEL       20
#define BDIV_SORENSEN      21
#define BDIV_SQUARED_CHISQ 22
#define BDIV_SQUARED_CHORD 23
#define BDIV_WAVE_HEDGES   24

// Baselga's partition of beta diversity: one component of one family.
#define BDIV_SORENSEN_TURNOVER   25
#define BDIV_SORENSEN_NESTEDNESS 26
#define BDIV_JACCARD_TURNOVER    27
#define BDIV_JACCARD_NESTEDNESS  28
#define BDIV_BRAY_BALANCED       29
#define BDIV_BRAY_GRADIENT       30
#define BDIV_SOERGEL_BALANCED    31
#define BDIV_SOERGEL_GRADIENT    32

// The vegdist() metrics (plans/vegdist_metrics.md). Written from vegan's
// documentation, the papers and vegan's outputs only.
#define BDIV_KULCZYNSKI   33
#define BDIV_MOUNTFORD    34
#define BDIV_RAUP_CRICK   35
#define BDIV_BINOMIAL     36
#define BDIV_CHAO_JACCARD 37
#define BDIV_CAO          38
#define BDIV_MAHALANOBIS  39
#define BDIV_ALT_GOWER    40

// The families of beta_partition() and beta_partition_multi().
#define PART_SORENSEN 1
#define PART_JACCARD  2
#define PART_BRAY     3
#define PART_RUZICKA  4

static int     n_samples;
static int     n_otus;
static int     n_dist;
static int     n_pairs;
static int    *pos_vec;
static int    *otu_vec;
static double *val_vec;
static double *clr_vec;
static int    *pairs_vec;
static double *dist_vec;
static SEXP   *sexp_extra;

// For the partition kernels: which component goes to dist_vec (PART_OUT()),
// and the triangles of the other two when all three are wanted.
static int     part_mode;
static double *dist_vec2;
static double *dist_vec3;


/*
 * The FOREACH_PAIR macro iterates through all combinations of 
 * samples. Skips unwanted pairings and pairings not assigned to 
 * the current thread. Ensures all threads process the same 
 * number of pairs. The code should assign to `distance`.
 * 
 * The FOREACH_OTU macro iterates through all OTU abundances for 
 * a given pair of samples, assigning the values to `x` and `y`.
 * 
 * Implemented as macros to avoid the overhead of a function
 * call or the messiness of duplicated code.
 */
#define FOREACH_PAIR(expression)                               \
  do {                                                         \
    int thread_i  = ((worker_t *)arg)->i;                      \
    int n_threads = ((worker_t *)arg)->n;                      \
    int dist_idx  = 0;                                         \
                                                               \
    if (pairs_vec == NULL) { /* All vs All */                  \
                                                               \
      for (int sam_i = 0; sam_i < n_samples - 1; sam_i++) {    \
        for (int sam_j = sam_i + 1; sam_j < n_samples; sam_j++) {\
          if (dist_idx % n_threads == thread_i) {              \
                                                               \
            double distance = 0;                               \
            int    out_idx  = dist_idx;                        \
                                                               \
            expression;                                        \
                                                               \
            (void)out_idx;                                     \
            dist_vec[dist_idx] = distance;                     \
          }                                                    \
          dist_idx++;                                          \
        }                                                      \
      }                                                        \
                                                               \
    } else { /* Specific Pairs of Samples */                   \
                                                               \
      int pair_idx = thread_i;                                 \
      for (; pair_idx < n_pairs; pair_idx += n_threads) {      \
                                                               \
        dist_idx = pairs_vec[pair_idx]; /* 1-based */          \
                                                               \
        int sam_i          = 0;                                \
        int sam_j          = dist_idx;                         \
        int pairs_in_block = n_samples - 1;                    \
                                                               \
        while (sam_j > pairs_in_block) {                       \
          sam_i++;                                             \
          sam_j -= pairs_in_block;                             \
          pairs_in_block--;                                    \
        }                                                      \
                                                               \
        sam_j += sam_i;                                        \
                                                               \
        double distance = 0;                                   \
        int    out_idx  = dist_idx - 1;                        \
                                                               \
        expression;                                            \
                                                               \
        (void)out_idx;                                         \
        dist_vec[dist_idx - 1] = distance;                     \
      }                                                        \
    }                                                          \
  } while (0)


#define FOREACH_OTU(expression)                                \
  do {                                                         \
    int    *i      = otu_vec + pos_vec[sam_i];                 \
    int    *j      = otu_vec + pos_vec[sam_j];                 \
    int    *i_end  = otu_vec + pos_vec[sam_i + 1];             \
    int    *j_end  = otu_vec + pos_vec[sam_j + 1];             \
    double *val_i  = val_vec + pos_vec[sam_i];                 \
    double *val_j  = val_vec + pos_vec[sam_j];                 \
    double  x_zero = clr_vec ? clr_vec[sam_i] : 0;             \
    double  y_zero = clr_vec ? clr_vec[sam_j] : 0;             \
    int     otu = 0, n_ops = 0;                                \
    double  x, y;                                              \
    while (1) {                                                \
      if (i != i_end && j != j_end) {                          \
        if (*i == *j) {                                        \
          otu = *i;                                            \
          x   = *val_i; i++; val_i++;                          \
          y   = *val_j; j++; val_j++;                          \
        }                                                      \
        else if (*i < *j) {                                    \
          otu = *i;                                            \
          x   = *val_i; i++; val_i++;                          \
          y   = y_zero;                                        \
        }                                                      \
        else {                                                 \
          otu = *j;                                            \
          x   = x_zero;                                        \
          y   = *val_j; j++; val_j++;                          \
        }                                                      \
      }                                                        \
      else if (i != i_end) {                                   \
        otu = *i;                                              \
        x   = *val_i; i++; val_i++;                            \
        y   = y_zero;                                          \
      }                                                        \
      else if (j != j_end) {                                   \
        otu = *j;                                              \
        x   = x_zero;                                          \
        y   = *val_j; j++; val_j++;                            \
      }                                                        \
      else {                                                   \
        break;                                                 \
      }                                                        \
      n_ops++;                                                 \
      expression;                                              \
    }                                                          \
    if (clr_vec) { /* Double zeros */                          \
      x = x_zero, y = y_zero;                                  \
      while (n_ops < n_otus) {                                 \
        n_ops++;                                               \
        expression;                                            \
      }                                                        \
    }                                                          \
    (void)otu;                                                 \
  } while (0)


#define WITH_ABJ(expression)                                   \
  do {                                                         \
    int *i     = otu_vec + pos_vec[sam_i];                     \
    int *j     = otu_vec + pos_vec[sam_j];                     \
    int *i_end = otu_vec + pos_vec[sam_i + 1];                 \
    int *j_end = otu_vec + pos_vec[sam_j + 1];                 \
    double A = 0, B = 0, J = 0;                                \
    while (i != i_end && j != j_end) {                         \
      if      (*i == *j) { A++; B++; J++; i++; j++; }          \
      else if (*i < *j)  { A++; i++; }                         \
      else               { B++; j++; }                         \
    }                                                          \
    A += i_end - i;                                            \
    B += j_end - j;                                            \
                                                               \
    expression;                                                \
                                                               \
  } while (0)








//======================================================
// Bhattacharyya
// -log(sum(sqrt(x * y)))
//======================================================
static void *bhattacharyya(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(distance += sqrt(x * y));
  
    distance = -1 * log(distance);
  );
  
  return NULL;
}



//======================================================
// Dice-Sorensen; Bray-Curtis
// sum(abs(x-y)) / sum(x+y)
//======================================================
static void *bray(void *arg) {
  FOREACH_PAIR(
  
    double diffs = 0;
    double sums  = 0;
    
    FOREACH_OTU(
      sums  += x + y;
      diffs += fabs(x - y);
    );
  
    distance = diffs / sums;
  );
  
  return NULL;
}



//======================================================
// Canberra
// nz = (x+y) > 0; x = x[nz]; y = y[nz]
// sum(abs(x-y) / (x + y)) / sum(nz)
//======================================================
static void *canberra(void *arg) {
  FOREACH_PAIR(
    FOREACH_OTU(distance += fabs(x - y) / (x + y));
  );
  
  return NULL;
}


//======================================================
// Chebyshev
// max(abs(x - y))
//======================================================
static void *chebyshev(void *arg) {
  FOREACH_PAIR(
    FOREACH_OTU(
      double d = fabs(x - y);
      if (d > distance) distance = d;
    );
  );
  
  return NULL;
}


//======================================================
// Clark
// sqrt(sum((abs(x - y) / (x + y)) ^ 2))
//======================================================
static void *clark(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double d = (x - y) / (x + y);
      distance += d * d;
    );
  
    distance = sqrt(distance);
  );
  
  return NULL;
}


//======================================================
// Divergence
// 2 * sum((x-y)^2 / (x+y)^2)
//======================================================
static void *divergence(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double diff = x - y;
      double sum  = x + y;
      distance += (diff * diff) / (sum * sum);
    );
  
    distance = 2 * distance;
  );
  
  return NULL;
}


//======================================================
// Euclidean
// sqrt(sum((x-y)^2))
//======================================================
static void *euclidean(void *arg) {
  FOREACH_PAIR(
    
    FOREACH_OTU(
      double d = x - y;
      distance += d * d;
    );
  
    distance = sqrt(distance);
  );
  
  return NULL;
}



//======================================================
// Gower
// sum(abs(x-y) / r) / n
//======================================================

static double *gower_range_vec;

static void *gower(void *arg) {
  
  FOREACH_PAIR(

    if (clr_vec) {

      // FOREACH_OTU doesn't track `otu` for double zeros, which are
      // non-zero after CLR. Walk every OTU so each gets its own range.
      int    *i     = otu_vec + pos_vec[sam_i];
      int    *j     = otu_vec + pos_vec[sam_j];
      int    *i_end = otu_vec + pos_vec[sam_i + 1];
      int    *j_end = otu_vec + pos_vec[sam_j + 1];
      double *val_i = val_vec + pos_vec[sam_i];
      double *val_j = val_vec + pos_vec[sam_j];

      for (int otu = 0; otu < n_otus; otu++) {

        double x = clr_vec[sam_i];
        double y = clr_vec[sam_j];
        if (i != i_end && *i == otu) { x = *val_i; i++; val_i++; }
        if (j != j_end && *j == otu) { y = *val_j; j++; val_j++; }

        double range = gower_range_vec[otu];
        if (range) distance += fabs(x - y) / range;
      }
    }

    else {
      FOREACH_OTU(

        double range = gower_range_vec[otu];

        if (range) {
          distance += fabs(x - y) / range;
        }
      );
    }

    distance /= n_otus;
  );
  
  return NULL;
}

static pthread_func_t gower_setup(void) {
  
  gower_range_vec = (double*) safe_malloc(n_otus * sizeof(double));
  double *min_vec = (double*) safe_malloc(n_otus * sizeof(double));
  double *max_vec = (double*) safe_malloc(n_otus * sizeof(double));


  // Search for min/max values across all samples, including the
  // implicit values of absent OTUs. Those are 0, except after CLR
  // normalization, where they are clr_vec[sam].
  for (int sam = 0; sam < n_samples; sam++) {

    double zero = clr_vec ? clr_vec[sam] : 0;
    int    i    = pos_vec[sam];
    int    end  = pos_vec[sam + 1];

    for (int otu = 0; otu < n_otus; otu++) {

      double val = zero;
      if (i < end && otu_vec[i] == otu) { val = val_vec[i]; i++; }

      if (sam == 0 || val < min_vec[otu]) min_vec[otu] = val;
      if (sam == 0 || val > max_vec[otu]) max_vec[otu] = val;
    }
  }

  // Assign final range values
  for (int otu = 0; otu < n_otus; otu++) {
    gower_range_vec[otu] = max_vec[otu] - min_vec[otu];
  }

  free_one(min_vec);
  free_one(max_vec);
  
  return gower;
}


//======================================================
// Hamming
// sum(xor(x, y))
//======================================================
static void *hamming(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = A + B - 2 * J);
  );
  
  return NULL;
}


//======================================================
// Horn
// z <- sum(x^2) / sum(x)^2 + sum(y^2) / sum(y)^2
// 1 - ((2 * sum(x * y)) / (z * sum(x) * sum(y)))
//======================================================
static void *horn(void *arg) {
  
  FOREACH_PAIR(
    
    double sum_x  = 0;
    double sum_y  = 0;
    double sum_x2 = 0;
    double sum_y2 = 0;
    
    FOREACH_OTU(
      distance += x * y;
      sum_x    += x;
      sum_y    += y;
      sum_x2   += x * x;
      sum_y2   += y * y;
    );
    
    sum_x2 /= sum_x * sum_x;
    sum_y2 /= sum_y * sum_y;
    
    distance = 1 - (2 * distance) / ((sum_x2 + sum_y2) * sum_x * sum_y);
  );
  
  return NULL;
}


//======================================================
// Jaccard
// sum(xor(x, y)) / sum(x | y)
//======================================================
static void *jaccard(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = (A + B - 2 * J) / (A + B - J));
  );
  
  return NULL;
}


//======================================================
// Jensen-Shannon Divergence (JSD)
// sum(x * log(2*x / (x+y)), y * log(2*y / (x+y))) / 2
//======================================================
static void *jsd(void *arg) {
  
  FOREACH_PAIR(
    
    FOREACH_OTU(
      if (x) distance += x * log(2 * x / (x + y));
      if (y) distance += y * log(2 * y / (x + y));
    );
    
    distance /= 2;
  );
  
  return NULL;
}


//======================================================
// Lorentzian
// sum(log(1 + abs(x - y)))
//======================================================
static void *lorentzian(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(distance += log(1 + fabs(x - y)));
  );
  
  return NULL;
}


//======================================================
// Manhattan
// sum(abs(x-y))
//======================================================
static void *manhattan(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(distance += fabs(x - y));
  );
  
  return NULL;
}


//======================================================
// Minkowski
// sum(abs(x - y)^p) ^ (1/p)
//======================================================
static void *minkowski(void *arg) {
  
  double power     = asReal(*sexp_extra);
  double inv_power = 1 / power;
  
  FOREACH_PAIR(
    
    FOREACH_OTU(distance += pow(fabs(x - y), power));
  
    distance = pow(distance, inv_power);
  );
  
  return NULL;
}


//======================================================
// Morisita
// 
// simpson_x <- sum(x * (x - 1)) / (sum(x) * (sum(x) - 1))
// simpson_y <- sum(y * (y - 1)) / (sum(y) * (sum(y) - 1))
// 1 - ((2 * sum(x * y)) / ((simpson_x + simpson_y) * sum(x) * sum(y)))
//======================================================
static void *morisita(void *arg) {
  
  FOREACH_PAIR(
    
    double sum_x = 0;
    double sum_y = 0;
    double simpson_x = 0;
    double simpson_y = 0;
    
    FOREACH_OTU(
      distance += x * y;
      sum_x    += x; simpson_x += x * (x - 1);
      sum_y    += y; simpson_y += y * (y - 1);
    );
    
    simpson_x /= sum_x * (sum_x - 1);
    simpson_y /= sum_y * (sum_y - 1);
    
    distance = 1 - (2 * distance) / ((simpson_x + simpson_y) * sum_x * sum_y);
  );
  
  return NULL;
}




//======================================================
// Motyka
// sum(pmax(x, y)) / sum(x, y)
//======================================================
static void *motyka(void *arg) {
  
  FOREACH_PAIR(
    
    double sums = 0;
  
    FOREACH_OTU(
      distance += (x > y) ? x : y;
      sums     += x + y;
    );
    
    distance /= sums;
  );
  
  return NULL;
}


//======================================================
// Dice-Sorensen
// 1 - 2 * sum(x & y) / sum(x>0, y>0)
//======================================================
static void *sorensen(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = 1 - (2 * J) / (A + B));
  );
  
  return NULL;
}



//======================================================
// Ochiai
// sum((x & y)) / sqrt(sum(x > 0) * sum(y > 0))
//======================================================
static void *ochiai(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(distance = 1 - J / sqrt(A * B));
  );
  
  return NULL;
}
  
  
//======================================================
// Soergel
// 1 - sum(pmin(x, y)) / sum(pmax(x, y))
//======================================================
static void *soergel(void *arg) {
  
  FOREACH_PAIR(
    
    double min_sum = 0;
    double max_sum = 0;
  
    FOREACH_OTU(
      if (x < y) { min_sum += x; max_sum += y; } 
      else       { min_sum += y; max_sum += x; }
    );
    
    distance = 1 - (min_sum / max_sum);
  );
  
  return NULL;
}


//======================================================
// Squared Ch-Squared
// sum((x - y) ^ 2 / (x + y))
//======================================================
static void *squared_chisq(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      double d  = (x - y);
      distance += (d * d) / (x + y);
    );
  );
  
  return NULL;
}


//======================================================
// Squared Chord
// sum((sqrt(x) - sqrt(y)) ^ 2)
//======================================================
static void *squared_chord(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      double d = sqrt(x) - sqrt(y);
      distance += d * d;
    );
  );
  
  return NULL;
}


//======================================================
// Wave Hedges
// sum(abs(x - y) / pmax(x, y))
//======================================================
static void *wave_hedges(void *arg) {
  
  FOREACH_PAIR(
    FOREACH_OTU(
      if (x > y) { distance += (x - y) / x; }
      else       { distance += (y - x) / y; }
    );
  );

  return NULL;
}




//======================================================
// The vegdist() metrics ecodive lacked. Formulas from
// vegan's vegdist() documentation and the cited papers;
// details (Mountford's map, Raup-Crick's tail, Cao's
// zeros, Chao's a2 = 0 and cap) settled from vegan's
// outputs, as plans/vegdist_metrics.md records.
//======================================================


//======================================================
// Kulczynski (quantitative; Faith et al. 1987)
// 1 - (sum(pmin(x, y)) / sum(x) + sum(pmin(x, y)) / sum(y)) / 2
//======================================================
static void *kulczynski(void *arg) {
  FOREACH_PAIR(
    double min_sum = 0; double sum_x = 0; double sum_y = 0;
    FOREACH_OTU(
      min_sum += (x < y) ? x : y;
      sum_x   += x;
      sum_y   += y;
    );
    distance = 1 - (min_sum / sum_x + min_sum / sum_y) / 2;
  );
  return NULL;
}


//======================================================
// Binomial deviance (Anderson & Millar 2004)
// n = x + y over features with n > 0;
// sum((x log(x/n) + y log(y/n) - n log(1/2)) / n), 0 log 0 = 0
//======================================================
static void *binomial(void *arg) {
  FOREACH_PAIR(
    double total = 0;
    FOREACH_OTU(
      double n = x + y;
      total += n;
      if (n > 0) {
        double t = n * M_LN2;
        if (x > 0) t += x * log(x / n);
        if (y > 0) t += y * log(y / n);
        distance += t / n;
      }
    );
    if (total == 0) distance = R_NaN; // Two empty samples, as vegan.
  );
  return NULL;
}


//======================================================
// Cao's CYd (Cao et al. 1997), natural logs
// Over the S features present in either sample, each
// value raised to at least 0.1: x' = max(x, 0.1),
// y' = max(y, 0.1), n' = x' + y';
// mean(log(n'/2) - (x' log(y') + y' log(x')) / n')
//======================================================
static void *cao(void *arg) {
  FOREACH_PAIR(
    double n_feat = 0;
    FOREACH_OTU(
      if (x != 0 || y != 0) {
        double xp = (x < 0.1) ? 0.1 : x;
        double yp = (y < 0.1) ? 0.1 : y;
        double n  = xp + yp;
        distance += log(n / 2) - (xp * log(yp) + yp * log(xp)) / n;
        n_feat++;
      }
    );
    distance /= n_feat;
  );
  return NULL;
}


//======================================================
// Alternative Gower (Anderson et al. 2006)
// sum(abs(x - y)) / sum(x != 0 | y != 0)
// With clr (no zeros), every feature counts.
//======================================================
static void *alt_gower(void *arg) {
  FOREACH_PAIR(
    double nz = 0;
    FOREACH_OTU(
      distance += fabs(x - y);
      nz++;
    );
    distance /= nz;
  );
  return NULL;
}


//======================================================
// Chao's abundance-based Jaccard (Chao et al. 2005)
// 1 - U V / (U + V - U V), with
//   U = C_x / N_x + (N_y - 1) / N_y * f1 / (2 f2) * S1_x / N_x
// C_x: individuals in x of the shared features; N_x: all
// individuals in x; f1, f2: shared features with 1 and 2
// individuals in y (f2 = 0 counts as 1); S1_x: individuals
// in x of the shared features with 1 individual in y. V the
// same with x and y swapped. U and V are capped at 1.
// No shared features: 1 (NaN if both samples are empty).
//======================================================
static void *chao_jaccard(void *arg) {
  FOREACH_PAIR(
    double n_x = 0; double n_y = 0; double c_x = 0; double c_y = 0;
    double shared = 0; double s1_x = 0; double s1_y = 0;
    double f1_x = 0; double f2_x = 0; double f1_y = 0; double f2_y = 0;
    FOREACH_OTU(
      n_x += x;
      n_y += y;
      if (x > 0 && y > 0) {
        shared++;
        c_x += x;
        c_y += y;
        if      (y == 1) { f1_y++; s1_x += x; }
        else if (y == 2) { f2_y++;            }
        if      (x == 1) { f1_x++; s1_y += y; }
        else if (x == 2) { f2_x++;            }
      }
    );
    if (shared == 0) {
      distance = (n_x > 0 || n_y > 0) ? 1 : R_NaN;
    }
    else {
      double u = c_x / n_x + (n_y - 1) / n_y * f1_y / (2 * (f2_y > 0 ? f2_y : 1)) * s1_x / n_x;
      double v = c_y / n_y + (n_x - 1) / n_x * f1_x / (2 * (f2_x > 0 ? f2_x : 1)) * s1_y / n_y;
      if (u > 1) u = 1;
      if (v > 1) v = 1;
      distance = 1 - u * v / (u + v - u * v);
    }
  );
  return NULL;
}


//======================================================
// Mountford (1962)
// M is the positive root of exp(aM) + exp(bM) = 1 +
// exp((a + b - j) M), with a and b the samples' richness
// and j the shared features. Divided by exp((a + b - j) M)
// and with p = b - j, q = a - j, t = a + b - j, that is
//   g(M) = expm1(-pM) + expm1(-qM) - expm1(-tM) = 0,
// which can't overflow. g(0) = 0, g'(0) = j > 0, and g
// falls to -1, so there is one positive root. The result
// is log(2) - M: log(2) for j = 0, and 0 when j equals a
// or b; NaN when a sample is empty, as vegan. Solved by
// Newton's method, safeguarded by bisection, to the last
// bits; vegan's own solution agrees to about 1e-10 while
// a and b are at most about 50.
//======================================================
static double mountford_root(double a, double b, double j) {

  double p = b - j, q = a - j, t = a + b - j;

  // The start: the root of g's second-order Taylor expansion.
  double m  = 2 * j / (t * t - p * p - q * q);
  double lo = 0, hi = m;

  // Bracket the root: g(lo) > 0 > g(hi).
  while (expm1(-p * hi) + expm1(-q * hi) - expm1(-t * hi) > 0) {
    lo = hi;
    hi *= 2;
  }

  for (int iter = 0; iter < 200; iter++) {
    double ep = exp(-p * m), eq = exp(-q * m), et = exp(-t * m);
    double g  = expm1(-p * m) + expm1(-q * m) - expm1(-t * m);
    double dg = -p * ep - q * eq + t * et;

    if (g > 0) lo = m; else if (g < 0) hi = m; else return m;

    double next = m - g / dg;
    if (!(next > lo && next < hi)) next = (lo + hi) / 2;

    if (fabs(next - m) <= 4 * DBL_EPSILON * next) return next;
    m = next;
    if (hi - lo <= 4 * DBL_EPSILON * m) return m;
  }

  return m; // # nocov
}

static void *mountford(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(
      if      (A == 0 || B == 0)               distance = R_NaN;
      else if (J == 0)                         distance = M_LN2;
      else if (J == A || J == B)               distance = 0;
      else    distance = M_LN2 - mountford_root(A, B, J);
    );
  );
  return NULL;
}


//======================================================
// Raup-Crick, analytic (Raup & Crick 1979)
// P(X >= j) for X hypergeometric: a features of the
// table's S drawn for one sample, b for the other, and X
// shared. S is the number of features in the table,
// empty ones included. The smaller tail is summed from its
// far end inwards, and the other found as its complement,
// so that small probabilities keep their precision.
//======================================================
static double *raup_lfact; // log(k!) for k = 0 .. n_otus

static double raup_lpmf(int x, int a, int b, int s) {
  return
    raup_lfact[a] - raup_lfact[x] - raup_lfact[a - x] +
    raup_lfact[s - a] - raup_lfact[b - x] - raup_lfact[s - a - b + x] -
    raup_lfact[s] + raup_lfact[b] + raup_lfact[s - b];
}

static void *raup_crick(void *arg) {
  int s = n_otus;
  FOREACH_PAIR(
    WITH_ABJ(
      int na = (int) A; int nb = (int) B; int nj = (int) J;
      int x_lo = na + nb - s; if (x_lo < 0) x_lo = 0;
      int x_hi = (na < nb) ? na : nb;

      if (nj <= x_lo) {
        distance = 1;
      }
      else {
        // The mean, a * b / s: the upper tail from j is the
        // smaller one when j is above it.
        double mean = (double) na * nb / s;
        double tail = 0;
        if (nj > mean) {
          for (int x = x_hi; x >= nj; x--) tail += exp(raup_lpmf(x, na, nb, s));
          distance = tail;
        }
        else {
          for (int x = x_lo; x < nj; x++) tail += exp(raup_lpmf(x, na, nb, s));
          distance = 1 - tail;
        }
        if (distance > 1) distance = 1;
        if (distance < 0) distance = 0;
      }
    );
  );
  return NULL;
}

static pthread_func_t raup_crick_setup(void) {
  raup_lfact = (double*) safe_malloc((size_t) (n_otus + 1) * sizeof(double));
  for (int k = 0; k <= n_otus; k++) raup_lfact[k] = lgammafn(k + 1.0);
  return raup_crick;
}


//======================================================
// Mahalanobis (Mahalanobis 1936; Mardia et al. 1979)
// sqrt((x - y)' S^-1 (x - y)), S the features' sample
// covariance (n - 1). Without S^-1: the centred n x D
// table X = Q R (Householder), and Y = sqrt(n - 1) X R^-1,
// whose Euclidean distances are the Mahalanobis distances.
// The setup runs serially; the rows of Y and the distances
// in parallel, each in a fixed order.
//======================================================
static double *maha_x;   // n_samples x n_otus, row-major: the centred table
static double *maha_r;   // n_otus x n_otus, row-major: R (upper triangle)
static double *maha_y;   // n_samples x n_otus, row-major: the whitened table

static void *maha_rows(void *arg) {
  int    thread_i  = ((worker_t *)arg)->i;
  int    n_threads = ((worker_t *)arg)->n;
  int    d         = n_otus;
  double scale     = sqrt((double) n_samples - 1);
  for (int s = thread_i; s < n_samples; s += n_threads) {
    const double *x = maha_x + (size_t) s * d;
    double       *y = maha_y + (size_t) s * d;
    // y R = x, solved left to right.
    for (int k = 0; k < d; k++) {
      double v = x[k];
      for (int i = 0; i < k; i++) v -= y[i] * maha_r[(size_t) i * d + k];
      y[k] = v / maha_r[(size_t) k * d + k];
    }
    for (int k = 0; k < d; k++) y[k] *= scale;
  }
  return NULL;
}

static void *mahalanobis(void *arg) {
  int d = n_otus;
  FOREACH_PAIR(
    const double *yi = maha_y + (size_t) sam_i * d;
    const double *yj = maha_y + (size_t) sam_j * d;
    for (int k = 0; k < d; k++) {
      double diff = yi[k] - yj[k];
      distance += diff * diff;
    }
    distance = sqrt(distance);
  );
  return NULL;
}

static pthread_func_t mahalanobis_setup(int n_threads) {

  int n = n_samples, d = n_otus;

  if (d >= n) {
    free_all();
    error(
      "Mahalanobis distances need more samples than features and no "
      "collinear features; reduce the features first (e.g. pca()).");
  }

  size_t nd = (size_t) n * d;
  maha_x = (double*) safe_malloc(nd * sizeof(double));
  maha_y = (double*) safe_malloc(nd * sizeof(double));
  maha_r = (double*) safe_malloc((size_t) d * d * sizeof(double));
  double *a = maha_y; // column-major n x d work copy for the QR

  // The dense table (absent features are 0, or the sample's
  // value for zero after CLR), then centred by column.
  for (size_t i = 0; i < nd; i++) maha_x[i] = 0;
  for (int s = 0; s < n; s++) {
    double zero = clr_vec ? clr_vec[s] : 0;
    double *row = maha_x + (size_t) s * d;
    if (zero != 0) for (int k = 0; k < d; k++) row[k] = zero;
    for (int e = pos_vec[s]; e < pos_vec[s + 1]; e++) row[otu_vec[e]] = val_vec[e];
  }
  for (int k = 0; k < d; k++) {
    double mean = 0;
    for (int s = 0; s < n; s++) mean += maha_x[(size_t) s * d + k];
    mean /= n;
    for (int s = 0; s < n; s++) {
      maha_x[(size_t) s * d + k] -= mean;
      a[s + (size_t) k * n] = maha_x[(size_t) s * d + k];
    }
  }

  // Householder QR of the centred table; R's diagonal is
  // checked against a relative tolerance.
  double max_norm = 0;
  for (int k = 0; k < d; k++) {
    double nrm = 0;
    for (int s = 0; s < n; s++) nrm += a[s + (size_t) k * n] * a[s + (size_t) k * n];
    nrm = sqrt(nrm);
    if (nrm > max_norm) max_norm = nrm;
  }

  int bad = (max_norm == 0);
  for (int k = 0; k < d && !bad; k++) {
    double *col = a + (size_t) k * n;
    double  nrm = 0;
    for (int s = k; s < n; s++) nrm += col[s] * col[s];
    nrm = sqrt(nrm);
    if (!(nrm > 1e-9 * max_norm)) { bad = 1; break; }
    double alpha = (col[k] > 0) ? -nrm : nrm;
    double v0    = col[k] - alpha;
    // v = (v0, col[k+1..n-1]); H = I - 2 v v' / (v'v)
    double vtv = v0 * v0;
    for (int s = k + 1; s < n; s++) vtv += col[s] * col[s];
    for (int c = k + 1; c < d; c++) {
      double *cc  = a + (size_t) c * n;
      double  dot = v0 * cc[k];
      for (int s = k + 1; s < n; s++) dot += col[s] * cc[s];
      double f = 2 * dot / vtv;
      cc[k] -= f * v0;
      for (int s = k + 1; s < n; s++) cc[s] -= f * col[s];
    }
    col[k] = alpha;
    for (int c = k; c < d; c++) maha_r[(size_t) k * d + c] = a[k + (size_t) c * n];
  }

  if (bad) {
    free_all();
    error(
      "Mahalanobis distances need more samples than features and no "
      "collinear features; reduce the features first (e.g. pca()).");
  }

  run_parallel(maha_rows, n_threads, n);
  return mahalanobis;
}




//======================================================
// Baselga's partition of beta diversity
//
// Incidence (Baselga 2010, 2012), with a shared features, and b and c
// only in the first and second sample:
//   Sorensen  total (b+c)/(2a+b+c), turnover min(b,c) / (a + min(b,c))
//   Jaccard   total (b+c)/(a+b+c),  turnover 2 min(b,c) / (a + 2 min(b,c))
// Abundance (Baselga 2013, 2017), with A = sum(min(x, y)),
// B = sum(max(x - y, 0)) and C = sum(max(y - x, 0)):
//   Bray-Curtis  total (B+C)/(2A+B+C), balanced min(B,C) / (A + min(B,C))
//   Ruzicka      total (B+C)/(A+B+C),  balanced 2 min(B,C) / (A + 2 min(B,C))
// The nestedness (abundance gradient) component is total - turnover
// (balanced), as abdiv 0.2.0 (MIT, Kyle Bittinger) computes it. Each total
// is computed exactly as sorensen(), jaccard(), bray() and soergel()
// compute it, so that it is identical to theirs.
//
// part_mode selects the component that goes to dist_vec: 0 turnover
// (balanced), 1 nestedness (gradient), 2 total; 3 writes all three, to
// dist_vec, dist_vec2 and dist_vec3.
//======================================================

#define PART_OUT(tu, to)                                       \
  do {                                                         \
    double tu_ = (tu), to_ = (to), ne_ = to_ - tu_;            \
    if (part_mode == 3) {                                      \
      distance           = tu_;                                \
      dist_vec2[out_idx] = ne_;                                \
      dist_vec3[out_idx] = to_;                                \
    }                                                          \
    else if (part_mode == 0) distance = tu_;                   \
    else if (part_mode == 1) distance = ne_;                   \
    else                     distance = to_;                   \
  } while (0)


static void *sorensen_part(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(
      double mn = (A < B ? A : B) - J;
      PART_OUT(mn / (J + mn), 1 - (2 * J) / (A + B));
    );
  );
  return NULL;
}

static void *jaccard_part(void *arg) {
  FOREACH_PAIR(
    WITH_ABJ(
      double mn = (A < B ? A : B) - J;
      PART_OUT(2 * mn / (J + 2 * mn), (A + B - 2 * J) / (A + B - J));
    );
  );
  return NULL;
}

static void *bray_part(void *arg) {
  FOREACH_PAIR(
    double diffs = 0; double sums = 0;
    double a = 0; double b = 0; double c = 0;
    FOREACH_OTU(
      sums  += x + y;
      diffs += fabs(x - y);
      if (x < y) { a += x; c += y - x; }
      else       { a += y; b += x - y; }
    );
    double mn = b < c ? b : c;
    PART_OUT(mn / (a + mn), diffs / sums);
  );
  return NULL;
}

static void *soergel_part(void *arg) {
  FOREACH_PAIR(
    double min_sum = 0; double max_sum = 0;
    double b = 0; double c = 0;
    FOREACH_OTU(
      if (x < y) { min_sum += x; max_sum += y; c += y - x; }
      else       { min_sum += y; max_sum += x; b += x - y; }
    );
    double mn = b < c ? b : c;
    PART_OUT(2 * mn / (min_sum + 2 * mn), 1 - (min_sum / max_sum));
  );
  return NULL;
}


// For the multiple-site values: the smaller and larger of the pair's two
// unshared parts, (b, c) or (B, C), to dist_vec and dist_vec2.
static int part_incidence;

static inline void part_bc(int sam_i, int sam_j, double *b_out, double *c_out) {
  if (part_incidence) {
    WITH_ABJ(*b_out = A - J; *c_out = B - J;);
  }
  else {
    double b = 0, c = 0;
    FOREACH_OTU(
      if (x < y) c += y - x;
      else       b += x - y;
    );
    *b_out = b; *c_out = c;
  }
}

// The sums over the pairs (s, t > s) of sites[] of the smaller and larger
// unshared parts, row by row: row s's sums in order of t, then the rows'
// sums in order of s. The same order whether a row is summed here or by
// multi_rows() on another thread.
static inline void row_sums(const int *sites, int n, int s, double *mn, double *mx) {
  double sum_min = 0, sum_max = 0;
  for (int t = s + 1; t < n; t++) {
    double b, c;
    part_bc(sites[s], sites[t], &b, &c);
    sum_min += b < c ? b : c;
    sum_max += b < c ? c : b;
  }
  *mn = sum_min;
  *mx = sum_max;
}

// All the sites: each row is a task, its two sums to dist_vec[s] and
// dist_vec2[s] (n_samples doubles each).
static int *multi_all;

static void *multi_rows(void *arg) {
  int thread_i  = ((worker_t *)arg)->i;
  int n_threads = ((worker_t *)arg)->n;
  for (int s = thread_i; s < n_samples - 1; s += n_threads)
    row_sums(multi_all, n_samples, s, dist_vec + s, dist_vec2 + s);
  return NULL;
}


// The three multiple-site values (turnover, nestedness, total) of family
// `index` from: sum_min and sum_max, the sums over pairs of sites of the
// smaller and larger unshared parts; and shared, the multiple-site shared
// part: sum(S_i) - S_T for incidence (Baselga 2010, 2012), and
// sum_k (sum_i x_ik - max_i x_ik) for abundance (Baselga 2017).
static void multi_values(int index, double shared, double sum_min, double sum_max, double *out) {
  double tu, to;
  switch (index) {
    case PART_SORENSEN:
    case PART_BRAY:
      tu = sum_min / (shared + sum_min);
      to = (sum_min + sum_max) / (2 * shared + sum_min + sum_max);
      break;
    default: // PART_JACCARD, PART_RUZICKA
      tu = 2 * sum_min / (shared + 2 * sum_min);
      to = (sum_min + sum_max) / (shared + sum_min + sum_max);
      break;
  }
  out[0] = tu;
  out[1] = to - tu;
  out[2] = to;
}


// The multiple-site shared part of the `n_sites` sites in `sites` (in
// order), with `fmax` a scratch vector of n_otus doubles. The two sums run
// in a fixed order: sites, then features.
static double multi_shared(const int *sites, int n_sites, double *fmax) {
  double total = 0, maxes = 0;
  for (int k = 0; k < n_otus; k++) fmax[k] = 0;
  for (int s = 0; s < n_sites; s++) {
    int sam = sites[s];
    for (int p = pos_vec[sam]; p < pos_vec[sam + 1]; p++) {
      double v = part_incidence ? 1 : val_vec[p];
      total += v;
      if (v > fmax[otu_vec[p]]) fmax[otu_vec[p]] = v;
    }
  }
  for (int k = 0; k < n_otus; k++) maxes += fmax[k];
  return total - maxes;
}


// Resampled multiple-site values (betapart's beta.sample()): draw r picks
// multi_n_sites of the sites without replacement, with PCG stream
// (seed, r) and Floyd's algorithm, and writes its three values to row r of
// the n_draws x 3 matrix multi_out. Each draw is one task.
static int     multi_index;
static int     multi_n_sites;
static int     multi_n_draws;
static int     multi_seed;
static double *multi_out;
static double *multi_scratch;  // per worker: n_otus doubles, n_samples ints

static int cmp_int(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

static void *multi_draws(void *arg) {
  int thread_i  = ((worker_t *)arg)->i;
  int n_threads = ((worker_t *)arg)->n;

  size_t  stride = (size_t) n_otus + (size_t) n_samples;
  double *fmax   = multi_scratch + (size_t) thread_i * stride;
  int    *sites  = (int *)(fmax + n_otus);
  int     k      = multi_n_sites;

  for (int r = thread_i; r < multi_n_draws; r += n_threads) {

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (uint64_t) (uint32_t) multi_seed, (uint64_t) r);

    // Floyd's algorithm: k distinct sites of n_samples.
    int n = 0;
    for (int j = n_samples - k; j < n_samples; j++) {
      int t = (int) pcg32_bounded(&rng, (uint32_t) j + 1u), seen = 0;
      for (int s = 0; s < n; s++) if (sites[s] == t) { seen = 1; break; }
      sites[n++] = seen ? j : t;
    }
    qsort(sites, (size_t) k, sizeof(int), cmp_int);

    double sum_min = 0, sum_max = 0;
    for (int s = 0; s < k - 1; s++) {
      double mn, mx;
      row_sums(sites, k, s, &mn, &mx);
      sum_min += mn;
      sum_max += mx;
    }

    double vals[3];
    multi_values(multi_index, multi_shared(sites, k, fmax), sum_min, sum_max, vals);
    for (int v = 0; v < 3; v++)
      multi_out[r + (size_t) v * multi_n_draws] = vals[v];
  }

  return NULL;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_beta_div(
    SEXP sexp_algorithm,   SEXP sexp_otu_mtx,   
    SEXP sexp_margin,      SEXP sexp_norm, 
    SEXP sexp_pairs_vec,   SEXP sexp_n_threads, 
    SEXP sexp_pseudocount, SEXP sexp_extra_args ) {
  
  int norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int n_threads   = asInteger(sexp_n_threads);
  sexp_extra      = &sexp_extra_args;
  init_n_ptrs(10);
  
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  if (norm) normalize(em, norm, n_threads, pseudocount);
  
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  otu_vec   = em->otu_vec;
  val_vec   = em->val_vec;
  clr_vec   = em->clr_vec;
  
  
  // function to run
  // void * (*bdiv_func)(void *) = NULL;
  pthread_func_t bdiv_func = NULL;
  
  switch (asInteger(sexp_algorithm)) {
    case BDIV_BHATTACHARYYA: bdiv_func = bhattacharyya; break;
    case BDIV_BRAY:          bdiv_func = bray;          break;
    case BDIV_CANBERRA:      bdiv_func = canberra;      break;
    case BDIV_CHEBYSHEV:     bdiv_func = chebyshev;     break;
    case BDIV_CLARK:         bdiv_func = clark;         break;
    case BDIV_DIVERGENCE:    bdiv_func = divergence;    break;
    case BDIV_EUCLIDEAN:     bdiv_func = euclidean;     break;
    case BDIV_GOWER:         bdiv_func = gower_setup(); break;
    case BDIV_HAMMING:       bdiv_func = hamming;       break;
    case BDIV_HORN:          bdiv_func = horn;          break;
    case BDIV_JACCARD:       bdiv_func = jaccard;       break;
    case BDIV_JSD:           bdiv_func = jsd;           break;
    case BDIV_LORENTZIAN:    bdiv_func = lorentzian;    break;
    case BDIV_MANHATTAN:     bdiv_func = manhattan;     break;
    case BDIV_MINKOWSKI:     bdiv_func = minkowski;     break;
    case BDIV_MORISITA:      bdiv_func = morisita;      break;
    case BDIV_MOTYKA:        bdiv_func = motyka;        break;
    case BDIV_OCHIAI:        bdiv_func = ochiai;        break;
    case BDIV_SOERGEL:       bdiv_func = soergel;       break;
    case BDIV_SORENSEN:      bdiv_func = sorensen;      break;
    case BDIV_SQUARED_CHISQ: bdiv_func = squared_chisq; break;
    case BDIV_SQUARED_CHORD: bdiv_func = squared_chord; break;
    case BDIV_WAVE_HEDGES:   bdiv_func = wave_hedges;   break;
    case BDIV_SORENSEN_TURNOVER:   part_mode = 0; bdiv_func = sorensen_part; break;
    case BDIV_SORENSEN_NESTEDNESS: part_mode = 1; bdiv_func = sorensen_part; break;
    case BDIV_JACCARD_TURNOVER:    part_mode = 0; bdiv_func = jaccard_part;  break;
    case BDIV_JACCARD_NESTEDNESS:  part_mode = 1; bdiv_func = jaccard_part;  break;
    case BDIV_BRAY_BALANCED:       part_mode = 0; bdiv_func = bray_part;     break;
    case BDIV_BRAY_GRADIENT:       part_mode = 1; bdiv_func = bray_part;     break;
    case BDIV_SOERGEL_BALANCED:    part_mode = 0; bdiv_func = soergel_part;  break;
    case BDIV_SOERGEL_GRADIENT:    part_mode = 1; bdiv_func = soergel_part;  break;
    case BDIV_KULCZYNSKI:    bdiv_func = kulczynski;          break;
    case BDIV_MOUNTFORD:     bdiv_func = mountford;           break;
    case BDIV_RAUP_CRICK:    bdiv_func = raup_crick_setup();  break;
    case BDIV_BINOMIAL:      bdiv_func = binomial;            break;
    case BDIV_CHAO_JACCARD:  bdiv_func = chao_jaccard;        break;
    case BDIV_CAO:           bdiv_func = cao;                 break;
    case BDIV_MAHALANOBIS:   bdiv_func = mahalanobis_setup(n_threads); break;
    case BDIV_ALT_GOWER:     bdiv_func = alt_gower;           break;
  }
  
  if (bdiv_func == NULL) { // # nocov start
    error("Invalid beta diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Create the dist object to return
  n_dist                = n_samples * (n_samples - 1) / 2;
  SEXP sexp_result_dist = PROTECT(allocVector(REALSXP, n_dist));
  dist_vec              = REAL(sexp_result_dist);
  
  SEXP sexp_dist_class = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(sexp_dist_class, 0, mkChar("ecodive_dist"));
  SET_STRING_ELT(sexp_dist_class, 1, mkChar("dist"));
  SEXP sexp_size_val   = PROTECT(ScalarInteger(n_samples));
  SEXP sexp_diag_val   = PROTECT(ScalarLogical(0));
  SEXP sexp_upper_val  = PROTECT(ScalarLogical(0));
  
  setAttrib(sexp_result_dist, R_ClassSymbol,     sexp_dist_class);
  setAttrib(sexp_result_dist, install("Size"),   sexp_size_val);
  setAttrib(sexp_result_dist, install("Diag"),   sexp_diag_val);
  setAttrib(sexp_result_dist, install("Upper"),  sexp_upper_val);
  setAttrib(sexp_result_dist, install("Labels"), em->sexp_sample_names);
  
  
  // Avoid allocating pairs_vec for common all-vs-all case
  if (isNull(sexp_pairs_vec)) {
    
    pairs_vec = NULL;
    n_pairs   = n_dist;
    
  } else {
    
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
    
    for (int i = 0; i < n_dist; i++)
      dist_vec[i] = NA_REAL;
    
    if (n_pairs == 0) {
      free_all();
      UNPROTECT(5);
      return sexp_result_dist;
    }
  }
  
  
  run_parallel(bdiv_func, n_threads, n_pairs);
  
  free_all();
  UNPROTECT(5);
  return sexp_result_dist;
}



//======================================================
// A dist object of n_samples samples, labelled as `em`'s, with every
// distance NA when only some pairs will be computed. Unprotected.
//======================================================
static SEXP part_dist(ecomatrix_t *em, int na) {

  SEXP d = PROTECT(allocVector(REALSXP, n_dist));
  if (na) for (int i = 0; i < n_dist; i++) REAL(d)[i] = NA_REAL;

  SEXP cls = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(cls, 0, mkChar("ecodive_dist"));
  SET_STRING_ELT(cls, 1, mkChar("dist"));
  setAttrib(d, R_ClassSymbol,     cls);
  setAttrib(d, install("Size"),   PROTECT(ScalarInteger(n_samples)));
  setAttrib(d, install("Diag"),   PROTECT(ScalarLogical(0)));
  setAttrib(d, install("Upper"),  PROTECT(ScalarLogical(0)));
  setAttrib(d, install("Labels"), em->sexp_sample_names);

  UNPROTECT(5);
  return d;
}


// The table of beta_partition() and beta_partition_multi(): read, then
// normalized, and its arrays set as this file's globals.
static ecomatrix_t *part_table(SEXP sexp_otu_mtx, SEXP sexp_margin,
                               SEXP sexp_norm, int n_threads) {
  int norm = asInteger(sexp_norm);
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  if (norm) normalize(em, norm, n_threads, 0);
  n_samples = em->n_samples;
  n_otus    = em->n_otus;
  pos_vec   = em->pos_vec;
  otu_vec   = em->otu_vec;
  val_vec   = em->val_vec;
  clr_vec   = em->clr_vec;
  n_dist    = n_samples * (n_samples - 1) / 2;
  return em;
}


//======================================================
// beta_partition(): the turnover (balanced), nestedness (gradient) and
// total dissimilarity matrices of family `index`, in one pass over the
// pairs. Returns a list of the three dist objects.
//======================================================
SEXP C_beta_partition(
    SEXP sexp_index,     SEXP sexp_otu_mtx,  SEXP sexp_margin,
    SEXP sexp_norm,      SEXP sexp_pairs_vec, SEXP sexp_n_threads ) {

  int n_threads = asInteger(sexp_n_threads);
  init_n_ptrs(16);
  ecomatrix_t *em = part_table(sexp_otu_mtx, sexp_margin, sexp_norm, n_threads);

  pthread_func_t func = NULL;
  switch (asInteger(sexp_index)) {
    case PART_SORENSEN: func = sorensen_part; break;
    case PART_JACCARD:  func = jaccard_part;  break;
    case PART_BRAY:     func = bray_part;     break;
    case PART_RUZICKA:  func = soergel_part;  break;
  }
  if (func == NULL) { // # nocov start
    free_all();
    error("Invalid beta partition index.");
  } // # nocov end

  int some = !isNull(sexp_pairs_vec);

  SEXP res = PROTECT(allocVector(VECSXP, 3));
  for (int k = 0; k < 3; k++) SET_VECTOR_ELT(res, k, part_dist(em, some));
  dist_vec  = REAL(VECTOR_ELT(res, 0));
  dist_vec2 = REAL(VECTOR_ELT(res, 1));
  dist_vec3 = REAL(VECTOR_ELT(res, 2));
  part_mode = 3;

  if (some) {
    pairs_vec = INTEGER(sexp_pairs_vec);
    n_pairs   = LENGTH(sexp_pairs_vec);
  } else {
    pairs_vec = NULL;
    n_pairs   = n_dist;
  }

  if (n_pairs > 0) run_parallel(func, n_threads, n_pairs);

  free_all();
  UNPROTECT(1);
  return res;
}


//======================================================
// beta_partition_multi(): the multiple-site turnover (balanced),
// nestedness (gradient) and total dissimilarity of family `index`. With
// n_sites = 0, of all the sites: a 1 x 3 matrix. Otherwise, of n_draws
// random sets of n_sites sites: an n_draws x 3 matrix.
//======================================================
SEXP C_beta_multi(
    SEXP sexp_index,   SEXP sexp_otu_mtx, SEXP sexp_margin,
    SEXP sexp_norm,    SEXP sexp_n_sites, SEXP sexp_n_draws,
    SEXP sexp_seed,    SEXP sexp_n_threads ) {

  int n_threads  = asInteger(sexp_n_threads);
  multi_index    = asInteger(sexp_index);
  multi_n_sites  = asInteger(sexp_n_sites);
  multi_n_draws  = multi_n_sites ? asInteger(sexp_n_draws) : 1;
  multi_seed     = asInteger(sexp_seed);
  part_incidence = multi_index == PART_SORENSEN || multi_index == PART_JACCARD;

  init_n_ptrs(16);
  part_table(sexp_otu_mtx, sexp_margin, sexp_norm, n_threads);

  SEXP res  = PROTECT(allocMatrix(REALSXP, multi_n_draws, 3));
  multi_out = REAL(res);

  if (multi_n_sites == 0) {

    // All sites: the rows' sums in parallel, then summed in row order.
    double *mins  = (double *) safe_malloc((size_t) n_samples * sizeof(double));
    double *maxs  = (double *) safe_malloc((size_t) n_samples * sizeof(double));
    int    *sites = (int *)    safe_malloc((size_t) n_samples * sizeof(int));
    double *fmax  = (double *) safe_malloc((size_t) n_otus    * sizeof(double));
    for (int i = 0; i < n_samples; i++) sites[i] = i;
    dist_vec  = mins;
    dist_vec2 = maxs;
    multi_all = sites;
    run_parallel(multi_rows, n_threads, n_samples - 1);

    double sum_min = 0, sum_max = 0;
    for (int i = 0; i < n_samples - 1; i++) { sum_min += mins[i]; sum_max += maxs[i]; }

    multi_values(multi_index, multi_shared(sites, n_samples, fmax), sum_min, sum_max, multi_out);
  }

  else {

    // One scratch block per worker.
    int    n_workers = n_threads < 1 ? 1 : n_threads;
    size_t stride    = (size_t) n_otus + (size_t) n_samples;
    multi_scratch    = (double *) safe_malloc((size_t) n_workers * stride * sizeof(double));
    run_parallel(multi_draws, n_workers, multi_n_draws);
  }

  free_all();
  UNPROTECT(1);
  return res;
}
