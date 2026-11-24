// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The metrics from "scikit-bio ports" below (dominance, gini, goods_coverage,
// heip, hill, mcintosh_evenness, pielou, renyi, simpson_evenness, strong,
// tsallis) and chao1()'s bias correction are C ports of the formulas and
// edge-case branches of skbio.diversity.alpha (scikit-bio 0.7.4,
// skbio/diversity/alpha/_base.py, _gini.py and _chao1.py), under this
// licence:
//
// Copyright (c) 2013--, scikit-bio development team.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// * Neither the names scikit-bio, skbio, or biocore nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ecodive.h"
#include <Rmath.h> // lgammafn

#define ADIV_ACE          1
#define ADIV_BERGER       2
#define ADIV_BRILLOUIN    3
#define ADIV_CHAO1        4
#define ADIV_FAITH        5
#define ADIV_FISHER       6
#define ADIV_INV_SIMPSON  7
#define ADIV_MARGALEF     8
#define ADIV_MCINTOSH     9
#define ADIV_MENHINICK   10
#define ADIV_OBSERVED    11
#define ADIV_SHANNON     12
#define ADIV_SIMPSON     13
#define ADIV_SQUARES     14
#define ADIV_DOMINANCE   15
#define ADIV_GINI        16
#define ADIV_GOODS       17
#define ADIV_HEIP        18
#define ADIV_HILL        19
#define ADIV_MCINTOSH_E  20
#define ADIV_PIELOU      21
#define ADIV_RENYI       22
#define ADIV_SIMPSON_E   23
#define ADIV_STRONG      24
#define ADIV_TSALLIS     25

static int     n_samples;
static int    *pos_vec;
static int    *otu_vec;
static double *val_vec;
static SEXP   *sexp_extra;
static double *result_vec;
static int     n_otus;

/*
 * FOREACH_SAMPLE iterates over all samples, ensuring that each 
 * is assigned to only a single thread. Provides `sample`, `nnz`, 
 * and `result` initialized to 0; expects `result` at the end.
 */
#define FOREACH_SAMPLE(expression)                             \
  do {                                                         \
    int sample    = ((worker_t *)arg)->i;                      \
    int n_threads = ((worker_t *)arg)->n;                      \
                                                               \
    for (; sample < n_samples; sample += n_threads) {          \
      double *val_begin = val_vec + pos_vec[sample];           \
      double *val_end   = val_vec + pos_vec[sample + 1];       \
      int     nnz       = val_end - val_begin;                 \
      double result     = 0;                                   \
      if (nnz) {                                               \
        expression;                                            \
      }                                                        \
      else {                                                   \
        result = NA_REAL;                                      \
      }                                                        \
      result_vec[sample] = result;                             \
    }                                                          \
  } while (0)

/*
 * FOREACH_VAL iterates over each value for the current sample.
 * Sets val to point to each value in val_vec in turn.
 */
#define FOREACH_VAL(expression)                                \
  do {                                                         \
    double *val = val_begin;                                   \
    for (; val != val_end; val++) {                            \
      expression;                                              \
    }                                                          \
  } while (0)



//======================================================
// Abundance-based Coverage Estimator (Chao & Lee 1992)
//======================================================

static int     ace_cutoff;
static double *ace_rare_nnz_k_mtx;

static void *ace(void *arg) {
  
  int     thread_i       = ((worker_t *)arg)->i;
  int     cutoff         = ace_cutoff;
  double *rare_nnz_k_vec = ace_rare_nnz_k_mtx + (thread_i * cutoff);
  
  FOREACH_SAMPLE(
    
    double abund_nnz = 0;
    double rare_sum  = 0;
    double rare_nnz  = 0;
    
    memset(rare_nnz_k_vec, 0, cutoff * sizeof(double));
    
    FOREACH_VAL(
      int x_int = (int)(ceil(*val));
      if (x_int < cutoff) {
        rare_sum += x_int;
        rare_nnz++;
        rare_nnz_k_vec[x_int]++;
      }
      else {
        abund_nnz++;
      }
    ); // FOREACH_VAL
    
    for (int k = 1; k < cutoff; k++)
      result += k * (k - 1) * rare_nnz_k_vec[k];
    
    double p = (1 - (rare_nnz_k_vec[1] / rare_sum));
    
    result = result * rare_nnz/(p * rare_sum * (rare_sum - 1)) - 1;
    if (result < 0) result = 0;
    result = abund_nnz + rare_nnz/p + result * rare_nnz_k_vec[1]/p;
    
  ); // FOREACH_SAMPLE
  
  return NULL;
}

static pthread_func_t ace_setup(int n_threads) {
  
  ace_cutoff         = asInteger(*sexp_extra) + 1;
  otu_vec            = maybe_free_one(otu_vec);
  ace_rare_nnz_k_mtx = safe_malloc(n_threads * ace_cutoff * sizeof(double));
  
  return ace;
}
  
  
//======================================================
// Berger-Parker (Berger & Parker 1970)
// max(x) / sum(x)
//======================================================
static void *berger(void *arg) {
  FOREACH_SAMPLE(
    FOREACH_VAL(
      if (*val > result) result = *val;
    );
  );
  
  return NULL;
}
  
  
//======================================================
// Brillouin (Brillouin 1956)
// (log(sum(x)!) - sum(log(x!))) / sum(x)
// note: lgamma(x + 1) == log(x!). R's lgammafn() keeps no global state (C's
// lgamma() writes signgam, a data race across threads) and never warns
// for x >= 1.
//======================================================
static void *brillouin(void *arg) {
  FOREACH_SAMPLE(
    
    double depth = 0;
  
    FOREACH_VAL(
      depth  += *val;
      result += lgammafn(*val + 1);
    );
    
    result = (lgammafn(depth + 1) - result) / depth;
  );
  
  return NULL;
}
  


//======================================================
// Chao1 (Chao 1984, 1987)
// bias-corrected: sum(x>0) + F1 * (F1 - 1) / (2 * (F2 + 1))
// classic:        sum(x>0) + F1^2 / (2 * F2), when F1 > 0 and F2 > 0
// As scikit-bio's chao1(), the classic form falls back on the
// bias-corrected one when there are no singletons or no doubletons.
//======================================================
static int chao1_bias_corrected;

static void *chao1(void *arg) {
  
  int bias_corrected = chao1_bias_corrected;
  
  FOREACH_SAMPLE(
    
    double ones = 0;
    double twos = 0;
    
    FOREACH_VAL(
      if      (*val == 1) ones++;
      else if (*val == 2) twos++;
    );
    
    if (!bias_corrected && ones && twos)
      result = nnz + (ones * ones) / (twos * 2);
    else
      result = nnz + ones * (ones - 1) / (2 * (twos + 1));
  );
  
  return NULL;
}

static pthread_func_t chao1_setup(void) {
  chao1_bias_corrected = asLogical(*sexp_extra);
  return chao1;
}



//======================================================
// Faith's Phylogenetic Diversity
// Find the edges that an otu passes through.
//======================================================

static ecotree_t *faith_et;
static char      *faith_has_edge_mtx;

static void *faith(void *arg) {
  
  int     thread_i     = ((worker_t *)arg)->i;
  int     n_edges      = faith_et->n_edges;
  node_t *node_vec     = faith_et->node_vec;
  char   *has_edge_vec = faith_has_edge_mtx + (thread_i * n_edges);
  
  FOREACH_SAMPLE(
    
    memset(has_edge_vec, 0, n_edges * sizeof(char));
    
    int *otu_begin = otu_vec + pos_vec[sample];
    int *otu_end   = otu_vec + pos_vec[sample + 1];
    
    for (int *otu = otu_begin; otu != otu_end; otu++) {
      
      int node_i = *otu;    // start at OTU tip/leaf in tree
      while (node_i > -1) { // traverse until we hit the tree's root
        
        node_t *node   = node_vec + node_i;
        char *has_edge = has_edge_vec + node->edge;
        
        if (*has_edge) break; // already traversed
        *has_edge = 1;
        
        result += node->length;
        node_i  = node->parent;
      }
    }
    
  );
  
  return NULL;
}

static pthread_func_t faith_setup(int n_threads) {
  
  faith_et           = new_ecotree(*sexp_extra);
  faith_has_edge_mtx = safe_malloc(n_threads * faith_et->n_edges * sizeof(char));
  
  return faith;
}


//======================================================
// Fisher's diversity index (Fisher 1943)
// otus = fisher * log(1 + depth/fisher)
//======================================================
static void *fisher(void *arg) {
  
  int    digits = asInteger(*sexp_extra);
  double mult   = pow(10, digits);
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    FOREACH_VAL(depth += *val);
    
    if (depth == nnz) {
      result = R_PosInf;  // All singletons -> infinite loop
    }
    else {
      
      double alpha;
      double lo = 2;
      double hi = 16;
      
      // Sometimes the result will be less than 2 or greater than 16.
      while (lo * log(1 + depth/lo) > nnz) { hi = lo; lo /= 2; }
      while (hi * log(1 + depth/hi) < nnz) { lo = hi; hi *= 2; }
      
      // Check if range has converged to same value after rounding.
      while (round(hi * mult) != round(lo * mult)) {
        
        // This loop's guess for the alpha term.
        alpha = (lo + hi) / 2;
        
        // Update the range we need to examine.
        if (alpha * log(1 + depth/alpha) > nnz) { hi = alpha; }
        else                                    { lo = alpha; }
        
      }
      
      result = round(hi * mult) / mult;
    }
  );
  
  return NULL;
}


//======================================================
// Inverse Simpson
// p <- x / sum(x)
// 1 / sum(p ** 2)
//======================================================
static void *inv_simpson(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * *val);
  
    result = 1 / result;
  );
  
  return NULL;
}


//======================================================
// Margalef (Margalef 1958)
// (sum(x > 0) - 1) / log(sum(x))
//======================================================
static void *margalef(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    
    FOREACH_VAL(depth += *val);
    
    result = (nnz - 1) / log(depth);
  );
  
  return NULL;
}


//======================================================
// McIntosh (McIntosh 1967)
// (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))
//======================================================
static void *mcintosh(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    
    FOREACH_VAL(
      depth  += *val;
      result += *val * *val;
    );
    
    result = (depth - sqrt(result)) / (depth - sqrt(depth));
  );
  
  return NULL;
}


//======================================================
// Menhinick (Menhinick 1964)
// sum(x > 0) / sqrt(sum(x))
//======================================================
static void *menhinick(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth = 0;
  
    FOREACH_VAL(depth += *val);
    
    result = nnz / sqrt(depth);
  );
  
  return NULL;
}


//======================================================
// Observed Features
// sum(x > 0)
//======================================================
static void *observed(void *arg) {
  
  FOREACH_SAMPLE(result = nnz);
  
  return NULL;
}


//======================================================
// Shannon (Shannon 1948)
// p <- x / sum(x)
// -sum(p * log(p))
//======================================================
static void *shannon(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * log(*val));
  
    result *= -1;
  );
  
  return NULL;
}


//======================================================
// Simpson (Simpson 1949)
// p <- x / sum(x)
// 1 - sum(p ** 2)
//======================================================
static void *simpson(void *arg) {
  
  FOREACH_SAMPLE(
    
    FOREACH_VAL(result += *val * *val);
  
    result = 1 - result;
  );
  
  return NULL;
}


//======================================================
// Squares Estimator (Alroy 2018)
// N = sum(x)      # sampling depth
// S = sum(x > 0)  # number of non-zero OTUs
// F1 = sum(x == 1) # singletons
// ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S)) + S
//======================================================
static void *squares(void *arg) {
  
  FOREACH_SAMPLE(
    
    double depth      = 0;
    double singletons = 0;
    
    FOREACH_VAL(
      depth  += *val;
      result += *val * *val;
      if (*val == 1) singletons++;
    );
    
    double denominator = (depth * depth) - (singletons * nnz);
    
    if (denominator == 0) {
      result = R_PosInf; // All singletons
    } else {
      result *= (singletons * singletons);
      result /= denominator;
      result += nnz;
    }
  );
  
  return NULL;
}




//======================================================
// scikit-bio ports (see the licence at the top of this file).
// Proportion-based metrics are called with norm = 'percent', so `val`
// holds p_i = x_i / sum(x); goods_coverage and mcintosh_evenness with
// norm = 'none', so `val` holds counts. strong and gini also take
// norm = 'none' and divide cumulative sums by the total, as scikit-bio.
//======================================================


// Simpson's Dominance: sum(p^2)
static void *dominance(void *arg) {
  FOREACH_SAMPLE(
    FOREACH_VAL(result += *val * *val);
  );
  return NULL;
}


// Good's Coverage (Good 1953): 1 - F1 / N
static void *goods(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    double ones  = 0;
    FOREACH_VAL(
      depth += *val;
      if (*val == 1) ones++;
    );
    result = 1 - ones / depth;
  );
  return NULL;
}


// Shannon entropy of the current sample's proportions.
#define SHANNON_H(H)                                           \
  double H = 0;                                                \
  FOREACH_VAL(H += -*val * log(*val));


// Heip's Evenness (Heip 1974): (exp(H) - 1) / (S - 1); 1 when S = 1
static void *heip(void *arg) {
  FOREACH_SAMPLE(
    if (nnz == 1) {
      result = 1;
    } else {
      SHANNON_H(H);
      result = (exp(H) - 1) / (nnz - 1);
    }
  );
  return NULL;
}


// Pielou's Evenness (Pielou 1966): H / log(S); 1 when S = 1
static void *pielou(void *arg) {
  FOREACH_SAMPLE(
    if (nnz == 1) {
      result = 1;
    } else {
      SHANNON_H(H);
      result = H / log(nnz);
    }
  );
  return NULL;
}


// Simpson's Evenness: 1 / (S * sum(p^2))
static void *simpson_e(void *arg) {
  FOREACH_SAMPLE(
    double D = 0;
    FOREACH_VAL(D += *val * *val);
    result = 1 / (nnz * D);
  );
  return NULL;
}


// McIntosh's Evenness (Heip & Engels 1974):
// sqrt(sum(x^2)) / sqrt((N - S + 1)^2 + S - 1)
static void *mcintosh_e(void *arg) {
  FOREACH_SAMPLE(
    double depth = 0;
    double sumsq = 0;
    FOREACH_VAL(
      depth += *val;
      sumsq += *val * *val;
    );
    double a = depth - nnz + 1;
    result = sqrt(sumsq) / sqrt(a * a + nnz - 1);
  );
  return NULL;
}


// Hill numbers, Renyi and Tsallis entropies of order q.
static double adiv_q;

// Hill (1973): sum(p^q)^(1/(1-q)); exp(H) at q = 1; 1/max(p) at q = Inf
static void *hill(void *arg) {
  double q = adiv_q;
  FOREACH_SAMPLE(
    if (q == 1) {
      SHANNON_H(H);
      result = exp(H);
    } else if (q == R_PosInf) {
      FOREACH_VAL(if (*val > result) result = *val);
      result = 1 / result;
    } else {
      FOREACH_VAL(result += pow(*val, q));
      result = pow(result, 1 / (1 - q));
    }
  );
  return NULL;
}

// Renyi (1961): log(sum(p^q)) / (1 - q); log(S) at 0, H at 1,
// -log(max(p)) at Inf; 0 when S = 1
static void *renyi(void *arg) {
  double q = adiv_q;
  FOREACH_SAMPLE(
    if (nnz == 1) {
      result = 0;
    } else if (q == 0) {
      result = log(nnz);
    } else if (q == 1) {
      SHANNON_H(H);
      result = H;
    } else if (q == R_PosInf) {
      FOREACH_VAL(if (*val > result) result = *val);
      result = -log(result);
    } else {
      FOREACH_VAL(result += pow(*val, q));
      result = log(result) / (1 - q);
    }
  );
  return NULL;
}

// Tsallis (1988): (1 - sum(p^q)) / (q - 1); H at 1, 0 at Inf; 0 when S = 1
static void *tsallis(void *arg) {
  double q = adiv_q;
  FOREACH_SAMPLE(
    if (nnz == 1 || q == R_PosInf) {
      result = 0;
    } else if (q == 1) {
      SHANNON_H(H);
      result = H;
    } else {
      FOREACH_VAL(result += pow(*val, q));
      result = (1 - result) / (q - 1);
    }
  );
  return NULL;
}

static pthread_func_t q_setup(pthread_func_t func) {
  adiv_q = asReal(*sexp_extra);
  return func;
}


// Strong's Dominance and the Gini index sort each sample's values into a
// buffer of their own thread's.
static double *sort_buf_mtx;
static int     sort_buf_len;

static int cmp_increasing(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static int cmp_decreasing(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) - (x > y);
}

// Strong (2002): max_i(c_i / N - i / S), c_i the cumulative sums of the
// values sorted in decreasing order.
static void *strong(void *arg) {
  
  double *buf = sort_buf_mtx + (size_t)((worker_t *)arg)->i * sort_buf_len;
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    FOREACH_VAL(depth += *val);
    
    memcpy(buf, val_begin, nnz * sizeof(double));
    qsort(buf, nnz, sizeof(double), cmp_decreasing);
    
    double cumsum = 0;
    result = R_NegInf;
    for (int i = 0; i < nnz; i++) {
      cumsum  += buf[i];
      double d = cumsum / depth - (double)(i + 1) / nnz;
      if (d > result) result = d;
    }
  );
  
  return NULL;
}

// Gini index, by rectangles: max(0, 1 - 2B), B = mean(L_j) over all K
// features of the table, L the Lorenz curve (cumulative shares of the
// values sorted increasing). Absent features sort first and add nothing.
static void *gini(void *arg) {
  
  double *buf = sort_buf_mtx + (size_t)((worker_t *)arg)->i * sort_buf_len;
  double  dx  = 1 / (double)n_otus;
  
  FOREACH_SAMPLE(
    
    double depth = 0;
    FOREACH_VAL(depth += *val);
    
    memcpy(buf, val_begin, nnz * sizeof(double));
    qsort(buf, nnz, sizeof(double), cmp_increasing);
    
    double cumsum = 0;
    double B      = 0;
    for (int i = 0; i < nnz; i++) {
      cumsum += buf[i];
      B      += cumsum / depth;
    }
    B     *= dx;
    result = 1 - 2 * B;
    if (result < 0) result = 0;
  );
  
  return NULL;
}

static pthread_func_t sort_setup(int n_threads, pthread_func_t func) {
  
  sort_buf_len = 1;
  for (int i = 0; i < n_samples; i++) {
    int nnz = pos_vec[i + 1] - pos_vec[i];
    if (nnz > sort_buf_len) sort_buf_len = nnz;
  }
  
  sort_buf_mtx = safe_malloc((size_t)n_threads * sort_buf_len * sizeof(double));
  
  return func;
}



//======================================================
// R interface. Distributes work across threads.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_algorithm, SEXP sexp_otu_mtx, 
    SEXP sexp_margin,    SEXP sexp_norm, 
    SEXP sexp_n_threads, SEXP sexp_extra_args ) {
  
  init_n_ptrs(10);
  
  int norm       = asInteger(sexp_norm);
  int n_threads  = asInteger(sexp_n_threads);
  sexp_extra     = &sexp_extra_args;
  
  ecomatrix_t *em = new_ecomatrix(sexp_otu_mtx, sexp_margin);
  if (norm) normalize(em, norm, n_threads, 0);
  
  n_samples = em->n_samples;
  pos_vec   = em->pos_vec;
  otu_vec   = em->otu_vec;
  val_vec   = em->val_vec;
  n_otus    = em->n_otus;
  
  
  // function to run
  // void * (*adiv_func)(void *) = NULL;
  pthread_func_t adiv_func = NULL;
  switch (asInteger(sexp_algorithm)) {
    case ADIV_ACE:         adiv_func = ace_setup(n_threads);   break;
    case ADIV_BERGER:      adiv_func = berger;                 break;
    case ADIV_BRILLOUIN:   adiv_func = brillouin;              break;
    case ADIV_CHAO1:       adiv_func = chao1_setup();          break;
    case ADIV_FAITH:       adiv_func = faith_setup(n_threads); break;
    case ADIV_FISHER:      adiv_func = fisher;                 break;
    case ADIV_INV_SIMPSON: adiv_func = inv_simpson;            break;
    case ADIV_MARGALEF:    adiv_func = margalef;               break;
    case ADIV_MCINTOSH:    adiv_func = mcintosh;               break;
    case ADIV_MENHINICK:   adiv_func = menhinick;              break;
    case ADIV_OBSERVED:    adiv_func = observed;               break;
    case ADIV_SHANNON:     adiv_func = shannon;                break;
    case ADIV_SIMPSON:     adiv_func = simpson;                break;
    case ADIV_SQUARES:     adiv_func = squares;                break;
    case ADIV_DOMINANCE:   adiv_func = dominance;              break;
    case ADIV_GINI:        adiv_func = sort_setup(n_threads, gini);   break;
    case ADIV_GOODS:       adiv_func = goods;                  break;
    case ADIV_HEIP:        adiv_func = heip;                   break;
    case ADIV_HILL:        adiv_func = q_setup(hill);          break;
    case ADIV_MCINTOSH_E:  adiv_func = mcintosh_e;             break;
    case ADIV_PIELOU:      adiv_func = pielou;                 break;
    case ADIV_RENYI:       adiv_func = q_setup(renyi);         break;
    case ADIV_SIMPSON_E:   adiv_func = simpson_e;              break;
    case ADIV_STRONG:      adiv_func = sort_setup(n_threads, strong); break;
    case ADIV_TSALLIS:     adiv_func = q_setup(tsallis);       break;
  }
  
  if (adiv_func == NULL) { // # nocov start
    error("Invalid alpha diversity algorithm.");
    return R_NilValue;
  } // # nocov end
  
  
  // Create the diversity vector to return
  SEXP sexp_result_vec = PROTECT(allocVector(REALSXP, n_samples));
  result_vec           = REAL(sexp_result_vec);
  setAttrib(sexp_result_vec, R_NamesSymbol, em->sexp_sample_names);
  
  
  run_parallel(adiv_func, n_threads, n_samples);
  
  free_all();
  UNPROTECT(1);
  return sexp_result_vec;
}
