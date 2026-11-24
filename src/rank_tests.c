// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Rank tests shared by ancom() (ancom.c) and lefse() (lefse.c):
//
//   rank_sort_vg     Sorts (value, group) pairs by value.
//   rank_kruskal_h   Kruskal-Wallis H with average ranks and the tie
//                    correction, ties by exact equality, in SciPy's order
//                    of operations (so H is SciPy's to the bit; R's
//                    kruskal.test() agrees to rounding). NaN when every
//                    value is tied.
//   rank_sum_ties    Wilcoxon's rank sum of group 1 among two groups, with
//                    average ranks, and the ties' sum of t^3 - t.

#include "ecodive.h"


// Sorts a[0 .. n - 1] by v: quicksort (median of three), insertion sort
// for short runs. Ties end up in any order, which changes nothing: tied
// values share their average rank, and rank sums of half-integers are
// exact in any order.
void rank_sort_vg (vg_t *a, int n) {
  while (n > 16) {
    int    mid = n / 2;
    vg_t   t;
    if (a[mid].v   < a[0].v)   { t = a[mid];   a[mid]   = a[0];   a[0]   = t; }
    if (a[n - 1].v < a[0].v)   { t = a[n - 1]; a[n - 1] = a[0];   a[0]   = t; }
    if (a[n - 1].v < a[mid].v) { t = a[n - 1]; a[n - 1] = a[mid]; a[mid] = t; }
    double pivot = a[mid].v;
    int i = 0, j = n - 1;
    for (;;) {
      while (a[i].v < pivot) i++;
      while (pivot < a[j].v) j--;
      if (i >= j) break;
      t = a[i]; a[i] = a[j]; a[j] = t;
      i++; j--;
    }
    // Recurse on the smaller part, loop on the larger.
    if (j + 1 < n - j - 1) { rank_sort_vg(a, j + 1); a += j + 1; n -= j + 1; }
    else                   { rank_sort_vg(a + j + 1, n - j - 1); n = j + 1; }
  }
  for (int i = 1; i < n; i++) {
    vg_t x = a[i];
    int  j = i - 1;
    while (j >= 0 && x.v < a[j].v) { a[j + 1] = a[j]; j--; }
    a[j + 1] = x;
  }
}


// Kruskal-Wallis H of the n pairs vg (sorted here, in place), with groups
// 0 .. k - 1 of sizes ng; rsum holds k doubles.
double rank_kruskal_h (vg_t *vg, int n, const int *ng, int k, double *rsum) {

  rank_sort_vg(vg, n);

  // Average ranks of tie runs, each group's rank sum, and sum(t^3 - t),
  // as SciPy's rankdata(method = 'average') and tie correction.
  for (int g = 0; g < k; g++) rsum[g] = 0;
  double T = 0;
  for (int a = 0; a < n; ) {
    int b = a + 1;
    while (b < n && vg[b].v == vg[a].v) b++;
    double r = (a + 1 + b) / 2.0;          // mean of ranks a + 1 .. b
    for (int c = a; c < b; c++) rsum[vg[c].g] += r;
    double t = b - a;
    T += t * t * t - t;
    a = b;
  }

  double N = n;
  double ties = 1 - T / (N * N * N - N);
  if (ties == 0) return R_NaN;             // every value tied

  double ssbn = 0;
  for (int g = 0; g < k; g++) ssbn += (rsum[g] * rsum[g]) / ng[g];

  double h = 12.0 / (N * (N + 1)) * ssbn - 3 * (N + 1);
  return h / ties;
}


// The rank sum of the pairs with g == 1 among the n pairs vg (sorted here,
// in place), with average ranks for ties; writes sum(t^3 - t) to `ties`.
double rank_sum_ties (vg_t *vg, int n, double *ties) {

  rank_sort_vg(vg, n);

  double r1 = 0, T = 0;
  for (int a = 0; a < n; ) {
    int b = a + 1;
    while (b < n && vg[b].v == vg[a].v) b++;
    double r = (a + 1 + b) / 2.0;
    for (int c = a; c < b; c++) if (vg[c].g == 1) r1 += r;
    double t = b - a;
    if (t > 1) T += t * t * t - t;
    a = b;
  }
  *ties = T;
  return r1;
}
