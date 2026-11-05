// Moran's I and Geary's C over a distance-weighted COMPLETE graph, computed
// without materializing the graph.
//
// WHY THIS EXISTS
//
// create_listw_from_dist() (R/utils.R) builds a listw in which every spot
// neighbours every other, as an R list of named vectors. Measured: 1.1 MB at 120
// spots, 1216.8 MB at 4,096 -- 1,100x the memory for 34x the spots, because each
// neighbour vector carries character names, and membership is computed with an
// O(n) string comparison inside an O(n) loop. Extrapolated to 50,000 cells it is
// roughly 181 GB, so SThet()'s DEFAULT weighting simply cannot run on
// Xenium/CosMx data.
//
// The graph itself is never needed. Both statistics reduce to three quantities:
//
//   s_i = sum_j w_ij                  (per-spot weight row sum; gene-independent)
//   S0  = sum_i s_i                   (Szero)
//   q_g = z_g' W z_g                  (per gene)
//
// because W is symmetric with a zero diagonal, so
//
//   sum_ij w_ij (z_i - z_j)^2 = 2 * sum_i z_i^2 s_i  -  2 * q
//
// giving Geary's numerator from the same q as Moran's. So one pass over the
// (implicit) distance matrix yields every gene's statistic, in O(n) memory
// instead of O(n^2).
//
// MATCHING spdep EXACTLY
//
// spdep::moran is  I = (n/S0) * sum(z * Wz) / sum(z^2), with z = x - mean(x).
// spdep::geary is  C = (n1/(2*S0)) * sum_ij w_ij (z_i-z_j)^2 / sum(z^2), and it
// standardizes via scale(x) rather than merely centering -- but both C and K are
// invariant to that scaling (numerator and zz each pick up the same sd^2 factor),
// so centered input reproduces both.
//
// Weights are w_ij = dmax / d_ij for i != j, where dmax is the largest pairwise
// distance. That is exactly what the R version computes: it forms
// as.matrix(dist(coords)), divides by max(), sets the diagonal aside, then
// inverts. Coincident spots give d_ij == 0 and therefore an infinite weight; that
// is preserved rather than "fixed", because the R version does the same and
// changing it would move published numbers.

#include <Rcpp.h>
#include <cmath>
#include <vector>

// [[Rcpp::export(rng = false)]]
Rcpp::List computeMoranGeary(const Rcpp::NumericMatrix& coords,
                             const Rcpp::NumericMatrix& zt,
                             bool want_moran,
                             bool want_geary) {
  const int n = coords.nrow();
  const int ndim = coords.ncol();
  const int G = zt.ncol();

  if (zt.nrow() != n) {
    Rcpp::stop("computeMoranGeary: zt must have one row per spot (n = %d), got %d.",
               n, zt.nrow());
  }
  if (n < 2) {
    Rcpp::stop("computeMoranGeary: need at least two spots.");
  }

  // Local copies with predictable layout. coords is column-major from R, so
  // reading coordinate k of spot i as cx[k*n + i] is contiguous in i.
  const double* cx = &coords[0];
  const double* z = &zt[0];  // zt is n x G column-major: z[g*n + j] is spot j, gene g

  // ---- pass 1: largest pairwise distance -----------------------------------
  // Needed before any weight can be formed, since weights are normalized by it.
  // Gene-independent and cheap relative to pass 2 (no inner gene loop).
  double dmax = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double acc = 0.0;
      for (int k = 0; k < ndim; ++k) {
        const double diff = cx[k * n + i] - cx[k * n + j];
        acc += diff * diff;
      }
      const double d = std::sqrt(acc);
      if (d > dmax) dmax = d;
    }
  }

  // ---- pass 2: accumulate s, S0 and q_g ------------------------------------
  // q is accumulated in place rather than storing Wz, so memory is O(n + G)
  // rather than O(n * G) -- 80 MB saved at n = 50,000 with 200 genes.
  // long double accumulators throughout, to match R. Base sum() and spdep's C
  // code both accumulate in LDOUBLE, and at n = 1000 a plain double S0 already
  // diverged from Szero(listw) by 1.1e-12 relative -- which then propagated into
  // every gene's statistic, since both I and C are scaled by 1/S0.
  std::vector<double> w(n);
  std::vector<long double> s(n, 0.0L);
  std::vector<long double> q(G, 0.0L);
  std::vector<long double> zsq_s(G, 0.0L);

  // S0 is accumulated in ONE continuous long-double running total, in exactly the
  // order spdep visits the weights -- Szero(listw) is sum(unlist(listw$weights)),
  // i.e. row 1's weights in increasing j, then row 2's, and so on. Summing
  // per-row totals and then adding those totals gives a different grouping and
  // diverged from Szero by 1.1e-12 relative at n = 1000, which propagated
  // straight into every gene's statistic since both I and C scale by 1/S0.
  long double S0l = 0.0L;

  for (int i = 0; i < n; ++i) {
    long double s_i = 0.0L;
    for (int j = 0; j < n; ++j) {
      if (j == i) { w[j] = 0.0; continue; }
      double acc = 0.0;
      for (int k = 0; k < ndim; ++k) {
        const double diff = cx[k * n + i] - cx[k * n + j];
        acc += diff * diff;
      }
      const double d = std::sqrt(acc) / dmax;
      const double wij = 1.0 / d;   // d == 0 for coincident spots -> Inf, as in R
      w[j] = wij;
      s_i += wij;
      S0l += wij;
    }
    s[i] = s_i;
    const double s_i_d = (double)s_i;

    for (int g = 0; g < G; ++g) {
      const double* zg = z + (std::size_t)g * n;
      long double dot = 0.0L;
      for (int j = 0; j < n; ++j) dot += (long double)w[j] * zg[j];
      q[g] += (long double)zg[i] * dot;
      zsq_s[g] += (long double)zg[i] * zg[i] * s_i_d;
    }
  }

  const double S0 = (double)S0l;

  // ---- assemble per-gene statistics ---------------------------------------
  Rcpp::NumericVector moran(want_moran ? G : 0);
  Rcpp::NumericVector geary(want_geary ? G : 0);
  Rcpp::NumericVector kurt(G);

  const double dn = (double)n;
  for (int g = 0; g < G; ++g) {
    const double* zg = z + (std::size_t)g * n;
    long double zzl = 0.0L, z4l = 0.0L;
    for (int j = 0; j < n; ++j) {
      const long double zj2 = (long double)zg[j] * zg[j];
      zzl += zj2;
      z4l += zj2 * zj2;
    }
    const double zz = (double)zzl;
    const double z4 = (double)z4l;
    kurt[g] = (dn * z4) / (zz * zz);

    if (want_moran) {
      moran[g] = (dn / S0) * ((double)q[g] / zz);
    }
    if (want_geary) {
      // sum_ij w_ij (z_i - z_j)^2 == 2*sum_i z_i^2 s_i - 2*q, for symmetric W
      // with zero diagonal.
      const double num = (double)(2.0L * zsq_s[g] - 2.0L * q[g]);
      geary[g] = ((dn - 1.0) / (2.0 * S0)) * (num / zz);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("moran") = moran,
    Rcpp::Named("geary") = geary,
    Rcpp::Named("K") = kurt,
    Rcpp::Named("S0") = S0,
    Rcpp::Named("dmax") = dmax);
}
