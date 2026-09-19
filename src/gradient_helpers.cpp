// Distance summarization for STgradient's per-sample distance step, computed
// without materializing the n x n distance matrix.
//
// WHY THIS EXISTS
//
// stgradient_dists_one_sample() (R/STgradient.R) built
// as.matrix(stats::dist(coords)) -- a full n x n double matrix, ~128 MB at
// 4,096 spots and ~20 GB at 50,000 cells, so STgradient could not run on
// Xenium/CosMx data at all. Only three things were read from that matrix:
//
//   * min_sample, the minimum pairwise distance over all spots;
//   * for each reference spot, how many reference spots fall in
//     [min_sample*thr_low, min_sample*thr_high] (the "enough neighbours" test);
//   * for each non-reference spot, the min or mean distance to the kept
//     reference spots.
//
// All three are computed here in O(n) memory. Time remains O(n^2) for the
// pairwise minimum, which is inherent to a distance-based neighbourhood with no
// spatial index; the win is that the matrix is never allocated.
//
// BIT-FIDELITY
//
// The euclidean distance is sqrt(sum over dimensions of (a-b)^2), the same
// formula and accumulation order stats::dist uses, so values are bit-identical.
// The mean reproduces R's real_mean() (summary.c) exactly: a two-pass
// corrected mean (sum, divide, add sum-of-residuals/n), not a one-pass sum.
// That is what makes distsumm='avg' bit-identical -- a one-pass long-double
// sum disagrees by ~3.6e-15 relative.

#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace {

inline double euclid(const Rcpp::NumericMatrix& coords, int a, int b) {
  double s = 0.0;
  const int ndim = coords.ncol();
  for (int k = 0; k < ndim; ++k) {
    const double diff = coords(a, k) - coords(b, k);
    s += diff * diff;
  }
  return std::sqrt(s);
}

}  // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::List computeGradientDists(const Rcpp::NumericMatrix& coords,
                                const Rcpp::IntegerVector& ref_idx,
                                const Rcpp::IntegerVector& nonref_idx,
                                double thr_low,
                                double thr_high,
                                int min_nb,
                                bool want_mean) {
  const int n = coords.nrow();
  const int nref = ref_idx.size();
  const int nnonref = nonref_idx.size();

  for (int i = 0; i < nref; ++i) {
    if (ref_idx[i] < 0 || ref_idx[i] >= n) {
      Rcpp::stop("computeGradientDists: ref_idx out of bounds (%d of %d).",
                 ref_idx[i] + 1, n);
    }
  }
  for (int i = 0; i < nnonref; ++i) {
    if (nonref_idx[i] < 0 || nonref_idx[i] >= n) {
      Rcpp::stop("computeGradientDists: nonref_idx out of bounds (%d of %d).",
                 nonref_idx[i] + 1, n);
    }
  }

  // Minimum pairwise distance over all spots, matching
  // min(dist_tmp[lower.tri(dist_tmp)]). R_PosInf when there is no pair, which
  // makes the neighbour window empty -- the same outcome the R path had.
  double min_sample = R_PosInf;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double d = euclid(coords, i, j);
      if (d < min_sample) min_sample = d;
    }
  }

  // Keep reference spots with at least min_nb reference neighbours inside the
  // tolerance window. Preserves ref_idx order, which is the order the R
  // version's colnames()/colSums() produced.
  const double lo = min_sample * thr_low;
  const double hi = min_sample * thr_high;
  std::vector<int> kept;
  kept.reserve(nref);
  for (int j = 0; j < nref; ++j) {
    int count = 0;
    for (int k = 0; k < nref; ++k) {
      const double d = euclid(coords, ref_idx[j], ref_idx[k]);
      if (d >= lo && d <= hi) ++count;
    }
    if (count >= min_nb) kept.push_back(j);
  }

  Rcpp::NumericVector dist2ref(nnonref, NA_REAL);
  if (nnonref > 1 && !kept.empty()) {
    const std::size_t nk = kept.size();
    const double dn = (double) nk;
    for (int i = 0; i < nnonref; ++i) {
      if (want_mean) {
        // R's mean() (real_mean() in R's summary.c) is a TWO-PASS corrected
        // mean, not sum/n: it sums, divides, then adds the sum of residuals
        // divided by n. Reproducing that here is what keeps distsumm='avg'
        // bit-identical; a one-pass long-double sum does not agree (measured
        // 3.6e-15 relative disagreement on the fixture).
        double s = 0.0;
        for (std::size_t m = 0; m < nk; ++m) {
          s += euclid(coords, nonref_idx[i], ref_idx[kept[m]]);
        }
        s /= dn;
        double t = 0.0;
        for (std::size_t m = 0; m < nk; ++m) {
          t += euclid(coords, nonref_idx[i], ref_idx[kept[m]]) - s;
        }
        s += t / dn;
        dist2ref[i] = s;
      } else {
        double mn = R_PosInf;
        for (std::size_t m = 0; m < nk; ++m) {
          const double d = euclid(coords, nonref_idx[i], ref_idx[kept[m]]);
          if (d < mn) mn = d;
        }
        dist2ref[i] = mn;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("dist2ref") = dist2ref,
    Rcpp::Named("n_kept") = (int) kept.size(),
    Rcpp::Named("min_sample") = min_sample);
}
