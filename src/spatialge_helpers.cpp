//


#include <Rcpp.h>
// #include <iostream>
using namespace Rcpp;

// computeSubsampleSums
// Calculates sum of distances among spots or cells randomly selected
// Arguments:
// coords is a matrix with x locations in the first column and y locations in the second
// n_subsample is the number of cells that are positive or above a thresold
// n_samples is the number of permutations ot perform to estimate the distribution
// Code by Dr. Alex "The Lab Warrior" Soupir, slightly modified
// [[Rcpp::export]]
NumericVector computeSubsampleSums(
    NumericMatrix coords,
    int n_subsample,
    int n_samples) {

  NumericVector x_coords = coords.column(0);
  NumericVector y_coords = coords.column(1);
  int N = coords.nrow();

  // vector for ith cell
  IntegerVector all_indices(N);
  std::iota(all_indices.begin(), all_indices.end(), 0);

  NumericVector sums(n_samples);

  for (int s = 0; s < n_samples; ++s) {
    // dont select same cell twice for permutations
    IntegerVector sample_indices = Rcpp::sample(all_indices, n_subsample, false);
    
    double total_dist = 0;
    int ns = sample_indices.size();

    // sum i,j only, not j,i/ upper triangle
    for (int i = 0; i < ns; ++i) {
      int idx_i = sample_indices[i];
      double xi = x_coords[idx_i], yi = y_coords[idx_i];

      for (int j = i + 1; j < ns; ++j) {
        int idx_j = sample_indices[j];
        double dx = xi - x_coords[idx_j];
        double dy = yi - y_coords[idx_j];

        total_dist += std::sqrt((dx * dx) + (dy * dy));
        //std::cout << total_dist << std::endl;
      }
    }

    sums[s] = total_dist;
  }
  return sums;
}

