// Copyright (c) 2021 Seurat authors
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
// to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or 
// substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.


#include <RcppEigen.h>
#include <progress.hpp>
//#include <cmath>
//#include <unordered_map>
//#include <fstream>
//#include <string>
//#include <Rinternals.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

/* use this if you know the row means */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu,
                            bool display_progress){
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], 2);
    }
    colSum += pow(mu[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* standardize matrix rows using given mean and standard deviation,
   clip values larger than vmax to vmax,
   then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax,
                              bool display_progress){
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating feature variances of standardized and clipped values" << std::endl;
  }
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

