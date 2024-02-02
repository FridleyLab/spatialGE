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
#include <Rcpp.h>
using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SparseRowVar2
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat, NumericVector mu
//, bool display_progress
);
RcppExport SEXP _SeuratMod_SparseRowVar2(SEXP matSEXP, SEXP muSEXP 
//, SEXP display_progressSEXP
) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    //Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(SparseRowVar2(mat, mu 
    //, display_progress
    ));
    return rcpp_result_gen;
END_RCPP
}

// SparseRowVarStd
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sd, double vmax
//, bool display_progress
);
RcppExport SEXP _SeuratMod_SparseRowVarStd(SEXP matSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP vmaxSEXP
//, SEXP display_progressSEXP
) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type vmax(vmaxSEXP);
    //Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(SparseRowVarStd(mat, mu, sd, vmax
    //, display_progress
    ));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
	//{"_SeuratMod_SparseRowVar2", (DL_FUNC) &_SeuratMod_SparseRowVar2, 3},
	//{"_SeuratMod_SparseRowVarStd", (DL_FUNC) &_SeuratMod_SparseRowVarStd, 5},
	{"_SeuratMod_SparseRowVar2", (DL_FUNC) &_SeuratMod_SparseRowVar2, 2},
	{"_SeuratMod_SparseRowVarStd", (DL_FUNC) &_SeuratMod_SparseRowVarStd, 4},
	{NULL, NULL, 0}
};

RcppExport void R_init_anRpackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

