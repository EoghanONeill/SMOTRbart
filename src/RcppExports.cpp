// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// phi_app
arma::mat phi_app(arma::mat X_stand, arma::mat anc, double tau);
RcppExport SEXP _SMOTRbart_phi_app(SEXP X_standSEXP, SEXP ancSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_stand(X_standSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_app(X_stand, anc, tau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMOTRbart_phi_app", (DL_FUNC) &_SMOTRbart_phi_app, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMOTRbart(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
