// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// validReport
bool validReport(Rcpp::LogicalMatrix report, int n_drugs, int n_events);
RcppExport SEXP _SRSim_validReport(SEXP reportSEXP, SEXP n_drugsSEXP, SEXP n_eventsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type report(reportSEXP);
    Rcpp::traits::input_parameter< int >::type n_drugs(n_drugsSEXP);
    Rcpp::traits::input_parameter< int >::type n_events(n_eventsSEXP);
    rcpp_result_gen = Rcpp::wrap(validReport(report, n_drugs, n_events));
    return rcpp_result_gen;
END_RCPP
}
// validCorrelation
bool validCorrelation(double rho, double p_x, double p_y);
RcppExport SEXP _SRSim_validCorrelation(SEXP rhoSEXP, SEXP p_xSEXP, SEXP p_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< double >::type p_y(p_ySEXP);
    rcpp_result_gen = Rcpp::wrap(validCorrelation(rho, p_x, p_y));
    return rcpp_result_gen;
END_RCPP
}
// sampleRandomRows
Rcpp::IntegerMatrix sampleRandomRows(Rcpp::IntegerMatrix mat, int nrows, int n_samples);
RcppExport SEXP _SRSim_sampleRandomRows(SEXP matSEXP, SEXP nrowsSEXP, SEXP n_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleRandomRows(mat, nrows, n_samples));
    return rcpp_result_gen;
END_RCPP
}
// returnRandomDrugEventPairsRcpp
Rcpp::IntegerMatrix returnRandomDrugEventPairsRcpp(Rcpp::NumericVector prob_drugs, Rcpp::NumericVector prob_events, int n_wanted_pairs, double rho);
RcppExport SEXP _SRSim_returnRandomDrugEventPairsRcpp(SEXP prob_drugsSEXP, SEXP prob_eventsSEXP, SEXP n_wanted_pairsSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_drugs(prob_drugsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_events(prob_eventsSEXP);
    Rcpp::traits::input_parameter< int >::type n_wanted_pairs(n_wanted_pairsSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(returnRandomDrugEventPairsRcpp(prob_drugs, prob_events, n_wanted_pairs, rho));
    return rcpp_result_gen;
END_RCPP
}
// returnRandomPairsRcpp
Rcpp::IntegerMatrix returnRandomPairsRcpp(Rcpp::NumericVector margprob, int n_wanted_pairs, double rho);
RcppExport SEXP _SRSim_returnRandomPairsRcpp(SEXP margprobSEXP, SEXP n_wanted_pairsSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type margprob(margprobSEXP);
    Rcpp::traits::input_parameter< int >::type n_wanted_pairs(n_wanted_pairsSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(returnRandomPairsRcpp(margprob, n_wanted_pairs, rho));
    return rcpp_result_gen;
END_RCPP
}
// generateCorrelationMatrixRcpp
Rcpp::NumericMatrix generateCorrelationMatrixRcpp(Rcpp::NumericVector prob_drugs, Rcpp::NumericVector prob_events, int n_correlated_drugs, double rho_drugs, int n_correlated_events, double rho_events, int n_correlated_pairs, double rho);
RcppExport SEXP _SRSim_generateCorrelationMatrixRcpp(SEXP prob_drugsSEXP, SEXP prob_eventsSEXP, SEXP n_correlated_drugsSEXP, SEXP rho_drugsSEXP, SEXP n_correlated_eventsSEXP, SEXP rho_eventsSEXP, SEXP n_correlated_pairsSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_drugs(prob_drugsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_events(prob_eventsSEXP);
    Rcpp::traits::input_parameter< int >::type n_correlated_drugs(n_correlated_drugsSEXP);
    Rcpp::traits::input_parameter< double >::type rho_drugs(rho_drugsSEXP);
    Rcpp::traits::input_parameter< int >::type n_correlated_events(n_correlated_eventsSEXP);
    Rcpp::traits::input_parameter< double >::type rho_events(rho_eventsSEXP);
    Rcpp::traits::input_parameter< int >::type n_correlated_pairs(n_correlated_pairsSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(generateCorrelationMatrixRcpp(prob_drugs, prob_events, n_correlated_drugs, rho_drugs, n_correlated_events, rho_events, n_correlated_pairs, rho));
    return rcpp_result_gen;
END_RCPP
}
// corr2cov
Rcpp::NumericMatrix corr2cov(Rcpp::NumericMatrix corrmat, Rcpp::NumericVector margprob);
RcppExport SEXP _SRSim_corr2cov(SEXP corrmatSEXP, SEXP margprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type corrmat(corrmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type margprob(margprobSEXP);
    rcpp_result_gen = Rcpp::wrap(corr2cov(corrmat, margprob));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat L, int ncols);
RcppExport SEXP _SRSim_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP LSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, L, ncols));
    return rcpp_result_gen;
END_RCPP
}
// threshold
Rcpp::LogicalMatrix threshold(arma::mat raw);
RcppExport SEXP _SRSim_threshold(SEXP rawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type raw(rawSEXP);
    rcpp_result_gen = Rcpp::wrap(threshold(raw));
    return rcpp_result_gen;
END_RCPP
}
// generateReports
Rcpp::LogicalMatrix generateReports(int n_reports, Rcpp::NumericVector means, arma::mat L, bool valid_reports, int n_drugs, int n_events, bool verbose);
RcppExport SEXP _SRSim_generateReports(SEXP n_reportsSEXP, SEXP meansSEXP, SEXP LSEXP, SEXP valid_reportsSEXP, SEXP n_drugsSEXP, SEXP n_eventsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_reports(n_reportsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type means(meansSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< bool >::type valid_reports(valid_reportsSEXP);
    Rcpp::traits::input_parameter< int >::type n_drugs(n_drugsSEXP);
    Rcpp::traits::input_parameter< int >::type n_events(n_eventsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(generateReports(n_reports, means, L, valid_reports, n_drugs, n_events, verbose));
    return rcpp_result_gen;
END_RCPP
}
// create2x2TablesRcpp
Rcpp::DataFrame create2x2TablesRcpp(Rcpp::IntegerMatrix reports, Rcpp::NumericMatrix corrmat, Rcpp::NumericVector prob_drugs, Rcpp::NumericVector prob_events);
RcppExport SEXP _SRSim_create2x2TablesRcpp(SEXP reportsSEXP, SEXP corrmatSEXP, SEXP prob_drugsSEXP, SEXP prob_eventsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type reports(reportsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type corrmat(corrmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_drugs(prob_drugsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_events(prob_eventsSEXP);
    rcpp_result_gen = Rcpp::wrap(create2x2TablesRcpp(reports, corrmat, prob_drugs, prob_events));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SRSim_validReport", (DL_FUNC) &_SRSim_validReport, 3},
    {"_SRSim_validCorrelation", (DL_FUNC) &_SRSim_validCorrelation, 3},
    {"_SRSim_sampleRandomRows", (DL_FUNC) &_SRSim_sampleRandomRows, 3},
    {"_SRSim_returnRandomDrugEventPairsRcpp", (DL_FUNC) &_SRSim_returnRandomDrugEventPairsRcpp, 4},
    {"_SRSim_returnRandomPairsRcpp", (DL_FUNC) &_SRSim_returnRandomPairsRcpp, 3},
    {"_SRSim_generateCorrelationMatrixRcpp", (DL_FUNC) &_SRSim_generateCorrelationMatrixRcpp, 8},
    {"_SRSim_corr2cov", (DL_FUNC) &_SRSim_corr2cov, 2},
    {"_SRSim_mvrnormArma", (DL_FUNC) &_SRSim_mvrnormArma, 4},
    {"_SRSim_threshold", (DL_FUNC) &_SRSim_threshold, 1},
    {"_SRSim_generateReports", (DL_FUNC) &_SRSim_generateReports, 7},
    {"_SRSim_create2x2TablesRcpp", (DL_FUNC) &_SRSim_create2x2TablesRcpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SRSim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
