// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ysim
List ysim(int nsim, int nlv, IntegerVector root, IntegerVector leaf, Nullable<IntegerVector> ulv, Nullable<IntegerVector> vlv, IntegerVector nclass, int nroot, int nleaf, int nedge, List ncat, IntegerVector cstr_leaf, List pi, List tau, List rho, bool print_class);
RcppExport SEXP _catlvm_ysim(SEXP nsimSEXP, SEXP nlvSEXP, SEXP rootSEXP, SEXP leafSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclassSEXP, SEXP nrootSEXP, SEXP nleafSEXP, SEXP nedgeSEXP, SEXP ncatSEXP, SEXP cstr_leafSEXP, SEXP piSEXP, SEXP tauSEXP, SEXP rhoSEXP, SEXP print_classSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< List >::type pi(piSEXP);
    Rcpp::traits::input_parameter< List >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< List >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type print_class(print_classSEXP);
    rcpp_result_gen = Rcpp::wrap(ysim(nsim, nlv, root, leaf, ulv, vlv, nclass, nroot, nleaf, nedge, ncat, cstr_leaf, pi, tau, rho, print_class));
    return rcpp_result_gen;
END_RCPP
}
// calcfreq
List calcfreq(IntegerVector mis, IntegerVector nrep, int nmis, IntegerVector freq, IntegerVector xobs, int nc, int N, double tol, int max_iter);
RcppExport SEXP _catlvm_calcfreq(SEXP misSEXP, SEXP nrepSEXP, SEXP nmisSEXP, SEXP freqSEXP, SEXP xobsSEXP, SEXP ncSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type mis(misSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type nmis(nmisSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(calcfreq(mis, nrep, nmis, freq, xobs, nc, N, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// par_gnr
List par_gnr(int nobs, IntegerVector nvar, List ncat, int nroot, int nedge, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass, IntegerVector nclass_leaf, LogicalVector init, List init_param);
RcppExport SEXP _catlvm_par_gnr(SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP initSEXP, SEXP init_paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< List >::type init_param(init_paramSEXP);
    rcpp_result_gen = Rcpp::wrap(par_gnr(nobs, nvar, ncat, nroot, nedge, nleaf_unique, root, ulv, vlv, nclass, nclass_leaf, init, init_param));
    return rcpp_result_gen;
END_RCPP
}
// emFit
List emFit(IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nedge, int nleaf, int nleaf_unique, IntegerVector root, IntegerVector tree_index, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf, LogicalVector init, List init_param, int max_iter, double tol, bool verbose, int periter);
RcppExport SEXP _catlvm_emFit(SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleafSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP tree_indexSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP initSEXP, SEXP init_paramSEXP, SEXP max_iterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP periterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tree_index(tree_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< List >::type init_param(init_paramSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type periter(periterSEXP);
    rcpp_result_gen = Rcpp::wrap(emFit(y, nobs, nvar, ncat, nlv, nroot, nedge, nleaf, nleaf_unique, root, tree_index, ulv, vlv, leaf, cstr_leaf, nclass, nclass_leaf, init, init_param, max_iter, tol, verbose, periter));
    return rcpp_result_gen;
END_RCPP
}
// floglik
double floglik(NumericVector param, IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nedge, int nleaf, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf);
RcppExport SEXP _catlvm_floglik(SEXP paramSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleafSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(floglik(param, y, nobs, nvar, ncat, nlv, nroot, nedge, nleaf, nleaf_unique, root, ulv, vlv, leaf, cstr_leaf, nclass, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}
// calclli
NumericVector calclli(NumericVector param, IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nedge, int nleaf, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf);
RcppExport SEXP _catlvm_calclli(SEXP paramSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleafSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(calclli(param, y, nobs, nvar, ncat, nlv, nroot, nedge, nleaf, nleaf_unique, root, ulv, vlv, leaf, cstr_leaf, nclass, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}
// calcPost
List calcPost(List param, IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nedge, int nleaf, int nleaf_unique, IntegerVector root, IntegerVector tree_index, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf);
RcppExport SEXP _catlvm_calcPost(SEXP paramSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleafSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP tree_indexSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tree_index(tree_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPost(param, y, nobs, nvar, ncat, nlv, nroot, nedge, nleaf, nleaf_unique, root, tree_index, ulv, vlv, leaf, cstr_leaf, nclass, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}
// logit2log
List logit2log(NumericVector param, int nobs, List ncat, int nroot, int nedge, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass, IntegerVector nclass_leaf);
RcppExport SEXP _catlvm_logit2log(SEXP paramSEXP, SEXP nobsSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(logit2log(param, nobs, ncat, nroot, nedge, nleaf_unique, root, ulv, vlv, nclass, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}
// splitSE
List splitSE(NumericVector se, List ncat, int nroot, int nedge, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass, IntegerVector nclass_leaf);
RcppExport SEXP _catlvm_splitSE(SEXP seSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nedgeSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type se(seSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(splitSE(se, ncat, nroot, nedge, nleaf_unique, root, ulv, vlv, nclass, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_catlvm_ysim", (DL_FUNC) &_catlvm_ysim, 16},
    {"_catlvm_calcfreq", (DL_FUNC) &_catlvm_calcfreq, 9},
    {"_catlvm_par_gnr", (DL_FUNC) &_catlvm_par_gnr, 13},
    {"_catlvm_emFit", (DL_FUNC) &_catlvm_emFit, 23},
    {"_catlvm_floglik", (DL_FUNC) &_catlvm_floglik, 17},
    {"_catlvm_calclli", (DL_FUNC) &_catlvm_calclli, 17},
    {"_catlvm_calcPost", (DL_FUNC) &_catlvm_calcPost, 18},
    {"_catlvm_logit2log", (DL_FUNC) &_catlvm_logit2log, 11},
    {"_catlvm_splitSE", (DL_FUNC) &_catlvm_splitSE, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_catlvm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
