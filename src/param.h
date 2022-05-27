#ifndef PARAM_H
#define PARAM_H

#include <Rcpp.h>
Rcpp::NumericVector pi_gnr(int nk);
Rcpp::NumericMatrix tau_gnr(int nk, int nl);
Rcpp::NumericVector rho_gnr(int nk, Rcpp::IntegerVector ncat);
Rcpp::List par_gnr(
      int nobs, Rcpp::IntegerVector nvar, Rcpp::List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
      Rcpp::IntegerVector nclass, Rcpp::IntegerVector nclass_u,
      Rcpp::IntegerVector nclass_v, Rcpp::IntegerVector nclass_leaf,
      Rcpp::LogicalVector init, Rcpp::List init_param
);

#endif
