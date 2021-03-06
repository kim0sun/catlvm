#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Rcpp.h>
Rcpp::NumericVector logit_pi(Rcpp::NumericVector pi, int nclass);
Rcpp::NumericMatrix logit_tau(Rcpp::NumericMatrix tau, int nk, int nl);
Rcpp::NumericVector logit_rho(Rcpp::NumericVector rho, Rcpp::IntegerVector ncat, int nclass);
Rcpp::NumericVector logistic_pi(double *lpi, int nclass);
Rcpp::NumericMatrix logistic_tau(double *ltau, int nk, int nl);
Rcpp::NumericVector logistic_rho(double *lrho, int nclass, Rcpp::IntegerVector ncat);
Rcpp::List logit2log(
      Rcpp::NumericVector logit, Rcpp::List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
      Rcpp::IntegerVector nclass_root, Rcpp::IntegerVector nclass_leaf,
      Rcpp::IntegerVector nclass_u, Rcpp::IntegerVector nclass_v
);
Rcpp::List splitSE(
      Rcpp::NumericVector se, Rcpp::List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
      Rcpp::IntegerVector nclass, Rcpp::IntegerVector nclass_u,
      Rcpp::IntegerVector nclass_v, Rcpp::IntegerVector nclass_leaf
);

#endif
