#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <Rcpp.h>
int sample1(int n, double *prob);
Rcpp::NumericVector logit_pi(Rcpp::NumericVector pi, int nclass);
Rcpp::NumericMatrix logit_tau(Rcpp::NumericMatrix tau, int nk, int nl);
Rcpp::NumericVector logistic_pi(double *lpi, int nclass);
Rcpp::NumericMatrix logistic_tau(double *ltau, int nk, int nl);
Rcpp::NumericVector logistic_rho(double *lrho, int nclass, Rcpp::IntegerVector ncat);
Rcpp::NumericVector elogdiri(Rcpp::NumericVector a);
Rcpp::NumericVector plogdiri(Rcpp::NumericVector a, Rcpp::NumericVector b);
Rcpp::NumericVector pi_gnr(int nk);
Rcpp::NumericMatrix tau_gnr(int nk, int nl);
void prev_gnr(double *prev, double *upr, double *tau, int nobs, int nk, int nl);
Rcpp::NumericVector rho_gnr(int nk, Rcpp::IntegerVector ncat);
Rcpp::IntegerVector root_gnr(int nobs, int nk, Rcpp::Nullable<Rcpp::NumericVector> prob = R_NilValue);
Rcpp::IntegerVector cls_gnr(int nobs, int nk, int nl, Rcpp::IntegerVector v, Rcpp::Nullable<Rcpp::NumericMatrix> prob = R_NilValue);
Rcpp::IntegerMatrix y_gnr(int nobs, int nk, Rcpp::IntegerVector ncat, Rcpp::IntegerVector cls, Rcpp::Nullable<Rcpp::NumericVector> prob = R_NilValue);


#endif
