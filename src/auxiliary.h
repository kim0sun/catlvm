#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <Rcpp.h>
int sample1(int n, double *prob);
Rcpp::NumericVector elogdiri(Rcpp::NumericVector a);
Rcpp::NumericVector plogdiri(Rcpp::NumericVector a, Rcpp::NumericVector b);

#endif
