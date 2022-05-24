#ifndef UPWARD_H
#define UPWARD_H

#include <Rcpp.h>
void upInit(
   int *y, const double *ptr_rho, double *lambda,
   int nk, int nobs, int nvar, Rcpp::IntegerVector ncat
);

void upRec(
   double *lambda, double *jlambda, double *llambda,
   const double *tau, int nobs, int nk, int nl
);

void upRec2(
   double *lambda, double *llambda, const double *tau,
   int nobs, int nk, int nl
);

#endif
