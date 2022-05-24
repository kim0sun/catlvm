#ifndef MSTEP_H
#define MSTEP_H

#include <Rcpp.h>

// M-step
void updatePi(double *pi, double *post, int nobs, int nclass);
void cumTau(
      double *joint, const double *ntau,
      int nobs, int nk, int nl
);
void updateTau(
      double *tau, double *ntau,
      int nk, int nl
);
void cumRho(
      const double *denom, const double *numer,
      int *y, int nobs, int nvar, Rcpp::IntegerVector ncat,
      int nk, double *post, const double *old_rho
);
void updateRho(
      double *rho, double *numer, double *denom,
      int nobs, int nclass, int nvar, Rcpp::IntegerVector ncat
);

#endif
