#ifndef DOWNWARD_H
#define DOWNWARD_H

#include <Rcpp.h>

void dnInit(
   double *alpha, double *lambda, double *pi,
   double *post, double *ll, int nobs, int nclass
);

double calcll(
   double *lambda, double *pi,
   int nobs, int nclass
);

void calclli(
   double *lambda, double *pi,
   double *ll, int nobs, int nclass
);

void dnRec(
   double *alpha, double *ualpha,
   double *lambda, double *ulambda, double *jlambda,
   int nobs, int nk, int nl, double *tau,
   double *post, double *joint, double *ll
);

#endif
