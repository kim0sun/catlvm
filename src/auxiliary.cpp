#include <Rcpp.h>
using namespace Rcpp;

int sample1(int n, double *prob) {
   double ran = R::runif(0, 1);
   double cp = 0;
   for (int i = 0; i < n; i ++) {
      cp += exp(prob[i]);
      if (ran < cp) return i;
   }
   return n - 1;
}

NumericVector elogdiri(NumericVector a) {
   return digamma(a) - R::digamma(sum(a));
}

NumericVector plogdiri(
      NumericVector a, NumericVector b
) {
   return lgammaf(sum(a)) - lgamma(a) +
      sum( (a - 1) * elogdiri(b) );;
}
