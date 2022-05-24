#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List calcfreq(
      IntegerVector mis, IntegerVector nrep,
      int nmis, IntegerVector freq,
      IntegerVector xobs, int nc, int N,
      double tol, int max_iter
) {
   const int *mis_begin = mis.begin();
   NumericVector theta(nc, 1 / (double) nc);
   NumericVector x(nc);
   double beta = 0;
   double loglik = 0;

   double diff = R_PosInf;
   int iter = 0;
   while (diff > tol || iter < max_iter) {
      iter ++;
      loglik = 0;
      for (int c = 0; c < nc; c ++) x[c] = xobs[c];
      int *mis_ = (int*) mis_begin;
      for (int i = 0; i < nmis; i ++) {
         beta = 0;
         for (int j = 0; j < nrep[i]; j ++)
            beta += theta[mis_[j]];
         loglik += freq[i] * log(beta);
         for (int j = 0; j < nrep[i]; j ++)
            x[mis_[j]] += freq[i] * theta[mis_[j]] / beta;
         mis_ += nrep[i];
      }
      diff = max(abs(x / N - theta));
      theta = x / N;
   }

   List res;
   res["freq"] = x;
   res["theta"] = theta;
   res["loglik"] = loglik;

   return res;
}

