#include <Rcpp.h>
using namespace Rcpp;

void dnInit(
      double *alpha, double *lambda, double *pi,
      double *post, double *ll, int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         alpha[k] = pi[k];
         post[k] = alpha[k] + lambda[k];
         lik += exp(post[k]);
      }
      ll[i] = log(lik);

      for (int k = 0; k < nclass; k ++) {
         post[k] -= ll[i];
      }

      alpha += nclass;
      lambda += nclass;
      post += nclass;
   }
}

double calcll(
      double *lambda, double *pi,
      int nobs, int nclass
) {
   double ll = 0;
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         lik += exp(pi[k] + lambda[k]);
      }
      ll += log(lik);

      lambda += nclass;
   }
   return ll;
}

void calclli(
      double *lambda, double *pi,
      double *ll, int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         lik += exp(pi[k] + lambda[k]);
      }
      ll[i] += log(lik);

      lambda += nclass;
   }
}

void dnRec(
      double *alpha, double *ualpha,
      double *lambda, double *ulambda, double *jlambda,
      int nobs, int nk, int nl, double *tau,
      double *post, double *joint, double *ll
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         double svl = 0;
         for (int l = 0; l < nl; l ++) {
            val = tau[k + l * nk] + ualpha[l] + ulambda[l] - jlambda[l];
            joint[k + l * nk] = val + lambda[k] - ll[i];
            svl += exp(val);
         }
         alpha[k] = log(svl);
         post[k] = alpha[k] + lambda[k] - ll[i];
      }

      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      lambda += nk; ulambda += nl; jlambda += nl;
   }
}
