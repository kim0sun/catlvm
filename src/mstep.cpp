#include <Rcpp.h>
using namespace Rcpp;

void updatePi(
      double *pi, double *post, int nobs, int nclass
) {
   NumericVector npi(nclass);
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++) {
         npi[k] += exp(post[k]);
      }
      post += nclass;
   }
   for (int k = 0; k < nclass; k ++) {
      pi[k] = log(npi[k]) - log(sum(npi));
   }
}

void cumTau(
      double *joint, const double *ntau,
      int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *ntau_ = (double *) ntau;
      for (int l = 0; l < nl; l ++) {
         for (int k = 0; k < nk; k ++) {
            ntau_[k] += exp(joint[k]);
         }
         joint += nk;
         ntau_ += nk;
      }
   }
}

void updateTau(
      double *tau, double *ntau,
      int nk, int nl
) {
   for (int l = 0; l < nl; l ++) {
      double stau = 0;
      for (int k = 0; k < nk; k ++) {
         stau += ntau[k];
      }
      for (int k = 0; k < nk; k ++) {
         tau[k] = log(ntau[k]) - log(stau);
      }
      ntau += nk;
      tau  += nk;
   }
}


void cumRho(
      const double *denom, const double *numer,
      int *y, int nobs, int nvar, IntegerVector ncat,
      int nk, double *post, const double *old_rho
) {
   for (int i = 0; i < nobs; i ++) {
      double *dnm = (double *) denom;
      double *nmr = (double *) numer;
      double *rho = (double *) old_rho;
      for (int k = 0; k < nk; k ++) {
         for (int m = 0; m < nvar; m ++) {
            dnm[m] += exp(post[k]);
            if (y[m] > 0) nmr[y[m] - 1] += exp(post[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  nmr[r] += exp(post[k] + rho[r]);
               }
            }
            nmr += ncat[m];
            rho += ncat[m];
         }
         dnm += nvar;
      }
      post += nk;
      y += nvar;
   }
}

void updateRho(
      double *rho, double *numer, double *denom,
      int nobs, int nclass, int nvar, IntegerVector ncat
) {
   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < nvar; m ++) {
         for (int r = 0; r < ncat[m]; r ++) {
            rho[r] = log(numer[r]) - log(denom[m]);
         }
         rho   += ncat[m];
         numer += ncat[m];
      }
      denom += nvar;
   }
}

