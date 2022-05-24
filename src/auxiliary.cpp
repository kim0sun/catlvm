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

NumericVector logit_pi(
      NumericVector pi, int nclass
) {
   NumericVector lpi(nclass - 1);

   for (int k = 0; k < nclass - 1; k ++)
      lpi[k] = pi[k] - pi[nclass - 1];

   return lpi;
}

NumericMatrix logit_tau(
      NumericMatrix tau, int nk, int nl
) {
   NumericMatrix ltau(nk - 1, nl);
   double *tau_ = tau.begin();
   double *ltau_ = ltau.begin();

   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk - 1; k ++) {
         ltau_[k] = tau_[k] - tau_[nk - 1];
      }
      ltau_ += nk - 1;
      tau_  += nk;
   }

   return ltau;
}


NumericVector logit_rho(
      NumericVector rho, IntegerVector ncat, int nclass
) {
   NumericVector lrho(nclass * sum(ncat - 1));
   double *rho_ = rho.begin();
   double *lrho_ = lrho.begin();

   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         for (int r = 0; r < ncat[m] - 1; r ++) {
            lrho_[r] = rho_[r] - rho_[ncat[m] - 1];
         }
         lrho_ += ncat[m] - 1;
         rho_  += ncat[m];
      }
   }
   return lrho;
}

NumericVector logistic_pi(double *lpi, int nclass) {
   NumericVector pi(nclass);
   double denom = 1;
   for (int i = 0; i < nclass - 1; i ++) {
      pi[i] = lpi[i];
      denom += exp(lpi[i]);
   }
   for (int i = 0; i < nclass; i ++)
      pi[i] -= log(denom);

   return pi;
}

NumericMatrix logistic_tau(double *ltau, int nk, int nl) {
   NumericMatrix tau(nk, nl);
   double *tau_ = tau.begin();
   double denom = 0;
   for (int l = 0; l < nl; l ++) {
      denom = 1;
      for (int k = 0; k < nk - 1; k ++) {
         tau_[k] = ltau[k];
         denom += exp(ltau[k]);
      }
      for (int k = 0; k < nk; k ++) {
         tau_[k] -= log(denom);
      }
      tau_ += nk;
      ltau += nk - 1;
   }

   return tau;
}

NumericVector logistic_rho(double *lrho, int nclass, IntegerVector ncat) {
   NumericVector rho(nclass * sum(ncat));
   double *rho_ = rho.begin();
   double denom = 0;
   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         denom = 1;
         for (int r = 0; r < ncat[m] - 1; r ++) {
            rho_[r] = lrho[r];
            denom += exp(lrho[r]);
         }
         for (int r = 0; r < ncat[m]; r ++) {
            rho_[r] -= log(denom);
         }
         rho_ += ncat[m];
         lrho += ncat[m] - 1;
      }
   }

   return rho;
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

NumericVector pi_gnr(int nk) {
   NumericVector pi(nk);
   pi.fill(-log(nk));
   return pi;
}

NumericMatrix tau_gnr(int nk, int nl) {
   NumericMatrix tau(nk, nl);
   double *ptau = tau.begin();

   for (int l = 0; l < nl; l ++) {
      NumericVector p = runif(nk, 0, 1);
      for (int k = 0; k < nk; k ++) {
         ptau[k] = log(p[k]) - log(sum(p));
      }
      ptau += nk;
   }

   return tau;
}

void prev_gnr(
      double *prev, double *upr, double *tau,
      int nobs, int nk, int nl
) {
   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk; k ++)
         prev[k] += exp(tau[k] + upr[l]);
      tau += nk;
   }
}


NumericVector rho_gnr(int nk, IntegerVector ncat) {
   NumericVector pv(nk * sum(ncat));
   double *ppv = pv.begin();

   for (int k = 0; k < nk; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         NumericVector p = runif(ncat[m], 0, 1);
         for (int r = 0; r < ncat[m]; r ++) {
            ppv[r] = log(p[r] / sum(p));
         }
         ppv += ncat[m];
      }
   }
   return pv;
}

IntegerVector root_gnr(
      int nobs, int nk,
      Nullable<NumericVector> prob = R_NilValue
) {
   NumericVector pi;
   if (prob.isNull()) {
      pi = pi_gnr(nk);
   } else {
      pi = as<NumericVector>(prob);
   }

   IntegerVector cls(nobs);
   for (int i = 0; i < nobs; i ++) {
      cls[i] = sample1(nk, pi.begin());
   }

   return cls;
}

IntegerVector cls_gnr(
      int nobs, int nk, int nl, IntegerVector v,
      Nullable<NumericMatrix> prob = R_NilValue
) {
   NumericMatrix tau;
   if (prob.isNull()) {
      tau = tau_gnr(nk, nl);
   } else {
      tau = as<NumericMatrix>(prob);
   }

   IntegerVector cls(nobs);
   double *ptau = tau.begin();
   for (int i = 0; i < nobs; i ++) {
      cls[i] = sample1(nk, ptau + v[i] * nk);
   }

   return cls;
}

IntegerMatrix y_gnr(
      int nobs, int nk, IntegerVector ncat,
      IntegerVector cls,
      Nullable<NumericVector> prob = R_NilValue
) {
   NumericVector rho;
   if (prob.isNull()) {
      rho = rho_gnr(nk, ncat);
   } else {
      rho = as<NumericVector>(prob);
   }

   int nvar = ncat.length();
   IntegerMatrix y(nvar, nobs);
   for (int i = 0; i < nobs; i ++) {
      double *pos = rho.begin() + cls[i] * sum(ncat);
      for (int m = 0; m < nvar; m ++) {
         y[i * nvar + m] = sample1(ncat[m], pos) + 1;
         pos += ncat[m];
      }
   }

   return y;
}
