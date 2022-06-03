#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logit_pi(NumericVector pi, int nclass) {
   NumericVector lpi(nclass - 1);

   for (int k = 0; k < nclass - 1; k ++)
      lpi[k] = pi[k] - pi[nclass - 1];

   return lpi;
}

// [[Rcpp::export]]
NumericMatrix logit_tau(NumericMatrix tau, int nk, int nl) {
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

// [[Rcpp::export]]
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


// [[Rcpp::export]]
List logit2log(
      NumericVector param, int nobs, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *param_ = param.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi = logistic_pi(param_, nclass[root[r]]);
      lst_pi[r] = pi;
      param_ += nclass[root[r]] - 1;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau = logistic_tau(param_, nclass_u[d], nclass_v[d]);
      lst_tau[d] = tau;
      param_ += nclass_u[d] * (nclass_v[d] - 1);
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho = logistic_rho(param_, nclass_leaf[v], ncatv);
      lst_rho[v] = rho;
      param_ += nclass_leaf[v] * sum(ncatv - 1);
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}

// [[Rcpp::export]]
List splitSE(
      NumericVector se, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass, IntegerVector nclass_u,
      IntegerVector nclass_v, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *se_ = se.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi(nclass[root[r]] - 1);
      for (int k = 0; k < nclass[root[r]] - 1; k ++) {
         pi[k] = se_[k];
      }
      se_ += nclass[root[r]] - 1;
      lst_pi[r] = pi;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau(nclass_u[d] - 1, nclass_v[d]);
      for (int l = 0; l < nclass_v[d]; l ++) {
         for (int k = 0; k < nclass_u[d] - 1; k ++) {
            tau(k, l) = se_[k];
         }
         se_ += nclass_u[d] - 1;
      }
      lst_tau[d] = tau;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho(nclass_leaf[v] * sum(ncatv - 1));
      double *rho_ = rho.begin();
      for (int k = 0; k < nclass_leaf[v]; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            for (int r = 0; r < ncatv[m] - 1; r ++) {
               rho_[r] = se_[r];
            }
            se_ += ncatv[m] - 1;
            rho_ += ncatv[m] - 1;
         }
      }
      lst_rho[v] = rho;
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}
