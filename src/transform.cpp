#include <Rcpp.h>
using namespace Rcpp;


void logit_pi(double *logit, NumericVector pi, int nclass) {
   for (int k = 0; k < nclass - 1; k ++)
      logit[k] = pi[k] - pi[nclass - 1];
}

void logit_tau(double *logit, NumericMatrix tau, int nk, int nl) {
   double *tau_ = tau.begin();

   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk - 1; k ++) {
         logit[k] = tau_[k] - tau_[nk - 1];
      }
      logit += nk - 1;
      tau_  += nk;
   }
}

void logit_rho(
   double *logit, NumericVector rho, IntegerVector ncat, int nclass
) {
   double *rho_ = rho.begin();

   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         for (int r = 0; r < ncat[m] - 1; r ++) {
            logit[r] = rho_[r] - rho_[ncat[m] - 1];
         }
         logit += ncat[m] - 1;
         rho_  += ncat[m];
      }
   }
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
NumericVector log2logit(
   List param, int nparam, List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   IntegerVector nclass_root, IntegerVector nclass_leaf,
   IntegerVector nclass_u, IntegerVector nclass_v
) {
   List lst_pi = param["pi"];
   List lst_tau = param["tau"];
   List lst_rho = param["rho"];

   NumericVector logit(nparam);
   double *logit_ = logit.begin();

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass_root[r];
      logit_pi(logit_, lst_pi[r], nk);
      logit_ += nk - 1;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d];
      int nl = nclass_v[d];
      logit_tau(logit_, lst_tau[d], nk, nl);
      logit_ += nl * (nk - 1);
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      int nk = nclass_leaf[v];
      logit_rho(logit_, lst_rho[v], ncat[v], nk);
      logit_ += nk * sum(ncatv - 1);
   }

   return logit;
}

// [[Rcpp::export]]
List logit2log(
      NumericVector param, int nobs, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass_root, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *param_ = param.begin();

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass_root[r];
      lst_pi[r] = logistic_pi(param_, nk);
      param_ += nk;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d];
      int nl = nclass_v[d];
      lst_tau[d] = logistic_tau(param_, nk, nl);
      param_ += nl * (nk - 1);
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      int nk = nclass_leaf[v];
      lst_rho[v] = logistic_rho(param_, nk, ncatv);
      param_ += nk * sum(ncatv - 1);
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}


// [[Rcpp::export]]
List splitlogit(
      NumericVector logit, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass_root, IntegerVector nclass_u,
      IntegerVector nclass_v, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *logit_ = logit.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi(nclass_root[r] - 1);
      for (int k = 0; k < nclass_root[r] - 1; k ++) {
         pi[k] = logit_[k];
      }
      logit_ += nclass_root[r] - 1;
      lst_pi[r] = pi;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau(nclass_u[d] - 1, nclass_v[d]);
      for (int l = 0; l < nclass_v[d]; l ++) {
         for (int k = 0; k < nclass_u[d] - 1; k ++) {
            tau(k, l) = logit_[k];
         }
         logit_ += nclass_u[d] - 1;
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
               rho_[r] = logit_[r];
            }
            logit_ += ncatv[m] - 1;
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


// [[Rcpp::export]]
List splitparam(
      NumericVector param, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass_root, IntegerVector nclass_u,
      IntegerVector nclass_v, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *param_ = param.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi(nclass_root[r]);
      for (int k = 0; k < nclass_root[r]; k ++) {
         pi[k] = param_[k];
      }
      param_ += nclass_root[r];
      lst_pi[r] = pi;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau(nclass_u[d], nclass_v[d]);
      for (int l = 0; l < nclass_v[d]; l ++) {
         for (int k = 0; k < nclass_u[d]; k ++) {
            tau(k, l) = param_[k];
         }
         param_ += nclass_u[d] - 1;
      }
      lst_tau[d] = tau;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho(nclass_leaf[v] * sum(ncatv));
      double *rho_ = rho.begin();
      for (int k = 0; k < nclass_leaf[v]; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            for (int r = 0; r < ncatv[m]; r ++) {
               rho_[r] = param_[r];
            }
            param_ += ncatv[m];
            rho_ += ncatv[m];
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
