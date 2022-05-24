#include <Rcpp.h>
#include "param.h"
#include "auxiliary.h"
using namespace Rcpp;

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
List ysim(
      int nsim, int nlv, IntegerVector root, IntegerVector leaf,
      Nullable<IntegerVector> ulv, Nullable<IntegerVector> vlv,
      IntegerVector nclass, int nroot, int nleaf, int nedge, List ncat,
      IntegerVector cstr_leaf, List pi, List tau, List rho, bool print_class
) {
   List cls(nlv);

   for (int r = 0; r < nroot; r ++)
      cls[root[r]] = root_gnr(nsim, nclass[root[r]], pi[r]);

   if (ulv.isNotNull() && vlv.isNotNull()) {
      IntegerVector u = as<IntegerVector>(ulv);
      IntegerVector v = as<IntegerVector>(vlv);
      for (int d = 0; d < nedge; d ++) {
         cls[u[d]] = cls_gnr(nsim, nclass[u[d]], nclass[v[d]],
                             cls[v[d]], tau[d]);
      }
   }

   List y(nleaf);
   for (int v = 0; v < nleaf; v ++)
      y[v] = y_gnr(nsim, nclass[leaf[v]], ncat[cstr_leaf[v]],
                   cls[leaf[v]], rho[cstr_leaf[v]]);

   if (print_class) {
      List ret;
      ret["y"] = y;
      ret["class"] = cls;
      return ret;
   }

   return y;
}
